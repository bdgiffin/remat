import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import minimize


THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
sys.path.append(str(REPO_ROOT / "install" / "package"))

import REMAT

# Matplotlib style rules from AGENTS.md
plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)

# -----------------------------------------------------------------------------
# Inversion configuration (single high-fidelity profile)
# -----------------------------------------------------------------------------

# Geometry and mesh
WIDTH = 15.0
HEIGHT = 3.0
NX = int(WIDTH * 10)
NY = int(HEIGHT * 10)

# Time integration
DT = 4.0e-3
N_STEPS = 1000
N_SUB_STEPS = 1
INTEGRATOR_TYPE = "fixed_visco"

# Excitation and sensors
IMPACT_VELOCITY = 1.50
IMPACT_WINDOW_WIDTH = 1.3
N_SENSORS = 5
SENSOR_DISTRIBUTION_WIDTH = 4.7

# Material model constants (fixed)
DENSITY = 1.0
YOUNGS_MODULUS = 5.0
POISSONS_RATIO = 0.28
SHEAR_MODULUS_MAXWELL = 2.0
MASS_DAMPING_FACTOR = 0.0
OVERFLOW_LIMIT = 10.0
MAT_OVERFLOW_LIMIT = 10.0

# Unknown bounds
STIFFNESS_MIN = 1.0
STIFFNESS_MAX = 20.0
TAU_MIN = 0.1
TAU_MAX = 2.0

# Initial guesses
INIT_STIFFNESS = 10.0
INIT_TAU = 1.0

# Regularization (applied to both unknown fields)
REG_L2_STIFFNESS = 7.5e-4
REG_TV_STIFFNESS = 3.5e-4
REG_L2_TAU = 6.0e-3
REG_TV_TAU = 2.5e-3
REG_TV_EPS = 1.0e-6

# Optimization
MAX_ITERS = 60
LBFGSB_GTOL = 1.0e-8
LBFGSB_MAXLS = 40

# Material discovery (post-hoc)
K_MIN = 1
K_MAX = 2
KMEANS_MAX_ITERS = 10
KMEANS_RESTARTS = 5
MODEL_SELECTION_PENALTY = 2.0
RANDOM_SEED = 7

# Output
OUTPUT_DIR = THIS_DIR / "inverse_outputs_elementwise"
OUTPUT_JSON = OUTPUT_DIR / "inverse_result.json"
OUTPUT_PLOT = OUTPUT_DIR / "inverse_summary.svg"

PLOT_COLOR_PRIMARY = "#2b738eff"
PLOT_COLOR_SECONDARY = "#f9826bff"

# -----------------------------------------------------------------------------
# Spatial callbacks consumed by REMAT
# -----------------------------------------------------------------------------

_ACTIVE_STIFFNESS = np.ones(1, dtype=np.double)
_ACTIVE_TAU = np.ones(1, dtype=np.double)
_ELEM_NX = 1
_ELEM_NY = 1
_ELEM_DX = 1.0
_ELEM_DY = 1.0
_ELEM_WIDTH = 1.0
_ELEM_HEIGHT = 1.0


def _element_id_from_xy(x, y):
    x_clamped = min(max(float(x), 0.0), _ELEM_WIDTH - 1.0e-12)
    y_clamped = min(max(float(y), 0.0), _ELEM_HEIGHT - 1.0e-12)
    ix = int(np.clip(np.floor(x_clamped / _ELEM_DX), 0, _ELEM_NX - 1))
    iy = int(np.clip(np.floor(y_clamped / _ELEM_DY), 0, _ELEM_NY - 1))
    return iy * _ELEM_NX + ix


def elementwise_stiffness_scaling(x, y):
    return float(_ACTIVE_STIFFNESS[_element_id_from_xy(x, y)])


def elementwise_relaxation_time(x, y):
    return float(_ACTIVE_TAU[_element_id_from_xy(x, y)])


# -----------------------------------------------------------------------------
# Problem setup
# -----------------------------------------------------------------------------


def make_structured_quad_problem():
    xs = np.linspace(0.0, WIDTH, NX + 1)
    ys = np.linspace(0.0, HEIGHT, NY + 1)

    num_nodes = (NX + 1) * (NY + 1)
    coordinates = np.zeros((num_nodes, 2), dtype=np.double)
    velocities = np.zeros((num_nodes, 2), dtype=np.double)
    fixity = np.zeros((num_nodes, 2), dtype=np.bool_)

    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            nid = j * (NX + 1) + i
            coordinates[nid, 0] = x
            coordinates[nid, 1] = y

    num_elems = NX * NY
    connectivity = np.zeros((num_elems, 4), dtype=np.int32)
    e = 0
    for j in range(NY):
        for i in range(NX):
            n0 = j * (NX + 1) + i
            n1 = n0 + 1
            n3 = (j + 1) * (NX + 1) + i
            n2 = n3 + 1
            connectivity[e, :] = [n0, n1, n2, n3]
            e += 1

    eps = 1.0e-12
    node_x = coordinates[:, 0]
    node_y = coordinates[:, 1]

    bottom_nodes = np.where(np.abs(node_y - 0.0) < eps)[0]
    top_nodes = np.where(np.abs(node_y - HEIGHT) < eps)[0]
    fixity[bottom_nodes, :] = True

    xmid = 0.5 * WIDTH
    impact_span = min(IMPACT_WINDOW_WIDTH, WIDTH)
    impact_half_width = 0.5 * impact_span
    impact_lo = xmid - impact_half_width
    impact_hi = xmid + impact_half_width
    impact_nodes = top_nodes[(node_x[top_nodes] >= impact_lo - eps) & (node_x[top_nodes] <= impact_hi + eps)]
    if impact_nodes.size == 0:
        impact_nodes = top_nodes[np.array([int(np.argmin(np.abs(node_x[top_nodes] - xmid)))], dtype=np.int32)]
    velocities[impact_nodes, 1] = -abs(IMPACT_VELOCITY)

    top_interior = top_nodes[(node_x[top_nodes] > eps) & (node_x[top_nodes] < WIDTH - eps)]
    sensor_span = min(SENSOR_DISTRIBUTION_WIDTH, WIDTH)
    sx_lo = xmid - 0.5 * sensor_span
    sx_hi = xmid + 0.5 * sensor_span
    sensor_candidates = top_interior[(node_x[top_interior] >= sx_lo - eps) & (node_x[top_interior] <= sx_hi + eps)]

    nsensors = min(N_SENSORS, sensor_candidates.size)
    sensor_pick = np.unique(np.round(np.linspace(0, sensor_candidates.size - 1, nsensors)).astype(int))
    sensor_nodes = sensor_candidates[sensor_pick]

    elem_centers = np.mean(coordinates[connectivity, :], axis=1)

    return {
        "coordinates": coordinates,
        "velocities": velocities,
        "fixity": fixity,
        "connectivity": connectivity,
        "sensor_nodes": sensor_nodes.astype(np.int32),
        "elem_centers": elem_centers,
        "nx": NX,
        "ny": NY,
        "width": WIDTH,
        "height": HEIGHT,
    }


# -----------------------------------------------------------------------------
# Forward / adjoint solve
# -----------------------------------------------------------------------------


def configure_run(problem, stiffness_elem, tau_elem):
    global _ACTIVE_STIFFNESS, _ACTIVE_TAU, _ELEM_NX, _ELEM_NY, _ELEM_DX, _ELEM_DY, _ELEM_WIDTH, _ELEM_HEIGHT

    _ACTIVE_STIFFNESS = np.asarray(stiffness_elem, dtype=np.double).copy()
    _ACTIVE_TAU = np.asarray(tau_elem, dtype=np.double).copy()

    _ELEM_NX = int(problem["nx"])
    _ELEM_NY = int(problem["ny"])
    _ELEM_WIDTH = float(problem["width"])
    _ELEM_HEIGHT = float(problem["height"])
    _ELEM_DX = _ELEM_WIDTH / float(_ELEM_NX)
    _ELEM_DY = _ELEM_HEIGHT / float(_ELEM_NY)

    REMAT.API.set_integrator_type(INTEGRATOR_TYPE.encode("utf-8"))
    REMAT.API.define_parameter(b"body_force_x", 0.0)
    REMAT.API.define_parameter(b"body_force_y", 0.0)
    REMAT.API.define_parameter(b"mass_damping_factor", MASS_DAMPING_FACTOR)
    REMAT.API.define_parameter(b"contact_stiffness", 0.0)
    REMAT.API.define_parameter(b"overflow_limit", OVERFLOW_LIMIT)
    REMAT.API.define_parameter(b"mat_overflow_limit", MAT_OVERFLOW_LIMIT)
    REMAT.API.define_parameter(b"adjoint_material_objective_weight", 0.0)

    REMAT.API.define_parameter(b"density", DENSITY)
    REMAT.API.define_parameter(b"youngs_modulus", YOUNGS_MODULUS)
    REMAT.API.define_parameter(b"poissons_ratio", POISSONS_RATIO)
    REMAT.API.define_parameter(b"shear_modulus_Maxwell_element", SHEAR_MODULUS_MAXWELL)

    # Fallback scalar tau for backward compatibility in the material constructor.
    REMAT.API.define_parameter(b"relaxation_time", float(np.mean(_ACTIVE_TAU)))

    truss_connectivity = np.zeros((0, 2), dtype=np.int32)
    REMAT.create_geometry(
        problem["coordinates"],
        problem["velocities"],
        problem["fixity"],
        problem["connectivity"],
        [],
        truss_connectivity,
    )
    REMAT.define_variable_properties(elementwise_stiffness_scaling)
    REMAT.define_variable_relaxation_time(elementwise_relaxation_time)
    REMAT.API.initialize()


def run_forward_or_adjoint(problem, stiffness_elem, tau_elem, observed_history=None, compute_gradients=False):
    configure_run(problem, stiffness_elem, tau_elem)

    sensor_nodes = problem["sensor_nodes"]
    nsensors = sensor_nodes.size
    nelem = problem["connectivity"].shape[0]

    sensor_history = np.zeros((N_STEPS, nsensors), dtype=np.double)
    data_loss = 0.0

    for k in range(N_STEPS):
        REMAT.API.update_state(DT, N_SUB_STEPS, REMAT.PASS_FORWARD)
        vy = np.asarray(REMAT.get_field(b"node", "velocity_Y"), dtype=np.double)[sensor_nodes]
        sensor_history[k, :] = vy
        if observed_history is not None:
            residual = vy - observed_history[k, :]
            data_loss += 0.5 * float(np.dot(residual, residual))

    grad_stiff = np.zeros(nelem, dtype=np.double)
    grad_tau = np.zeros(nelem, dtype=np.double)

    if compute_gradients:
        if observed_history is None:
            raise ValueError("observed_history is required when compute_gradients=True")

        REMAT.clear_adjoint_state()
        for rev in range(N_STEPS):
            k = N_STEPS - 1 - rev
            residual = sensor_history[k, :] - observed_history[k, :]
            seed_xy = np.zeros((nsensors, 2), dtype=np.double)
            seed_xy[:, 1] = residual
            REMAT.add_nodal_velocity_adjoint_seed(sensor_nodes, seed_xy)
            REMAT.API.update_state(DT, N_SUB_STEPS, REMAT.PASS_BACKWARD_ADJOINT)

        grad_stiff = np.asarray(REMAT.get_field(b"element", "dparam_stiffness_scaling_factor"), dtype=np.double)
        grad_tau = np.asarray(REMAT.get_field(b"element", "dparam_relaxation_time"), dtype=np.double)

    return {
        "data_loss": float(data_loss),
        "sensor_history": sensor_history,
        "grad_stiff": grad_stiff,
        "grad_tau": grad_tau,
    }


# -----------------------------------------------------------------------------
# Regularization
# -----------------------------------------------------------------------------


def build_edge_pairs(nx, ny):
    i_idx = []
    j_idx = []

    for j in range(ny):
        row_start = j * nx
        for i in range(nx - 1):
            a = row_start + i
            b = a + 1
            i_idx.append(a)
            j_idx.append(b)

    for j in range(ny - 1):
        row_start = j * nx
        next_row_start = (j + 1) * nx
        for i in range(nx):
            a = row_start + i
            b = next_row_start + i
            i_idx.append(a)
            j_idx.append(b)

    return np.asarray(i_idx, dtype=np.int32), np.asarray(j_idx, dtype=np.int32)


def regularization_loss_and_grad(field, i_idx, j_idx, l2_weight, tv_weight, tv_eps):
    diff = field[j_idx] - field[i_idx]

    l2_loss = 0.5 * l2_weight * float(np.sum(diff * diff))
    tv_norm = np.sqrt(diff * diff + tv_eps * tv_eps)
    tv_loss = tv_weight * float(np.sum(tv_norm))

    edge_grad = l2_weight * diff + tv_weight * (diff / tv_norm)
    grad = np.zeros_like(field)
    np.add.at(grad, i_idx, -edge_grad)
    np.add.at(grad, j_idx, +edge_grad)

    return l2_loss + tv_loss, grad


# -----------------------------------------------------------------------------
# Synthetic truth and post-hoc material discovery
# -----------------------------------------------------------------------------


def make_true_fields(problem):
    y = problem["elem_centers"][:, 1]

    # Simpler synthetic target: two horizontal stiffness regions.
    labels = np.zeros(y.size, dtype=np.int32)
    labels[y >= 0.5 * HEIGHT] = 1

    stiffness_by_material = np.array([15.0, 7.7], dtype=np.double)

    stiffness_true = stiffness_by_material[labels]
    # Keep tau consistent with fixed inversion bounds (TAU_MIN == TAU_MAX == 0.1).
    tau_true = np.full(y.size, TAU_MIN, dtype=np.double)

    return stiffness_true, tau_true, labels


def _kmeans_pp_init(features, k, rng):
    n, d = features.shape
    centers = np.zeros((k, d), dtype=np.double)

    first = int(rng.integers(0, n))
    centers[0] = features[first]
    min_sq = np.sum((features - centers[0]) ** 2, axis=1)

    for c in range(1, k):
        denom = float(np.sum(min_sq))
        if denom <= 0.0:
            centers[c] = features[int(rng.integers(0, n))]
        else:
            probs = min_sq / denom
            idx = int(rng.choice(n, p=probs))
            centers[c] = features[idx]
        sq = np.sum((features - centers[c]) ** 2, axis=1)
        min_sq = np.minimum(min_sq, sq)

    return centers


def run_kmeans(features, k, rng_seed):
    rng = np.random.default_rng(rng_seed)
    n = features.shape[0]

    centers = _kmeans_pp_init(features, k, rng)
    labels = np.zeros(n, dtype=np.int32)

    for _ in range(KMEANS_MAX_ITERS):
        sq_dist = np.sum((features[:, None, :] - centers[None, :, :]) ** 2, axis=2)
        new_labels = np.argmin(sq_dist, axis=1).astype(np.int32)

        if np.array_equal(new_labels, labels):
            break
        labels = new_labels

        for c in range(k):
            members = features[labels == c]
            if members.size == 0:
                farthest_idx = int(np.argmax(np.min(sq_dist, axis=1)))
                centers[c] = features[farthest_idx]
            else:
                centers[c] = np.mean(members, axis=0)

    sq_dist = np.sum((features[:, None, :] - centers[None, :, :]) ** 2, axis=2)
    inertia = float(np.sum(sq_dist[np.arange(n), labels]))
    return centers, labels, inertia


def discover_materials(stiffness, tau):
    features = np.column_stack([np.log(stiffness), np.log(tau)])
    n, d = features.shape

    best = None
    for k in range(K_MIN, K_MAX + 1):
        best_k = None
        for r in range(KMEANS_RESTARTS):
            centers, labels, inertia = run_kmeans(features, k, RANDOM_SEED + 100 * k + r)
            if best_k is None or inertia < best_k["inertia"]:
                best_k = {
                    "centers": centers,
                    "labels": labels,
                    "inertia": inertia,
                }

        variance = best_k["inertia"] / max(float(n * d), 1.0)
        score = float(n * d * np.log(variance + 1.0e-12) + MODEL_SELECTION_PENALTY * k * np.log(max(n, 2)))

        candidate = {
            "k": k,
            "score": score,
            "inertia": best_k["inertia"],
            "centers_log": best_k["centers"],
            "labels": best_k["labels"],
        }
        if best is None or candidate["score"] < best["score"]:
            best = candidate

    centers_log = best["centers_log"]
    centers_phys = np.column_stack([np.exp(centers_log[:, 0]), np.exp(centers_log[:, 1])])

    return {
        "k": int(best["k"]),
        "score": float(best["score"]),
        "inertia": float(best["inertia"]),
        "centers_log": centers_log,
        "centers": centers_phys,
        "labels": best["labels"],
    }


# -----------------------------------------------------------------------------
# Inversion driver
# -----------------------------------------------------------------------------


def invert(problem, observed_history, init_stiffness, init_tau):
    if minimize is None:
        raise ImportError("SciPy is required for L-BFGS-B inversion. Install with: pip install scipy")

    nelem = problem["connectivity"].shape[0]
    i_idx, j_idx = build_edge_pairs(problem["nx"], problem["ny"])

    x0 = np.concatenate([init_stiffness, init_tau])
    bounds = [(STIFFNESS_MIN, STIFFNESS_MAX)] * nelem + [(TAU_MIN, TAU_MAX)] * nelem

    history = []
    cache = {
        "x": None,
        "total_loss": None,
        "data_loss": None,
        "reg_loss": None,
        "grad": None,
    }

    def objective_with_grad(x):
        x = np.asarray(x, dtype=np.double)
        stiffness = x[:nelem]
        tau = x[nelem:]

        run = run_forward_or_adjoint(problem, stiffness, tau, observed_history=observed_history, compute_gradients=True)

        reg_s_loss, reg_s_grad = regularization_loss_and_grad(
            stiffness, i_idx, j_idx, REG_L2_STIFFNESS, REG_TV_STIFFNESS, REG_TV_EPS
        )
        reg_t_loss, reg_t_grad = regularization_loss_and_grad(
            tau, i_idx, j_idx, REG_L2_TAU, REG_TV_TAU, REG_TV_EPS
        )

        reg_loss = reg_s_loss + reg_t_loss
        total_loss = run["data_loss"] + reg_loss

        grad_stiff = run["grad_stiff"] + reg_s_grad
        grad_tau = run["grad_tau"] + reg_t_grad
        grad = np.concatenate([grad_stiff, grad_tau])

        cache["x"] = x.copy()
        cache["total_loss"] = float(total_loss)
        cache["data_loss"] = float(run["data_loss"])
        cache["reg_loss"] = float(reg_loss)
        cache["grad"] = grad.copy()

        history.append(
            {
                "eval": len(history),
                "loss_total": float(total_loss),
                "loss_data": float(run["data_loss"]),
                "loss_reg": float(reg_loss),
                "grad_inf_stiff": float(np.max(np.abs(grad_stiff))),
                "grad_inf_tau": float(np.max(np.abs(grad_tau))),
            }
        )

        return cache["total_loss"], grad

    iter_counter = {"k": 0}

    def callback(xk):
        if cache["x"] is None or not np.array_equal(np.asarray(xk, dtype=np.double), cache["x"]):
            objective_with_grad(xk)
        k = iter_counter["k"]
        print(
            f"iter={k:03d}  total={cache['total_loss']:.6e}  data={cache['data_loss']:.6e}  "
            f"reg={cache['reg_loss']:.6e}"
        )
        iter_counter["k"] += 1

    # Log initial point.
    objective_with_grad(x0)
    print(
        f"iter=000  total={cache['total_loss']:.6e}  data={cache['data_loss']:.6e}  "
        f"reg={cache['reg_loss']:.6e}"
    )

    result = minimize(
        objective_with_grad,
        x0=x0,
        method="L-BFGS-B",
        jac=True,
        bounds=bounds,
        callback=callback,
        options={
            "maxiter": MAX_ITERS,
            "gtol": LBFGSB_GTOL,
            "maxls": LBFGSB_MAXLS,
        },
    )

    return result, history


# -----------------------------------------------------------------------------
# Visualization/output
# -----------------------------------------------------------------------------


def _reshape(problem, flat_field):
    return np.asarray(flat_field, dtype=np.double).reshape(problem["ny"], problem["nx"])


def save_summary_plot(
    problem,
    true_stiff,
    true_tau,
    init_stiff,
    init_tau,
    rec_stiff,
    rec_tau,
    true_labels,
    discovered_labels,
    history,
):
    fig, axes = plt.subplots(2, 4, figsize=(18, 8))

    ext = [0.0, problem["width"], 0.0, problem["height"]]

    smin = min(np.min(true_stiff), np.min(init_stiff), np.min(rec_stiff))
    smax = max(np.max(true_stiff), np.max(init_stiff), np.max(rec_stiff))
    tmin = min(np.min(true_tau), np.min(init_tau), np.min(rec_tau))
    tmax = max(np.max(true_tau), np.max(init_tau), np.max(rec_tau))

    panels = [
        (axes[0, 0], _reshape(problem, true_stiff), "True stiffness", smin, smax, "viridis"),
        (axes[0, 1], _reshape(problem, init_stiff), "Init stiffness", smin, smax, "viridis"),
        (axes[0, 2], _reshape(problem, rec_stiff), "Recovered stiffness", smin, smax, "viridis"),
        (axes[1, 0], _reshape(problem, true_tau), "True tau", tmin, tmax, "magma"),
        (axes[1, 1], _reshape(problem, init_tau), "Init tau", tmin, tmax, "magma"),
        (axes[1, 2], _reshape(problem, rec_tau), "Recovered tau", tmin, tmax, "magma"),
    ]

    for ax, field, title, vmin, vmax, cmap in panels:
        im = ax.imshow(field, origin="lower", extent=ext, aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize="medium")
        ax.set_xlabel("x", fontsize="large")
        ax.set_ylabel("y", fontsize="large")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    ax_true_labels = axes[0, 3]
    im0 = ax_true_labels.imshow(_reshape(problem, true_labels), origin="lower", extent=ext, aspect="auto", cmap="tab20")
    ax_true_labels.set_title("True material regions", fontsize="medium")
    ax_true_labels.set_xlabel("x", fontsize="large")
    ax_true_labels.set_ylabel("y", fontsize="large")
    fig.colorbar(im0, ax=ax_true_labels, fraction=0.046, pad=0.04)

    ax_disc_labels = axes[1, 3]
    im1 = ax_disc_labels.imshow(_reshape(problem, discovered_labels), origin="lower", extent=ext, aspect="auto", cmap="tab20")
    ax_disc_labels.set_title("Discovered materials", fontsize="medium")
    ax_disc_labels.set_xlabel("x", fontsize="large")
    ax_disc_labels.set_ylabel("y", fontsize="large")
    fig.colorbar(im1, ax=ax_disc_labels, fraction=0.046, pad=0.04)

    # Overlay loss history on discovered-material panel as a compact inset.
    inset = ax_disc_labels.inset_axes([0.05, 0.05, 0.55, 0.35])
    loss_vals = [row["loss_total"] for row in history]
    inset.plot(np.arange(len(loss_vals)), loss_vals, color=PLOT_COLOR_PRIMARY, linewidth=1.6)
    inset.set_title("Loss", fontsize="medium")
    inset.tick_params(labelsize=7)
    inset.grid(True, alpha=0.25)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUTPUT_PLOT, dpi=200, metadata={"Title": "Dissipative wave inverse summary"})
    fig.clf()
    plt.close(fig)


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    problem = make_structured_quad_problem()
    nelem = problem["connectivity"].shape[0]

    true_stiff, true_tau, true_labels = make_true_fields(problem)
    init_stiff = np.full(nelem, INIT_STIFFNESS, dtype=np.double)
    init_tau = np.full(nelem, INIT_TAU, dtype=np.double)

    print("Generating synthetic observations...")
    observed = run_forward_or_adjoint(problem, true_stiff, true_tau, observed_history=None, compute_gradients=False)
    observed_history = observed["sensor_history"]

    print("Running element-wise inversion (stiffness + tau)...")
    result, history = invert(problem, observed_history, init_stiff, init_tau)

    rec_stiff = np.asarray(result.x[:nelem], dtype=np.double)
    rec_tau = np.asarray(result.x[nelem:], dtype=np.double)

    discovered = discover_materials(rec_stiff, rec_tau)

    final_eval = run_forward_or_adjoint(problem, rec_stiff, rec_tau, observed_history=observed_history, compute_gradients=False)

    payload = {
        "config": {
            "mesh": {"nx": NX, "ny": NY, "width": WIDTH, "height": HEIGHT},
            "time": {"dt": DT, "n_steps": N_STEPS, "n_sub_steps": N_SUB_STEPS},
            "bounds": {
                "stiffness": [STIFFNESS_MIN, STIFFNESS_MAX],
                "tau": [TAU_MIN, TAU_MAX],
            },
            "regularization": {
                "l2_stiffness": REG_L2_STIFFNESS,
                "tv_stiffness": REG_TV_STIFFNESS,
                "l2_tau": REG_L2_TAU,
                "tv_tau": REG_TV_TAU,
                "tv_eps": REG_TV_EPS,
            },
            "optimizer": {
                "method": "L-BFGS-B",
                "max_iters": MAX_ITERS,
                "gtol": LBFGSB_GTOL,
                "maxls": LBFGSB_MAXLS,
            },
            "material_discovery": {
                "k_min": K_MIN,
                "k_max": K_MAX,
                "kmeans_max_iters": KMEANS_MAX_ITERS,
                "kmeans_restarts": KMEANS_RESTARTS,
                "model_selection_penalty": MODEL_SELECTION_PENALTY,
                "seed": RANDOM_SEED,
            },
        },
        "result": {
            "status": int(result.status),
            "message": str(result.message),
            "objective_final": float(result.fun),
            "data_loss_final": float(final_eval["data_loss"]),
        },
        "fields": {
            "true_stiffness": true_stiff.tolist(),
            "true_tau": true_tau.tolist(),
            "init_stiffness": init_stiff.tolist(),
            "init_tau": init_tau.tolist(),
            "recovered_stiffness": rec_stiff.tolist(),
            "recovered_tau": rec_tau.tolist(),
        },
        "material_discovery": {
            "k": int(discovered["k"]),
            "score": float(discovered["score"]),
            "inertia": float(discovered["inertia"]),
            "centers": discovered["centers"].tolist(),
            "labels": discovered["labels"].tolist(),
            "true_region_labels": true_labels.tolist(),
        },
        "history": history,
    }

    with OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    save_summary_plot(
        problem,
        true_stiff,
        true_tau,
        init_stiff,
        init_tau,
        rec_stiff,
        rec_tau,
        true_labels,
        discovered["labels"],
        history,
    )

    print("\nRun complete")
    print(f"  status             : {result.status}")
    print(f"  message            : {result.message}")
    print(f"  recovered materials: K = {discovered['k']}")
    print(f"  final total loss   : {result.fun:.6e}")
    print(f"  final data loss    : {final_eval['data_loss']:.6e}")
    print(f"  wrote json         : {OUTPUT_JSON}")
    print(f"  wrote plot         : {OUTPUT_PLOT}")


if __name__ == "__main__":
    main()
