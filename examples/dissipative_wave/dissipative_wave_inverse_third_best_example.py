import contextlib
import json
import os
import sys
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

import dissipative_wave_inverse as inv
import REMAT


THIS_DIR = Path(__file__).resolve().parent

# -----------------------------------------------------------------------------
# Element-wise third-best example:
#   element-wise stiffness + three region-wise relaxation times
# -----------------------------------------------------------------------------

WIDTH = 8.0
HEIGHT = 2.0
NX = 28
NY = 8

DT = 4.0e-3
N_STEPS = 420
N_SUB_STEPS = 1
INTEGRATOR_TYPE = "fixed_visco"

IMPACT_VELOCITY = 1.50
IMPACT_WINDOW_WIDTH = 0.8
IMPACT_CENTERS = [0.8, 2.0, 3.2, 4.4, 5.6, 7.2]

STIFFNESS_MIN = 1.0
STIFFNESS_MAX = 20.0
INIT_BOTTOM_STIFFNESS = 14.0
INIT_TOP_STIFFNESS = 8.0

# Keep tau comfortably above the very dissipative tau < 0.1 regime.
TAU_MIN = 0.2
TAU_MAX = 1.2
INIT_REGION_TAU = np.array([0.5, 0.5, 0.5], dtype=np.double)

TRUE_BOTTOM_STIFFNESS = 15.0
TRUE_TOP_STIFFNESS = 7.7
TRUE_ROUND_STIFFNESS = 12.0

TRUE_BOTTOM_TAU = 0.55
TRUE_TOP_TAU = 0.35
TRUE_ROUND_TAU = 0.80

ROUND_CENTER = (5.4, 1.45)
ROUND_RADIUS = 0.45

REG_L2_STIFFNESS = 5.0e-4
REG_TV_STIFFNESS = 5.0e-4
REG_TV_EPS = 1.0e-5

MAX_ITERS = 200
LBFGSB_GTOL = 1.0e-10
LBFGSB_MAXLS = 60
MAX_OBJECTIVE_EVALS = 90

OVERFLOW_LIMIT = 10.0
MAT_OVERFLOW_LIMIT = 10.0

OUTPUT_DIR = THIS_DIR / "inverse_outputs_third_best_example"
OUTPUT_JSON = OUTPUT_DIR / "third_best_example_result.json"
OUTPUT_PLOT = OUTPUT_DIR / "third_best_example_summary.svg"

PLOT_COLOR_PRIMARY = "#2b738eff"
PLOT_COLOR_SECONDARY = "#f9826bff"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)


class EarlyStop(Exception):
    pass


@contextlib.contextmanager
def suppress_solver_output():
    sys.stdout.flush()
    sys.stderr.flush()
    stdout_fd = os.dup(1)
    stderr_fd = os.dup(2)
    with open(os.devnull, "w", encoding="utf-8") as devnull:
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            sys.stdout.flush()
            sys.stderr.flush()
            os.dup2(stdout_fd, 1)
            os.dup2(stderr_fd, 2)
            os.close(stdout_fd)
            os.close(stderr_fd)


def configure_inverse_backend():
    inv.WIDTH = WIDTH
    inv.HEIGHT = HEIGHT
    inv.NX = NX
    inv.NY = NY
    inv.DT = DT
    inv.N_STEPS = N_STEPS
    inv.N_SUB_STEPS = N_SUB_STEPS
    inv.INTEGRATOR_TYPE = INTEGRATOR_TYPE
    inv.IMPACT_VELOCITY = IMPACT_VELOCITY
    inv.IMPACT_WINDOW_WIDTH = IMPACT_WINDOW_WIDTH
    inv.STIFFNESS_MIN = STIFFNESS_MIN
    inv.STIFFNESS_MAX = STIFFNESS_MAX
    inv.TAU_MIN = TAU_MIN
    inv.TAU_MAX = TAU_MAX
    inv.OVERFLOW_LIMIT = OVERFLOW_LIMIT
    inv.MAT_OVERFLOW_LIMIT = MAT_OVERFLOW_LIMIT


def configure_impact(problem, center_x):
    coordinates = problem["coordinates"]
    velocities = problem["velocities"]
    node_x = coordinates[:, 0]
    node_y = coordinates[:, 1]
    eps = 1.0e-12

    velocities[:, :] = 0.0
    top_nodes = np.where(np.abs(node_y - HEIGHT) < eps)[0]
    half_width = 0.5 * min(IMPACT_WINDOW_WIDTH, WIDTH)
    impact_nodes = top_nodes[
        (node_x[top_nodes] >= center_x - half_width - eps)
        & (node_x[top_nodes] <= center_x + half_width + eps)
    ]
    if impact_nodes.size == 0:
        nearest = int(np.argmin(np.abs(node_x[top_nodes] - center_x)))
        impact_nodes = top_nodes[np.array([nearest], dtype=np.int32)]

    velocities[impact_nodes, 1] = -abs(IMPACT_VELOCITY)
    problem["impact_center"] = float(center_x)
    problem["impact_nodes"] = impact_nodes.astype(np.int32)


def use_dense_velocity_sensors(problem):
    coordinates = problem["coordinates"]
    node_y = coordinates[:, 1]
    problem["sensor_nodes"] = np.where(node_y > 1.0e-12)[0].astype(np.int32)


def make_experiment(center_x):
    problem = inv.make_structured_quad_problem()
    configure_impact(problem, center_x)
    use_dense_velocity_sensors(problem)
    return problem


def make_labels(problem):
    xy = np.asarray(problem["elem_centers"], dtype=np.double)
    y = xy[:, 1]

    labels = np.zeros(y.size, dtype=np.int32)
    labels[y >= 0.5 * HEIGHT] = 1

    dx = xy[:, 0] - ROUND_CENTER[0]
    dy = xy[:, 1] - ROUND_CENTER[1]
    round_mask = (dx * dx + dy * dy <= ROUND_RADIUS * ROUND_RADIUS) & (y >= 0.5 * HEIGHT)
    labels[round_mask] = 2
    return labels


def region_stiffness_values():
    return np.array(
        [TRUE_BOTTOM_STIFFNESS, TRUE_TOP_STIFFNESS, TRUE_ROUND_STIFFNESS],
        dtype=np.double,
    )


def region_tau_values():
    return np.array(
        [TRUE_BOTTOM_TAU, TRUE_TOP_TAU, TRUE_ROUND_TAU],
        dtype=np.double,
    )


def expand_region_values(values, labels):
    return np.asarray(values, dtype=np.double)[labels]


def make_initial_stiffness(labels):
    return np.where(labels == 0, INIT_BOTTOM_STIFFNESS, INIT_TOP_STIFFNESS).astype(np.double)


def pack_controls(stiffness_elem, region_tau):
    return np.concatenate(
        [
            np.asarray(stiffness_elem, dtype=np.double),
            np.asarray(region_tau, dtype=np.double),
        ]
    )


def unpack_controls(x, labels):
    x = np.asarray(x, dtype=np.double)
    nelem = labels.size
    stiffness_elem = x[:nelem]
    region_tau = x[nelem : nelem + 3]
    tau_elem = expand_region_values(region_tau, labels)
    return stiffness_elem, region_tau, tau_elem


def run_velocity_history(problem, stiffness_elem, tau_elem, observed_history=None, compute_gradients=False):
    with suppress_solver_output():
        inv.configure_run(problem, stiffness_elem, tau_elem)

        sensor_nodes = problem["sensor_nodes"]
        nsensors = sensor_nodes.size
        nelem = problem["connectivity"].shape[0]

        history = np.zeros((N_STEPS, nsensors, 2), dtype=np.double)
        data_loss = 0.0

        for k in range(N_STEPS):
            REMAT.API.update_state(DT, N_SUB_STEPS, REMAT.PASS_FORWARD)
            history[k, :, 0] = np.asarray(REMAT.get_field(b"node", "velocity_X"), dtype=np.double)[sensor_nodes]
            history[k, :, 1] = np.asarray(REMAT.get_field(b"node", "velocity_Y"), dtype=np.double)[sensor_nodes]
            if observed_history is not None:
                residual = history[k, :, :] - observed_history[k, :, :]
                data_loss += 0.5 * float(np.sum(residual * residual))

        grad_stiff = np.zeros(nelem, dtype=np.double)
        grad_tau = np.zeros(nelem, dtype=np.double)
        if compute_gradients:
            if observed_history is None:
                raise ValueError("observed_history is required when compute_gradients=True")

            REMAT.clear_adjoint_state()
            for rev in range(N_STEPS):
                k = N_STEPS - 1 - rev
                seed_xy = history[k, :, :] - observed_history[k, :, :]
                REMAT.add_nodal_velocity_adjoint_seed(sensor_nodes, np.asarray(seed_xy, dtype=np.double))
                REMAT.API.update_state(DT, N_SUB_STEPS, REMAT.PASS_BACKWARD_ADJOINT)

            grad_stiff = np.asarray(
                REMAT.get_field(b"element", "dparam_stiffness_scaling_factor"),
                dtype=np.double,
            )
            grad_tau = np.asarray(
                REMAT.get_field(b"element", "dparam_relaxation_time"),
                dtype=np.double,
            )

    return {
        "sensor_history": history,
        "data_loss": float(data_loss),
        "grad_stiff": grad_stiff,
        "grad_tau": grad_tau,
    }


def generate_observations(experiments, true_stiffness, true_tau):
    observations = []
    norm_sq = 0.0
    for problem in experiments:
        run = run_velocity_history(problem, true_stiffness, true_tau)
        history = run["sensor_history"]
        observations.append(history)
        norm_sq += float(np.sum(history * history))
    return observations, max(norm_sq, 1.0)


def region_tau_gradient(grad_tau_elem, labels):
    grad = np.zeros(3, dtype=np.double)
    for label in range(3):
        grad[label] = float(np.sum(grad_tau_elem[labels == label]))
    return grad


def invert_elementwise_with_region_tau(experiments, observations, obs_norm_sq, labels, init_stiffness, init_region_tau):
    nelem = labels.size
    i_idx, j_idx = inv.build_edge_pairs(NX, NY)
    x0 = pack_controls(init_stiffness, init_region_tau)
    bounds = [(STIFFNESS_MIN, STIFFNESS_MAX)] * nelem + [(TAU_MIN, TAU_MAX)] * 3

    state = {
        "x": x0.copy(),
        "total_loss": None,
        "data_loss": None,
        "reg_stiffness_loss": None,
        "grad_inf_stiffness": None,
        "grad_inf_tau": None,
        "status": "running",
    }
    history = []

    def objective_with_grad(x):
        stiffness_elem, region_tau, tau_elem = unpack_controls(x, labels)
        total_data_loss = 0.0
        grad_data_stiff = np.zeros(nelem, dtype=np.double)
        grad_data_region_tau = np.zeros(3, dtype=np.double)

        for problem, observed in zip(experiments, observations):
            run = run_velocity_history(
                problem,
                stiffness_elem,
                tau_elem,
                observed_history=observed,
                compute_gradients=True,
            )
            total_data_loss += run["data_loss"]
            grad_data_stiff += run["grad_stiff"]
            grad_data_region_tau += region_tau_gradient(run["grad_tau"], labels)

        data_loss = total_data_loss / obs_norm_sq
        grad_data_stiff /= obs_norm_sq
        grad_data_region_tau /= obs_norm_sq

        reg_s_loss_raw, reg_s_grad_raw = inv.regularization_loss_and_grad(
            stiffness_elem,
            i_idx,
            j_idx,
            REG_L2_STIFFNESS,
            REG_TV_STIFFNESS,
            REG_TV_EPS,
        )
        reg_s_loss = reg_s_loss_raw / nelem
        reg_s_grad = reg_s_grad_raw / nelem

        total_loss = data_loss + reg_s_loss
        grad_stiff = grad_data_stiff + reg_s_grad
        grad = np.concatenate([grad_stiff, grad_data_region_tau])

        state["x"] = np.asarray(x, dtype=np.double).copy()
        state["total_loss"] = float(total_loss)
        state["data_loss"] = float(data_loss)
        state["reg_stiffness_loss"] = float(reg_s_loss)
        state["grad_inf_stiffness"] = float(np.max(np.abs(grad_stiff)))
        state["grad_inf_tau"] = float(np.max(np.abs(grad_data_region_tau)))

        history.append(
            {
                "eval": len(history),
                "loss_total": float(total_loss),
                "loss_data": float(data_loss),
                "loss_reg_stiffness": float(reg_s_loss),
                "grad_inf_stiffness": float(np.max(np.abs(grad_stiff))),
                "grad_inf_tau": float(np.max(np.abs(grad_data_region_tau))),
                "region_tau": region_tau.tolist(),
            }
        )

        if len(history) % 10 == 0:
            print(
                f"  eval={len(history):03d} total={total_loss:.6e} "
                f"data={data_loss:.6e} reg_s={reg_s_loss:.6e} "
                f"tau=[{region_tau[0]:.4f}, {region_tau[1]:.4f}, {region_tau[2]:.4f}]"
            )

        if len(history) >= MAX_OBJECTIVE_EVALS:
            state["status"] = "max_objective_evals"
            raise EarlyStop()

        return float(total_loss), grad

    t0 = time.time()
    try:
        result = minimize(
            objective_with_grad,
            x0=x0,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds,
            options={
                "maxiter": MAX_ITERS,
                "gtol": LBFGSB_GTOL,
                "maxls": LBFGSB_MAXLS,
                "ftol": 1.0e-13,
            },
        )
        state["x"] = np.asarray(result.x, dtype=np.double)
        state["status"] = f"scipy_status_{int(result.status)}"
        message = str(result.message)
        success = bool(result.success)
    except EarlyStop:
        message = state["status"]
        success = True

    return {
        "x": state["x"],
        "history": history,
        "status": state["status"],
        "message": message,
        "success": success,
        "elapsed_seconds": float(time.time() - t0),
        "n_objective_evals": len(history),
        "final_total_loss": state["total_loss"],
        "final_data_loss": state["data_loss"],
        "final_reg_stiffness_loss": state["reg_stiffness_loss"],
        "final_grad_inf_stiffness": state["grad_inf_stiffness"],
        "final_grad_inf_tau": state["grad_inf_tau"],
    }


def evaluate_recovered_histories(experiments, observations, recovered_stiffness, recovered_tau, obs_norm_sq):
    recovered_histories = []
    data_loss = 0.0
    for problem, observed in zip(experiments, observations):
        run = run_velocity_history(problem, recovered_stiffness, recovered_tau, observed_history=observed)
        recovered_histories.append(run["sensor_history"])
        data_loss += run["data_loss"]
    return recovered_histories, float(data_loss / obs_norm_sq)


def summarize_regions(labels, recovered_stiffness, recovered_region_tau):
    true_stiffness = region_stiffness_values()
    true_tau = region_tau_values()
    names = ["bottom", "top", "round"]
    rows = []
    for label, name in zip([0, 1, 2], names):
        values = recovered_stiffness[labels == label]
        rows.append(
            {
                "name": name,
                "n_elements": int(values.size),
                "true_stiffness": float(true_stiffness[label]),
                "mean_recovered_stiffness": float(np.mean(values)),
                "mean_absolute_stiffness_error": float(np.mean(np.abs(values - true_stiffness[label]))),
                "max_absolute_stiffness_error": float(np.max(np.abs(values - true_stiffness[label]))),
                "true_tau": float(true_tau[label]),
                "recovered_tau": float(recovered_region_tau[label]),
                "absolute_tau_error": float(abs(recovered_region_tau[label] - true_tau[label])),
            }
        )
    return rows


def save_summary_plot(
    problem,
    labels,
    true_stiffness,
    init_stiffness,
    recovered_stiffness,
    true_tau,
    recovered_tau,
    observed_history,
    recovered_history,
):
    time_axis = DT * np.arange(N_STEPS, dtype=np.double)
    extent = [0.0, WIDTH, 0.0, HEIGHT]

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0))
    field_panels = [
        (axes[0, 0], true_stiffness, "True stiffness", "viridis"),
        (axes[0, 1], init_stiffness, "Initial stiffness", "viridis"),
        (axes[0, 2], recovered_stiffness, "Recovered stiffness", "viridis"),
        (axes[1, 0], true_tau, "True tau", "magma"),
        (axes[1, 1], recovered_tau, "Recovered tau", "magma"),
    ]

    stiffness_vmin = float(min(np.min(true_stiffness), np.min(init_stiffness), np.min(recovered_stiffness)))
    stiffness_vmax = float(max(np.max(true_stiffness), np.max(init_stiffness), np.max(recovered_stiffness)))
    tau_vmin = float(min(np.min(true_tau), np.min(recovered_tau)))
    tau_vmax = float(max(np.max(true_tau), np.max(recovered_tau)))

    for ax, field, title, cmap in field_panels:
        if "tau" in title.lower():
            vmin = tau_vmin
            vmax = tau_vmax
        else:
            vmin = stiffness_vmin
            vmax = stiffness_vmax
        im = ax.imshow(
            np.asarray(field, dtype=np.double).reshape(NY, NX),
            origin="lower",
            extent=extent,
            aspect="equal",
            interpolation="nearest",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(title, fontsize="medium")
        ax.set_xlabel("x", fontsize="large")
        ax.set_ylabel("y", fontsize="large")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    ax = axes[1, 2]
    sensor_nodes = problem["sensor_nodes"]
    coordinates = problem["coordinates"]
    sensor_xy = coordinates[sensor_nodes]
    target = np.array([problem["impact_center"], HEIGHT], dtype=np.double)
    sensor_idx = int(np.argmin(np.linalg.norm(sensor_xy - target[None, :], axis=1)))
    ax.plot(
        time_axis,
        observed_history[:, sensor_idx, 1],
        color=PLOT_COLOR_PRIMARY,
        linewidth=1.6,
        label="observed",
    )
    ax.plot(
        time_axis,
        recovered_history[:, sensor_idx, 1],
        color=PLOT_COLOR_SECONDARY,
        linewidth=1.6,
        linestyle="--",
        label="recovered",
    )
    ax.set_title("Top sensor velocity", fontsize="medium")
    ax.set_xlabel("time (s)", fontsize="large")
    ax.set_ylabel("velocity y", fontsize="large")
    ax.set_xlim(0.0, N_STEPS * DT)
    ax.legend(fontsize="medium")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUTPUT_PLOT, dpi=200, metadata={"Title": "Element-wise dissipative wave third-best example"})
    fig.clf()
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    configure_inverse_backend()

    experiments = [make_experiment(center) for center in IMPACT_CENTERS]
    labels = make_labels(experiments[0])
    nelem = labels.size

    true_region_stiffness = region_stiffness_values()
    true_region_tau = region_tau_values()
    true_stiffness = expand_region_values(true_region_stiffness, labels)
    true_tau = expand_region_values(true_region_tau, labels)

    init_stiffness = make_initial_stiffness(labels)
    init_region_tau = INIT_REGION_TAU.copy()
    init_tau = expand_region_values(init_region_tau, labels)

    print("Generating synthetic observations for stiffness + region-wise tau target...")
    observations, obs_norm_sq = generate_observations(experiments, true_stiffness, true_tau)

    print("Running element-wise stiffness plus region-wise tau inversion...")
    result = invert_elementwise_with_region_tau(
        experiments,
        observations,
        obs_norm_sq,
        labels,
        init_stiffness,
        init_region_tau,
    )
    recovered_stiffness, recovered_region_tau, recovered_tau = unpack_controls(result["x"], labels)

    recovered_histories, normalized_data_loss = evaluate_recovered_histories(
        experiments,
        observations,
        recovered_stiffness,
        recovered_tau,
        obs_norm_sq,
    )

    stiffness_error = recovered_stiffness - true_stiffness
    tau_error = recovered_region_tau - true_region_tau
    relative_stiffness_error = float(np.linalg.norm(stiffness_error) / np.linalg.norm(true_stiffness))
    relative_tau_error = float(np.linalg.norm(tau_error) / np.linalg.norm(true_region_tau))
    stiffness_rmse = float(np.sqrt(np.mean(stiffness_error * stiffness_error)))
    stiffness_mae = float(np.mean(np.abs(stiffness_error)))
    region_summary = summarize_regions(labels, recovered_stiffness, recovered_region_tau)

    payload = {
        "config": {
            "mesh": {"nx": NX, "ny": NY, "width": WIDTH, "height": HEIGHT},
            "time": {"dt": DT, "n_steps": N_STEPS, "n_sub_steps": N_SUB_STEPS},
            "impact": {
                "velocity": IMPACT_VELOCITY,
                "window_width": IMPACT_WINDOW_WIDTH,
                "centers": IMPACT_CENTERS,
            },
            "sensors": {
                "type": "dense non-bottom nodal velocity sensors",
                "components": ["velocity_x", "velocity_y"],
                "n_sensors_per_impact": int(experiments[0]["sensor_nodes"].size),
            },
            "truth": {
                "region_order": ["bottom", "top", "round"],
                "region_stiffness": true_region_stiffness.tolist(),
                "region_tau": true_region_tau.tolist(),
                "round_center": list(ROUND_CENTER),
                "round_radius": ROUND_RADIUS,
                "round_region_element_count": int(np.sum(labels == 2)),
            },
            "unknowns": {
                "stiffness_type": "element-wise stiffness scaling factor",
                "n_stiffness_unknowns": int(nelem),
                "tau_type": "one relaxation time per material region",
                "n_tau_unknowns": 3,
                "total_unknowns": int(nelem + 3),
                "initial_bottom_stiffness": INIT_BOTTOM_STIFFNESS,
                "initial_top_stiffness": INIT_TOP_STIFFNESS,
                "initial_region_tau": init_region_tau.tolist(),
                "stiffness_bounds": [STIFFNESS_MIN, STIFFNESS_MAX],
                "tau_bounds": [TAU_MIN, TAU_MAX],
            },
            "regularization": {
                "l2_stiffness": REG_L2_STIFFNESS,
                "tv_stiffness": REG_TV_STIFFNESS,
                "tv_eps": REG_TV_EPS,
                "tau_regularization": "none; tau has only region-wise bounds",
            },
            "safety": {
                "overflow_limit": OVERFLOW_LIMIT,
                "mat_overflow_limit": MAT_OVERFLOW_LIMIT,
                "tau_min_is_above_strong_dissipation_regime": True,
            },
            "optimizer": {
                "method": "L-BFGS-B with adjoint gradient",
                "max_iters": MAX_ITERS,
                "max_objective_evals": MAX_OBJECTIVE_EVALS,
            },
        },
        "result": {
            "status": result["status"],
            "message": result["message"],
            "success": result["success"],
            "n_objective_evals": result["n_objective_evals"],
            "elapsed_seconds": result["elapsed_seconds"],
            "normalized_data_loss": normalized_data_loss,
            "relative_stiffness_error": relative_stiffness_error,
            "stiffness_rmse": stiffness_rmse,
            "stiffness_mae": stiffness_mae,
            "max_absolute_stiffness_error": float(np.max(np.abs(stiffness_error))),
            "relative_tau_error": relative_tau_error,
            "recovered_region_tau": recovered_region_tau.tolist(),
            "absolute_tau_error": np.abs(tau_error).tolist(),
            "region_summary": region_summary,
        },
        "fields": {
            "labels": labels.tolist(),
            "true_stiffness": true_stiffness.tolist(),
            "initial_stiffness": init_stiffness.tolist(),
            "recovered_stiffness": recovered_stiffness.tolist(),
            "true_tau": true_tau.tolist(),
            "initial_tau": init_tau.tolist(),
            "recovered_tau": recovered_tau.tolist(),
        },
        "history": result["history"],
    }

    with OUTPUT_JSON.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    save_summary_plot(
        experiments[0],
        labels,
        true_stiffness,
        init_stiffness,
        recovered_stiffness,
        true_tau,
        recovered_tau,
        observations[0],
        recovered_histories[0],
    )

    print("\nRun complete")
    print(f"  stiffness unknowns    : {nelem} element-wise values")
    print("  tau unknowns          : 3 region-wise values")
    print(f"  normalized data loss  : {normalized_data_loss:.6e}")
    print(f"  stiffness rel. error  : {relative_stiffness_error:.6e}")
    print(f"  stiffness RMSE        : {stiffness_rmse:.6e}")
    print(f"  stiffness MAE         : {stiffness_mae:.6e}")
    print(f"  tau rel. error        : {relative_tau_error:.6e}")
    for row in region_summary:
        print(
            f"  {row['name']:>6s} region      : "
            f"stiff true={row['true_stiffness']:.6f}, "
            f"stiff rec={row['mean_recovered_stiffness']:.6f}, "
            f"tau true={row['true_tau']:.6f}, "
            f"tau rec={row['recovered_tau']:.6f}"
        )
    print(f"  tau bounds            : [{TAU_MIN:.3f}, {TAU_MAX:.3f}]")
    print(f"  wrote json            : {OUTPUT_JSON}")
    print(f"  wrote plot            : {OUTPUT_PLOT}")


if __name__ == "__main__":
    main()
