import argparse
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
# Element-wise fifth-best example:
#   same benchmark as the fourth example, but tau is an independent
#   element-wise field rather than a region-wise control.
# -----------------------------------------------------------------------------

WIDTH = 12.0
DEPTH = 2.0
NX = 72
NZ = 12

DT = 4.0e-3
N_STEPS = 520
N_SUB_STEPS = 1
INTEGRATOR_TYPE = "fixed_visco"

IMPACT_VELOCITY = 1.5
IMPACT_WINDOW_WIDTH = 0.8
IMPACT_CENTERS = [0.8, 2.1, 3.4, 4.7, 6.0, 7.3, 8.6, 9.9, 11.2]
SURFACE_SENSOR_COUNT = 17

STIFFNESS_MIN = 1.0
STIFFNESS_MAX = 20.0
INIT_BOTTOM_STIFFNESS = 14.4
INIT_TOP_STIFFNESS = 11.0
INIT_CROSS_STIFFNESS = 11.0

# Keep tau comfortably above the very dissipative tau < 0.1 regime.
TAU_MIN = 0.2
TAU_MAX = 1.2
INIT_BOTTOM_TAU = 0.50
INIT_TOP_TAU = 0.40
INIT_CROSS_TAU = 0.40

TRUE_BOTTOM_STIFFNESS = 15.0
TRUE_TOP_STIFFNESS = 7.7
TRUE_CROSS_STIFFNESS = 12.0

TRUE_BOTTOM_TAU = 0.55
TRUE_TOP_TAU = 0.35
TRUE_CROSS_TAU = 0.80

CROSS_CENTER = (8.4, 0.55)
CROSS_ARM_HALF_LENGTH = 0.48
CROSS_ARM_HALF_WIDTH = 0.20

REG_L2_STIFFNESS = 2.0e-4
REG_TV_STIFFNESS = 2.0e-3
REG_L2_TAU = 1.0e-3
REG_TV_TAU = 4.0e-3
REG_TV_EPS = 1.0e-5

MAX_ITERS = 260
LBFGSB_GTOL = 1.0e-10
LBFGSB_MAXLS = 60
MAX_OBJECTIVE_EVALS = 240

OVERFLOW_LIMIT = 10.0
MAT_OVERFLOW_LIMIT = 10.0

OUTPUT_DIR = THIS_DIR / "inverse_outputs_fifth_best_example"
OUTPUT_JSON = OUTPUT_DIR / "fifth_best_example_result.json"
OUTPUT_PLOT = OUTPUT_DIR / "fifth_best_example_summary.svg"
OUTPUT_VELOCITY_PLOT = OUTPUT_DIR / "fifth_best_example_velocity.svg"

PLOT_COLOR_PRIMARY = "#2b738eff"
PLOT_COLOR_SECONDARY = "#f9826bff"
STIFFNESS_CMAP = "viridis"
TAU_CMAP = "magma"

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
    inv.HEIGHT = DEPTH
    inv.NX = NX
    inv.NY = NZ
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
    node_vertical = coordinates[:, 1]
    eps = 1.0e-12

    velocities[:, :] = 0.0
    surface_nodes = np.where(np.abs(node_vertical - DEPTH) < eps)[0]
    half_width = 0.5 * min(IMPACT_WINDOW_WIDTH, WIDTH)
    impact_nodes = surface_nodes[
        (node_x[surface_nodes] >= center_x - half_width - eps)
        & (node_x[surface_nodes] <= center_x + half_width + eps)
    ]
    if impact_nodes.size == 0:
        nearest = int(np.argmin(np.abs(node_x[surface_nodes] - center_x)))
        impact_nodes = surface_nodes[np.array([nearest], dtype=np.int32)]

    velocities[impact_nodes, 1] = -abs(IMPACT_VELOCITY)
    problem["impact_center"] = float(center_x)
    problem["impact_nodes"] = impact_nodes.astype(np.int32)


def surface_sensor_targets():
    if SURFACE_SENSOR_COUNT == len(IMPACT_CENTERS):
        return np.asarray(IMPACT_CENTERS, dtype=np.double)
    return np.linspace(min(IMPACT_CENTERS), max(IMPACT_CENTERS), SURFACE_SENSOR_COUNT, dtype=np.double)


def use_sparse_surface_velocity_sensors(problem):
    coordinates = problem["coordinates"]
    node_x = coordinates[:, 0]
    node_vertical = coordinates[:, 1]
    eps = 1.0e-12
    surface_nodes = np.where(np.abs(node_vertical - DEPTH) < eps)[0]
    targets = surface_sensor_targets()
    sensor_nodes = []
    for target_x in targets:
        nearest = int(np.argmin(np.abs(node_x[surface_nodes] - target_x)))
        sensor_nodes.append(int(surface_nodes[nearest]))

    sensor_nodes = np.unique(np.asarray(sensor_nodes, dtype=np.int32))
    problem["sensor_nodes"] = sensor_nodes
    problem["sensor_x_targets"] = targets
    problem["sensor_x_positions"] = node_x[sensor_nodes].astype(np.double)
    problem["sensor_z_positions"] = np.zeros(sensor_nodes.size, dtype=np.double)


def make_experiment(center_x):
    problem = inv.make_structured_quad_problem()
    configure_impact(problem, center_x)
    use_sparse_surface_velocity_sensors(problem)
    return problem


def make_labels(problem):
    centers = np.asarray(problem["elem_centers"], dtype=np.double)
    elem_z = DEPTH - centers[:, 1]

    labels = np.zeros(elem_z.size, dtype=np.int32)
    labels[elem_z <= 0.5 * DEPTH] = 1

    dx = np.abs(centers[:, 0] - CROSS_CENTER[0])
    dz = np.abs(elem_z - CROSS_CENTER[1])
    cross_mask = (
        ((dx <= CROSS_ARM_HALF_WIDTH) & (dz <= CROSS_ARM_HALF_LENGTH))
        | ((dx <= CROSS_ARM_HALF_LENGTH) & (dz <= CROSS_ARM_HALF_WIDTH))
    ) & (elem_z <= 0.5 * DEPTH)
    labels[cross_mask] = 2
    return labels


def region_stiffness_values():
    return np.array(
        [TRUE_BOTTOM_STIFFNESS, TRUE_TOP_STIFFNESS, TRUE_CROSS_STIFFNESS],
        dtype=np.double,
    )


def region_tau_values():
    return np.array(
        [TRUE_BOTTOM_TAU, TRUE_TOP_TAU, TRUE_CROSS_TAU],
        dtype=np.double,
    )


def initial_stiffness_values():
    return np.array(
        [INIT_BOTTOM_STIFFNESS, INIT_TOP_STIFFNESS, INIT_CROSS_STIFFNESS],
        dtype=np.double,
    )


def initial_tau_values():
    return np.array(
        [INIT_BOTTOM_TAU, INIT_TOP_TAU, INIT_CROSS_TAU],
        dtype=np.double,
    )


def expand_region_values(values, labels):
    return np.asarray(values, dtype=np.double)[labels]


def make_initial_stiffness(labels):
    return expand_region_values(initial_stiffness_values(), labels)


def make_initial_tau(labels):
    return expand_region_values(initial_tau_values(), labels)


def pack_controls(stiffness_elem, tau_elem):
    return np.concatenate(
        [
            np.asarray(stiffness_elem, dtype=np.double),
            np.asarray(tau_elem, dtype=np.double),
        ]
    )


def unpack_controls(x, nelem):
    x = np.asarray(x, dtype=np.double)
    stiffness_elem = x[:nelem]
    tau_elem = x[nelem : 2 * nelem]
    return stiffness_elem, tau_elem


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
            history[k, :, 1] = -np.asarray(REMAT.get_field(b"node", "velocity_Y"), dtype=np.double)[sensor_nodes]
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
                residual_xz = history[k, :, :] - observed_history[k, :, :]
                backend_seed = np.empty_like(residual_xz)
                backend_seed[:, 0] = residual_xz[:, 0]
                backend_seed[:, 1] = -residual_xz[:, 1]
                REMAT.add_nodal_velocity_adjoint_seed(sensor_nodes, np.asarray(backend_seed, dtype=np.double))
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


def invert_elementwise_with_full_tau(experiments, observations, obs_norm_sq, init_stiffness, init_tau):
    nelem = init_stiffness.size
    i_idx, j_idx = inv.build_edge_pairs(NX, NZ)
    x0 = pack_controls(init_stiffness, init_tau)
    bounds = [(STIFFNESS_MIN, STIFFNESS_MAX)] * nelem + [(TAU_MIN, TAU_MAX)] * nelem

    state = {
        "x": x0.copy(),
        "total_loss": None,
        "data_loss": None,
        "reg_stiffness_loss": None,
        "reg_tau_loss": None,
        "grad_inf_stiffness": None,
        "grad_inf_tau": None,
        "status": "running",
    }
    history = []

    def objective_with_grad(x):
        stiffness_elem, tau_elem = unpack_controls(x, nelem)
        total_data_loss = 0.0
        grad_data_stiff = np.zeros(nelem, dtype=np.double)
        grad_data_tau = np.zeros(nelem, dtype=np.double)

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
            grad_data_tau += run["grad_tau"]

        data_loss = total_data_loss / obs_norm_sq
        grad_data_stiff /= obs_norm_sq
        grad_data_tau /= obs_norm_sq

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
        reg_t_loss_raw, reg_t_grad_raw = inv.regularization_loss_and_grad(
            tau_elem,
            i_idx,
            j_idx,
            REG_L2_TAU,
            REG_TV_TAU,
            REG_TV_EPS,
        )
        reg_t_loss = reg_t_loss_raw / nelem
        reg_t_grad = reg_t_grad_raw / nelem

        total_loss = data_loss + reg_s_loss + reg_t_loss
        grad_stiff = grad_data_stiff + reg_s_grad
        grad_tau = grad_data_tau + reg_t_grad
        grad = np.concatenate([grad_stiff, grad_tau])

        state["x"] = np.asarray(x, dtype=np.double).copy()
        state["total_loss"] = float(total_loss)
        state["data_loss"] = float(data_loss)
        state["reg_stiffness_loss"] = float(reg_s_loss)
        state["reg_tau_loss"] = float(reg_t_loss)
        state["grad_inf_stiffness"] = float(np.max(np.abs(grad_stiff)))
        state["grad_inf_tau"] = float(np.max(np.abs(grad_tau)))

        history.append(
            {
                "eval": len(history),
                "loss_total": float(total_loss),
                "loss_data": float(data_loss),
                "loss_reg_stiffness": float(reg_s_loss),
                "loss_reg_tau": float(reg_t_loss),
                "grad_inf_stiffness": float(np.max(np.abs(grad_stiff))),
                "grad_inf_tau": float(np.max(np.abs(grad_tau))),
                "tau_min": float(np.min(tau_elem)),
                "tau_mean": float(np.mean(tau_elem)),
                "tau_max": float(np.max(tau_elem)),
            }
        )

        if len(history) % 10 == 0:
            print(
                f"  eval={len(history):03d} total={total_loss:.6e} "
                f"data={data_loss:.6e} reg_s={reg_s_loss:.6e} reg_tau={reg_t_loss:.6e} "
                f"tau_range=[{np.min(tau_elem):.4f}, {np.max(tau_elem):.4f}]"
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
        "final_reg_tau_loss": state["reg_tau_loss"],
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


def summarize_regions(labels, recovered_stiffness, recovered_tau):
    true_stiffness = region_stiffness_values()
    true_tau = region_tau_values()
    names = ["bottom", "top", "cross"]
    rows = []
    for label, name in zip([0, 1, 2], names):
        mask = labels == label
        stiffness_values = recovered_stiffness[mask]
        tau_values = recovered_tau[mask]
        rows.append(
            {
                "name": name,
                "n_elements": int(stiffness_values.size),
                "true_stiffness": float(true_stiffness[label]),
                "mean_recovered_stiffness": float(np.mean(stiffness_values)),
                "mean_absolute_stiffness_error": float(np.mean(np.abs(stiffness_values - true_stiffness[label]))),
                "max_absolute_stiffness_error": float(np.max(np.abs(stiffness_values - true_stiffness[label]))),
                "true_tau": float(true_tau[label]),
                "mean_recovered_tau": float(np.mean(tau_values)),
                "mean_absolute_tau_error": float(np.mean(np.abs(tau_values - true_tau[label]))),
                "max_absolute_tau_error": float(np.max(np.abs(tau_values - true_tau[label]))),
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
    init_tau,
    recovered_tau,
):
    extent = [0.0, WIDTH, DEPTH, 0.0]

    fig, axes = plt.subplots(3, 2, figsize=(9.4, 5.80))
    stiffness_vmin = float(min(np.min(true_stiffness), np.min(init_stiffness), np.min(recovered_stiffness)))
    stiffness_vmax = float(max(np.max(true_stiffness), np.max(init_stiffness), np.max(recovered_stiffness)))
    tau_vmin = float(min(np.min(true_tau), np.min(init_tau), np.min(recovered_tau)))
    tau_vmax = float(max(np.max(true_tau), np.max(init_tau), np.max(recovered_tau)))

    field_panels = [
        (axes[0, 0], true_stiffness, "True stiffness", STIFFNESS_CMAP, stiffness_vmin, stiffness_vmax),
        (axes[1, 0], init_stiffness, "Initial stiffness", STIFFNESS_CMAP, stiffness_vmin, stiffness_vmax),
        (axes[2, 0], recovered_stiffness, "Recovered stiffness", STIFFNESS_CMAP, stiffness_vmin, stiffness_vmax),
        (axes[0, 1], true_tau, "True relaxation time", TAU_CMAP, tau_vmin, tau_vmax),
        (axes[1, 1], init_tau, "Initial relaxation time", TAU_CMAP, tau_vmin, tau_vmax),
        (axes[2, 1], recovered_tau, "Recovered relaxation time", TAU_CMAP, tau_vmin, tau_vmax),
    ]

    stiffness_image = None
    tau_image = None
    for ax, field, title, cmap, vmin, vmax in field_panels:
        im = ax.imshow(
            np.flipud(np.asarray(field, dtype=np.double).reshape(NZ, NX)),
            origin="upper",
            extent=extent,
            aspect="equal",
            interpolation="nearest",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
        ax.set_title(title, fontsize="medium")
        ax.set_xlabel("x", fontsize="large")
        ax.set_ylabel("z", fontsize="large")
        if "stiffness" in title.lower():
            stiffness_image = im
        else:
            tau_image = im

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=[0.03, 0.03, .93, .97], w_pad=5.0, h_pad=0)
    
    left_positions = [ax.get_position() for ax in axes[:, 0]]
    right_positions = [ax.get_position() for ax in axes[:, 1]]
    left_x1 = max(position.x1 for position in left_positions)
    right_x0 = min(position.x0 for position in right_positions)
    right_x1 = max(position.x1 for position in right_positions)
    panel_bottom = min(position.bounds[1] for position in left_positions + right_positions)
    panel_top = max(position.bounds[1] + position.bounds[3] for position in left_positions + right_positions)

    colorbar_width = 0.012
    colorbar_pad = 0.012
    stiffness_cbar_x = left_x1 + colorbar_pad
    tau_cbar_x = right_x1 + colorbar_pad
    stiffness_cax = fig.add_axes([stiffness_cbar_x, panel_bottom, colorbar_width, panel_top - panel_bottom])
    tau_cax = fig.add_axes([tau_cbar_x, panel_bottom, colorbar_width, panel_top - panel_bottom])
    fig.colorbar(stiffness_image, cax=stiffness_cax, orientation="vertical")
    fig.colorbar(tau_image, cax=tau_cax, orientation="vertical")
    divider_x = 0.5 * (stiffness_cbar_x + colorbar_width + right_x0)
    fig.add_artist(
        plt.Line2D(
            [divider_x, divider_x],
            [panel_bottom - 0.02, panel_top + 0.02],
            transform=fig.transFigure,
            color="black",
            linewidth=2.4,
        )
    )
    fig.savefig(OUTPUT_PLOT, dpi=200, metadata={"Title": "Element-wise dissipative wave fifth-best example"})
    fig.clf()
    plt.close(fig)


def save_velocity_plot(problem, observed_history, recovered_history):
    time_axis = DT * np.arange(N_STEPS, dtype=np.double)
    fig, ax = plt.subplots(figsize=(6.2, 3.4))

    sensor_nodes = problem["sensor_nodes"]
    coordinates = problem["coordinates"]
    sensor_coords = coordinates[sensor_nodes]
    target = np.array([problem["impact_center"], DEPTH], dtype=np.double)
    sensor_idx = int(np.argmin(np.linalg.norm(sensor_coords - target[None, :], axis=1)))
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
    ax.set_title("Surface sensor velocity", fontsize="medium")
    ax.set_xlabel("time (s)", fontsize="large")
    ax.set_ylabel("velocity z", fontsize="large")
    ax.set_xlim(0.0, N_STEPS * DT)
    ax.legend(fontsize="medium")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(
        OUTPUT_VELOCITY_PLOT,
        dpi=200,
        metadata={"Title": "Element-wise dissipative wave fifth-best velocity history"},
    )
    fig.clf()
    plt.close(fig)


def load_saved_fields():
    if not OUTPUT_JSON.exists():
        raise FileNotFoundError(f"Saved result JSON does not exist: {OUTPUT_JSON}")

    with OUTPUT_JSON.open(encoding="utf-8") as stream:
        payload = json.load(stream)

    fields = payload["fields"]
    return {
        "labels": np.asarray(fields["labels"], dtype=np.int32),
        "true_stiffness": np.asarray(fields["true_stiffness"], dtype=np.double),
        "initial_stiffness": np.asarray(fields["initial_stiffness"], dtype=np.double),
        "recovered_stiffness": np.asarray(fields["recovered_stiffness"], dtype=np.double),
        "true_tau": np.asarray(fields["true_tau"], dtype=np.double),
        "initial_tau": np.asarray(fields["initial_tau"], dtype=np.double),
        "recovered_tau": np.asarray(fields["recovered_tau"], dtype=np.double),
    }


def regenerate_plots_from_saved_result(include_velocity):
    configure_inverse_backend()
    fields = load_saved_fields()
    experiments = [make_experiment(center) for center in IMPACT_CENTERS]

    save_summary_plot(
        experiments[0],
        fields["labels"],
        fields["true_stiffness"],
        fields["initial_stiffness"],
        fields["recovered_stiffness"],
        fields["true_tau"],
        fields["initial_tau"],
        fields["recovered_tau"],
    )

    if include_velocity:
        print("Regenerating velocity plot from saved fields; no optimization will run.")
        observations, obs_norm_sq = generate_observations(
            experiments,
            fields["true_stiffness"],
            fields["true_tau"],
        )
        recovered_histories, _ = evaluate_recovered_histories(
            experiments,
            observations,
            fields["recovered_stiffness"],
            fields["recovered_tau"],
            obs_norm_sq,
        )
        save_velocity_plot(experiments[0], observations[0], recovered_histories[0])
        print(f"  wrote velocity plot   : {OUTPUT_VELOCITY_PLOT}")

    print(f"  wrote plot            : {OUTPUT_PLOT}")


def apply_runtime_overrides(args):
    global SURFACE_SENSOR_COUNT
    global REG_L2_STIFFNESS, REG_TV_STIFFNESS, REG_L2_TAU, REG_TV_TAU
    global MAX_ITERS, MAX_OBJECTIVE_EVALS
    global OUTPUT_JSON, OUTPUT_PLOT, OUTPUT_VELOCITY_PLOT

    if args.surface_sensors is not None:
        SURFACE_SENSOR_COUNT = int(args.surface_sensors)
    if args.reg_l2_stiffness is not None:
        REG_L2_STIFFNESS = float(args.reg_l2_stiffness)
    if args.reg_tv_stiffness is not None:
        REG_TV_STIFFNESS = float(args.reg_tv_stiffness)
    if args.reg_l2_tau is not None:
        REG_L2_TAU = float(args.reg_l2_tau)
    if args.reg_tv_tau is not None:
        REG_TV_TAU = float(args.reg_tv_tau)
    if args.max_iters is not None:
        MAX_ITERS = int(args.max_iters)
    if args.max_objective_evals is not None:
        MAX_OBJECTIVE_EVALS = int(args.max_objective_evals)

    if args.output_suffix:
        suffix = args.output_suffix
        if not suffix.startswith("_"):
            suffix = f"_{suffix}"
        OUTPUT_JSON = OUTPUT_DIR / f"fifth_best_example_result{suffix}.json"
        OUTPUT_PLOT = OUTPUT_DIR / f"fifth_best_example_summary{suffix}.svg"
        OUTPUT_VELOCITY_PLOT = OUTPUT_DIR / f"fifth_best_example_velocity{suffix}.svg"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plot-only",
        action="store_true",
        help="Regenerate both SVG plots from the saved JSON and skip the optimizer.",
    )
    parser.add_argument(
        "--field-plot-only",
        action="store_true",
        help="Regenerate only the field summary SVG from the saved JSON.",
    )
    parser.add_argument("--surface-sensors", type=int, default=None)
    parser.add_argument("--reg-l2-stiffness", type=float, default=None)
    parser.add_argument("--reg-tv-stiffness", type=float, default=None)
    parser.add_argument("--reg-l2-tau", type=float, default=None)
    parser.add_argument("--reg-tv-tau", type=float, default=None)
    parser.add_argument("--max-iters", type=int, default=None)
    parser.add_argument("--max-objective-evals", type=int, default=None)
    parser.add_argument("--output-suffix", type=str, default="")
    return parser.parse_args()


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
    init_tau = make_initial_tau(labels)

    print("Generating synthetic observations for stiffness + full tau target...")
    observations, obs_norm_sq = generate_observations(experiments, true_stiffness, true_tau)

    print("Running element-wise stiffness plus independent element-wise tau inversion...")
    result = invert_elementwise_with_full_tau(
        experiments,
        observations,
        obs_norm_sq,
        init_stiffness,
        init_tau,
    )
    recovered_stiffness, recovered_tau = unpack_controls(result["x"], nelem)

    recovered_histories, normalized_data_loss = evaluate_recovered_histories(
        experiments,
        observations,
        recovered_stiffness,
        recovered_tau,
        obs_norm_sq,
    )

    stiffness_error = recovered_stiffness - true_stiffness
    tau_error = recovered_tau - true_tau
    relative_stiffness_error = float(np.linalg.norm(stiffness_error) / np.linalg.norm(true_stiffness))
    relative_tau_error = float(np.linalg.norm(tau_error) / np.linalg.norm(true_tau))
    stiffness_rmse = float(np.sqrt(np.mean(stiffness_error * stiffness_error)))
    stiffness_mae = float(np.mean(np.abs(stiffness_error)))
    tau_rmse = float(np.sqrt(np.mean(tau_error * tau_error)))
    tau_mae = float(np.mean(np.abs(tau_error)))
    recovery_score = float(relative_stiffness_error + 2.0 * relative_tau_error + 10.0 * normalized_data_loss)
    region_summary = summarize_regions(labels, recovered_stiffness, recovered_tau)

    payload = {
        "config": {
            "mesh": {
                "nx": NX,
                "nz": NZ,
                "width": WIDTH,
                "depth": DEPTH,
                "coordinate_convention": "z is depth from the surface",
            },
            "time": {"dt": DT, "n_steps": N_STEPS, "n_sub_steps": N_SUB_STEPS},
            "impact": {
                "velocity": IMPACT_VELOCITY,
                "window_width": IMPACT_WINDOW_WIDTH,
                "centers": IMPACT_CENTERS,
            },
            "sensors": {
                "type": "sparse surface nodal velocity sensors",
                "components": ["velocity_x", "velocity_z"],
                "n_sensors_per_impact": int(experiments[0]["sensor_nodes"].size),
                "surface_only": True,
                "x_targets": experiments[0]["sensor_x_targets"].tolist(),
                "x_positions": experiments[0]["sensor_x_positions"].tolist(),
                "z_positions": experiments[0]["sensor_z_positions"].tolist(),
            },
            "truth": {
                "region_order": ["bottom", "top", "cross"],
                "region_stiffness": true_region_stiffness.tolist(),
                "region_tau": true_region_tau.tolist(),
                "cross_center_xz": list(CROSS_CENTER),
                "cross_arm_half_length": CROSS_ARM_HALF_LENGTH,
                "cross_arm_half_width": CROSS_ARM_HALF_WIDTH,
                "cross_region_element_count": int(np.sum(labels == 2)),
            },
            "unknowns": {
                "stiffness_type": "element-wise stiffness scaling factor",
                "n_stiffness_unknowns": int(nelem),
                "tau_type": "element-wise relaxation time independent of stiffness",
                "n_tau_unknowns": int(nelem),
                "total_unknowns": int(2 * nelem),
                "initial_region_stiffness": initial_stiffness_values().tolist(),
                "initial_region_tau": initial_tau_values().tolist(),
                "stiffness_bounds": [STIFFNESS_MIN, STIFFNESS_MAX],
                "tau_bounds": [TAU_MIN, TAU_MAX],
                "tau_parameterization_note": "labels are used only for truth generation and post-hoc summaries, not as tau controls",
            },
            "regularization": {
                "l2_stiffness": REG_L2_STIFFNESS,
                "tv_stiffness": REG_TV_STIFFNESS,
                "l2_tau": REG_L2_TAU,
                "tv_tau": REG_TV_TAU,
                "tv_eps": REG_TV_EPS,
                "tau_regularization": "same mesh-neighbor L2+smoothed-TV form as stiffness, with tau-specific weights",
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
            "recovery_score": recovery_score,
            "relative_stiffness_error": relative_stiffness_error,
            "stiffness_rmse": stiffness_rmse,
            "stiffness_mae": stiffness_mae,
            "max_absolute_stiffness_error": float(np.max(np.abs(stiffness_error))),
            "relative_tau_error": relative_tau_error,
            "tau_rmse": tau_rmse,
            "tau_mae": tau_mae,
            "max_absolute_tau_error": float(np.max(np.abs(tau_error))),
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
        init_tau,
        recovered_tau,
    )
    save_velocity_plot(
        experiments[0],
        observations[0],
        recovered_histories[0],
    )

    print("\nRun complete")
    print(f"  surface sensors       : {experiments[0]['sensor_nodes'].size}")
    print(f"  stiffness unknowns    : {nelem} element-wise values")
    print(f"  tau unknowns          : {nelem} independent element-wise values")
    print(f"  normalized data loss  : {normalized_data_loss:.6e}")
    print(f"  recovery score        : {recovery_score:.6e}")
    print(f"  stiffness rel. error  : {relative_stiffness_error:.6e}")
    print(f"  stiffness RMSE        : {stiffness_rmse:.6e}")
    print(f"  stiffness MAE         : {stiffness_mae:.6e}")
    print(f"  tau rel. error        : {relative_tau_error:.6e}")
    print(f"  tau RMSE              : {tau_rmse:.6e}")
    print(f"  tau MAE               : {tau_mae:.6e}")
    for row in region_summary:
        print(
            f"  {row['name']:>6s} region      : "
            f"stiff true={row['true_stiffness']:.6f}, "
            f"stiff rec={row['mean_recovered_stiffness']:.6f}, "
            f"tau true={row['true_tau']:.6f}, "
            f"tau rec mean={row['mean_recovered_tau']:.6f}"
        )
    print(f"  tau bounds            : [{TAU_MIN:.3f}, {TAU_MAX:.3f}]")
    print(f"  wrote json            : {OUTPUT_JSON}")
    print(f"  wrote plot            : {OUTPUT_PLOT}")
    print(f"  wrote velocity plot   : {OUTPUT_VELOCITY_PLOT}")


if __name__ == "__main__":
    args = parse_args()
    apply_runtime_overrides(args)
    if args.field_plot_only:
        regenerate_plots_from_saved_result(include_velocity=False)
    elif args.plot_only:
        regenerate_plots_from_saved_result(include_velocity=True)
    else:
        main()
