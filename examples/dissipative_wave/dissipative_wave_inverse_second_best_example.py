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
# Element-wise second-best example: wider domain + circular third material
# -----------------------------------------------------------------------------

WIDTH = 8.0
HEIGHT = 2.0
NX = 28
NY = 8

DT = 4.0e-3
N_STEPS = 340
N_SUB_STEPS = 1
INTEGRATOR_TYPE = "fixed_visco"

IMPACT_VELOCITY = 1.50
IMPACT_WINDOW_WIDTH = 0.8
IMPACT_CENTERS = [0.8, 2.0, 3.2, 4.4, 5.6, 7.2]

STIFFNESS_MIN = 1.0
STIFFNESS_MAX = 20.0
INIT_BOTTOM_STIFFNESS = 14.0
INIT_TOP_STIFFNESS = 8.0
FIXED_TAU = 0.1

TRUE_BOTTOM_STIFFNESS = 15.0
TRUE_TOP_STIFFNESS = 7.7
TRUE_ROUND_STIFFNESS = 12.0
ROUND_CENTER = (5.4, 1.45)
ROUND_RADIUS = 0.45

REG_L2_STIFFNESS = 5.0e-4
REG_TV_STIFFNESS = 5.0e-4
REG_TV_EPS = 1.0e-5

MAX_ITERS = 180
LBFGSB_GTOL = 1.0e-10
LBFGSB_MAXLS = 60
MAX_OBJECTIVE_EVALS = 80

OUTPUT_DIR = THIS_DIR / "inverse_outputs_second_best_example"
OUTPUT_JSON = OUTPUT_DIR / "second_best_example_result.json"
OUTPUT_PLOT = OUTPUT_DIR / "second_best_example_summary.svg"

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
    inv.TAU_MIN = FIXED_TAU
    inv.TAU_MAX = FIXED_TAU


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


def make_true_fields(problem):
    xy = np.asarray(problem["elem_centers"], dtype=np.double)
    y = xy[:, 1]

    labels = np.zeros(y.size, dtype=np.int32)
    labels[y >= 0.5 * HEIGHT] = 1

    dx = xy[:, 0] - ROUND_CENTER[0]
    dy = xy[:, 1] - ROUND_CENTER[1]
    round_mask = (dx * dx + dy * dy <= ROUND_RADIUS * ROUND_RADIUS) & (y >= 0.5 * HEIGHT)
    labels[round_mask] = 2

    stiffness_by_label = np.array(
        [TRUE_BOTTOM_STIFFNESS, TRUE_TOP_STIFFNESS, TRUE_ROUND_STIFFNESS],
        dtype=np.double,
    )
    return stiffness_by_label[labels], labels


def make_initial_stiffness(labels):
    return np.where(labels == 0, INIT_BOTTOM_STIFFNESS, INIT_TOP_STIFFNESS).astype(np.double)


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

    return {
        "sensor_history": history,
        "data_loss": float(data_loss),
        "grad_stiff": grad_stiff,
    }


def generate_observations(experiments, true_stiffness, tau_elem):
    observations = []
    norm_sq = 0.0
    for problem in experiments:
        run = run_velocity_history(problem, true_stiffness, tau_elem)
        history = run["sensor_history"]
        observations.append(history)
        norm_sq += float(np.sum(history * history))
    return observations, max(norm_sq, 1.0)


def invert_elementwise(experiments, observations, obs_norm_sq, tau_elem, init_stiffness):
    nelem = experiments[0]["connectivity"].shape[0]
    i_idx, j_idx = inv.build_edge_pairs(NX, NY)

    state = {
        "x": init_stiffness.copy(),
        "total_loss": None,
        "data_loss": None,
        "reg_loss": None,
        "grad_inf": None,
        "status": "running",
    }
    history = []

    def objective_with_grad(x):
        x = np.asarray(x, dtype=np.double)
        total_data_loss = 0.0
        grad_data = np.zeros(nelem, dtype=np.double)

        for problem, observed in zip(experiments, observations):
            run = run_velocity_history(
                problem,
                x,
                tau_elem,
                observed_history=observed,
                compute_gradients=True,
            )
            total_data_loss += run["data_loss"]
            grad_data += run["grad_stiff"]

        data_loss = total_data_loss / obs_norm_sq
        grad_data /= obs_norm_sq

        reg_loss_raw, grad_reg_raw = inv.regularization_loss_and_grad(
            x,
            i_idx,
            j_idx,
            REG_L2_STIFFNESS,
            REG_TV_STIFFNESS,
            REG_TV_EPS,
        )
        reg_loss = reg_loss_raw / nelem
        grad_reg = grad_reg_raw / nelem

        total_loss = data_loss + reg_loss
        grad = grad_data + grad_reg

        state["x"] = x.copy()
        state["total_loss"] = float(total_loss)
        state["data_loss"] = float(data_loss)
        state["reg_loss"] = float(reg_loss)
        state["grad_inf"] = float(np.max(np.abs(grad)))

        history.append(
            {
                "eval": len(history),
                "loss_total": float(total_loss),
                "loss_data": float(data_loss),
                "loss_reg": float(reg_loss),
                "grad_inf": float(np.max(np.abs(grad))),
            }
        )

        if len(history) % 10 == 0:
            print(
                f"  eval={len(history):03d} total={total_loss:.6e} "
                f"data={data_loss:.6e} reg={reg_loss:.6e}"
            )

        if len(history) >= MAX_OBJECTIVE_EVALS:
            state["status"] = "max_objective_evals"
            raise EarlyStop()

        return float(total_loss), grad

    t0 = time.time()
    try:
        result = minimize(
            objective_with_grad,
            x0=init_stiffness,
            method="L-BFGS-B",
            jac=True,
            bounds=[(STIFFNESS_MIN, STIFFNESS_MAX)] * nelem,
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
        "final_reg_loss": state["reg_loss"],
        "final_grad_inf": state["grad_inf"],
    }


def evaluate_recovered_histories(experiments, observations, recovered_stiffness, tau_elem, obs_norm_sq):
    recovered_histories = []
    data_loss = 0.0
    for problem, observed in zip(experiments, observations):
        run = run_velocity_history(problem, recovered_stiffness, tau_elem, observed_history=observed)
        recovered_histories.append(run["sensor_history"])
        data_loss += run["data_loss"]
    return recovered_histories, float(data_loss / obs_norm_sq)


def summarize_regions(labels, recovered_stiffness):
    truth = [TRUE_BOTTOM_STIFFNESS, TRUE_TOP_STIFFNESS, TRUE_ROUND_STIFFNESS]
    names = ["bottom", "top", "round"]
    rows = []
    for label, name, true_value in zip([0, 1, 2], names, truth):
        values = recovered_stiffness[labels == label]
        rows.append(
            {
                "name": name,
                "n_elements": int(values.size),
                "true_stiffness": float(true_value),
                "mean_recovered_stiffness": float(np.mean(values)),
                "mean_absolute_error": float(np.mean(np.abs(values - true_value))),
                "max_absolute_error": float(np.max(np.abs(values - true_value))),
            }
        )
    return rows


def save_summary_plot(problem, labels, true_stiffness, init_stiffness, recovered_stiffness, observed_history, recovered_history):
    time_axis = DT * np.arange(N_STEPS, dtype=np.double)
    extent = [0.0, WIDTH, 0.0, HEIGHT]

    fig, axes = plt.subplots(2, 3, figsize=(13.0, 7.0))
    field_panels = [
        (axes[0, 0], true_stiffness, "True stiffness"),
        (axes[0, 1], init_stiffness, "Initial stiffness"),
        (axes[0, 2], recovered_stiffness, "Recovered stiffness"),
        (axes[1, 0], labels, "True material regions"),
        (axes[1, 1], np.abs(recovered_stiffness - true_stiffness), "Absolute stiffness error"),
    ]

    stiffness_vmin = float(min(np.min(true_stiffness), np.min(init_stiffness), np.min(recovered_stiffness)))
    stiffness_vmax = float(max(np.max(true_stiffness), np.max(init_stiffness), np.max(recovered_stiffness)))
    for ax, field, title in field_panels:
        is_error = "error" in title.lower()
        is_labels = "regions" in title.lower()
        if is_labels:
            cmap = "tab20"
            vmin = 0
            vmax = 2
        else:
            cmap = "viridis"
            vmin = 0.0 if is_error else stiffness_vmin
            vmax = float(np.max(field)) if is_error else stiffness_vmax
        im = ax.imshow(
            np.asarray(field, dtype=np.double).reshape(NY, NX),
            origin="lower",
            extent=extent,
            aspect="auto",
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
    fig.savefig(OUTPUT_PLOT, dpi=200, metadata={"Title": "Element-wise dissipative wave second-best example"})
    fig.clf()
    plt.close(fig)


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    configure_inverse_backend()

    experiments = [make_experiment(center) for center in IMPACT_CENTERS]
    nelem = experiments[0]["connectivity"].shape[0]
    true_stiffness, labels = make_true_fields(experiments[0])
    init_stiffness = make_initial_stiffness(labels)
    tau_elem = np.full(nelem, FIXED_TAU, dtype=np.double)

    print("Generating synthetic observations for the wider three-material target...")
    observations, obs_norm_sq = generate_observations(experiments, true_stiffness, tau_elem)

    print("Running element-wise stiffness inversion...")
    result = invert_elementwise(experiments, observations, obs_norm_sq, tau_elem, init_stiffness)
    recovered_stiffness = np.asarray(result["x"], dtype=np.double)

    recovered_histories, normalized_data_loss = evaluate_recovered_histories(
        experiments,
        observations,
        recovered_stiffness,
        tau_elem,
        obs_norm_sq,
    )

    error = recovered_stiffness - true_stiffness
    abs_error = np.abs(error)
    relative_field_error = float(np.linalg.norm(error) / np.linalg.norm(true_stiffness))
    rmse = float(np.sqrt(np.mean(error * error)))
    mae = float(np.mean(abs_error))
    region_summary = summarize_regions(labels, recovered_stiffness)

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
                "bottom_stiffness": TRUE_BOTTOM_STIFFNESS,
                "top_stiffness": TRUE_TOP_STIFFNESS,
                "round_stiffness": TRUE_ROUND_STIFFNESS,
                "round_center": list(ROUND_CENTER),
                "round_radius": ROUND_RADIUS,
                "round_region_element_count": int(np.sum(labels == 2)),
            },
            "unknowns": {
                "type": "element-wise stiffness scaling factor",
                "n_unknowns": int(nelem),
                "initial_bottom_stiffness": INIT_BOTTOM_STIFFNESS,
                "initial_top_stiffness": INIT_TOP_STIFFNESS,
                "initial_round_region_uses_top_value": True,
                "bounds": [STIFFNESS_MIN, STIFFNESS_MAX],
                "fixed_tau": FIXED_TAU,
            },
            "regularization": {
                "l2_stiffness": REG_L2_STIFFNESS,
                "tv_stiffness": REG_TV_STIFFNESS,
                "tv_eps": REG_TV_EPS,
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
            "relative_field_error": relative_field_error,
            "rmse": rmse,
            "mae": mae,
            "max_absolute_error": float(np.max(abs_error)),
            "recovered_min": float(np.min(recovered_stiffness)),
            "recovered_max": float(np.max(recovered_stiffness)),
            "region_summary": region_summary,
        },
        "fields": {
            "labels": labels.tolist(),
            "true_stiffness": true_stiffness.tolist(),
            "initial_stiffness": init_stiffness.tolist(),
            "recovered_stiffness": recovered_stiffness.tolist(),
            "tau": tau_elem.tolist(),
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
        observations[0],
        recovered_histories[0],
    )

    print("\nRun complete")
    print(f"  unknowns              : {nelem} element-wise stiffness values")
    print(f"  round material elems  : {int(np.sum(labels == 2))}")
    print(f"  normalized data loss  : {normalized_data_loss:.6e}")
    print(f"  relative field error  : {relative_field_error:.6e}")
    print(f"  stiffness RMSE        : {rmse:.6e}")
    print(f"  stiffness MAE         : {mae:.6e}")
    for row in region_summary:
        print(
            f"  {row['name']:>6s} region mean : true={row['true_stiffness']:.6f}, "
            f"recovered={row['mean_recovered_stiffness']:.6f}, "
            f"MAE={row['mean_absolute_error']:.6e}"
        )
    print(f"  fixed tau             : {FIXED_TAU:.6f}")
    print(f"  wrote json            : {OUTPUT_JSON}")
    print(f"  wrote plot            : {OUTPUT_PLOT}")


if __name__ == "__main__":
    main()
