import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
sys.path.append(str(REPO_ROOT / "install" / "package"))

import REMAT


_ACTIVE_LAYER_COEFFS = np.ones(1, dtype=np.double)
_LAYER_BOUNDS = np.array([0.0, 1.0], dtype=np.double)


def layered_stiffness_scaling(_, y):
    idx = np.searchsorted(_LAYER_BOUNDS, y, side="right") - 1
    idx = int(np.clip(idx, 0, _ACTIVE_LAYER_COEFFS.size - 1))
    return float(_ACTIVE_LAYER_COEFFS[idx])


def make_structured_quad_problem(
    nx,
    ny,
    width,
    height,
    impact_velocity,
    impact_window_width,
    nsensors,
    n_layers,
    sensor_distribution_width=None,
):
    xs = np.linspace(0.0, width, nx + 1)
    ys = np.linspace(0.0, height, ny + 1)

    num_nodes = (nx + 1) * (ny + 1)
    coordinates = np.zeros((num_nodes, 2), dtype=np.double)
    velocities = np.zeros((num_nodes, 2), dtype=np.double)
    fixity = np.zeros((num_nodes, 2), dtype=np.bool_)

    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            n = j * (nx + 1) + i
            coordinates[n, 0] = x
            coordinates[n, 1] = y

    num_elems = nx * ny
    connectivity = np.zeros((num_elems, 4), dtype=np.int32)
    e = 0
    for j in range(ny):
        for i in range(nx):
            n0 = j * (nx + 1) + i
            n1 = n0 + 1
            n3 = (j + 1) * (nx + 1) + i
            n2 = n3 + 1
            connectivity[e, :] = [n0, n1, n2, n3]
            e += 1

    eps = 1.0e-12
    node_x = coordinates[:, 0]
    node_y = coordinates[:, 1]
    bottom_nodes = np.where(np.abs(node_y - 0.0) < eps)[0]
    top_nodes = np.where(np.abs(node_y - height) < eps)[0]

    fixity[bottom_nodes, :] = True
    # Free side boundaries: do not constrain left/right edges.

    xmid = 0.5 * width
    impact_window_width = float(impact_window_width)
    if impact_window_width <= 0.0:
        raise ValueError(f"impact_window_width must be > 0, got {impact_window_width}")

    impact_span = min(impact_window_width, width)
    impact_half_width = 0.5 * impact_span
    impact_x_lo = xmid - impact_half_width
    impact_x_hi = xmid + impact_half_width
    impact_nodes = top_nodes[
        (node_x[top_nodes] >= impact_x_lo - eps) & (node_x[top_nodes] <= impact_x_hi + eps)
    ]
    if impact_nodes.size == 0:
        nearest = int(np.argmin(np.abs(node_x[top_nodes] - xmid)))
        impact_nodes = top_nodes[np.array([nearest], dtype=np.int64)]
    velocities[impact_nodes, 1] = -abs(impact_velocity)

    top_interior = top_nodes[(node_x[top_nodes] > eps) & (node_x[top_nodes] < width - eps)]
    if top_interior.size == 0:
        raise ValueError("No interior top-surface nodes found for sensors.")

    if sensor_distribution_width is None:
        sensor_distribution_width = width
    sensor_distribution_width = float(sensor_distribution_width)
    if sensor_distribution_width <= 0.0:
        raise ValueError(f"sensor_distribution_width must be > 0, got {sensor_distribution_width}")

    sensor_span = min(sensor_distribution_width, width)
    sensor_x_lo = xmid - 0.5 * sensor_span
    sensor_x_hi = xmid + 0.5 * sensor_span
    sensor_candidates = top_interior[
        (node_x[top_interior] >= sensor_x_lo - eps) & (node_x[top_interior] <= sensor_x_hi + eps)
    ]
    if sensor_candidates.size == 0:
        raise ValueError("No top-surface nodes found inside sensor distribution width.")

    if nsensors > sensor_candidates.size:
        nsensors = sensor_candidates.size
    sensor_pick = np.unique(np.round(np.linspace(0, sensor_candidates.size - 1, nsensors)).astype(int))
    sensor_nodes = sensor_candidates[sensor_pick]

    layer_bounds = np.linspace(0.0, height, n_layers + 1)
    elem_center_y = np.mean(coordinates[connectivity, 1], axis=1)
    elem_layer_ids = np.searchsorted(layer_bounds, elem_center_y, side="right") - 1
    elem_layer_ids = np.clip(elem_layer_ids, 0, n_layers - 1).astype(np.int32)

    return {
        "coordinates": coordinates,
        "velocities": velocities,
        "fixity": fixity,
        "connectivity": connectivity,
        "sensor_nodes": sensor_nodes.astype(np.int32),
        "impact_nodes": impact_nodes.astype(np.int32),
        "impact_center_x": float(xmid),
        "impact_half_width": float(impact_half_width),
        "impact_window_width": float(impact_span),
        "impact_x_min": float(impact_x_lo),
        "impact_x_max": float(impact_x_hi),
        "sensor_distribution_width": float(sensor_span),
        "sensor_distribution_x_min": float(sensor_x_lo),
        "sensor_distribution_x_max": float(sensor_x_hi),
        "layer_bounds": layer_bounds,
        "elem_layer_ids": elem_layer_ids,
        "width": width,
        "height": height,
    }


def configure_run(problem, layer_coeffs, tau, args):
    global _ACTIVE_LAYER_COEFFS, _LAYER_BOUNDS
    _ACTIVE_LAYER_COEFFS = np.asarray(layer_coeffs, dtype=np.double).copy()
    _LAYER_BOUNDS = np.asarray(problem["layer_bounds"], dtype=np.double).copy()
    adjoint_debug_dump = bool(getattr(args, "adjoint_debug_dump", False))
    adjoint_debug_threshold = float(getattr(args, "adjoint_debug_threshold", 0.0))
    adjoint_debug_max_rows = int(getattr(args, "adjoint_debug_max_rows", 200000))
    adjoint_debug_stride = int(getattr(args, "adjoint_debug_stride", 1))
    dual_overflow_warn = bool(getattr(args, "dual_overflow_warn", False))
    dual_overflow_warn_fraction = float(getattr(args, "dual_overflow_warn_fraction", 0.95))
    dual_overflow_warn_limit = int(getattr(args, "dual_overflow_warn_limit", 20))

    REMAT.API.set_integrator_type(args.integrator_type.encode("utf-8"))
    REMAT.API.define_parameter(b"body_force_x", 0.0)
    REMAT.API.define_parameter(b"body_force_y", 0.0)
    REMAT.API.define_parameter(b"mass_damping_factor", args.mass_damping_factor)
    REMAT.API.define_parameter(b"contact_stiffness", 0.0)
    REMAT.API.define_parameter(b"overflow_limit", float(args.overflow_limit))
    REMAT.API.define_parameter(b"mat_overflow_limit", float(args.mat_overflow_limit))
    REMAT.API.define_parameter(b"adjoint_material_objective_weight", 0.0)
    REMAT.API.define_parameter(b"adjoint_debug_dump_enable", 1.0 if adjoint_debug_dump else 0.0)
    REMAT.API.define_parameter(b"adjoint_debug_dump_threshold", float(adjoint_debug_threshold))
    REMAT.API.define_parameter(b"adjoint_debug_dump_max_rows", float(adjoint_debug_max_rows))
    REMAT.API.define_parameter(b"adjoint_debug_dump_stride", float(adjoint_debug_stride))
    REMAT.API.define_parameter(b"dual_overflow_warn_enable", 1.0 if dual_overflow_warn else 0.0)
    REMAT.API.define_parameter(b"dual_overflow_warn_fraction", float(dual_overflow_warn_fraction))
    REMAT.API.define_parameter(b"dual_overflow_warn_limit", float(dual_overflow_warn_limit))

    REMAT.API.define_parameter(b"density", args.density)
    REMAT.API.define_parameter(b"youngs_modulus", args.youngs_modulus)
    REMAT.API.define_parameter(b"poissons_ratio", args.poissons_ratio)
    REMAT.API.define_parameter(b"relaxation_time", tau)
    REMAT.API.define_parameter(b"shear_modulus_Maxwell_element", args.shear_modulus_maxwell)

    truss_connectivity = np.zeros((0, 2), dtype=np.int32)
    REMAT.create_geometry(
        problem["coordinates"],
        problem["velocities"],
        problem["fixity"],
        problem["connectivity"],
        [],
        truss_connectivity,
    )
    REMAT.define_variable_properties(layered_stiffness_scaling)
    REMAT.API.initialize()


def _sum_by_layer(values, elem_layer_ids, n_layers):
    out = np.zeros(n_layers, dtype=np.double)
    for lid in range(n_layers):
        out[lid] = float(np.sum(values[elem_layer_ids == lid]))
    return out


def _l2_by_layer(values, elem_layer_ids, n_layers):
    out = np.zeros(n_layers, dtype=np.double)
    for lid in range(n_layers):
        val = values[elem_layer_ids == lid]
        out[lid] = float(np.linalg.norm(val))
    return out


def _maxabs_by_layer(values, elem_layer_ids, n_layers):
    out = np.zeros(n_layers, dtype=np.double)
    for lid in range(n_layers):
        val = values[elem_layer_ids == lid]
        out[lid] = float(np.max(np.abs(val))) if val.size else 0.0
    return out


def _node_field_or_zeros(problem, field_name):
    values = REMAT.get_field(b"node", field_name)
    if values is None:
        return np.zeros(problem["coordinates"].shape[0], dtype=np.double)
    return np.asarray(values, dtype=np.double).reshape(-1)


def _adjoint_node_field(problem, primary_name, dual_fallback_name):
    values = REMAT.get_field(b"node", primary_name)
    if values is not None:
        return np.asarray(values, dtype=np.double).reshape(-1)
    return _node_field_or_zeros(problem, dual_fallback_name)


def _capture_adjoint_snapshot(problem, n_layers):
    sensor_nodes = problem["sensor_nodes"]
    elem_layer_ids = problem["elem_layer_ids"]

    dual_vy = _node_field_or_zeros(problem, "dual_velocity_Y")
    dual_uy = _node_field_or_zeros(problem, "dual_displacement_Y")
    adjoint_vy = _adjoint_node_field(problem, "adjoint_velocity_Y", "dual_velocity_Y")
    adjoint_uy = _adjoint_node_field(problem, "adjoint_displacement_Y", "dual_displacement_Y")

    lambda_xx = REMAT.get_field(b"element", "lambda_q_xx")
    lambda_yy = REMAT.get_field(b"element", "lambda_q_yy")
    lambda_xy = REMAT.get_field(b"element", "lambda_q_xy")
    lambda_mag = np.sqrt(lambda_xx * lambda_xx + lambda_yy * lambda_yy + lambda_xy * lambda_xy)

    elem_grad_scaling = REMAT.get_field(b"element", "dparam_stiffness_scaling_factor")
    layer_grad = _sum_by_layer(elem_grad_scaling, elem_layer_ids, n_layers)
    layer_lambda_l2 = _l2_by_layer(lambda_mag, elem_layer_ids, n_layers)
    layer_lambda_maxabs = _maxabs_by_layer(lambda_mag, elem_layer_ids, n_layers)

    return {
        "grad_tau": float(REMAT.get_field(b"global", "dL_dparam_relaxation_time")[0]),
        "grad_layers": layer_grad,
        "adjoint_velocity_y_maxabs_all_nodes": float(np.max(np.abs(adjoint_vy))) if adjoint_vy.size else 0.0,
        "adjoint_velocity_y_maxabs_sensors": float(np.max(np.abs(adjoint_vy[sensor_nodes]))) if sensor_nodes.size else 0.0,
        "adjoint_displacement_y_maxabs_all_nodes": float(np.max(np.abs(adjoint_uy))) if adjoint_uy.size else 0.0,
        "dual_velocity_y_maxabs_all_nodes": float(np.max(np.abs(dual_vy))) if dual_vy.size else 0.0,
        "dual_velocity_y_maxabs_sensors": float(np.max(np.abs(dual_vy[sensor_nodes]))) if sensor_nodes.size else 0.0,
        "dual_displacement_y_maxabs_all_nodes": float(np.max(np.abs(dual_uy))) if dual_uy.size else 0.0,
        "lambda_q_l2_by_layer": layer_lambda_l2,
        "lambda_q_maxabs_by_layer": layer_lambda_maxabs,
    }


def _write_trace_csv(trace, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    nsteps = int(trace["obs_step"].size)
    n_layers = int(trace["grad_layers_cum"].shape[1])
    order = np.argsort(trace["obs_step"])
    header = [
        "obs_step",
        "obs_time",
        "residual_l2",
        "residual_rms",
        "dual_vy_sensor_max_pre",
        "dual_vy_sensor_max_post",
        "dual_vy_all_max_post",
        "dual_uy_all_max_post",
        "adjoint_vy_sensor_max_pre",
        "adjoint_vy_sensor_max_post",
        "adjoint_vy_all_max_post",
        "adjoint_uy_all_max_post",
        "grad_tau_cum",
        "grad_tau_delta",
    ]
    for lid in range(n_layers):
        header.append(f"grad_layer{lid}_cum")
    for lid in range(n_layers):
        header.append(f"grad_layer{lid}_delta")
    for lid in range(n_layers):
        header.append(f"lambda_layer{lid}_l2")
    for lid in range(n_layers):
        header.append(f"lambda_layer{lid}_maxabs")

    with path.open("w", encoding="utf-8") as f:
        f.write(",".join(header) + "\n")
        for idx in order:
            row = [
                str(int(trace["obs_step"][idx])),
                f"{float(trace['obs_time'][idx]):.12e}",
                f"{float(trace['residual_l2'][idx]):.12e}",
                f"{float(trace['residual_rms'][idx]):.12e}",
                f"{float(trace['dual_vy_sensor_max_pre'][idx]):.12e}",
                f"{float(trace['dual_vy_sensor_max_post'][idx]):.12e}",
                f"{float(trace['dual_vy_all_max_post'][idx]):.12e}",
                f"{float(trace['dual_uy_all_max_post'][idx]):.12e}",
                f"{float(trace['adjoint_vy_sensor_max_pre'][idx]):.12e}",
                f"{float(trace['adjoint_vy_sensor_max_post'][idx]):.12e}",
                f"{float(trace['adjoint_vy_all_max_post'][idx]):.12e}",
                f"{float(trace['adjoint_uy_all_max_post'][idx]):.12e}",
                f"{float(trace['grad_tau_cum'][idx]):.12e}",
                f"{float(trace['grad_tau_delta'][idx]):.12e}",
            ]
            for lid in range(n_layers):
                row.append(f"{float(trace['grad_layers_cum'][idx, lid]):.12e}")
            for lid in range(n_layers):
                row.append(f"{float(trace['grad_layers_delta'][idx, lid]):.12e}")
            for lid in range(n_layers):
                row.append(f"{float(trace['lambda_layers_l2'][idx, lid]):.12e}")
            for lid in range(n_layers):
                row.append(f"{float(trace['lambda_layers_maxabs'][idx, lid]):.12e}")
            f.write(",".join(row) + "\n")
    return path


def _plot_trace(trace, n_layers, label, output_file):
    obs_step = trace["obs_step"]
    order_t = np.argsort(obs_step)
    t = trace["obs_time"][order_t]

    rev = np.arange(obs_step.size)

    fig, axes = plt.subplots(2, 2, figsize=(12.8, 8.4))
    ax_res, ax_grad_inc, ax_grad_cum, ax_adj = axes.ravel()

    ax_res.plot(t, trace["residual_rms"][order_t], color="#b03a2e", linewidth=1.8)
    ax_res.set_title(f"Residual RMS vs Time ({label})")
    ax_res.set_xlabel("time")
    ax_res.set_ylabel("RMS residual at sensors")
    ax_res.grid(True, alpha=0.25)

    for lid in range(n_layers):
        ax_grad_inc.semilogy(
            t,
            np.abs(trace["grad_layers_delta"][order_t, lid]) + 1.0e-30,
            linewidth=1.7,
            label=f"|dL/dlayer{lid}| step contrib",
        )
    ax_grad_inc.semilogy(
        t,
        np.abs(trace["grad_tau_delta"][order_t]) + 1.0e-30,
        color="#111111",
        linewidth=1.3,
        linestyle="--",
        label="|dL/dtau| step contrib",
    )
    ax_grad_inc.set_title("Per-Time-Step Gradient Contribution Magnitude")
    ax_grad_inc.set_xlabel("time")
    ax_grad_inc.set_ylabel("absolute contribution")
    ax_grad_inc.grid(True, alpha=0.25)
    ax_grad_inc.legend(fontsize=8)

    for lid in range(n_layers):
        ax_grad_cum.plot(rev, trace["grad_layers_cum"][:, lid], linewidth=1.8, label=f"layer {lid}")
    ax_grad_cum.plot(rev, trace["grad_tau_cum"], color="#111111", linewidth=1.3, linestyle="--", label="tau")
    ax_grad_cum.set_title("Cumulative Gradient During Reverse Sweep")
    ax_grad_cum.set_xlabel("reverse-sweep index (0 = last observation)")
    ax_grad_cum.set_ylabel("cumulative gradient")
    ax_grad_cum.grid(True, alpha=0.25)
    ax_grad_cum.legend(fontsize=8)

    for lid in range(n_layers):
        ax_adj.semilogy(
            t,
            trace["lambda_layers_l2"][order_t, lid] + 1.0e-30,
            linewidth=1.8,
            label=f"layer {lid} ||lambda_q||_2",
        )
    sensor_global_adjoint = (
        trace["adjoint_vy_sensor_max_post"][order_t]
        if "adjoint_vy_sensor_max_post" in trace
        else trace["dual_vy_sensor_max_post"][order_t]
    )
    ax_adj.semilogy(t, sensor_global_adjoint + 1.0e-30, color="#111111", linewidth=1.2, linestyle="--",
                    label="max|v*| sensors")
    ax_adj.set_title("Adjoint State Magnitude vs Time")
    ax_adj.set_xlabel("time")
    ax_adj.set_ylabel("magnitude")
    ax_adj.grid(True, alpha=0.25)
    ax_adj.legend(fontsize=8)

    fig.tight_layout()
    output_file = Path(output_file)
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_file, dpi=180)
    plt.close(fig)
    return output_file


def _save_diagnostics_bundle(run, problem, args, label, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    trace = run["adjoint_trace"]
    csv_path = _write_trace_csv(trace, output_dir / f"{label}_adjoint_trace.csv")
    fig_path = _plot_trace(trace, args.n_layers, label, output_dir / f"{label}_adjoint_trace.png")

    endpoint = {
        "label": label,
        "loss": float(run["loss"]),
        "grad_tau": float(run["grad_tau"]),
        "grad_layers": run["grad_layers"].tolist(),
        "terminal_before_seeding": {
            k: (v.tolist() if isinstance(v, np.ndarray) else float(v))
            for k, v in trace["endpoint_terminal"].items()
        },
        "initial_time_after_full_reverse": {
            k: (v.tolist() if isinstance(v, np.ndarray) else float(v))
            for k, v in trace["endpoint_initial"].items()
        },
    }

    endpoint_path = output_dir / f"{label}_adjoint_endpoints.json"
    with endpoint_path.open("w", encoding="utf-8") as f:
        json.dump(endpoint, f, indent=2)

    return csv_path, fig_path, endpoint_path


def _fd_report(problem, layer_coeffs, tau, args, observed_history):
    adj = run_forward_or_adjoint(
        problem,
        layer_coeffs,
        tau,
        args,
        observed_history=observed_history,
        compute_gradients=True,
    )
    fd_layers, fd_tau = finite_difference_gradients(problem, layer_coeffs, tau, args, observed_history, args.fd_step)
    eps = 1.0e-14
    rel_layer = np.abs(adj["grad_layers"] - fd_layers) / np.maximum.reduce(
        [np.abs(adj["grad_layers"]), np.abs(fd_layers), np.full_like(fd_layers, eps)]
    )
    rel_tau = abs(adj["grad_tau"] - fd_tau) / max(abs(adj["grad_tau"]), abs(fd_tau), eps)
    return {
        "adjoint_layers": adj["grad_layers"],
        "adjoint_tau": float(adj["grad_tau"]),
        "fd_layers": fd_layers,
        "fd_tau": float(fd_tau),
        "rel_error_layers": rel_layer,
        "rel_error_tau": float(rel_tau),
    }


def run_forward_or_adjoint(
    problem,
    layer_coeffs,
    tau,
    args,
    observed_history=None,
    compute_gradients=False,
    collect_adjoint_trace=False,
):
    configure_run(problem, layer_coeffs, tau, args)

    sensor_nodes = problem["sensor_nodes"]
    nsensors = sensor_nodes.size
    sensor_history = np.zeros((args.nsteps, nsensors), dtype=np.double)
    loss = 0.0

    for k in range(args.nsteps):
        REMAT.API.update_state(args.dt, args.nsub_steps, REMAT.PASS_FORWARD)
        vy = REMAT.get_field(b"node", "velocity_Y")[sensor_nodes]
        sensor_history[k, :] = vy
        if observed_history is not None:
            residual = vy - observed_history[k, :]
            loss += 0.5 * float(np.dot(residual, residual))

    grad_tau = 0.0
    grad_layers = np.zeros_like(layer_coeffs)
    adjoint_trace = None
    if compute_gradients:
        if observed_history is None:
            raise ValueError("compute_gradients=True requires observed_history.")

        n_layers = int(layer_coeffs.size)
        elem_layer_ids = problem["elem_layer_ids"]
        REMAT.clear_adjoint_state()
        grad_tau_prev = 0.0
        grad_layers_prev = np.zeros(n_layers, dtype=np.double)

        if collect_adjoint_trace:
            adjoint_trace = {
                "obs_step": np.zeros(args.nsteps, dtype=np.int32),
                "obs_time": np.zeros(args.nsteps, dtype=np.double),
                "residual_l2": np.zeros(args.nsteps, dtype=np.double),
                "residual_rms": np.zeros(args.nsteps, dtype=np.double),
                "dual_vy_sensor_max_pre": np.zeros(args.nsteps, dtype=np.double),
                "dual_vy_sensor_max_post": np.zeros(args.nsteps, dtype=np.double),
                "dual_vy_all_max_post": np.zeros(args.nsteps, dtype=np.double),
                "dual_uy_all_max_post": np.zeros(args.nsteps, dtype=np.double),
                "adjoint_vy_sensor_max_pre": np.zeros(args.nsteps, dtype=np.double),
                "adjoint_vy_sensor_max_post": np.zeros(args.nsteps, dtype=np.double),
                "adjoint_vy_all_max_post": np.zeros(args.nsteps, dtype=np.double),
                "adjoint_uy_all_max_post": np.zeros(args.nsteps, dtype=np.double),
                "grad_tau_cum": np.zeros(args.nsteps, dtype=np.double),
                "grad_tau_delta": np.zeros(args.nsteps, dtype=np.double),
                "grad_layers_cum": np.zeros((args.nsteps, n_layers), dtype=np.double),
                "grad_layers_delta": np.zeros((args.nsteps, n_layers), dtype=np.double),
                "lambda_layers_l2": np.zeros((args.nsteps, n_layers), dtype=np.double),
                "lambda_layers_maxabs": np.zeros((args.nsteps, n_layers), dtype=np.double),
                "endpoint_terminal": _capture_adjoint_snapshot(problem, n_layers),
                "endpoint_initial": None,
            }

        for rev in range(args.nsteps):
            k = args.nsteps - 1 - rev
            residual = sensor_history[k, :] - observed_history[k, :]
            seed_xy = np.zeros((nsensors, 2), dtype=np.double)
            seed_xy[:, 1] = residual
            REMAT.add_nodal_velocity_adjoint_seed(sensor_nodes, seed_xy)

            if collect_adjoint_trace:
                dual_vy_pre = _node_field_or_zeros(problem, "dual_velocity_Y")
                adjoint_vy_pre = _adjoint_node_field(problem, "adjoint_velocity_Y", "dual_velocity_Y")
                adjoint_trace["obs_step"][rev] = int(k)
                adjoint_trace["obs_time"][rev] = float(k * args.dt)
                adjoint_trace["residual_l2"][rev] = float(np.linalg.norm(residual))
                adjoint_trace["residual_rms"][rev] = float(np.sqrt(np.mean(residual * residual)))
                adjoint_trace["dual_vy_sensor_max_pre"][rev] = (
                    float(np.max(np.abs(dual_vy_pre[sensor_nodes]))) if sensor_nodes.size else 0.0
                )
                adjoint_trace["adjoint_vy_sensor_max_pre"][rev] = (
                    float(np.max(np.abs(adjoint_vy_pre[sensor_nodes]))) if sensor_nodes.size else 0.0
                )

            REMAT.API.update_state(args.dt, args.nsub_steps, REMAT.PASS_BACKWARD_ADJOINT)

            if collect_adjoint_trace:
                dual_vy = _node_field_or_zeros(problem, "dual_velocity_Y")
                dual_uy = _node_field_or_zeros(problem, "dual_displacement_Y")
                adjoint_vy = _adjoint_node_field(problem, "adjoint_velocity_Y", "dual_velocity_Y")
                adjoint_uy = _adjoint_node_field(problem, "adjoint_displacement_Y", "dual_displacement_Y")
                elem_grad_scaling = REMAT.get_field(b"element", "dparam_stiffness_scaling_factor")
                grad_tau_now = float(REMAT.get_field(b"global", "dL_dparam_relaxation_time")[0])
                grad_layers_now = _sum_by_layer(elem_grad_scaling, elem_layer_ids, n_layers)

                lambda_xx = REMAT.get_field(b"element", "lambda_q_xx")
                lambda_yy = REMAT.get_field(b"element", "lambda_q_yy")
                lambda_xy = REMAT.get_field(b"element", "lambda_q_xy")
                lambda_mag = np.sqrt(lambda_xx * lambda_xx + lambda_yy * lambda_yy + lambda_xy * lambda_xy)

                adjoint_trace["dual_vy_sensor_max_post"][rev] = (
                    float(np.max(np.abs(dual_vy[sensor_nodes]))) if sensor_nodes.size else 0.0
                )
                adjoint_trace["dual_vy_all_max_post"][rev] = float(np.max(np.abs(dual_vy))) if dual_vy.size else 0.0
                adjoint_trace["dual_uy_all_max_post"][rev] = float(np.max(np.abs(dual_uy))) if dual_uy.size else 0.0
                adjoint_trace["adjoint_vy_sensor_max_post"][rev] = (
                    float(np.max(np.abs(adjoint_vy[sensor_nodes]))) if sensor_nodes.size else 0.0
                )
                adjoint_trace["adjoint_vy_all_max_post"][rev] = (
                    float(np.max(np.abs(adjoint_vy))) if adjoint_vy.size else 0.0
                )
                adjoint_trace["adjoint_uy_all_max_post"][rev] = (
                    float(np.max(np.abs(adjoint_uy))) if adjoint_uy.size else 0.0
                )

                adjoint_trace["grad_tau_cum"][rev] = grad_tau_now
                adjoint_trace["grad_tau_delta"][rev] = grad_tau_now - grad_tau_prev
                adjoint_trace["grad_layers_cum"][rev, :] = grad_layers_now
                adjoint_trace["grad_layers_delta"][rev, :] = grad_layers_now - grad_layers_prev
                adjoint_trace["lambda_layers_l2"][rev, :] = _l2_by_layer(lambda_mag, elem_layer_ids, n_layers)
                adjoint_trace["lambda_layers_maxabs"][rev, :] = _maxabs_by_layer(lambda_mag, elem_layer_ids, n_layers)

                grad_tau_prev = grad_tau_now
                grad_layers_prev = grad_layers_now

        if collect_adjoint_trace:
            adjoint_trace["endpoint_initial"] = _capture_adjoint_snapshot(problem, n_layers)
            grad_tau = float(grad_tau_prev)
            grad_layers = grad_layers_prev.copy()
        else:
            grad_tau = float(REMAT.get_field(b"global", "dL_dparam_relaxation_time")[0])
            elem_grad_scaling = REMAT.get_field(b"element", "dparam_stiffness_scaling_factor")
            grad_layers = _sum_by_layer(elem_grad_scaling, elem_layer_ids, n_layers)

    result = {
        "loss": float(loss),
        "sensor_history": sensor_history,
        "grad_tau": float(grad_tau),
        "grad_layers": grad_layers,
    }
    if adjoint_trace is not None:
        result["adjoint_trace"] = adjoint_trace
    return result


def finite_difference_gradients(problem, layer_coeffs, tau, args, observed_history, h):
    fd_layers = np.zeros_like(layer_coeffs)
    for i in range(layer_coeffs.size):
        plus = layer_coeffs.copy()
        minus = layer_coeffs.copy()
        plus[i] += h
        minus[i] -= h
        loss_plus = run_forward_or_adjoint(problem, plus, tau, args, observed_history, False)["loss"]
        loss_minus = run_forward_or_adjoint(problem, minus, tau, args, observed_history, False)["loss"]
        fd_layers[i] = (loss_plus - loss_minus) / (2.0 * h)

    loss_plus_tau = run_forward_or_adjoint(problem, layer_coeffs, tau + h, args, observed_history, False)["loss"]
    loss_minus_tau = run_forward_or_adjoint(problem, layer_coeffs, tau - h, args, observed_history, False)["loss"]
    fd_tau = (loss_plus_tau - loss_minus_tau) / (2.0 * h)
    return fd_layers, float(fd_tau)


def step_profile(ax, coeffs, layer_bounds, label, color):
    y = layer_bounds
    values = np.zeros_like(y)
    values[:-1] = coeffs
    values[-1] = coeffs[-1]
    ax.step(values, y, where="post", label=label, color=color, linewidth=2.0)


def plot_problem_setup(problem, true_layers, true_tau, args, output_file):
    coords = problem["coordinates"]
    sensor_nodes = problem["sensor_nodes"]
    impact_nodes = problem["impact_nodes"]
    layer_bounds = problem["layer_bounds"]
    width = float(problem["width"])
    height = float(problem["height"])

    fig = plt.figure(figsize=(12.5, 5.3))
    grid = fig.add_gridspec(1, 2, width_ratios=[3.3, 1.7])
    ax = fig.add_subplot(grid[0, 0])
    ax_info = fig.add_subplot(grid[0, 1])

    cmin = float(np.min(true_layers))
    cmax = float(np.max(true_layers))
    cmap = plt.get_cmap("YlGnBu")
    coeff_span = max(cmax - cmin, 1.0e-12)
    for lid in range(true_layers.size):
        y0 = layer_bounds[lid]
        y1 = layer_bounds[lid + 1]
        shade = (float(true_layers[lid]) - cmin) / coeff_span
        ax.axhspan(y0, y1, facecolor=cmap(0.25 + 0.65 * shade), alpha=0.55, zorder=0)
        ax.text(
            0.02 * width,
            0.5 * (y0 + y1),
            f"layer {lid}: {true_layers[lid]:.3f}",
            fontsize=8.8,
            va="center",
            ha="left",
            color="#202020",
            bbox={"boxstyle": "round,pad=0.2", "facecolor": "white", "edgecolor": "none", "alpha": 0.65},
        )

    for yb in layer_bounds:
        ax.hlines(yb, 0.0, width, colors="#2f2f2f", linewidth=0.9, linestyles="--", alpha=0.5, zorder=1)

    # Fixed boundaries in this setup are bottom, left, and right.
    ax.plot([0.0, width], [0.0, 0.0], color="#111111", linewidth=2.6, label="fixed boundary", zorder=2)
    ax.plot([0.0, 0.0], [0.0, height], color="#111111", linewidth=2.6, zorder=2)
    ax.plot([width, width], [0.0, height], color="#111111", linewidth=2.6, zorder=2)

    sensor_xy = coords[sensor_nodes]
    ax.scatter(
        sensor_xy[:, 0],
        sensor_xy[:, 1],
        s=44,
        c="#d1495b",
        edgecolors="white",
        linewidths=0.7,
        label=f"sensors ({sensor_nodes.size})",
        zorder=6,
    )

    impact_xy = coords[impact_nodes]
    ax.scatter(
        impact_xy[:, 0],
        impact_xy[:, 1],
        s=36,
        marker="s",
        c="#205a8e",
        edgecolors="white",
        linewidths=0.6,
        label=f"impact patch ({impact_nodes.size})",
        zorder=7,
    )

    impact_x = float(np.mean(impact_xy[:, 0])) if impact_xy.size else float(problem["impact_center_x"])
    arrow_length = (0.12 + 0.08 * np.tanh(abs(args.impact_velocity))) * height
    y_start = height + 0.17 * height
    y_end = height + 0.17 * height - arrow_length
    ax.annotate(
        "",
        xy=(impact_x, y_end),
        xytext=(impact_x, y_start),
        arrowprops={"arrowstyle": "-|>", "color": "#205a8e", "lw": 2.4},
        annotation_clip=False,
    )
    ax.text(
        impact_x + 0.02 * width,
        height + 0.11 * height,
        f"impact Vy = {-abs(args.impact_velocity):.3f}",
        color="#205a8e",
        fontsize=10,
        ha="left",
        va="center",
    )

    ax.set_xlim(-0.03 * width, 1.03 * width)
    ax.set_ylim(-0.06 * height, 1.28 * height)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title("Problem Setup: Geometry, Sensors, True Layers, and Impact")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, alpha=0.2, linestyle=":")
    ax.legend(loc="upper right", fontsize=9, framealpha=0.95)

    norm = plt.Normalize(vmin=cmin, vmax=cmax if cmax > cmin else cmin + 1.0)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.045, pad=0.02)
    cbar.set_label("true stiffness scaling")

    ax_info.axis("off")
    info_lines = [
        "Inverse Problem Definition",
        "",
        "Recover from sensor traces:",
        f"  stiffness_scaling_factor (layers={args.n_layers})",
        "  relaxation_time (tau)",
        "",
        "Observation model:",
        f"  {sensor_nodes.size} sensors on top surface",
        f"  sensor band width = {problem['sensor_distribution_width']:.3f}",
        (
            f"  sensor x-range = "
            f"[{problem['sensor_distribution_x_min']:.3f}, {problem['sensor_distribution_x_max']:.3f}]"
        ),
        "  measured field: velocity_Y(t)",
        "",
        "Initial condition:",
        f"  top impact patch width = {problem['impact_window_width']:.3f}",
        f"  impact x-range = [{problem['impact_x_min']:.3f}, {problem['impact_x_max']:.3f}]",
        f"  impact center x = {problem['impact_center_x']:.3f}",
        f"  impact velocity_Y = {-abs(args.impact_velocity):.3f}",
        "",
        "Numerical setup:",
        f"  mesh = {args.nx} x {args.ny} quads",
        f"  nsteps = {args.nsteps}, dt = {args.dt:.4e}",
        f"  true tau = {true_tau:.4f}",
        f"  integrator = {args.integrator_type}",
    ]
    ax_info.text(
        0.0,
        1.0,
        "\n".join(info_lines),
        ha="left",
        va="top",
        fontsize=10.0,
        family="monospace",
    )

    fig.tight_layout()
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, format="svg")
    plt.close(fig)
    return output_path


def parse_layer_values(text, n_layers, default_value):
    if text is None:
        return np.full(n_layers, default_value, dtype=np.double)
    values = np.fromstring(text, sep=",", dtype=np.double)
    if values.size != n_layers:
        raise ValueError(f"Expected {n_layers} layer values, got {values.size}")
    return values


def pack_params(layers, tau):
    return np.concatenate([np.asarray(layers, dtype=np.double), np.array([float(tau)], dtype=np.double)])


def unpack_params(params, n_layers):
    params = np.asarray(params, dtype=np.double)
    return params[:n_layers].copy(), float(params[n_layers])


def optimize_with_lbfgsb(problem, observed_history, init_layers, init_tau, args):
    try:
        from scipy.optimize import minimize
    except ImportError as exc:
        raise ImportError(
            "L-BFGS-B optimizer requires SciPy. Install it with: pip install scipy"
        ) from exc

    x0 = pack_params(init_layers, init_tau)
    bounds = [(args.min_layer, args.max_layer)] * args.n_layers + [(args.min_tau, args.max_tau)]

    loss_history = []
    tau_history = []
    layer_history = []
    iteration_counter = {"value": 0}
    cache = {"z": None, "x": None, "loss": None, "grad_z": None, "grad_x": None}

    def objective_with_grad_x(x):
        layers, tau = unpack_params(x, args.n_layers)
        run = run_forward_or_adjoint(
            problem,
            layers,
            tau,
            args,
            observed_history=observed_history,
            compute_gradients=True,
        )
        grad = np.concatenate(
            [run["grad_layers"], np.array([run["grad_tau"]], dtype=np.double)]
        )
        return float(run["loss"]), grad

    # Build a static diagonal scaling from the initial gradient (layers only).
    # This is a simple, defensible preconditioning via variable reparameterization.
    loss0, grad0 = objective_with_grad_x(x0)
    grad_layers0 = np.abs(grad0[:args.n_layers])
    eps = 1.0e-14
    target = float(np.median(np.maximum(grad_layers0, eps))) if args.n_layers > 0 else 1.0
    layer_scale = target / np.maximum(grad_layers0, eps)
    layer_scale = np.clip(layer_scale, 0.3, 3.0)
    param_scale = np.ones_like(x0)
    param_scale[:args.n_layers] = layer_scale

    z0 = x0 / param_scale
    bounds_z = [(lo / param_scale[i], hi / param_scale[i]) for i, (lo, hi) in enumerate(bounds)]

    print("Gradient balancing (static layer scales):", np.array2string(layer_scale, precision=3))

    def objective_with_grad(z):
        z = np.asarray(z, dtype=np.double)
        x = z * param_scale
        loss, grad_x = objective_with_grad_x(x)
        grad_z = grad_x * param_scale
        cache["z"] = z.copy()
        cache["x"] = x.copy()
        cache["loss"] = float(loss)
        cache["grad_z"] = grad_z.copy()
        cache["grad_x"] = grad_x.copy()
        return cache["loss"], grad_z

    def log_iteration(z):
        z = np.asarray(z, dtype=np.double)
        if cache["z"] is None or not np.array_equal(z, cache["z"]):
            loss, _ = objective_with_grad(z)
        else:
            loss = float(cache["loss"])
        x = cache["x"]
        grad = cache["grad_x"]

        layers, tau = unpack_params(x, args.n_layers)
        grad_layers = grad[:args.n_layers]
        grad_tau = float(grad[args.n_layers])

        it = iteration_counter["value"]
        print(
            f"iter={it:03d}  loss={loss: .8e}  tau={tau: .8e}  "
            f"|grad_layers|_inf={np.max(np.abs(grad_layers)): .8e}  grad_tau={grad_tau: .8e}"
        )

        loss_history.append(loss)
        tau_history.append(tau)
        layer_history.append(layers.copy())
        iteration_counter["value"] += 1

    # Log the starting point and optionally run FD check at the same point.
    cache["z"] = z0.copy()
    cache["x"] = x0.copy()
    cache["loss"] = float(loss0)
    cache["grad_x"] = grad0.copy()
    cache["grad_z"] = grad0 * param_scale
    log_iteration(z0)

    if args.fd_check:
        layers0, tau0 = unpack_params(x0, args.n_layers)
        fd_layers, fd_tau = finite_difference_gradients(
            problem, layers0, tau0, args, observed_history, args.fd_step
        )
        eps = 1.0e-14
        rel_layer = np.abs(grad0[:args.n_layers] - fd_layers) / np.maximum.reduce(
            [np.abs(grad0[:args.n_layers]), np.abs(fd_layers), np.full_like(fd_layers, eps)]
        )
        rel_tau = abs(grad0[args.n_layers] - fd_tau) / max(abs(grad0[args.n_layers]), abs(fd_tau), eps)
        print("FD check (layer gradients):")
        for i in range(args.n_layers):
            print(
                f"  layer {i}: adj={grad0[i]: .8e}  fd={fd_layers[i]: .8e}  rel_err={rel_layer[i]: .3e}"
            )
        print(
            f"FD check (tau): adj={grad0[args.n_layers]: .8e}  fd={fd_tau: .8e}  rel_err={rel_tau: .3e}"
        )

    total_nit = 0
    total_nfev = 0
    current_z = z0.copy()
    result = None
    best_result = None
    max_runs = max(1, int(args.lbfgsb_restarts) + 1)

    for run_id in range(max_runs):
        remaining_iters = max(1, int(args.max_iters) - total_nit)
        result = minimize(
            objective_with_grad,
            x0=current_z,
            method="L-BFGS-B",
            jac=True,
            bounds=bounds_z,
            callback=log_iteration,
            options={
                "maxiter": remaining_iters,
                "gtol": float(args.lbfgsb_gtol),
                "maxls": int(args.lbfgsb_maxls),
            },
        )
        total_nit += int(result.nit)
        total_nfev += int(result.nfev)

        if best_result is None or float(result.fun) < float(best_result.fun):
            best_result = result

        if total_nit >= int(args.max_iters):
            break
        if run_id + 1 >= max_runs:
            break

        message = str(result.message)
        if "PROJECTED GRADIENT" not in message.upper():
            break

        z = np.asarray(result.x, dtype=np.double).copy()
        z_pert = z.copy()
        had_active_bound = False
        for i, (lo, hi) in enumerate(bounds_z):
            span = max(float(hi) - float(lo), 1.0e-12)
            eps = float(args.lbfgsb_boundary_perturb) * span
            tol = 1.0e-12 * max(1.0, abs(float(lo)), abs(float(hi)))
            if z[i] <= float(lo) + tol:
                z_pert[i] = min(float(hi), float(lo) + eps)
                had_active_bound = True
            elif z[i] >= float(hi) - tol:
                z_pert[i] = max(float(lo), float(hi) - eps)
                had_active_bound = True

        if not had_active_bound:
            break

        print(
            f"lbfgsb-restart run={run_id + 1} "
            f"total_nit={total_nit}  fun={float(result.fun): .8e}"
        )
        current_z = z_pert

    if result is not None and best_result is not None:
        best_result.nit = total_nit
        best_result.nfev = total_nfev
        best_result.x = np.asarray(best_result.x, dtype=np.double) * param_scale
    return best_result, loss_history, tau_history, layer_history


def main():
    parser = argparse.ArgumentParser(
        description="Inverse dissipative-wave example: recover layered stiffness and relaxation_time from top-surface velocity sensors."
    )
    parser.add_argument("--nx", type=int, default=48)
    parser.add_argument("--ny", type=int, default=18)
    parser.add_argument("--width", type=float, default=15.0)
    parser.add_argument("--height", type=float, default=3.0)
    parser.add_argument("--nsteps", type=int, default=80)
    parser.add_argument("--nsub-steps", type=int, default=6)
    parser.add_argument("--dt", type=float, default=2.5e-3)
    parser.add_argument("--n-layers", type=int, default=4)
    parser.add_argument("--n-sensors", type=int, default=10)
    parser.add_argument(
        "--sensor-distribution-width",
        type=float,
        default=None,
        help="Centered top-surface width over which sensors are distributed (default: use full domain width).",
    )
    parser.add_argument(
        "--impact-window-width",
        type=float,
        default=1.8,
        help="Centered top-surface impact width (absolute length units).",
    )
    parser.add_argument("--impact-velocity", type=float, default=0.35)
    parser.add_argument("--integrator-type", type=str, default="fixed_visco")

    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--youngs-modulus", type=float, default=5.0)
    parser.add_argument("--poissons-ratio", type=float, default=0.28)
    parser.add_argument("--shear-modulus-maxwell", type=float, default=2.0)
    parser.add_argument("--mass-damping-factor", type=float, default=0.0)
    parser.add_argument("--overflow-limit", type=int, default=50)
    parser.add_argument("--mat-overflow-limit", type=int, default=50)
    parser.add_argument(
        "--adjoint-debug-dump",
        action="store_true",
        default=False,
        help=(
            "Enable C++ per-quadrature debug dump for stiffness-scaling adjoint contributions "
            "(writes adjoint_stiffness_qp_debug.csv in the working directory)."
        ),
    )
    parser.add_argument(
        "--adjoint-debug-threshold",
        type=float,
        default=0.0,
        help="Absolute delta threshold for writing rows to adjoint_stiffness_qp_debug.csv.",
    )
    parser.add_argument(
        "--adjoint-debug-max-rows",
        type=int,
        default=200000,
        help="Maximum number of debug rows written per run for adjoint_stiffness_qp_debug.csv.",
    )
    parser.add_argument(
        "--adjoint-debug-stride",
        type=int,
        default=1,
        help="Element stride for adjoint stiffness debug dump (1 = every element).",
    )
    parser.add_argument(
        "--dual-overflow-warn",
        action="store_true",
        default=False,
        help=(
            "Enable warnings when dual/ancillary fixed-point states approach int32 range "
            "before overflow checkpoint storage."
        ),
    )
    parser.add_argument(
        "--dual-overflow-warn-fraction",
        type=float,
        default=0.95,
        help="Warning threshold as a fraction of int32 max mantissa.",
    )
    parser.add_argument(
        "--dual-overflow-warn-limit",
        type=int,
        default=20,
        help="Maximum number of dual-overflow warnings per run.",
    )

    parser.add_argument("--true-tau", type=float, default=0.12)
    parser.add_argument("--init-tau", type=float, default=0.22)
    parser.add_argument("--true-layers", type=str, default=None)
    parser.add_argument("--init-layers", type=str, default=None)

    parser.add_argument("--max-iters", type=int, default=20)
    parser.add_argument("--lr-layers", type=float, default=6.0e-4)
    parser.add_argument("--lr-tau", type=float, default=2.0e-4)
    parser.add_argument("--min-layer", type=float, default=0.2)
    parser.add_argument("--max-layer", type=float, default=2.0)
    parser.add_argument("--min-tau", type=float, default=1.0e-3)
    parser.add_argument("--max-tau", type=float, default=1.0)

    parser.add_argument(
        "--fd-check",
        action="store_true",
        default=False,
        help="Enable one-time finite-difference gradient check at iteration 0 (default: off).",
    )
    parser.add_argument("--fd-step", type=float, default=1.0e-4)
    parser.add_argument(
        "--lbfgsb-restarts",
        type=int,
        default=2,
        help="Number of additional L-BFGS-B restarts when convergence occurs at active bounds.",
    )
    parser.add_argument(
        "--lbfgsb-boundary-perturb",
        type=float,
        default=1.0e-3,
        help="Relative inward perturbation size (fraction of bound span) used for restart points.",
    )
    parser.add_argument("--lbfgsb-gtol", type=float, default=1.0e-8)
    parser.add_argument("--lbfgsb-maxls", type=int, default=40)
    parser.add_argument(
        "--adjoint-diagnostics",
        action="store_true",
        default=False,
        help="Run extra adjoint-state diagnostics (endpoint snapshots, per-step traces, and plots).",
    )
    parser.add_argument(
        "--adjoint-diagnostics-dir",
        type=str,
        default="adjoint_diagnostics",
        help="Output directory for adjoint diagnostics files.",
    )
    parser.add_argument(
        "--adjoint-diagnostics-fd-check",
        action="store_true",
        default=False,
        help="Run finite-difference checks for initial and recovered states during diagnostics.",
    )
    parser.add_argument("--plot-file", type=str, default="dissipative_wave_inverse_summary.png")
    parser.add_argument("--setup-plot-file", type=str, default="dissipative_wave_inverse_setup.svg")
    parser.add_argument("--show", action="store_true")
    args = parser.parse_args()

    true_layers = parse_layer_values(args.true_layers, args.n_layers, default_value=1.0)
    if args.true_layers is None:
        # Default heterogeneous profile from bottom to top.
        true_layers = np.array([0.85, 1.35, 0.70, 1.20], dtype=np.double)
        if args.n_layers != 4:
            true_layers = np.linspace(0.8, 1.2, args.n_layers, dtype=np.double)

    init_layers = parse_layer_values(args.init_layers, args.n_layers, default_value=1.0)

    problem = make_structured_quad_problem(
        args.nx,
        args.ny,
        args.width,
        args.height,
        args.impact_velocity,
        args.impact_window_width,
        args.n_sensors,
        args.n_layers,
        args.sensor_distribution_width,
    )
    setup_plot_path = plot_problem_setup(problem, true_layers, args.true_tau, args, args.setup_plot_file)
    print(f"Saved: {setup_plot_path}")

    observed = run_forward_or_adjoint(problem, true_layers, args.true_tau, args, observed_history=None, compute_gradients=False)
    observed_history = observed["sensor_history"]

    layers = init_layers.copy()
    tau = float(args.init_tau)
    initial_run = run_forward_or_adjoint(problem, layers, tau, args, observed_history=observed_history, compute_gradients=False)

    print("Optimizer: L-BFGS-B")
    opt_result, loss_history, tau_history, layer_history = optimize_with_lbfgsb(
        problem, observed_history, layers, tau, args
    )
    layers, tau = unpack_params(opt_result.x, args.n_layers)

    final_run = run_forward_or_adjoint(problem, layers, tau, args, observed_history=observed_history, compute_gradients=False)

    print("\nOptimization summary:")
    print(f"  status        = {opt_result.status}")
    print(f"  message       = {opt_result.message}")
    print(f"  true tau      = {args.true_tau:.6f}")
    print(f"  recovered tau = {tau:.6f}")
    print(f"  true layers   = {true_layers}")
    print(f"  recovered     = {layers}")

    if args.adjoint_diagnostics:
        diag_dir = Path(args.adjoint_diagnostics_dir)
        diag_dir.mkdir(parents=True, exist_ok=True)
        diagnostic_cases = [
            ("initial", init_layers.copy(), float(args.init_tau)),
            ("recovered", layers.copy(), float(tau)),
            ("truth", true_layers.copy(), float(args.true_tau)),
        ]

        print("\nAdjoint diagnostics:")
        for label, diag_layers, diag_tau in diagnostic_cases:
            diag_run = run_forward_or_adjoint(
                problem,
                diag_layers,
                diag_tau,
                args,
                observed_history=observed_history,
                compute_gradients=True,
                collect_adjoint_trace=True,
            )
            csv_path, fig_path, endpoint_path = _save_diagnostics_bundle(diag_run, problem, args, label, diag_dir)
            print(f"  [{label}] loss={diag_run['loss']:.6e}  grad_tau={diag_run['grad_tau']:.6e}")
            print(f"    Saved: {csv_path}")
            print(f"    Saved: {fig_path}")
            print(f"    Saved: {endpoint_path}")

            if args.adjoint_diagnostics_fd_check and label in ("initial", "recovered"):
                fd = _fd_report(problem, diag_layers, diag_tau, args, observed_history)
                fd_payload = {
                    "label": label,
                    "adjoint_layers": fd["adjoint_layers"].tolist(),
                    "adjoint_tau": float(fd["adjoint_tau"]),
                    "fd_layers": fd["fd_layers"].tolist(),
                    "fd_tau": float(fd["fd_tau"]),
                    "rel_error_layers": fd["rel_error_layers"].tolist(),
                    "rel_error_tau": float(fd["rel_error_tau"]),
                }
                fd_path = diag_dir / f"{label}_fd_check.json"
                with fd_path.open("w", encoding="utf-8") as f:
                    json.dump(fd_payload, f, indent=2)
                print(
                    f"    FD rel_err: layers={np.max(fd['rel_error_layers']):.3e}, "
                    f"tau={fd['rel_error_tau']:.3e}"
                )
                print(f"    Saved: {fd_path}")

    fig, axes = plt.subplots(2, 2, figsize=(12.0, 8.0))
    ax_loss, ax_layers, ax_trace, ax_profile = axes.ravel()

    ax_loss.plot(np.arange(len(loss_history)), loss_history, color="#2b738e", linewidth=2.0)
    ax_loss.set_title("Loss History")
    ax_loss.set_xlabel("iteration")
    ax_loss.set_ylabel("J")
    ax_loss.grid(True, alpha=0.25)

    idx = np.arange(args.n_layers)
    w = 0.28
    ax_layers.bar(idx - w, true_layers, width=w, label="true", color="#7a7a7a")
    ax_layers.bar(idx, init_layers, width=w, label="initial", color="#c9c9c9")
    ax_layers.bar(idx + w, layers, width=w, label="recovered", color="#f9826b")
    ax_layers.set_title("Layer Coefficients")
    ax_layers.set_xlabel("layer id (bottom to top)")
    ax_layers.set_ylabel("stiffness_scaling_factor")
    ax_layers.legend()
    ax_layers.grid(True, axis="y", alpha=0.25)

    sensor_id = problem["sensor_nodes"].size // 2
    t = args.dt * np.arange(args.nsteps)
    ax_trace.plot(t, observed_history[:, sensor_id], label="observed", color="#111111", linewidth=1.8)
    ax_trace.plot(t, initial_run["sensor_history"][:, sensor_id], label="initial", color="#c9c9c9", linewidth=1.5)
    ax_trace.plot(t, final_run["sensor_history"][:, sensor_id], label="recovered", color="#2b738e", linewidth=1.8)
    ax_trace.set_title("Sensor Velocity History (Vy)")
    ax_trace.set_xlabel("time")
    ax_trace.set_ylabel("velocity_Y")
    ax_trace.legend()
    ax_trace.grid(True, alpha=0.25)

    step_profile(ax_profile, true_layers, problem["layer_bounds"], "true", "#111111")
    step_profile(ax_profile, init_layers, problem["layer_bounds"], "initial", "#c9c9c9")
    step_profile(ax_profile, layers, problem["layer_bounds"], "recovered", "#f9826b")
    ax_profile.invert_yaxis()
    ax_profile.set_title("Layer Profile vs Depth")
    ax_profile.set_xlabel("stiffness_scaling_factor")
    ax_profile.set_ylabel("y")
    ax_profile.legend()
    ax_profile.grid(True, alpha=0.25)

    fig.tight_layout()
    fig.savefig(args.plot_file, dpi=170)
    print(f"Saved: {args.plot_file}")
    if args.show:
        plt.show()
    else:
        plt.close(fig)


if __name__ == "__main__":
    main()
