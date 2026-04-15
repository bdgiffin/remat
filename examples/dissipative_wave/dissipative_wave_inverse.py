import argparse
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


def make_structured_quad_problem(nx, ny, width, height, impact_velocity, source_window_fraction, nsensors, n_layers):
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
    left_nodes = np.where(np.abs(node_x - 0.0) < eps)[0]
    right_nodes = np.where(np.abs(node_x - width) < eps)[0]
    top_nodes = np.where(np.abs(node_y - height) < eps)[0]

    fixity[bottom_nodes, :] = True
    fixity[left_nodes, :] = True
    fixity[right_nodes, :] = True

    source_half_width = 0.5 * source_window_fraction * width
    xmid = 0.5 * width
    source_nodes = top_nodes[np.abs(node_x[top_nodes] - xmid) <= source_half_width]
    velocities[source_nodes, 1] = -abs(impact_velocity)

    top_interior = top_nodes[(node_x[top_nodes] > eps) & (node_x[top_nodes] < width - eps)]
    if top_interior.size == 0:
        raise ValueError("No interior top-surface nodes found for sensors.")
    if nsensors > top_interior.size:
        nsensors = top_interior.size
    sensor_pick = np.unique(np.round(np.linspace(0, top_interior.size - 1, nsensors)).astype(int))
    sensor_nodes = top_interior[sensor_pick]

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
        "source_nodes": source_nodes.astype(np.int32),
        "source_center_x": float(xmid),
        "source_half_width": float(source_half_width),
        "layer_bounds": layer_bounds,
        "elem_layer_ids": elem_layer_ids,
        "width": width,
        "height": height,
    }


def configure_run(problem, layer_coeffs, tau, args):
    global _ACTIVE_LAYER_COEFFS, _LAYER_BOUNDS
    _ACTIVE_LAYER_COEFFS = np.asarray(layer_coeffs, dtype=np.double).copy()
    _LAYER_BOUNDS = np.asarray(problem["layer_bounds"], dtype=np.double).copy()

    REMAT.API.set_integrator_type(args.integrator_type.encode("utf-8"))
    REMAT.API.define_parameter(b"body_force_x", 0.0)
    REMAT.API.define_parameter(b"body_force_y", 0.0)
    REMAT.API.define_parameter(b"mass_damping_factor", args.mass_damping_factor)
    REMAT.API.define_parameter(b"contact_stiffness", 0.0)
    REMAT.API.define_parameter(b"overflow_limit", float(args.overflow_limit))
    REMAT.API.define_parameter(b"mat_overflow_limit", float(args.mat_overflow_limit))
    REMAT.API.define_parameter(b"adjoint_material_objective_weight", 0.0)

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


def run_forward_or_adjoint(problem, layer_coeffs, tau, args, observed_history=None, compute_gradients=False):
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
    if compute_gradients:
        if observed_history is None:
            raise ValueError("compute_gradients=True requires observed_history.")

        REMAT.clear_adjoint_state()
        for rev in range(args.nsteps):
            k = args.nsteps - 1 - rev
            residual = sensor_history[k, :] - observed_history[k, :]
            seed_xy = np.zeros((nsensors, 2), dtype=np.double)
            seed_xy[:, 1] = residual
            REMAT.add_nodal_velocity_adjoint_seed(sensor_nodes, seed_xy)
            REMAT.API.update_state(args.dt, args.nsub_steps, REMAT.PASS_BACKWARD_ADJOINT)

        grad_tau = float(REMAT.get_field(b"global", "dL_dparam_relaxation_time")[0])
        elem_grad_scaling = REMAT.get_field(b"element", "dparam_stiffness_scaling_factor")
        for lid in range(layer_coeffs.size):
            grad_layers[lid] = float(np.sum(elem_grad_scaling[problem["elem_layer_ids"] == lid]))

    return {
        "loss": float(loss),
        "sensor_history": sensor_history,
        "grad_tau": float(grad_tau),
        "grad_layers": grad_layers,
    }


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
    source_nodes = problem["source_nodes"]
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

    source_xy = coords[source_nodes]
    ax.scatter(
        source_xy[:, 0],
        source_xy[:, 1],
        s=36,
        marker="s",
        c="#205a8e",
        edgecolors="white",
        linewidths=0.6,
        label=f"impact patch ({source_nodes.size})",
        zorder=7,
    )

    impact_x = float(np.mean(source_xy[:, 0])) if source_xy.size else float(problem["source_center_x"])
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
        "  measured field: velocity_Y(t)",
        "",
        "Initial condition:",
        f"  top impact patch width = {2.0 * problem['source_half_width']:.3f}",
        f"  impact center x = {problem['source_center_x']:.3f}",
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


def main():
    parser = argparse.ArgumentParser(
        description="Inverse dissipative-wave example: recover layered stiffness and relaxation_time from top-surface velocity sensors."
    )
    parser.add_argument("--nx", type=int, default=48)
    parser.add_argument("--ny", type=int, default=18)
    parser.add_argument("--width", type=float, default=10.0)
    parser.add_argument("--height", type=float, default=3.0)
    parser.add_argument("--nsteps", type=int, default=80)
    parser.add_argument("--nsub-steps", type=int, default=6)
    parser.add_argument("--dt", type=float, default=2.5e-3)
    parser.add_argument("--n-layers", type=int, default=4)
    parser.add_argument("--n-sensors", type=int, default=10)
    parser.add_argument("--source-window-fraction", type=float, default=0.12)
    parser.add_argument("--impact-velocity", type=float, default=0.35)
    parser.add_argument("--integrator-type", type=str, default="fixed_visco")

    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--youngs-modulus", type=float, default=5.0)
    parser.add_argument("--poissons-ratio", type=float, default=0.28)
    parser.add_argument("--shear-modulus-maxwell", type=float, default=2.0)
    parser.add_argument("--mass-damping-factor", type=float, default=0.0)
    parser.add_argument("--overflow-limit", type=int, default=50)
    parser.add_argument("--mat-overflow-limit", type=int, default=50)

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
        args.source_window_fraction,
        args.n_sensors,
        args.n_layers,
    )
    setup_plot_path = plot_problem_setup(problem, true_layers, args.true_tau, args, args.setup_plot_file)
    print(f"Saved: {setup_plot_path}")

    observed = run_forward_or_adjoint(problem, true_layers, args.true_tau, args, observed_history=None, compute_gradients=False)
    observed_history = observed["sensor_history"]

    layers = init_layers.copy()
    tau = float(args.init_tau)
    initial_run = run_forward_or_adjoint(problem, layers, tau, args, observed_history=observed_history, compute_gradients=False)

    loss_history = []
    tau_history = []
    layer_history = []

    for it in range(args.max_iters):
        run = run_forward_or_adjoint(problem, layers, tau, args, observed_history=observed_history, compute_gradients=True)
        loss = run["loss"]
        grad_layers = run["grad_layers"]
        grad_tau = run["grad_tau"]

        loss_history.append(loss)
        tau_history.append(tau)
        layer_history.append(layers.copy())

        print(
            f"iter={it:03d}  loss={loss: .8e}  tau={tau: .8e}  "
            f"|grad_layers|_inf={np.max(np.abs(grad_layers)): .8e}  grad_tau={grad_tau: .8e}"
        )

        if args.fd_check and it == 0:
            fd_layers, fd_tau = finite_difference_gradients(problem, layers, tau, args, observed_history, args.fd_step)
            eps = 1.0e-14
            rel_layer = np.abs(grad_layers - fd_layers) / np.maximum.reduce(
                [np.abs(grad_layers), np.abs(fd_layers), np.full_like(fd_layers, eps)]
            )
            rel_tau = abs(grad_tau - fd_tau) / max(abs(grad_tau), abs(fd_tau), eps)
            print("FD check (layer gradients):")
            for i in range(args.n_layers):
                print(
                    f"  layer {i}: adj={grad_layers[i]: .8e}  fd={fd_layers[i]: .8e}  rel_err={rel_layer[i]: .3e}"
                )
            print(f"FD check (tau): adj={grad_tau: .8e}  fd={fd_tau: .8e}  rel_err={rel_tau: .3e}")

        layers = np.clip(layers - args.lr_layers * grad_layers, args.min_layer, args.max_layer)
        tau = float(np.clip(tau - args.lr_tau * grad_tau, args.min_tau, args.max_tau))

    final_run = run_forward_or_adjoint(problem, layers, tau, args, observed_history=observed_history, compute_gradients=False)

    print("\nOptimization summary:")
    print(f"  true tau      = {args.true_tau:.6f}")
    print(f"  recovered tau = {tau:.6f}")
    print(f"  true layers   = {true_layers}")
    print(f"  recovered     = {layers}")

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
