from __future__ import annotations

import argparse
import contextlib
import os
import pickle
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.colors import Normalize

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
sys.path.append(str(REPO_ROOT / "install" / "package"))

import REMAT


@dataclass
class ProblemData:
    coordinates: np.ndarray
    velocities: np.ndarray
    fixity: np.ndarray
    connectivity: np.ndarray
    contacts: list[tuple[np.ndarray, np.ndarray]]
    truss_connectivity: np.ndarray


@contextlib.contextmanager
def suppress_solver_output(enabled: bool):
    if not enabled:
        yield
        return
    try:
        stdout_fd = sys.stdout.fileno()
        stderr_fd = sys.stderr.fileno()
    except (AttributeError, OSError):
        yield
        return

    saved_stdout = os.dup(stdout_fd)
    saved_stderr = os.dup(stderr_fd)
    try:
        with open(os.devnull, "w") as sink:
            os.dup2(sink.fileno(), stdout_fd)
            os.dup2(sink.fileno(), stderr_fd)
            yield
    finally:
        os.dup2(saved_stdout, stdout_fd)
        os.dup2(saved_stderr, stderr_fd)
        os.close(saved_stdout)
        os.close(saved_stderr)


def load_problem(problem_path: Path) -> ProblemData:
    with problem_path.open("rb") as handle:
        coordinates, velocities, fixity, connectivity, contacts, truss_connectivity = pickle.load(handle)

    normalized_contacts = []
    for node_ids, segment_connectivity in contacts:
        normalized_contacts.append(
            (
                np.ascontiguousarray(node_ids, dtype=np.int32),
                np.ascontiguousarray(segment_connectivity, dtype=np.int32),
            )
        )

    return ProblemData(
        coordinates=np.ascontiguousarray(coordinates, dtype=np.double),
        velocities=np.ascontiguousarray(velocities, dtype=np.double),
        fixity=np.ascontiguousarray(fixity, dtype=np.bool_),
        connectivity=np.ascontiguousarray(connectivity, dtype=np.int32),
        contacts=normalized_contacts,
        truss_connectivity=np.ascontiguousarray(truss_connectivity, dtype=np.int32),
    )


def define_parameters(args: argparse.Namespace, tau: float) -> None:
    REMAT.API.define_parameter(b"body_force_x", args.body_force_x)
    REMAT.API.define_parameter(b"body_force_y", args.body_force_y)
    REMAT.API.define_parameter(b"initial_velocity_x", 0.0)
    REMAT.API.define_parameter(b"initial_velocity_y", 0.0)
    REMAT.API.define_parameter(b"mass_damping_factor", args.mass_damping_factor)
    REMAT.API.define_parameter(b"contact_stiffness", args.contact_stiffness)

    REMAT.API.define_parameter(b"truss_density", args.truss_density)
    REMAT.API.define_parameter(b"truss_youngs_modulus", args.truss_youngs_modulus)
    REMAT.API.define_parameter(b"density", args.density)
    REMAT.API.define_parameter(b"youngs_modulus", args.youngs_modulus)
    REMAT.API.define_parameter(b"yield_stress", args.yield_stress)
    REMAT.API.define_parameter(b"relaxation_time", float(tau))
    REMAT.API.define_parameter(b"area", args.area)
    REMAT.API.define_parameter(b"poissons_ratio", args.poissons_ratio)
    REMAT.API.define_parameter(b"eps_fail", args.eps_fail)
    REMAT.API.define_parameter(b"mat_overflow_limit", args.mat_overflow_limit)


def capture_snapshot(num_nodes: int, num_truss: int) -> dict[str, np.ndarray | float]:
    coords = np.zeros((num_nodes, 2), dtype=np.double)
    REMAT.API.get_node_coords(coords, True)
    eqps = REMAT.get_field(b"truss", "equivalent_plastic_strain")
    if eqps is None:
        eqps_arr = np.zeros((num_truss,), dtype=np.double)
    else:
        eqps_arr = np.asarray(eqps, dtype=np.double).copy()
    return {
        "coords": coords,
        "eqps": eqps_arr,
        "time": float(REMAT.API.get_time()),
    }


def evaluate_forward_loss_only(*, tau: float, args: argparse.Namespace, problem: ProblemData) -> float:
    REMAT.API.set_integrator_type(args.integrator_type.encode("utf-8"))
    define_parameters(args, tau)

    with suppress_solver_output(enabled=not args.verbose_solver):
        REMAT.create_geometry(
            problem.coordinates,
            problem.velocities,
            problem.fixity,
            problem.connectivity,
            problem.contacts,
            problem.truss_connectivity,
        )
        REMAT.API.initialize()

        loss = 0.0
        for _ in range(args.nsteps):
            axial_force = REMAT.get_field(b"truss", "axial_force")
            if axial_force is None:
                raise RuntimeError("Expected truss field 'axial_force' in viscoplastic run.")
            sigma = np.asarray(axial_force, dtype=np.double) / args.area
            loss += 0.5 * float(np.sum((sigma * sigma) / args.truss_youngs_modulus))
            REMAT.API.update_state(+args.dt, args.nsub_steps)

    return float(loss)


def run_iteration(
    *,
    iteration_id: int,
    tau: float,
    args: argparse.Namespace,
    problem: ProblemData,
) -> dict[str, object]:
    REMAT.API.set_integrator_type(args.integrator_type.encode("utf-8"))
    define_parameters(args, tau)

    with suppress_solver_output(enabled=not args.verbose_solver):
        REMAT.create_geometry(
            problem.coordinates,
            problem.velocities,
            problem.fixity,
            problem.connectivity,
            problem.contacts,
            problem.truss_connectivity,
        )
        REMAT.API.initialize()

        num_nodes = REMAT.API.get_num_entities(b"node")
        num_truss = REMAT.API.get_num_entities(b"truss")
        snapshots: list[dict[str, np.ndarray | float]] = [capture_snapshot(num_nodes, num_truss)]

        loss = 0.0
        peak_eqps = 0.0
        mean_force_history = []
        for step in range(args.nsteps):
            axial_force = REMAT.get_field(b"truss", "axial_force")
            if axial_force is None:
                raise RuntimeError("Expected truss field 'axial_force' in viscoplastic run.")

            sigma = np.asarray(axial_force, dtype=np.double) / args.area
            loss += 0.5 * float(np.sum((sigma * sigma) / args.truss_youngs_modulus))
            mean_force_history.append(float(np.mean(axial_force)))

            eqps = REMAT.get_field(b"truss", "equivalent_plastic_strain")
            if eqps is not None and np.size(eqps) > 0:
                peak_eqps = max(peak_eqps, float(np.max(eqps)))

            REMAT.API.update_state(+args.dt, args.nsub_steps)

            if ((step + 1) % args.frame_stride == 0) or ((step + 1) == args.nsteps):
                snapshots.append(capture_snapshot(num_nodes, num_truss))

        for _ in range(args.nsteps):
            REMAT.API.update_state(-args.dt, args.nsub_steps)

    grad_field = REMAT.get_field(b"global", "dL_dparam_relaxation_time")
    if grad_field is None:
        raise RuntimeError("Expected global field 'dL_dparam_relaxation_time' in adjoint mode.")
    grad_tau = float(grad_field[0])

    return {
        "iter": iteration_id,
        "tau_before": float(tau),
        "loss": float(loss),
        "grad_tau": grad_tau,
        "peak_eqps": float(peak_eqps),
        "mean_force_history": np.asarray(mean_force_history, dtype=np.double),
        "frames": snapshots,
    }


def collect_iteration_data(args: argparse.Namespace, problem: ProblemData) -> list[dict[str, object]]:
    tau = args.tau0
    records: list[dict[str, object]] = []
    for it in range(args.max_iters):
        record = run_iteration(iteration_id=it + 1, tau=tau, args=args, problem=problem)
        tau_next = tau
        record["line_search_trials"] = 0
        record["line_search_accepted"] = False
        if args.optimize_tau:
            step_scale = args.lr_tau
            for trial in range(args.line_search_trials):
                tau_candidate = float(np.clip(tau - step_scale * record["grad_tau"], args.min_tau, args.max_tau))
                if abs(tau_candidate - tau) <= 1.0e-16:
                    break
                candidate_loss = evaluate_forward_loss_only(tau=tau_candidate, args=args, problem=problem)
                record["line_search_trials"] = trial + 1
                if candidate_loss <= record["loss"]:
                    tau_next = tau_candidate
                    record["line_search_accepted"] = True
                    break
                step_scale *= args.line_search_shrink
        record["tau_after"] = tau_next
        record["delta_tau"] = tau_next - tau
        records.append(record)

        print(
            f"iter={record['iter']:02d}  "
            f"J={record['loss']:.6e}  "
            f"tau={record['tau_before']:.6e}  "
            f"dJ/dtau={record['grad_tau']:.6e}  "
            f"tau_next={record['tau_after']:.6e}  "
            f"ls_trials={record['line_search_trials']}"
        )
        tau = tau_next

    return records


def compute_coordinate_bounds(records: list[dict[str, object]]) -> tuple[float, float, float, float]:
    all_coords = []
    for record in records:
        for frame in record["frames"]:
            all_coords.append(frame["coords"])
    stacked = np.concatenate(all_coords, axis=0)
    xmin = float(np.min(stacked[:, 0]))
    xmax = float(np.max(stacked[:, 0]))
    ymin = float(np.min(stacked[:, 1]))
    ymax = float(np.max(stacked[:, 1]))
    dx = xmax - xmin
    dy = ymax - ymin
    pad_x = max(0.35, 0.08 * dx)
    pad_y = max(0.35, 0.08 * dy)
    return xmin - pad_x, xmax + pad_x, ymin - pad_y, ymax + pad_y


def make_animation(records: list[dict[str, object]], problem: ProblemData, args: argparse.Namespace) -> None:
    if len(records) == 0:
        raise ValueError("No optimization records were produced.")

    timeline: list[tuple[str, int, int]] = []
    for i, record in enumerate(records):
        num_frames = len(record["frames"])
        for local_frame in range(num_frames):
            timeline.append(("simulate", i, local_frame))
        for _ in range(args.hold_frames):
            timeline.append(("update", i, num_frames - 1))

    xlim = compute_coordinate_bounds(records)
    iter_ids = np.arange(1, len(records) + 1)
    tau_hist = np.asarray([record["tau_before"] for record in records], dtype=np.double)
    loss_hist = np.asarray([record["loss"] for record in records], dtype=np.double)
    eqps_all = np.concatenate([frame["eqps"] for rec in records for frame in rec["frames"]])
    if args.truss_field_max > 0.0:
        eqps_max = args.truss_field_max
    else:
        eqps_max = max(1.0e-6, float(np.percentile(eqps_all, 99.0)))

    fig = plt.figure(figsize=(13.0, 7.2))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.9, 1.0], hspace=0.30, wspace=0.25)
    ax_scene = fig.add_subplot(grid[:, 0])
    ax_tau = fig.add_subplot(grid[0, 1])
    ax_loss = fig.add_subplot(grid[1, 1])

    first_coords = records[0]["frames"][0]["coords"]
    solid_polys = [first_coords[idx, :] for idx in problem.connectivity]
    truss_segments = [first_coords[idx, :] for idx in problem.truss_connectivity]
    truss_eqps = records[0]["frames"][0]["eqps"]

    solid_collection = PolyCollection(
        solid_polys,
        facecolor="#a7d8ff",
        edgecolor="#3f6f96",
        linewidth=0.6,
        alpha=0.95,
    )
    ax_scene.add_collection(solid_collection)

    truss_collection = LineCollection(
        truss_segments,
        cmap="magma",
        norm=Normalize(vmin=0.0, vmax=eqps_max),
        linewidths=2.3,
    )
    truss_collection.set_array(truss_eqps)
    ax_scene.add_collection(truss_collection)

    fixed_nodes = np.where(np.logical_and(problem.fixity[:, 0], problem.fixity[:, 1]))[0]
    if fixed_nodes.size > 0:
        ax_scene.scatter(
            problem.coordinates[fixed_nodes, 0],
            problem.coordinates[fixed_nodes, 1],
            color="#1f1f1f",
            marker="s",
            s=30.0,
            label="fixed support",
            zorder=4,
        )

    ax_scene.axhline(0.0, color="#888888", linewidth=1.0, linestyle="--", alpha=0.6)
    ax_scene.set_aspect("equal", adjustable="box")
    ax_scene.set_xlim(xlim[0], xlim[1])
    ax_scene.set_ylim(xlim[2], xlim[3])
    ax_scene.set_xlabel("x")
    ax_scene.set_ylabel("y")
    ax_scene.set_title("Impact Simulation (Ball + Truss)")
    ax_scene.legend(loc="upper right", fontsize="small")

    cbar = fig.colorbar(truss_collection, ax=ax_scene, fraction=0.040, pad=0.018)
    cbar.set_label("equivalent plastic strain")

    tau_min = float(np.min(tau_hist))
    tau_max = float(np.max(tau_hist))
    tau_pad = max(1.0e-5, 0.10 * (tau_max - tau_min if tau_max > tau_min else tau_max + 1.0e-5))
    ax_tau.plot(iter_ids, tau_hist, "--", color="#c9c9c9", linewidth=1.2, label="all iterations")
    (tau_line,) = ax_tau.plot([], [], color="#2b6cb0", linewidth=2.1, label="played")
    (tau_marker,) = ax_tau.plot([], [], "o", color="#1e4f84", markersize=7)
    ax_tau.set_xlim(0.8, len(records) + 0.2)
    ax_tau.set_ylim(tau_min - tau_pad, tau_max + tau_pad)
    ax_tau.set_xlabel("optimization iteration")
    ax_tau.set_ylabel("tau (relaxation_time)")
    ax_tau.set_title("Design Parameter")
    ax_tau.legend(loc="best", fontsize="small")

    loss_min = float(np.min(loss_hist))
    loss_max = float(np.max(loss_hist))
    loss_pad = max(1.0e-6, 0.12 * (loss_max - loss_min if loss_max > loss_min else loss_max + 1.0e-6))
    ax_loss.plot(iter_ids, loss_hist, "--", color="#d0d0d0", linewidth=1.2, label="all iterations")
    (loss_line,) = ax_loss.plot([], [], color="#cb4b16", linewidth=2.1, label="played")
    (loss_marker,) = ax_loss.plot([], [], "o", color="#a43a10", markersize=7)
    ax_loss.set_xlim(0.8, len(records) + 0.2)
    ax_loss.set_ylim(loss_min - loss_pad, loss_max + loss_pad)
    ax_loss.set_xlabel("optimization iteration")
    ax_loss.set_ylabel("objective J")
    ax_loss.set_title("Stress Objective History")
    ax_loss.legend(loc="best", fontsize="small")

    info_text = ax_loss.text(
        0.02,
        0.04,
        "",
        transform=ax_loss.transAxes,
        ha="left",
        va="bottom",
        fontsize="small",
        family="monospace",
    )
    title_text = fig.suptitle("", fontsize="large")

    def update(frame_id: int):
        phase, iter_idx, local_idx = timeline[frame_id]
        record = records[iter_idx]
        frame = record["frames"][local_idx]
        coords = frame["coords"]

        solid_collection.set_verts([coords[idx, :] for idx in problem.connectivity])
        truss_collection.set_segments([coords[idx, :] for idx in problem.truss_connectivity])
        truss_collection.set_array(frame["eqps"])

        count = iter_idx + 1
        tau_line.set_data(iter_ids[:count], tau_hist[:count])
        tau_marker.set_data([count], [tau_hist[iter_idx]])

        loss_line.set_data(iter_ids[:count], loss_hist[:count])
        loss_marker.set_data([count], [loss_hist[iter_idx]])

        phase_label = "impact simulation" if phase == "simulate" else "design update"
        info_text.set_text(
            "\n".join(
                [
                    f"phase      : {phase_label}",
                    f"iter       : {record['iter']}/{len(records)}",
                    f"time       : {float(frame['time']):.3f} s",
                    f"tau        : {record['tau_before']:.5e}",
                    f"dJ/dtau    : {record['grad_tau']:.5e}",
                    f"tau_next   : {record['tau_after']:.5e}",
                    f"ls accepted: {record['line_search_accepted']}",
                    f"peak eqps  : {record['peak_eqps']:.4f}",
                    f"J          : {record['loss']:.5e}",
                ]
            )
        )
        title_text.set_text(
            "Adjoint Viscoplastic Truss Optimization "
            f"(iter {record['iter']}, phase: {phase_label})"
        )

        return [solid_collection, truss_collection, tau_line, tau_marker, loss_line, loss_marker, info_text, title_text]

    animation = FuncAnimation(
        fig,
        update,
        frames=len(timeline),
        interval=int(1000.0 / max(1, args.fps)),
        blit=False,
        repeat=True,
        repeat_delay=1200,
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    suffix = output_path.suffix.lower()
    if suffix == ".gif":
        animation.save(str(output_path), writer="pillow", fps=args.fps, dpi=args.dpi)
    else:
        animation.save(str(output_path), fps=args.fps, dpi=args.dpi)

    print(f"Saved optimization animation to: {output_path}")
    if args.show:
        plt.show()
    else:
        plt.close(fig)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Adjoint-driven viscoplastic truss optimization animation. "
            "Each optimization iteration runs a forward ball-impact simulation, "
            "then a reverse adjoint sweep, then updates tau."
        )
    )
    parser.add_argument(
        "--problem",
        type=str,
        default=str(REPO_ROOT / "examples" / "viscoplastic_truss" / "problem.pkl"),
        help="Path to the pre-generated viscoplastic truss pickle problem.",
    )
    parser.add_argument("--output", type=str, default=str(THIS_DIR / "adjoint_viscoplastic_optimize_animation.gif"))
    parser.add_argument("--integrator-type", type=str, default="float_truss_viscoplastic_adjoint")

    parser.add_argument("--dt", type=float, default=5.0e-3)
    parser.add_argument("--nsteps", type=int, default=70)
    parser.add_argument("--nsub-steps", type=int, default=1)
    parser.add_argument("--frame-stride", type=int, default=2)

    parser.add_argument("--max-iters", type=int, default=10)
    parser.add_argument("--tau0", type=float, default=4.0e-2)
    parser.add_argument("--lr-tau", type=float, default=5.0e-6)
    parser.add_argument("--min-tau", type=float, default=1.0e-3)
    parser.add_argument("--max-tau", type=float, default=2.0e-1)
    parser.add_argument("--disable-optimize-tau", action="store_true")
    parser.add_argument("--line-search-trials", type=int, default=5)
    parser.add_argument("--line-search-shrink", type=float, default=0.5)

    parser.add_argument("--body-force-x", type=float, default=0.0)
    parser.add_argument("--body-force-y", type=float, default=0.0)
    parser.add_argument("--mass-damping-factor", type=float, default=0.0)
    parser.add_argument("--contact-stiffness", type=float, default=2.0e3)

    parser.add_argument("--truss-density", type=float, default=1.0)
    parser.add_argument("--truss-youngs-modulus", type=float, default=2000.0)
    parser.add_argument("--density", type=float, default=3.0)
    parser.add_argument("--youngs-modulus", type=float, default=2000.0)
    parser.add_argument("--yield-stress", type=float, default=30.0)
    parser.add_argument("--area", type=float, default=1.0)
    parser.add_argument("--poissons-ratio", type=float, default=0.28)
    parser.add_argument("--eps-fail", type=float, default=0.3)
    parser.add_argument("--mat-overflow-limit", type=float, default=50.0)

    parser.add_argument("--truss-field-max", type=float, default=0.0)
    parser.add_argument("--hold-frames", type=int, default=5)
    parser.add_argument("--fps", type=int, default=12)
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--show", action="store_true")
    parser.add_argument("--verbose-solver", action="store_true")

    args = parser.parse_args()
    args.optimize_tau = not args.disable_optimize_tau
    if args.frame_stride <= 0:
        raise ValueError("--frame-stride must be > 0")
    if args.line_search_trials <= 0:
        raise ValueError("--line-search-trials must be > 0")
    if not (0.0 < args.line_search_shrink < 1.0):
        raise ValueError("--line-search-shrink must be in (0, 1)")
    return args


def main() -> None:
    args = parse_args()
    problem = load_problem(Path(args.problem))
    print("Collecting viscoplastic optimization iterations...")
    records = collect_iteration_data(args, problem)
    print("Building optimization animation...")
    make_animation(records, problem, args)


if __name__ == "__main__":
    main()
