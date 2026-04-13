from math import *
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
sys.path.append(str(REPO_ROOT / "install" / "package"))
sys.path.append(str(THIS_DIR))

import REMAT
import adjoint_visco_optimize as avo

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)


def run_iteration(
    *,
    num_elements,
    dt,
    nsteps,
    nsub_steps,
    tau,
    youngs_modulus,
    impact_velocity,
    lumped_point_mass,
    mat_overflow_limit,
    integrator_type,
):
    REMAT.API.set_integrator_type(integrator_type)
    avo.define_parameters(tau, youngs_modulus, mat_overflow_limit)
    avo.create_chain_problem(num_elements, impact_velocity, lumped_point_mass)
    REMAT.API.initialize()

    loss = 0.0
    time_history = []
    mean_stress_history = []
    right_disp_history = []
    right_vel_history = []

    for _ in range(nsteps):
        stress = REMAT.get_field(b"truss", "axial_stress")
        node_disp_x = REMAT.get_field(b"node", "displacement_X")
        node_vel_x = REMAT.get_field(b"node", "velocity_X")

        loss += 0.5 * np.sum((stress * stress) / youngs_modulus)

        time_history.append(REMAT.API.get_time())
        mean_stress_history.append(float(np.mean(stress)))
        right_disp_history.append(float(node_disp_x[-1]))
        right_vel_history.append(float(node_vel_x[-1]))

        REMAT.API.update_state(+dt, nsub_steps)

    for _ in range(nsteps):
        REMAT.API.update_state(-dt, nsub_steps)

    grad_tau = float(REMAT.get_field(b"global", "dL_dparam_relaxation_time")[0])

    return {
        "loss": float(loss),
        "grad_tau": grad_tau,
        "time_history": np.asarray(time_history),
        "mean_stress_history": np.asarray(mean_stress_history),
        "right_disp_history": np.asarray(right_disp_history),
        "right_vel_history": np.asarray(right_vel_history),
    }


def collect_iteration_data(args):
    tau = args.tau0
    youngs_modulus = args.youngs_modulus

    records = []

    for it in range(args.max_iters):
        run = run_iteration(
            num_elements=args.num_elements,
            dt=args.dt,
            nsteps=args.nsteps,
            nsub_steps=args.nsub_steps,
            tau=tau,
            youngs_modulus=youngs_modulus,
            impact_velocity=args.impact_velocity,
            lumped_point_mass=args.point_mass,
            mat_overflow_limit=args.mat_overflow_limit,
            integrator_type=args.integrator_type.encode("utf-8"),
        )

        record = {
            "iter": it + 1,
            "tau": tau,
            "loss": run["loss"],
            "grad_tau": run["grad_tau"],
            "time_history": run["time_history"],
            "mean_stress_history": run["mean_stress_history"],
            "right_disp_history": run["right_disp_history"],
            "right_vel_history": run["right_vel_history"],
        }
        records.append(record)

        print(
            f"iter={record['iter']:03d}  "
            f"loss={record['loss']: .8e}  "
            f"tau={record['tau']: .8e}  grad_tau={record['grad_tau']: .8e}"
        )

        if args.optimize_tau:
            tau = max(args.min_tau, tau - args.lr_tau * record["grad_tau"])

    return records


def pad_bounds(ymin, ymax, frac=0.08):
    if not np.isfinite(ymin) or not np.isfinite(ymax):
        return -1.0, 1.0
    if ymax <= ymin:
        c = 0.5 * (ymin + ymax)
        return c - 1.0, c + 1.0
    pad = frac * (ymax - ymin)
    return ymin - pad, ymax + pad


def make_animation(records, args):
    num_iters = len(records)
    if num_iters == 0:
        raise ValueError("No records were generated for animation.")

    all_stress = np.concatenate([r["mean_stress_history"] for r in records])
    all_disp = np.concatenate([r["right_disp_history"] for r in records])
    all_loss = np.asarray([r["loss"] for r in records], dtype=float)
    all_time = np.concatenate([r["time_history"] for r in records])

    stress_ylim = pad_bounds(float(np.min(all_stress)), float(np.max(all_stress)))
    disp_ylim = pad_bounds(float(np.min(all_disp)), float(np.max(all_disp)))
    loss_ylim = pad_bounds(float(np.min(all_loss)), float(np.max(all_loss)), frac=0.12)
    time_xlim = (float(np.min(all_time)), float(np.max(all_time)))

    fig = plt.figure(figsize=(11.0, 8.0))
    grid = fig.add_gridspec(2, 2, hspace=0.38, wspace=0.28)

    ax_setup = fig.add_subplot(grid[0, 0])
    ax_disp = fig.add_subplot(grid[0, 1])
    ax_stress = fig.add_subplot(grid[1, 0])
    ax_loss = fig.add_subplot(grid[1, 1])

    # ------------------------------------------------------------------
    # Setup panel: show geometry + impact loading explicitly.
    # ------------------------------------------------------------------
    node_x = np.arange(args.num_elements + 1, dtype=float)
    node_y = np.zeros_like(node_x)
    ax_setup.plot(node_x, node_y, "o-", color="#2b738e", linewidth=2.0, markersize=5.0)
    ax_setup.scatter([node_x[0]], [0.0], color="#111111", s=70, marker="s", label="fixed support")

    if args.point_mass > 0.0:
        ax_setup.scatter(
            [node_x[-1]],
            [0.0],
            s=170,
            facecolors="none",
            edgecolors="#f9826b",
            linewidths=2.0,
            label="lumped mass",
        )

    ax_setup.annotate(
        "",
        xy=(node_x[-1] + 0.85, 0.0),
        xytext=(node_x[-1], 0.0),
        arrowprops=dict(arrowstyle="->", linewidth=2.6, color="#d94841"),
    )
    ax_setup.text(
        node_x[-1] + 0.88,
        0.05,
        rf"$v_{{\mathrm{{imp}}}}={args.impact_velocity:.2f}$",
        color="#d94841",
        fontsize="medium",
        ha="left",
    )

    ax_setup.set_title("Impact Loading Setup", fontsize="medium")
    ax_setup.set_xlabel("node index / x-coordinate", fontsize="large")
    ax_setup.set_yticks([])
    ax_setup.set_xlim(-0.4, node_x[-1] + 1.5)
    ax_setup.set_ylim(-0.35, 0.35)
    ax_setup.legend(loc="upper left", fontsize="small")

    # ------------------------------------------------------------------
    # Right-node displacement response panel.
    # ------------------------------------------------------------------
    (disp_line,) = ax_disp.plot([], [], color="#6f6f6f", linewidth=1.8)
    ax_disp.set_title("Right-Node Displacement Response", fontsize="medium")
    ax_disp.set_xlabel("time (s)", fontsize="large")
    ax_disp.set_ylabel(r"$u_{N,x}$", fontsize="large")
    ax_disp.set_xlim(time_xlim)
    ax_disp.set_ylim(disp_ylim)

    # ------------------------------------------------------------------
    # Mean stress response panel.
    # ------------------------------------------------------------------
    (stress_line,) = ax_stress.plot([], [], color="#2b738e", linewidth=1.8)
    # ax_stress.set_title("Mean Axial Stress Response", fontsize="medium")
    ax_stress.set_xlabel("time (s)", fontsize="large")
    ax_stress.set_ylabel(r"$\sigma$", fontsize="large")
    ax_stress.set_xlim(time_xlim)
    ax_stress.set_ylim(stress_ylim)

    # ------------------------------------------------------------------
    # Loss history panel.
    # ------------------------------------------------------------------
    iter_id = np.arange(1, num_iters + 1)
    ax_loss.plot(iter_id, all_loss, "--", color="#c8c8c8", linewidth=1.2, label="all iterations")
    (loss_line,) = ax_loss.plot([], [], color="#f9826b", linewidth=2.0, label="played in animation")
    (loss_marker,) = ax_loss.plot([], [], "o", color="#d94841", markersize=7)

    ax_loss.set_title("Optimization Loss by Iteration", fontsize="medium")
    ax_loss.set_xlabel("optimization iteration", fontsize="large")
    ax_loss.set_ylabel("loss / energy objective", fontsize="large")
    ax_loss.set_xlim(0.8, num_iters + 0.2)
    ax_loss.set_ylim(loss_ylim)
    ax_loss.legend(loc="upper right", fontsize="small")

    info_text = ax_loss.text(
        0.03,
        0.05,
        "",
        transform=ax_loss.transAxes,
        ha="left",
        va="bottom",
        fontsize="small",
    )

    title_text = fig.suptitle("", fontsize="large")

    def update(frame_id):
        record = records[frame_id]

        disp_line.set_data(record["time_history"], record["right_disp_history"])
        stress_line.set_data(record["time_history"], record["mean_stress_history"])

        xloss = iter_id[: frame_id + 1]
        yloss = all_loss[: frame_id + 1]
        loss_line.set_data(xloss, yloss)
        loss_marker.set_data([frame_id + 1], [record["loss"]])

        loss0 = all_loss[0]
        loss_drop = loss0 - record["loss"]
        loss_drop_pct = 100.0 * loss_drop / abs(loss0) if abs(loss0) > 1.0e-16 else 0.0

        info_color = "#1c7c31" if loss_drop >= 0.0 else "#b71c1c"
        info_text.set_color(info_color)
        info_text.set_text(
            "\n".join([
                f"iter {record['iter']}/{num_iters}",
                f"J = {record['loss']:.6e}",
                f"tau = {record['tau']:.4e}",
                f"dJ/dtau = {record['grad_tau']:.3e}",
                f"decrease vs iter 1 = {loss_drop_pct:.2f}%",
            ])
)
        title_text.set_text(
            f"Adjoint Viscoelastic Optimization Animation  "
            f"(iter = {record['iter']})"
        )

        return [disp_line, stress_line, loss_line, loss_marker, info_text, title_text]

    animation = FuncAnimation(
        fig,
        update,
        frames=num_iters,
        interval=500,
        blit=False,
        repeat=True,
        repeat_delay=1000,
    )

    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fps = 2  # 0.5 sec per optimization iteration frame
    suffix = output_path.suffix.lower()
    if suffix == ".gif":
        animation.save(str(output_path), writer="pillow", fps=fps, dpi=args.dpi)
    else:
        animation.save(str(output_path), fps=fps, dpi=args.dpi)

    print(f"Saved animation to: {output_path}")

    if args.show:
        plt.show()
    else:
        plt.close(fig)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Animate optimization iterations of the viscoelastic adjoint example. "
            "Each frame corresponds to one optimization iteration (0.5 s per frame)."
        )
    )

    parser.add_argument("--num-elements", type=int, default=8)
    parser.add_argument("--dt", type=float, default=1.0e-3)
    parser.add_argument("--nsteps", type=int, default=200)
    parser.add_argument("--nsub-steps", type=int, default=1)

    parser.add_argument("--tau0", type=float, default=0.30)
    parser.add_argument("--youngs-modulus", type=float, default=1.0)
    parser.add_argument("--max-iters", type=int, default=8)
    parser.add_argument("--lr-tau", type=float, default=1.0e-3)
    parser.add_argument("--min-tau", type=float, default=1.0e-4)
    parser.add_argument("--disable-optimize-tau", action="store_true")

    parser.add_argument("--impact-velocity", type=float, default=2.0)
    parser.add_argument("--point-mass", type=float, default=0.5)
    parser.add_argument("--mat-overflow-limit", type=float, default=1.0e6)

    parser.add_argument("--integrator-type", type=str, default="float_truss_visco_adjoint")
    parser.add_argument(
        "--output",
        type=str,
        default="adjoint_visco_optimize_animation.gif",
        help="Output animation path (.gif recommended for portability).",
    )
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--show", action="store_true")

    args = parser.parse_args()

    args.optimize_tau = not args.disable_optimize_tau
    if not args.optimize_tau:
        raise ValueError("Tau optimization must be enabled.")

    return args


def main():
    args = parse_args()

    print("Collecting optimization iteration data for animation...")
    records = collect_iteration_data(args)

    print("Building animation...")
    make_animation(records, args)


if __name__ == "__main__":
    main()
