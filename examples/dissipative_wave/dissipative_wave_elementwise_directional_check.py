import argparse
import contextlib
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import dissipative_wave_inverse as inv


THIS_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = THIS_DIR / "inverse_case_outputs" / "case3_adjoint_diag"
DEFAULT_JSON = DEFAULT_OUTPUT_DIR / "elementwise_s_tau_directional_check.json"
DEFAULT_PLOT = DEFAULT_OUTPUT_DIR / "elementwise_s_tau_directional_check.svg"
DEFAULT_PDF = DEFAULT_OUTPUT_DIR / "elementwise_s_tau_directional_check.pdf"

PLOT_COLOR_PRIMARY = "#2b738eff"
PLOT_COLOR_SECONDARY = "#f9826bff"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)


def parse_h_values(text):
    values = [float(part.strip()) for part in text.split(",") if part.strip()]
    if not values:
        raise ValueError("At least one finite-difference step is required.")
    return sorted(values, reverse=True)


def build_h_values(args):
    if args.h_values is not None:
        return parse_h_values(args.h_values)
    if args.h_min <= 0.0 or args.h_max <= 0.0:
        raise ValueError("h-min and h-max must be positive.")
    if args.h_max <= args.h_min:
        raise ValueError("h-max must be larger than h-min.")
    if args.h_count < 2:
        raise ValueError("h-count must be at least 2.")
    return np.geomspace(args.h_max, args.h_min, args.h_count).tolist()


@contextlib.contextmanager
def maybe_suppress_solver_output(enabled):
    if not enabled:
        yield
        return

    stdout_fd = os.dup(1)
    stderr_fd = os.dup(2)
    with open(os.devnull, "w", encoding="utf-8") as devnull:
        try:
            os.dup2(devnull.fileno(), 1)
            os.dup2(devnull.fileno(), 2)
            yield
        finally:
            os.dup2(stdout_fd, 1)
            os.dup2(stderr_fd, 2)
            os.close(stdout_fd)
            os.close(stderr_fd)


def configure_inverse_globals(args):
    inv.WIDTH = float(args.width)
    inv.HEIGHT = float(args.height)
    inv.NX = int(args.nx)
    inv.NY = int(args.ny)
    inv.DT = float(args.dt)
    inv.N_STEPS = int(args.nsteps)
    inv.N_SUB_STEPS = int(args.nsub_steps)
    inv.INTEGRATOR_TYPE = str(args.integrator_type)
    inv.IMPACT_VELOCITY = float(args.impact_velocity)
    inv.IMPACT_WINDOW_WIDTH = float(args.impact_window_width)
    inv.N_SENSORS = int(args.n_sensors)
    inv.SENSOR_DISTRIBUTION_WIDTH = float(args.sensor_distribution_width)
    inv.OVERFLOW_LIMIT = float(args.overflow_limit)
    inv.MAT_OVERFLOW_LIMIT = float(args.mat_overflow_limit)
    inv.STIFFNESS_MIN = float(args.stiffness_min)
    inv.STIFFNESS_MAX = float(args.stiffness_max)
    inv.TAU_MIN = float(args.tau_min)
    inv.TAU_MAX = float(args.tau_max)


def smooth_random_field(problem, seed, max_mode_x, max_mode_y):
    rng = np.random.default_rng(seed)
    xy = np.asarray(problem["elem_centers"], dtype=np.double)
    x = xy[:, 0] / float(problem["width"])
    y = xy[:, 1] / float(problem["height"])

    field = np.zeros(x.size, dtype=np.double)
    for kx in range(max_mode_x + 1):
        for ky in range(max_mode_y + 1):
            if kx == 0 and ky == 0:
                continue
            weight = 1.0 / (1.0 + kx * kx + ky * ky)
            phase_x = rng.uniform(0.0, 2.0 * np.pi)
            phase_y = rng.uniform(0.0, 2.0 * np.pi)
            coeff = rng.normal() * weight
            field += coeff * np.cos(np.pi * kx * x + phase_x) * np.cos(np.pi * ky * y + phase_y)

    field -= np.mean(field)
    max_abs = float(np.max(np.abs(field)))
    if max_abs <= 0.0:
        raise ValueError("Smooth random field has zero amplitude.")
    return field / max_abs


def regularization_loss_and_grad(problem, stiffness, tau, objective):
    nelem = problem["connectivity"].shape[0]
    if objective == "data":
        return 0.0, np.zeros(nelem, dtype=np.double), np.zeros(nelem, dtype=np.double)

    i_idx, j_idx = inv.build_edge_pairs(problem["nx"], problem["ny"])
    reg_s_loss, reg_s_grad = inv.regularization_loss_and_grad(
        stiffness,
        i_idx,
        j_idx,
        inv.REG_L2_STIFFNESS,
        inv.REG_TV_STIFFNESS,
        inv.REG_TV_EPS,
    )
    reg_t_loss, reg_t_grad = inv.regularization_loss_and_grad(
        tau,
        i_idx,
        j_idx,
        inv.REG_L2_TAU,
        inv.REG_TV_TAU,
        inv.REG_TV_EPS,
    )
    return reg_s_loss + reg_t_loss, reg_s_grad, reg_t_grad


def run_forward_or_adjoint(problem, stiffness, tau, observed_history, compute_gradients, suppress_solver_output):
    with maybe_suppress_solver_output(suppress_solver_output):
        return inv.run_forward_or_adjoint(
            problem,
            stiffness,
            tau,
            observed_history=observed_history,
            compute_gradients=compute_gradients,
        )


def objective_value(problem, stiffness, tau, observed_history, objective, suppress_solver_output):
    run = run_forward_or_adjoint(
        problem,
        stiffness,
        tau,
        observed_history,
        False,
        suppress_solver_output,
    )
    reg_loss, _, _ = regularization_loss_and_grad(problem, stiffness, tau, objective)
    return float(run["data_loss"] + reg_loss)


def objective_and_gradient(problem, stiffness, tau, observed_history, objective, suppress_solver_output):
    run = run_forward_or_adjoint(
        problem,
        stiffness,
        tau,
        observed_history,
        True,
        suppress_solver_output,
    )
    reg_loss, reg_s_grad, reg_t_grad = regularization_loss_and_grad(problem, stiffness, tau, objective)
    grad_stiff = np.asarray(run["grad_stiff"], dtype=np.double) + reg_s_grad
    grad_tau = np.asarray(run["grad_tau"], dtype=np.double) + reg_t_grad
    return float(run["data_loss"] + reg_loss), grad_stiff, grad_tau


def bounds_ok(stiffness, tau):
    return (
        np.all(stiffness >= inv.STIFFNESS_MIN)
        and np.all(stiffness <= inv.STIFFNESS_MAX)
        and np.all(tau >= inv.TAU_MIN)
        and np.all(tau <= inv.TAU_MAX)
    )


def check_direction(
    name,
    problem,
    base_stiffness,
    base_tau,
    observed_history,
    base_loss,
    grad_stiff,
    grad_tau,
    ds,
    dtau,
    h_values,
    objective,
    suppress_solver_output,
):
    adjoint_directional = float(np.dot(grad_stiff, ds) + np.dot(grad_tau, dtau))
    rows = []

    for h in h_values:
        print(f"  {name}: h={h:.3e}", flush=True)
        plus_s = base_stiffness + h * ds
        plus_t = base_tau + h * dtau
        minus_s = base_stiffness - h * ds
        minus_t = base_tau - h * dtau

        if not bounds_ok(plus_s, plus_t) or not bounds_ok(minus_s, minus_t):
            rows.append(
                {
                    "h": float(h),
                    "status": "skipped_out_of_bounds",
                }
            )
            continue

        loss_plus = objective_value(problem, plus_s, plus_t, observed_history, objective, suppress_solver_output)
        loss_minus = objective_value(problem, minus_s, minus_t, observed_history, objective, suppress_solver_output)
        fd_directional = (loss_plus - loss_minus) / (2.0 * h)
        rel_error = abs(fd_directional - adjoint_directional) / max(
            abs(fd_directional), abs(adjoint_directional), 1.0e-14
        )
        first_remainder = abs(loss_plus - base_loss)
        second_remainder = abs(loss_plus - base_loss - h * adjoint_directional)

        rows.append(
            {
                "h": float(h),
                "status": "ok",
                "loss_plus": float(loss_plus),
                "loss_minus": float(loss_minus),
                "adjoint_directional": adjoint_directional,
                "fd_directional": float(fd_directional),
                "relative_error": float(rel_error),
                "first_remainder": float(first_remainder),
                "second_remainder": float(second_remainder),
            }
        )

    ok_rows = [row for row in rows if row["status"] == "ok"]
    for prev, curr in zip(ok_rows, ok_rows[1:]):
        if prev["second_remainder"] > 0.0 and curr["second_remainder"] > 0.0:
            curr["second_remainder_order_from_previous"] = float(
                np.log(prev["second_remainder"] / curr["second_remainder"]) / np.log(prev["h"] / curr["h"])
            )

    active = ds if name == "s_only" else dtau
    return {
        "name": name,
        "adjoint_directional": adjoint_directional,
        "direction_stats": {
            "max_abs": float(np.max(np.abs(active))),
            "mean": float(np.mean(active)),
            "l2": float(np.linalg.norm(active)),
            "physical_units": "stiffness_scaling" if name == "s_only" else "relaxation_time",
        },
        "rows": rows,
    }


def save_plot(results, output_plot, output_pdf):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

    style = {
        "s_only": {"label": "s only", "color": PLOT_COLOR_PRIMARY, "linestyle": "-"},
        "tau_only": {"label": "tau only", "color": PLOT_COLOR_SECONDARY, "linestyle": "--"},
    }

    for result in results:
        rows = [row for row in result["rows"] if row["status"] == "ok"]
        if not rows:
            continue
        h = np.array([row["h"] for row in rows], dtype=np.double)
        rel = np.array([row["relative_error"] for row in rows], dtype=np.double)
        rem2 = np.array([row["second_remainder"] for row in rows], dtype=np.double)
        kwargs = style[result["name"]]

        axes[0].loglog(
            h,
            rel,
            color=kwargs["color"],
            linestyle=kwargs["linestyle"],
            linewidth=1.6,
            label=kwargs["label"],
        )
        best = min(rows, key=lambda row: row["relative_error"])
        axes[0].axvline(best["h"], color=kwargs["color"], linewidth=1.0, alpha=0.35)
        axes[1].loglog(
            h,
            rem2,
            color=kwargs["color"],
            linestyle=kwargs["linestyle"],
            linewidth=1.6,
            label=kwargs["label"],
        )
        if rem2[0] > 0.0:
            axes[1].loglog(
                h,
                rem2[0] * (h / h[0]) ** 2,
                color="0.35",
                linewidth=1.0,
                dashes=(6, 8),
                alpha=0.7,
                label=r"$O(h^2)$ reference" if result["name"] == "s_only" else None,
            )

    axes[0].set_xlabel("dimensionless FD step h", fontsize="large")
    axes[0].set_ylabel("relative directional error", fontsize="large")
    axes[0].set_title("Directional FD", fontsize="medium")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(fontsize="medium")

    axes[1].set_xlabel("dimensionless FD step h", fontsize="large")
    axes[1].set_ylabel("second Taylor remainder", fontsize="large")
    axes[1].set_title("Taylor check", fontsize="medium")
    axes[1].grid(True, alpha=0.25)
    axes[1].legend(fontsize="medium")

    fig.tight_layout()
    fig.savefig(output_plot, dpi=200, metadata={"Title": "Elementwise s and tau directional checks"})
    fig.savefig(output_pdf, dpi=200, metadata={"Title": "Elementwise s and tau directional checks"})
    fig.clf()
    plt.close(fig)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Elementwise directional finite-difference checks for smooth s-only and tau-only fields."
    )
    parser.add_argument("--nx", type=int, default=inv.NX)
    parser.add_argument("--ny", type=int, default=inv.NY)
    parser.add_argument("--width", type=float, default=inv.WIDTH)
    parser.add_argument("--height", type=float, default=inv.HEIGHT)
    parser.add_argument("--nsteps", type=int, default=inv.N_STEPS)
    parser.add_argument("--nsub-steps", type=int, default=inv.N_SUB_STEPS)
    parser.add_argument("--dt", type=float, default=inv.DT)
    parser.add_argument("--integrator-type", type=str, default=inv.INTEGRATOR_TYPE)
    parser.add_argument("--impact-velocity", type=float, default=inv.IMPACT_VELOCITY)
    parser.add_argument("--impact-window-width", type=float, default=inv.IMPACT_WINDOW_WIDTH)
    parser.add_argument("--n-sensors", type=int, default=inv.N_SENSORS)
    parser.add_argument("--sensor-distribution-width", type=float, default=inv.SENSOR_DISTRIBUTION_WIDTH)
    parser.add_argument("--overflow-limit", type=float, default=inv.OVERFLOW_LIMIT)
    parser.add_argument("--mat-overflow-limit", type=float, default=inv.MAT_OVERFLOW_LIMIT)
    parser.add_argument("--stiffness-min", type=float, default=inv.STIFFNESS_MIN)
    parser.add_argument("--stiffness-max", type=float, default=inv.STIFFNESS_MAX)
    parser.add_argument("--tau-min", type=float, default=inv.TAU_MIN)
    parser.add_argument("--tau-max", type=float, default=inv.TAU_MAX)
    parser.add_argument("--base-stiffness", type=float, default=inv.INIT_STIFFNESS)
    parser.add_argument("--base-tau", type=float, default=inv.INIT_TAU)
    parser.add_argument("--objective", choices=("data", "total"), default="data")
    parser.add_argument("--h-values", type=str, default=None)
    parser.add_argument("--h-min", type=float, default=1.0e-6)
    parser.add_argument("--h-max", type=float, default=4.0)
    parser.add_argument("--h-count", type=int, default=33)
    parser.add_argument(
        "--stiffness-step-fraction",
        type=float,
        default=0.05,
        help="At h=1, max |delta s| is this fraction of the stiffness bound width.",
    )
    parser.add_argument(
        "--tau-target-step-fraction",
        type=float,
        default=0.05,
        help="When tau is aligned to the s-only dip, max |delta tau| at that h is this fraction of the tau range.",
    )
    parser.add_argument("--stiffness-direction-scale", type=float, default=None)
    parser.add_argument("--tau-direction-scale", type=float, default=None)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--max-mode-x", type=int, default=4)
    parser.add_argument("--max-mode-y", type=int, default=3)
    parser.add_argument("--verbose-remat", action="store_true")
    parser.add_argument("--output-json", type=str, default=str(DEFAULT_JSON))
    parser.add_argument("--output-plot", type=str, default=str(DEFAULT_PLOT))
    parser.add_argument("--output-pdf", type=str, default=str(DEFAULT_PDF))
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    configure_inverse_globals(args)

    problem = inv.make_structured_quad_problem()
    nelem = problem["connectivity"].shape[0]
    h_values = build_h_values(args)

    true_stiffness, true_tau, _ = inv.make_true_fields(problem)
    base_tau_value = float(args.base_tau)
    base_stiffness = np.full(nelem, float(args.base_stiffness), dtype=np.double)
    base_tau = np.full(nelem, base_tau_value, dtype=np.double)
    if not bounds_ok(base_stiffness, base_tau):
        raise ValueError("Base fields must be inside the configured stiffness/tau bounds.")

    print("Generating synthetic observations...", flush=True)
    observed = run_forward_or_adjoint(
        problem,
        true_stiffness,
        true_tau,
        None,
        False,
        not args.verbose_remat,
    )["sensor_history"]

    print("Computing elementwise adjoint gradient at the base point...", flush=True)
    base_loss, grad_stiff, grad_tau = objective_and_gradient(
        problem,
        base_stiffness,
        base_tau,
        observed,
        args.objective,
        not args.verbose_remat,
    )

    shared_phi = smooth_random_field(problem, args.seed, args.max_mode_x, args.max_mode_y)

    stiffness_scale = (
        float(args.stiffness_direction_scale)
        if args.stiffness_direction_scale is not None
        else float(args.stiffness_step_fraction) * (inv.STIFFNESS_MAX - inv.STIFFNESS_MIN)
    )
    zero = np.zeros(nelem, dtype=np.double)

    print(f"Using shared smooth spatial direction phi from seed={args.seed}.", flush=True)
    print(f"  stiffness scale at h=1: {stiffness_scale:.6e}", flush=True)

    print("Running s_only directional FD sweep...", flush=True)
    s_result = check_direction(
        "s_only",
        problem,
        base_stiffness,
        base_tau,
        observed,
        base_loss,
        grad_stiff,
        grad_tau,
        stiffness_scale * shared_phi,
        zero,
        h_values,
        args.objective,
        not args.verbose_remat,
    )
    s_ok_rows = [row for row in s_result["rows"] if row["status"] == "ok"]
    if not s_ok_rows:
        raise RuntimeError("The s_only check produced no in-bounds finite-difference rows.")
    s_best = min(s_ok_rows, key=lambda row: row["relative_error"])
    s_best_h = float(s_best["h"])

    if args.tau_direction_scale is None:
        tau_scale = float(args.tau_target_step_fraction) * (inv.TAU_MAX - inv.TAU_MIN) / s_best_h
        tau_scale_mode = "aligned_to_s_only_best_h"
    else:
        tau_scale = float(args.tau_direction_scale)
        tau_scale_mode = "manual"

    print(f"  s_only best h: {s_best_h:.6e}", flush=True)
    print(f"  tau scale at h=1: {tau_scale:.6e} ({tau_scale_mode})", flush=True)

    print("Running tau_only directional FD sweep...", flush=True)
    tau_result = check_direction(
        "tau_only",
        problem,
        base_stiffness,
        base_tau,
        observed,
        base_loss,
        grad_stiff,
        grad_tau,
        zero,
        tau_scale * shared_phi,
        h_values,
        args.objective,
        not args.verbose_remat,
    )
    results = [s_result, tau_result]

    for result in results:
        skipped = sum(1 for row in result["rows"] if row["status"] != "ok")
        if skipped:
            print(
                f"  {result['name']}: skipped {skipped} out-of-bounds large-step rows.",
                flush=True,
            )

    payload = {
        "config": {
            "mesh": {"nx": inv.NX, "ny": inv.NY, "width": inv.WIDTH, "height": inv.HEIGHT},
            "time": {"dt": inv.DT, "n_steps": inv.N_STEPS, "n_sub_steps": inv.N_SUB_STEPS},
            "integrator_type": inv.INTEGRATOR_TYPE,
            "objective": args.objective,
            "bounds": {
                "stiffness": [inv.STIFFNESS_MIN, inv.STIFFNESS_MAX],
                "tau": [inv.TAU_MIN, inv.TAU_MAX],
            },
            "base": {
                "stiffness": float(args.base_stiffness),
                "tau": float(base_tau_value),
            },
            "h_values": h_values,
            "smooth_direction": {
                "seed": int(args.seed),
                "max_mode_x": int(args.max_mode_x),
                "max_mode_y": int(args.max_mode_y),
                "same_spatial_phi_for_stiffness_and_tau": True,
                "raw_normalization": "max_abs_1",
                "stiffness_scale_at_h1": float(stiffness_scale),
                "tau_scale_at_h1": float(tau_scale),
                "stiffness_step_fraction": float(args.stiffness_step_fraction),
                "tau_scale_mode": tau_scale_mode,
                "tau_target_step_fraction_at_s_best_h": float(args.tau_target_step_fraction),
                "s_only_best_h_used_for_tau_alignment": s_best_h,
            },
        },
        "base_loss": float(base_loss),
        "gradient_norms": {
            "stiffness_l2": float(np.linalg.norm(grad_stiff)),
            "tau_l2": float(np.linalg.norm(grad_tau)),
            "stiffness_inf": float(np.max(np.abs(grad_stiff))),
            "tau_inf": float(np.max(np.abs(grad_tau))),
        },
        "checks": results,
    }

    output_json = Path(args.output_json)
    output_plot = Path(args.output_plot)
    output_pdf = Path(args.output_pdf)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_plot.parent.mkdir(parents=True, exist_ok=True)
    output_pdf.parent.mkdir(parents=True, exist_ok=True)

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    save_plot(results, output_plot, output_pdf)

    print("\nRun complete")
    print(f"  base loss : {base_loss:.6e}")
    for result in results:
        ok_rows = [row for row in result["rows"] if row["status"] == "ok"]
        if ok_rows:
            best = min(ok_rows, key=lambda row: row["relative_error"])
            print(
                f"  {result['name']}: best rel={best['relative_error']:.6e} "
                f"at h={best['h']:.3e}"
            )
    print(f"  wrote json: {output_json}")
    print(f"  wrote plot: {output_plot}")
    print(f"  wrote pdf : {output_pdf}")


if __name__ == "__main__":
    main()
