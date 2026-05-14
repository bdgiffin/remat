import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import dissipative_wave_elementwise_directional_check as check
import dissipative_wave_inverse as inv


THIS_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT_DIR = THIS_DIR / "inverse_case_outputs" / "case3_adjoint_diag"
DEFAULT_JSON = DEFAULT_OUTPUT_DIR / "elementwise_fd_stability_check.json"
DEFAULT_PLOT = DEFAULT_OUTPUT_DIR / "elementwise_fd_stability_check.svg"

PLOT_COLOR_PRIMARY = "#2b738eff"
PLOT_COLOR_SECONDARY = "#f9826bff"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)


def build_h_values(args):
    if args.h_values is not None:
        values = [float(part.strip()) for part in args.h_values.split(",") if part.strip()]
        if not values:
            raise ValueError("At least one h value is required.")
        return sorted(values, reverse=True)
    if args.h_min <= 0.0 or args.h_max <= 0.0:
        raise ValueError("h-min and h-max must be positive.")
    if args.h_max <= args.h_min:
        raise ValueError("h-max must be larger than h-min.")
    return np.geomspace(args.h_max, args.h_min, args.h_count).tolist()


def fd_directional_sweep(problem, base_stiffness, base_tau, observed, ds, dtau, h_values, objective, name, quiet):
    rows = []
    for h in h_values:
        print(f"  {name}: h={h:.3e}", flush=True)
        plus_s = base_stiffness + h * ds
        plus_t = base_tau + h * dtau
        minus_s = base_stiffness - h * ds
        minus_t = base_tau - h * dtau

        if not check.bounds_ok(plus_s, plus_t) or not check.bounds_ok(minus_s, minus_t):
            rows.append({"h": float(h), "status": "skipped_out_of_bounds"})
            continue

        loss_plus = check.objective_value(problem, plus_s, plus_t, observed, objective, quiet)
        loss_minus = check.objective_value(problem, minus_s, minus_t, observed, objective, quiet)
        fd_directional = (loss_plus - loss_minus) / (2.0 * h)

        active = ds if name == "s_only" else dtau
        rows.append(
            {
                "h": float(h),
                "status": "ok",
                "loss_plus": float(loss_plus),
                "loss_minus": float(loss_minus),
                "fd_directional": float(fd_directional),
                "max_physical_perturbation": float(h * np.max(np.abs(active))),
            }
        )
    return rows


def fd_self_consistency(rows, fd_signal_floor):
    ok_rows = [row for row in rows if row["status"] == "ok"]
    pairs = []
    for coarse, fine in zip(ok_rows, ok_rows[1:]):
        d0 = coarse["fd_directional"]
        d1 = fine["fd_directional"]
        denom = max(abs(d0), abs(d1))
        status = "ok" if denom > fd_signal_floor else "unresolved_fd_signal"
        pairs.append(
            {
                "h_mid": float(np.sqrt(coarse["h"] * fine["h"])),
                "h_coarse": float(coarse["h"]),
                "h_fine": float(fine["h"]),
                "fd_coarse": float(d0),
                "fd_fine": float(d1),
                "status": status,
                "relative_step_change": float(abs(d0 - d1) / max(denom, fd_signal_floor)),
                "max_physical_perturbation_mid": float(
                    np.sqrt(coarse["max_physical_perturbation"] * fine["max_physical_perturbation"])
                ),
            }
        )
    return pairs


def save_plot(results, output_plot):
    fig, ax = plt.subplots(1, 1, figsize=(6.2, 4.2))
    style = {
        "s_only": {"label": "s only", "color": PLOT_COLOR_PRIMARY, "linestyle": "-"},
        "tau_only": {"label": "tau only", "color": PLOT_COLOR_SECONDARY, "linestyle": "--"},
    }

    for result in results:
        pairs = [row for row in result["self_consistency"] if row["status"] == "ok"]
        if not pairs:
            continue
        h = np.array([row["h_mid"] for row in pairs], dtype=np.double)
        err = np.array([row["relative_step_change"] for row in pairs], dtype=np.double)
        kwargs = style[result["name"]]
        ax.loglog(h, err, color=kwargs["color"], linestyle=kwargs["linestyle"], linewidth=1.6, label=kwargs["label"])
        best = min(pairs, key=lambda row: row["relative_step_change"])
        ax.axvline(best["h_mid"], color=kwargs["color"], linewidth=1.0, alpha=0.35)

    ax.set_xlabel("relative perturbation fraction h", fontsize="large")
    ax.set_ylabel("relative FD step change", fontsize="large")
    ax.set_title("Centered FD Self-Consistency", fontsize="medium")
    ax.grid(True, alpha=0.25)
    ax.legend(fontsize="medium")

    fig.tight_layout()
    fig.savefig(output_plot, dpi=200, metadata={"Title": "Elementwise finite-difference stability check"})
    fig.clf()
    plt.close(fig)


def build_parser():
    parser = argparse.ArgumentParser(description="Centered FD self-consistency check for s-only and tau-only directions.")
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
    parser.add_argument("--h-min", type=float, default=1.0e-7)
    parser.add_argument("--h-max", type=float, default=4.0e-1)
    parser.add_argument("--h-count", type=int, default=41)
    parser.add_argument("--seed", type=int, default=9)
    parser.add_argument("--max-mode-x", type=int, default=4)
    parser.add_argument("--max-mode-y", type=int, default=3)
    parser.add_argument("--fd-signal-floor", type=float, default=1.0e-12)
    parser.add_argument("--verbose-remat", action="store_true")
    parser.add_argument("--output-json", type=str, default=str(DEFAULT_JSON))
    parser.add_argument("--output-plot", type=str, default=str(DEFAULT_PLOT))
    return parser


def main():
    args = build_parser().parse_args()
    check.configure_inverse_globals(args)
    h_values = build_h_values(args)

    problem = inv.make_structured_quad_problem()
    nelem = problem["connectivity"].shape[0]
    true_stiffness, true_tau, _ = inv.make_true_fields(problem)
    base_stiffness = np.full(nelem, float(args.base_stiffness), dtype=np.double)
    base_tau = np.full(nelem, float(args.base_tau), dtype=np.double)
    if not check.bounds_ok(base_stiffness, base_tau):
        raise ValueError("Base fields must be inside the configured stiffness/tau bounds.")

    quiet = not args.verbose_remat
    print("Generating synthetic observations...", flush=True)
    observed = check.run_forward_or_adjoint(problem, true_stiffness, true_tau, None, False, quiet)["sensor_history"]

    phi = check.smooth_random_field(problem, args.seed, args.max_mode_x, args.max_mode_y)
    zero = np.zeros(nelem, dtype=np.double)
    directions = [
        ("s_only", (inv.STIFFNESS_MAX - inv.STIFFNESS_MIN) * phi, zero),
        ("tau_only", zero, (inv.TAU_MAX - inv.TAU_MIN) * phi),
    ]

    results = []
    for name, ds, dtau in directions:
        print(f"Running {name} centered FD sweep...", flush=True)
        rows = fd_directional_sweep(problem, base_stiffness, base_tau, observed, ds, dtau, h_values, args.objective, name, quiet)
        consistency = fd_self_consistency(rows, args.fd_signal_floor)
        results.append(
            {
                "name": name,
                "direction_scale": {
                    "max_abs": float(np.max(np.abs(ds if name == "s_only" else dtau))),
                    "meaning": "h is the maximum perturbation as a fraction of this parameter range",
                },
                "fd_rows": rows,
                "self_consistency": consistency,
            }
        )

    output_json = Path(args.output_json)
    output_plot = Path(args.output_plot)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_plot.parent.mkdir(parents=True, exist_ok=True)

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
            "base": {"stiffness": float(args.base_stiffness), "tau": float(args.base_tau)},
            "h_values": h_values,
            "smooth_direction": {
                "seed": int(args.seed),
                "max_mode_x": int(args.max_mode_x),
                "max_mode_y": int(args.max_mode_y),
                "normalization": "max_abs_1",
                "same_spatial_phi_for_stiffness_and_tau": True,
            },
            "metric": "relative change between centered FD quotients at adjacent h values",
            "fd_signal_floor": float(args.fd_signal_floor),
        },
        "results": results,
    }

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)
    save_plot(results, output_plot)

    print("\nRun complete")
    for result in results:
        pairs = [row for row in result["self_consistency"] if row["status"] == "ok"]
        if pairs:
            best = min(pairs, key=lambda row: row["relative_step_change"])
            print(
                f"  {result['name']}: best self-change={best['relative_step_change']:.6e} "
                f"near h={best['h_mid']:.3e}"
            )
    print(f"  wrote json: {output_json}")
    print(f"  wrote plot: {output_plot}")


if __name__ == "__main__":
    main()
