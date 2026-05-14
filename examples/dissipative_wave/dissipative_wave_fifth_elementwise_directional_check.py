import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

import dissipative_wave_inverse as inv
import dissipative_wave_inverse_fifth_best_example as fifth


DEFAULT_JSON = fifth.OUTPUT_DIR / "elementwise_s_tau_directional_check.json"
DEFAULT_PLOT = fifth.OUTPUT_DIR / "elementwise_s_tau_directional_check.svg"

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


def smooth_random_field(problem, seed, max_mode_x, max_mode_z):
    rng = np.random.default_rng(seed)
    centers = np.asarray(problem["elem_centers"], dtype=np.double)
    x = centers[:, 0] / float(fifth.WIDTH)
    z = (fifth.DEPTH - centers[:, 1]) / float(fifth.DEPTH)

    field = np.zeros(x.size, dtype=np.double)
    for kx in range(max_mode_x + 1):
        for kz in range(max_mode_z + 1):
            if kx == 0 and kz == 0:
                continue
            weight = 1.0 / (1.0 + kx * kx + kz * kz)
            phase_x = rng.uniform(0.0, 2.0 * np.pi)
            phase_z = rng.uniform(0.0, 2.0 * np.pi)
            coeff = rng.normal() * weight
            field += coeff * np.cos(np.pi * kx * x + phase_x) * np.cos(np.pi * kz * z + phase_z)

    field -= np.mean(field)
    max_abs = float(np.max(np.abs(field)))
    if max_abs <= 0.0:
        raise ValueError("Smooth random field has zero amplitude.")
    return field / max_abs


def make_fifth_fields():
    fifth.configure_inverse_backend()
    experiments = [fifth.make_experiment(center) for center in fifth.IMPACT_CENTERS]
    labels = fifth.make_labels(experiments[0])
    true_stiffness = fifth.expand_region_values(fifth.region_stiffness_values(), labels)
    true_tau = fifth.expand_region_values(fifth.region_tau_values(), labels)
    initial_stiffness = fifth.make_initial_stiffness(labels)
    initial_tau = fifth.make_initial_tau(labels)
    return experiments, labels, true_stiffness, true_tau, initial_stiffness, initial_tau


def load_recovered_fields():
    fields = fifth.load_saved_fields()
    return fields["recovered_stiffness"], fields["recovered_tau"]


def regularization_loss_and_grad(stiffness, tau, objective):
    nelem = stiffness.size
    if objective == "data":
        return 0.0, np.zeros(nelem, dtype=np.double), np.zeros(nelem, dtype=np.double)

    i_idx, j_idx = inv.build_edge_pairs(fifth.NX, fifth.NZ)
    reg_s_loss, reg_s_grad = inv.regularization_loss_and_grad(
        stiffness,
        i_idx,
        j_idx,
        fifth.REG_L2_STIFFNESS,
        fifth.REG_TV_STIFFNESS,
        fifth.REG_TV_EPS,
    )
    reg_t_loss, reg_t_grad = inv.regularization_loss_and_grad(
        tau,
        i_idx,
        j_idx,
        fifth.REG_L2_TAU,
        fifth.REG_TV_TAU,
        fifth.REG_TV_EPS,
    )
    return (
        float((reg_s_loss + reg_t_loss) / nelem),
        reg_s_grad / nelem,
        reg_t_grad / nelem,
    )


def objective_value(experiments, observations, obs_norm_sq, stiffness, tau, objective):
    data_loss = 0.0
    for problem, observed in zip(experiments, observations):
        run = fifth.run_velocity_history(
            problem,
            stiffness,
            tau,
            observed_history=observed,
            compute_gradients=False,
        )
        data_loss += run["data_loss"]

    reg_loss, _, _ = regularization_loss_and_grad(stiffness, tau, objective)
    return float(data_loss / obs_norm_sq + reg_loss)


def objective_and_gradient(experiments, observations, obs_norm_sq, stiffness, tau, objective):
    nelem = stiffness.size
    data_loss = 0.0
    grad_stiff = np.zeros(nelem, dtype=np.double)
    grad_tau = np.zeros(nelem, dtype=np.double)

    for problem, observed in zip(experiments, observations):
        run = fifth.run_velocity_history(
            problem,
            stiffness,
            tau,
            observed_history=observed,
            compute_gradients=True,
        )
        data_loss += run["data_loss"]
        grad_stiff += run["grad_stiff"]
        grad_tau += run["grad_tau"]

    data_loss /= obs_norm_sq
    grad_stiff /= obs_norm_sq
    grad_tau /= obs_norm_sq

    reg_loss, reg_s_grad, reg_t_grad = regularization_loss_and_grad(stiffness, tau, objective)
    return (
        float(data_loss + reg_loss),
        grad_stiff + reg_s_grad,
        grad_tau + reg_t_grad,
    )


def bounds_ok(stiffness, tau):
    return (
        np.all(stiffness >= fifth.STIFFNESS_MIN)
        and np.all(stiffness <= fifth.STIFFNESS_MAX)
        and np.all(tau >= fifth.TAU_MIN)
        and np.all(tau <= fifth.TAU_MAX)
    )


def check_direction(
    name,
    experiments,
    observations,
    obs_norm_sq,
    base_stiffness,
    base_tau,
    base_loss,
    grad_stiff,
    grad_tau,
    ds,
    dtau,
    h_values,
    objective,
):
    adjoint_directional = float(np.dot(grad_stiff, ds) + np.dot(grad_tau, dtau))
    rows = []
    active = ds if name == "stiffness" else dtau

    for h in h_values:
        print(f"  {name}: h={h:.3e}", flush=True)
        plus_s = base_stiffness + h * ds
        plus_t = base_tau + h * dtau
        minus_s = base_stiffness - h * ds
        minus_t = base_tau - h * dtau

        if not bounds_ok(plus_s, plus_t) or not bounds_ok(minus_s, minus_t):
            rows.append({"h": float(h), "status": "skipped_out_of_bounds"})
            continue

        loss_plus = objective_value(experiments, observations, obs_norm_sq, plus_s, plus_t, objective)
        loss_minus = objective_value(experiments, observations, obs_norm_sq, minus_s, minus_t, objective)
        fd_directional = (loss_plus - loss_minus) / (2.0 * h)
        fd_taylor_remainder = abs(loss_plus - base_loss - h * fd_directional)
        adjoint_taylor_remainder = abs(loss_plus - base_loss - h * adjoint_directional)

        rows.append(
            {
                "h": float(h),
                "status": "ok",
                "loss_plus": float(loss_plus),
                "loss_minus": float(loss_minus),
                "adjoint_directional": adjoint_directional,
                "centered_fd_directional": float(fd_directional),
                "fd_taylor_remainder": float(fd_taylor_remainder),
                "adjoint_taylor_remainder": float(adjoint_taylor_remainder),
                "max_physical_perturbation": float(h * np.max(np.abs(active))),
            }
        )

    ok_rows = [row for row in rows if row["status"] == "ok"]
    for prev, curr in zip(ok_rows, ok_rows[1:]):
        if prev["fd_taylor_remainder"] > 0.0 and curr["fd_taylor_remainder"] > 0.0:
            curr["fd_taylor_remainder_order_from_previous"] = float(
                np.log(prev["fd_taylor_remainder"] / curr["fd_taylor_remainder"])
                / np.log(prev["h"] / curr["h"])
            )
        if prev["adjoint_taylor_remainder"] > 0.0 and curr["adjoint_taylor_remainder"] > 0.0:
            curr["adjoint_taylor_remainder_order_from_previous"] = float(
                np.log(prev["adjoint_taylor_remainder"] / curr["adjoint_taylor_remainder"])
                / np.log(prev["h"] / curr["h"])
            )

    return {
        "name": name,
        "adjoint_directional": adjoint_directional,
        "direction_stats": {
            "max_abs": float(np.max(np.abs(active))),
            "mean": float(np.mean(active)),
            "l2": float(np.linalg.norm(active)),
            "physical_units": "stiffness scaling factor" if name == "stiffness" else "relaxation time",
        },
        "rows": rows,
    }


def save_plot(results, output_plot):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True)

    style = {
        "stiffness": {"label": "stiffness", "color": PLOT_COLOR_PRIMARY, "linestyle": "-"},
        "relaxation_time": {"label": "relaxation time", "color": PLOT_COLOR_SECONDARY, "linestyle": "--"},
    }

    for result in results:
        rows = [row for row in result["rows"] if row["status"] == "ok"]
        if not rows:
            continue
        kwargs = style[result["name"]]

        fd_rows = [row for row in rows if row["fd_taylor_remainder"] > 0.0]
        if fd_rows:
            axes[0].loglog(
                [row["h"] for row in fd_rows],
                [row["fd_taylor_remainder"] for row in fd_rows],
                color=kwargs["color"],
                linestyle=kwargs["linestyle"],
                linewidth=1.6,
                label=kwargs["label"],
            )

        adjoint_rows = [row for row in rows if row["adjoint_taylor_remainder"] > 0.0]
        if adjoint_rows:
            axes[1].loglog(
                [row["h"] for row in adjoint_rows],
                [row["adjoint_taylor_remainder"] for row in adjoint_rows],
                color=kwargs["color"],
                linestyle=kwargs["linestyle"],
                linewidth=1.6,
                label=kwargs["label"],
            )

    reference_rows = [
        row
        for row in results[0]["rows"]
        if row["status"] == "ok" and row["fd_taylor_remainder"] > 0.0 and row["adjoint_taylor_remainder"] > 0.0
    ]
    if reference_rows:
        h_ref = np.array([row["h"] for row in reference_rows], dtype=np.double)
        fd_ref = float(reference_rows[0]["fd_taylor_remainder"])
        adj_ref = float(reference_rows[0]["adjoint_taylor_remainder"])
        if fd_ref > 0.0:
            axes[0].loglog(
                h_ref,
                fd_ref * (h_ref / h_ref[0]) ** 2,
                color="0.35",
                linewidth=1.0,
                dashes=(6, 8),
                alpha=0.7,
                label=r"$O(h^2)$ reference",
            )
        if adj_ref > 0.0:
            axes[1].loglog(
                h_ref,
                adj_ref * (h_ref / h_ref[0]) ** 2,
                color="0.35",
                linewidth=1.0,
                dashes=(6, 8),
                alpha=0.7,
                label=r"$O(h^2)$ reference",
            )

    axes[0].set_xlabel("dimensionless step h", fontsize="large")
    axes[0].set_ylabel("second Taylor remainder", fontsize="large")
    axes[0].set_title("FD Taylor remainder", fontsize="medium")
    axes[0].legend(fontsize="medium")

    axes[1].set_xlabel("dimensionless step h", fontsize="large")
    axes[1].set_title("Adjoint Taylor remainder", fontsize="medium")
    axes[1].legend(fontsize="medium")

    fig.tight_layout()
    fig.savefig(
        output_plot,
        dpi=200,
        metadata={"Title": "Fifth example elementwise stiffness and relaxation time Taylor checks"},
    )
    plt.close(fig)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Directional finite-difference checks for the fifth-best elementwise stiffness and relaxation-time example."
    )
    parser.add_argument("--base", choices=("initial", "truth", "recovered"), default="initial")
    parser.add_argument("--objective", choices=("data", "total"), default="data")
    parser.add_argument("--h-values", type=str, default=None)
    parser.add_argument("--h-min", type=float, default=1.0e-6)
    parser.add_argument("--h-max", type=float, default=4.0)
    parser.add_argument("--h-count", type=int, default=33)
    parser.add_argument(
        "--stiffness-step-fraction",
        type=float,
        default=0.05,
        help="At h=1, max |delta stiffness| is this fraction of the stiffness bound width.",
    )
    parser.add_argument(
        "--tau-step-fraction",
        type=float,
        default=0.05,
        help="At h=1, max |delta relaxation time| is this fraction of the relaxation-time bound width.",
    )
    parser.add_argument("--stiffness-direction-scale", type=float, default=None)
    parser.add_argument("--tau-direction-scale", type=float, default=None)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--max-mode-x", type=int, default=4)
    parser.add_argument("--max-mode-z", type=int, default=3)
    parser.add_argument("--output-json", type=str, default=str(DEFAULT_JSON))
    parser.add_argument("--output-plot", type=str, default=str(DEFAULT_PLOT))
    return parser


def main():
    args = build_parser().parse_args()
    h_values = build_h_values(args)

    experiments, labels, true_stiffness, true_tau, initial_stiffness, initial_tau = make_fifth_fields()
    nelem = labels.size

    if args.base == "initial":
        base_stiffness = initial_stiffness.copy()
        base_tau = initial_tau.copy()
    elif args.base == "truth":
        base_stiffness = true_stiffness.copy()
        base_tau = true_tau.copy()
    else:
        base_stiffness, base_tau = load_recovered_fields()

    if not bounds_ok(base_stiffness, base_tau):
        raise ValueError("Base fields must be inside the configured stiffness and relaxation-time bounds.")

    print("Generating fifth-example synthetic observations...", flush=True)
    observations, obs_norm_sq = fifth.generate_observations(experiments, true_stiffness, true_tau)

    print("Computing fifth-example elementwise adjoint gradient at the base point...", flush=True)
    base_loss, grad_stiff, grad_tau = objective_and_gradient(
        experiments,
        observations,
        obs_norm_sq,
        base_stiffness,
        base_tau,
        args.objective,
    )

    shared_phi = smooth_random_field(experiments[0], args.seed, args.max_mode_x, args.max_mode_z)
    stiffness_scale = (
        float(args.stiffness_direction_scale)
        if args.stiffness_direction_scale is not None
        else float(args.stiffness_step_fraction) * (fifth.STIFFNESS_MAX - fifth.STIFFNESS_MIN)
    )
    tau_scale = (
        float(args.tau_direction_scale)
        if args.tau_direction_scale is not None
        else float(args.tau_step_fraction) * (fifth.TAU_MAX - fifth.TAU_MIN)
    )
    relaxation_time_scale_mode = (
        "manual" if args.tau_direction_scale is not None else "same_fraction_of_relaxation_time_range"
    )

    zero = np.zeros(nelem, dtype=np.double)

    print(f"Using shared smooth spatial direction phi from seed={args.seed}.", flush=True)
    print(f"  stiffness scale at h=1: {stiffness_scale:.6e}", flush=True)
    print(f"  relaxation time scale at h=1: {tau_scale:.6e} ({relaxation_time_scale_mode})", flush=True)

    results = []
    print("Running stiffness directional FD sweep...", flush=True)
    results.append(
        check_direction(
            "stiffness",
            experiments,
            observations,
            obs_norm_sq,
            base_stiffness,
            base_tau,
            base_loss,
            grad_stiff,
            grad_tau,
            stiffness_scale * shared_phi,
            zero,
            h_values,
            args.objective,
        )
    )

    print("Running relaxation_time directional FD sweep...", flush=True)
    results.append(
        check_direction(
            "relaxation_time",
            experiments,
            observations,
            obs_norm_sq,
            base_stiffness,
            base_tau,
            base_loss,
            grad_stiff,
            grad_tau,
            zero,
            tau_scale * shared_phi,
            h_values,
            args.objective,
        )
    )

    for result in results:
        skipped = sum(1 for row in result["rows"] if row["status"] != "ok")
        if skipped:
            print(f"  {result['name']}: skipped {skipped} out-of-bounds large-step rows.", flush=True)

    payload = {
        "config": {
            "example": "fifth_best_elementwise_stiffness_and_relaxation_time",
            "mesh": {"nx": fifth.NX, "nz": fifth.NZ, "width": fifth.WIDTH, "depth": fifth.DEPTH},
            "time": {"dt": fifth.DT, "n_steps": fifth.N_STEPS, "n_sub_steps": fifth.N_SUB_STEPS},
            "integrator_type": fifth.INTEGRATOR_TYPE,
            "impact_centers": list(fifth.IMPACT_CENTERS),
            "sensors": {
                "type": "sparse surface nodal velocity sensors",
                "components": ["velocity_x", "velocity_z"],
                "n_sensors_per_impact": int(experiments[0]["sensor_nodes"].size),
            },
            "objective": args.objective,
            "base": args.base,
            "bounds": {
                "stiffness": [fifth.STIFFNESS_MIN, fifth.STIFFNESS_MAX],
                "relaxation_time": [fifth.TAU_MIN, fifth.TAU_MAX],
            },
            "h_values": h_values,
            "smooth_direction": {
                "seed": int(args.seed),
                "max_mode_x": int(args.max_mode_x),
                "max_mode_z": int(args.max_mode_z),
                "same_spatial_phi_for_stiffness_and_relaxation_time": True,
                "raw_normalization": "max_abs_1",
                "stiffness_scale_at_h1": float(stiffness_scale),
                "relaxation_time_scale_at_h1": float(tau_scale),
                "stiffness_step_fraction": float(args.stiffness_step_fraction),
                "relaxation_time_step_fraction": float(args.tau_step_fraction),
                "relaxation_time_scale_mode": relaxation_time_scale_mode,
                "relaxation_time_alignment_to_stiffness_best_h": False,
            },
            "unknowns": {
                "n_stiffness_unknowns": int(nelem),
                "n_relaxation_time_unknowns": int(nelem),
                "relaxation_time_parameterization": "independent element-wise relaxation time",
            },
        },
        "obs_norm_sq": float(obs_norm_sq),
        "base_loss": float(base_loss),
        "gradient_norms": {
            "stiffness_l2": float(np.linalg.norm(grad_stiff)),
            "relaxation_time_l2": float(np.linalg.norm(grad_tau)),
            "stiffness_inf": float(np.max(np.abs(grad_stiff))),
            "relaxation_time_inf": float(np.max(np.abs(grad_tau))),
        },
        "checks": results,
    }

    output_json = Path(args.output_json)
    output_plot = Path(args.output_plot)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_plot.parent.mkdir(parents=True, exist_ok=True)

    with output_json.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    save_plot(results, output_plot)

    print("\nRun complete")
    print(f"  base loss : {base_loss:.6e}")
    for result in results:
        ok_rows = [row for row in result["rows"] if row["status"] == "ok"]
        if ok_rows:
            positive_fd = [row for row in ok_rows if row["fd_taylor_remainder"] > 0.0]
            positive_adjoint = [row for row in ok_rows if row["adjoint_taylor_remainder"] > 0.0]
            if positive_fd:
                best_fd = min(positive_fd, key=lambda row: row["fd_taylor_remainder"])
                print(
                    f"  {result['name']}: min positive FD remainder={best_fd['fd_taylor_remainder']:.6e} "
                    f"at h={best_fd['h']:.3e}"
                )
            if positive_adjoint:
                best_adjoint = min(positive_adjoint, key=lambda row: row["adjoint_taylor_remainder"])
                print(
                    f"  {result['name']}: min positive adjoint remainder="
                    f"{best_adjoint['adjoint_taylor_remainder']:.6e} at h={best_adjoint['h']:.3e}"
                )
    print(f"  wrote json: {output_json}")
    print(f"  wrote plot: {output_plot}")


if __name__ == "__main__":
    main()
