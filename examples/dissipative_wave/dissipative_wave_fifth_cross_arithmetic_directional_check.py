import argparse
import json
from pathlib import Path

import numpy as np

import dissipative_wave_inverse as inv
import dissipative_wave_inverse_fifth_best_example as fifth


DEFAULT_JSON = fifth.OUTPUT_DIR / "cross_arithmetic_gradient_direction_check.json"


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


def configure_backend(integrator_type):
    fifth.configure_inverse_backend()
    inv.INTEGRATOR_TYPE = integrator_type
    inv.OVERFLOW_LIMIT = fifth.OVERFLOW_LIMIT
    inv.MAT_OVERFLOW_LIMIT = fifth.MAT_OVERFLOW_LIMIT


def make_fifth_fields():
    configure_backend("fixed_visco")
    experiments = [fifth.make_experiment(center) for center in fifth.IMPACT_CENTERS]
    labels = fifth.make_labels(experiments[0])
    true_stiffness = fifth.expand_region_values(fifth.region_stiffness_values(), labels)
    true_tau = fifth.expand_region_values(fifth.region_tau_values(), labels)
    initial_stiffness = fifth.make_initial_stiffness(labels)
    initial_tau = fifth.make_initial_tau(labels)
    return experiments, labels, true_stiffness, true_tau, initial_stiffness, initial_tau


def objective_and_gradient(integrator_type, experiments, observations, obs_norm_sq, stiffness, tau):
    configure_backend(integrator_type)
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

    return {
        "integrator": integrator_type,
        "loss": float(data_loss / obs_norm_sq),
        "grad_stiffness": grad_stiff / obs_norm_sq,
        "grad_relaxation_time": grad_tau / obs_norm_sq,
    }


def vector_metrics(fixed, floating, activity_tol):
    fixed = np.asarray(fixed, dtype=np.double)
    floating = np.asarray(floating, dtype=np.double)
    diff = fixed - floating
    fixed_norm = float(np.linalg.norm(fixed))
    floating_norm = float(np.linalg.norm(floating))
    diff_norm = float(np.linalg.norm(diff))
    denom = max(fixed_norm, floating_norm, 1.0e-300)
    active = np.maximum(np.abs(fixed), np.abs(floating)) >= activity_tol
    active_count = int(np.count_nonzero(active))
    if active_count:
        sign_agreement = float(np.mean(np.sign(fixed[active]) == np.sign(floating[active])))
    else:
        sign_agreement = None

    return {
        "n_components": int(fixed.size),
        "norm_fixed": fixed_norm,
        "norm_float": floating_norm,
        "norm_difference": diff_norm,
        "relative_l2_mismatch": float(diff_norm / denom),
        "cosine_similarity": float(np.dot(fixed, floating) / (fixed_norm * floating_norm + 1.0e-300)),
        "max_abs_fixed": float(np.max(np.abs(fixed))) if fixed.size else 0.0,
        "max_abs_float": float(np.max(np.abs(floating))) if floating.size else 0.0,
        "max_abs_difference": float(np.max(np.abs(diff))) if diff.size else 0.0,
        "activity_tol": float(activity_tol),
        "active_components": active_count,
        "active_sign_agreement_fraction": sign_agreement,
    }


def directional_projection(name, ds, dtau, fixed_run, float_run):
    fixed_value = float(
        np.dot(fixed_run["grad_stiffness"], ds) + np.dot(fixed_run["grad_relaxation_time"], dtau)
    )
    float_value = float(
        np.dot(float_run["grad_stiffness"], ds) + np.dot(float_run["grad_relaxation_time"], dtau)
    )
    return {
        "name": name,
        "fixed_directional": fixed_value,
        "float_directional": float_value,
        "absolute_difference": float(abs(fixed_value - float_value)),
        "relative_mismatch": float(abs(fixed_value - float_value) / max(abs(fixed_value), abs(float_value), 1.0e-300)),
        "max_abs_stiffness_direction": float(np.max(np.abs(ds))),
        "max_abs_relaxation_time_direction": float(np.max(np.abs(dtau))),
    }


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Cross-arithmetic fixed_visco vs float_visco adjoint-gradient direction "
            "check for the fifth-best dissipative-wave example."
        )
    )
    parser.add_argument("--output-json", type=str, default=str(DEFAULT_JSON))
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--max-mode-x", type=int, default=4)
    parser.add_argument("--max-mode-z", type=int, default=3)
    parser.add_argument("--stiffness-step-fraction", type=float, default=0.05)
    parser.add_argument("--tau-step-fraction", type=float, default=0.05)
    parser.add_argument("--activity-tol", type=float, default=1.0e-12)
    return parser


def main():
    args = build_parser().parse_args()
    experiments, labels, true_stiffness, true_tau, initial_stiffness, initial_tau = make_fifth_fields()
    nelem = labels.size
    sensor_history_shape = list((fifth.N_STEPS, int(experiments[0]["sensor_nodes"].size), 2))

    print("Generating fifth-example fixed-point synthetic observations...", flush=True)
    configure_backend("fixed_visco")
    observations, obs_norm_sq = fifth.generate_observations(experiments, true_stiffness, true_tau)

    print("Computing fixed_visco gradient at the fifth initial point...", flush=True)
    fixed_run = objective_and_gradient(
        "fixed_visco",
        experiments,
        observations,
        obs_norm_sq,
        initial_stiffness,
        initial_tau,
    )

    print("Computing float_visco gradient at the same fifth initial point...", flush=True)
    float_run = objective_and_gradient(
        "float_visco",
        experiments,
        observations,
        obs_norm_sq,
        initial_stiffness,
        initial_tau,
    )

    fixed_grad = np.concatenate([fixed_run["grad_stiffness"], fixed_run["grad_relaxation_time"]])
    float_grad = np.concatenate([float_run["grad_stiffness"], float_run["grad_relaxation_time"]])

    shared_phi = smooth_random_field(experiments[0], args.seed, args.max_mode_x, args.max_mode_z)
    zero = np.zeros(nelem, dtype=np.double)
    ds = args.stiffness_step_fraction * (fifth.STIFFNESS_MAX - fifth.STIFFNESS_MIN) * shared_phi
    dtau = args.tau_step_fraction * (fifth.TAU_MAX - fifth.TAU_MIN) * shared_phi

    comparisons = {
        "all_controls": vector_metrics(fixed_grad, float_grad, args.activity_tol),
        "stiffness": vector_metrics(
            fixed_run["grad_stiffness"],
            float_run["grad_stiffness"],
            args.activity_tol,
        ),
        "relaxation_time": vector_metrics(
            fixed_run["grad_relaxation_time"],
            float_run["grad_relaxation_time"],
            args.activity_tol,
        ),
    }
    directionals = [
        directional_projection("stiffness_only_shared_phi", ds, zero, fixed_run, float_run),
        directional_projection("relaxation_time_only_shared_phi", zero, dtau, fixed_run, float_run),
    ]

    payload = {
        "config": {
            "example": "fifth_best_elementwise_stiffness_and_relaxation_time",
            "mesh": {
                "nx": fifth.NX,
                "nz": fifth.NZ,
                "width": fifth.WIDTH,
                "depth": fifth.DEPTH,
                "n_elements": int(nelem),
            },
            "time": {
                "dt": fifth.DT,
                "n_steps": fifth.N_STEPS,
                "n_sub_steps": fifth.N_SUB_STEPS,
                "physical_horizon": float(fifth.DT * fifth.N_STEPS),
            },
            "observation_integrator": "fixed_visco",
            "gradient_integrators": ["fixed_visco", "float_visco"],
            "objective": {
                "type": "normalized data misfit",
                "regularization_included": False,
            },
            "impact": {
                "n_impacts": len(fifth.IMPACT_CENTERS),
                "centers": [float(center) for center in fifth.IMPACT_CENTERS],
                "window_width": fifth.IMPACT_WINDOW_WIDTH,
                "velocity": fifth.IMPACT_VELOCITY,
            },
            "sensors": {
                "type": "sparse top-surface nodal velocity sensors",
                "components": ["velocity_x", "velocity_z"],
                "n_sensors_per_impact": int(experiments[0]["sensor_nodes"].size),
                "x_targets": [float(x) for x in experiments[0]["sensor_x_targets"]],
                "x_positions": [float(x) for x in experiments[0]["sensor_x_positions"]],
                "history_shape_per_impact": sensor_history_shape,
                "total_scalar_samples": int(len(experiments) * np.prod(sensor_history_shape)),
            },
            "unknowns": {
                "control_order": "[s_1,...,s_ne,tau_1,...,tau_ne]",
                "n_stiffness_unknowns": int(nelem),
                "n_relaxation_time_unknowns": int(nelem),
                "total_unknowns": int(2 * nelem),
            },
            "bounds": {
                "stiffness": [fifth.STIFFNESS_MIN, fifth.STIFFNESS_MAX],
                "relaxation_time": [fifth.TAU_MIN, fifth.TAU_MAX],
            },
            "region_values": {
                "true_stiffness": fifth.region_stiffness_values().tolist(),
                "initial_stiffness": fifth.initial_stiffness_values().tolist(),
                "true_relaxation_time": fifth.region_tau_values().tolist(),
                "initial_relaxation_time": fifth.initial_tau_values().tolist(),
            },
            "finite_precision_regime": {
                "larger_relaxation_times_than_case3_diagnostic": True,
                "short_horizon_used_for_reversibility": True,
                "reason": (
                    "The comparison is made in the finite-precision regime used by the fifth "
                    "example so rematerialized backward states remain numerically meaningful."
                ),
            },
            "smooth_direction": {
                "seed": int(args.seed),
                "max_mode_x": int(args.max_mode_x),
                "max_mode_z": int(args.max_mode_z),
                "same_spatial_phi_for_stiffness_and_relaxation_time": True,
                "raw_normalization": "max_abs_1",
                "stiffness_scale_at_h1": float(np.max(np.abs(ds))),
                "relaxation_time_scale_at_h1": float(np.max(np.abs(dtau))),
                "stiffness_step_fraction": float(args.stiffness_step_fraction),
                "relaxation_time_step_fraction": float(args.tau_step_fraction),
            },
        },
        "obs_norm_sq": float(obs_norm_sq),
        "losses": {
            "fixed_visco": fixed_run["loss"],
            "float_visco": float_run["loss"],
        },
        "gradient_comparison": comparisons,
        "directional_projections": directionals,
    }

    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2)

    print("\nRun complete")
    print(f"  fixed loss          : {fixed_run['loss']:.6e}")
    print(f"  float loss          : {float_run['loss']:.6e}")
    print(f"  full-gradient cosine: {comparisons['all_controls']['cosine_similarity']:.8f}")
    print(f"  full-gradient rel L2: {comparisons['all_controls']['relative_l2_mismatch']:.6e}")
    for item in directionals:
        print(
            f"  {item['name']}: fixed={item['fixed_directional']:.6e}, "
            f"float={item['float_directional']:.6e}, rel={item['relative_mismatch']:.6e}"
        )
    print(f"  wrote json: {output_path}")


if __name__ == "__main__":
    main()
