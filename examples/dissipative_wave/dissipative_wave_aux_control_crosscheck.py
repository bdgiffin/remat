import argparse
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import dissipative_wave_inverse as inv


THIS_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    THIS_DIR
    / "inverse_case_outputs"
    / "case3_adjoint_diag"
    / "reduced_control_float_fixed_crosscheck.json"
)


def parse_layer_values(text, n_layers):
    values = np.fromstring(text, sep=",", dtype=np.double)
    if values.size != n_layers:
        raise ValueError(f"Expected {n_layers} layer values, got {values.size}.")
    return values


def parse_h_values(text):
    parts = [p.strip() for p in text.split(",") if p.strip()]
    if not parts:
        raise ValueError("At least one finite-difference step must be provided.")
    return [float(p) for p in parts]


def build_runtime_args(cli_args, integrator_type):
    if integrator_type == "fixed_visco":
        overflow_limit = cli_args.fixed_overflow_limit
        mat_overflow_limit = cli_args.fixed_mat_overflow_limit
    elif integrator_type == "float_visco":
        overflow_limit = cli_args.float_overflow_limit
        mat_overflow_limit = cli_args.float_mat_overflow_limit
    else:
        raise ValueError(f"Unsupported integrator type: {integrator_type}")

    return SimpleNamespace(
        nx=cli_args.nx,
        ny=cli_args.ny,
        width=cli_args.width,
        height=cli_args.height,
        nsteps=cli_args.nsteps,
        nsub_steps=cli_args.nsub_steps,
        dt=cli_args.dt,
        n_layers=cli_args.n_layers,
        n_sensors=cli_args.n_sensors,
        sensor_distribution_width=cli_args.sensor_distribution_width,
        impact_window_width=cli_args.impact_window_width,
        impact_velocity=cli_args.impact_velocity,
        integrator_type=integrator_type,
        density=cli_args.density,
        youngs_modulus=cli_args.youngs_modulus,
        poissons_ratio=cli_args.poissons_ratio,
        shear_modulus_maxwell=cli_args.shear_modulus_maxwell,
        mass_damping_factor=cli_args.mass_damping_factor,
        overflow_limit=overflow_limit,
        mat_overflow_limit=mat_overflow_limit,
        adjoint_debug_dump=False,
        adjoint_debug_threshold=0.0,
        adjoint_debug_max_rows=1000,
        adjoint_debug_stride=1,
        dual_overflow_warn=False,
        dual_overflow_warn_fraction=0.95,
        dual_overflow_warn_limit=20,
    )


def component_rel_error(adj, fd):
    denom = np.maximum.reduce([np.abs(adj), np.abs(fd), np.full_like(adj, 1.0e-14)])
    return np.abs(adj - fd) / denom


def active_mask(adj, fd, activity_tol):
    return np.maximum(np.abs(adj), np.abs(fd)) >= activity_tol


def directional_fd(problem, layers, tau, runtime_args, observed_history, direction, h):
    p0 = np.concatenate([np.asarray(layers, dtype=np.double), np.array([tau], dtype=np.double)])
    p_plus = p0 + h * direction
    p_minus = p0 - h * direction

    l_plus = inv.run_forward_or_adjoint(
        problem,
        p_plus[: runtime_args.n_layers],
        float(p_plus[runtime_args.n_layers]),
        runtime_args,
        observed_history=observed_history,
        compute_gradients=False,
    )["loss"]
    l_minus = inv.run_forward_or_adjoint(
        problem,
        p_minus[: runtime_args.n_layers],
        float(p_minus[runtime_args.n_layers]),
        runtime_args,
        observed_history=observed_history,
        compute_gradients=False,
    )["loss"]
    return float((l_plus - l_minus) / (2.0 * h))


def run_one_integrator(cli_args, integrator_type, true_layers, init_layers):
    runtime_args = build_runtime_args(cli_args, integrator_type)
    problem = inv.make_structured_quad_problem(
        runtime_args.nx,
        runtime_args.ny,
        runtime_args.width,
        runtime_args.height,
        runtime_args.impact_velocity,
        runtime_args.impact_window_width,
        runtime_args.n_sensors,
        runtime_args.n_layers,
        runtime_args.sensor_distribution_width,
    )

    observed = inv.run_forward_or_adjoint(
        problem,
        true_layers,
        cli_args.true_tau,
        runtime_args,
        observed_history=None,
        compute_gradients=False,
    )["sensor_history"]

    run = inv.run_forward_or_adjoint(
        problem,
        init_layers,
        cli_args.init_tau,
        runtime_args,
        observed_history=observed,
        compute_gradients=True,
    )
    g_adj = np.concatenate([run["grad_layers"], np.array([run["grad_tau"]], dtype=np.double)])

    rows = []
    for h in cli_args.h_values:
        fd_layers, fd_tau = inv.finite_difference_gradients(
            problem, init_layers, cli_args.init_tau, runtime_args, observed, h
        )
        g_fd = np.concatenate([fd_layers, np.array([fd_tau], dtype=np.double)])
        rel = component_rel_error(g_adj, g_fd)
        mask = active_mask(g_adj, g_fd, cli_args.activity_tol)
        max_rel_active = float(np.max(rel[mask])) if np.any(mask) else 0.0
        max_abs_error = float(np.max(np.abs(g_adj - g_fd)))
        rows.append(
            {
                "h": float(h),
                "fd_gradient": g_fd.tolist(),
                "rel_error": rel.tolist(),
                "active_components": mask.astype(int).tolist(),
                "max_rel_error_active": max_rel_active,
                "max_abs_error": max_abs_error,
            }
        )

    best = min(rows, key=lambda r: (r["max_rel_error_active"], r["max_abs_error"]))

    direction = np.array([0.7, -0.4, 0.2], dtype=np.double)
    direction = direction / np.linalg.norm(direction)
    dir_adj = float(np.dot(g_adj, direction))
    dir_fd = directional_fd(
        problem,
        init_layers,
        cli_args.init_tau,
        runtime_args,
        observed,
        direction,
        cli_args.directional_h,
    )
    dir_rel = abs(dir_adj - dir_fd) / max(abs(dir_adj), abs(dir_fd), 1.0e-14)

    return {
        "integrator": integrator_type,
        "loss": float(run["loss"]),
        "adjoint_gradient": g_adj.tolist(),
        "fd_sweep": rows,
        "best_fd_row": best,
        "directional_check": {
            "h": float(cli_args.directional_h),
            "direction": direction.tolist(),
            "adjoint_directional": dir_adj,
            "fd_directional": dir_fd,
            "rel_error": float(dir_rel),
        },
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reduced-stiffness/reduced-horizon control cross-check for "
            "fixed_visco vs float_visco adjoint gradients."
        )
    )
    parser.add_argument("--nx", type=int, default=60)
    parser.add_argument("--ny", type=int, default=12)
    parser.add_argument("--width", type=float, default=15.0)
    parser.add_argument("--height", type=float, default=3.0)
    parser.add_argument("--nsteps", type=int, default=450)
    parser.add_argument("--nsub-steps", type=int, default=1)
    parser.add_argument("--dt", type=float, default=2.0e-3)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--n-sensors", type=int, default=5)
    parser.add_argument("--sensor-distribution-width", type=float, default=13.5)
    parser.add_argument("--impact-window-width", type=float, default=1.35)
    parser.add_argument("--impact-velocity", type=float, default=0.5)

    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--youngs-modulus", type=float, default=3.0)
    parser.add_argument("--poissons-ratio", type=float, default=0.28)
    parser.add_argument("--shear-modulus-maxwell", type=float, default=1.0)
    parser.add_argument("--mass-damping-factor", type=float, default=0.0)

    parser.add_argument("--true-layers", type=str, default="6.0,3.5")
    parser.add_argument("--init-layers", type=str, default="4.8,4.2")
    parser.add_argument("--true-tau", type=float, default=0.22)
    parser.add_argument("--init-tau", type=float, default=0.30)

    parser.add_argument("--fixed-overflow-limit", type=int, default=200)
    parser.add_argument("--fixed-mat-overflow-limit", type=int, default=200)
    parser.add_argument("--float-overflow-limit", type=int, default=1000000)
    parser.add_argument("--float-mat-overflow-limit", type=int, default=1000000)

    parser.add_argument(
        "--h-values",
        type=str,
        default="1e-2,3e-3,1e-3,3e-4",
        help="Comma-separated finite-difference step sizes.",
    )
    parser.add_argument(
        "--directional-h",
        type=float,
        default=3.0e-3,
        help="Step size for directional-derivative check.",
    )
    parser.add_argument(
        "--min-fd-step",
        type=float,
        default=1.0e-5,
        help="Reject FD steps smaller than this threshold.",
    )
    parser.add_argument(
        "--activity-tol",
        type=float,
        default=1.0e-6,
        help=(
            "A gradient component is considered active when "
            "max(|adj|,|fd|) >= activity_tol."
        ),
    )
    parser.add_argument("--output-json", type=str, default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()

    h_values = parse_h_values(args.h_values)
    h_values = [h for h in h_values if h >= args.min_fd_step]
    if not h_values:
        raise ValueError(
            "No admissible FD step sizes remain after min-fd-step filtering. "
            "Increase --h-values or reduce --min-fd-step."
        )
    args.h_values = h_values
    if args.directional_h < args.min_fd_step:
        raise ValueError(
            f"directional-h={args.directional_h} is below min-fd-step={args.min_fd_step}."
        )

    true_layers = parse_layer_values(args.true_layers, args.n_layers)
    init_layers = parse_layer_values(args.init_layers, args.n_layers)

    fixed = run_one_integrator(args, "fixed_visco", true_layers, init_layers)
    floating = run_one_integrator(args, "float_visco", true_layers, init_layers)

    g_fixed = np.asarray(fixed["adjoint_gradient"], dtype=np.double)
    g_float = np.asarray(floating["adjoint_gradient"], dtype=np.double)
    cos = float(np.dot(g_fixed, g_float) / (np.linalg.norm(g_fixed) * np.linalg.norm(g_float) + 1.0e-30))
    sign_match = (np.sign(g_fixed) == np.sign(g_float)).astype(int).tolist()

    payload = {
        "config": {
            "nx": args.nx,
            "ny": args.ny,
            "width": args.width,
            "height": args.height,
            "nsteps": args.nsteps,
            "dt": args.dt,
            "n_layers": args.n_layers,
            "n_sensors": args.n_sensors,
            "impact_velocity": args.impact_velocity,
            "h_values": [float(h) for h in args.h_values],
            "directional_h": float(args.directional_h),
            "min_fd_step": float(args.min_fd_step),
            "activity_tol": float(args.activity_tol),
            "true_layers": true_layers.tolist(),
            "init_layers": init_layers.tolist(),
            "true_tau": float(args.true_tau),
            "init_tau": float(args.init_tau),
            "fixed_overflow_limit": int(args.fixed_overflow_limit),
            "fixed_mat_overflow_limit": int(args.fixed_mat_overflow_limit),
        },
        "fixed_visco": fixed,
        "float_visco": floating,
        "cross_integrator": {
            "adjoint_cosine_similarity": cos,
            "adjoint_sign_match_by_component": sign_match,
        },
    }

    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"Saved: {output_path}")
    print(
        "Cross-check summary: "
        f"cos(adj_fixed, adj_float)={payload['cross_integrator']['adjoint_cosine_similarity']:.6f}, "
        f"fixed best h={fixed['best_fd_row']['h']:.3e}, "
        f"float best h={floating['best_fd_row']['h']:.3e}"
    )


if __name__ == "__main__":
    main()
