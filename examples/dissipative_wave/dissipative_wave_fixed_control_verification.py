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
    / "reduced_control_fixedpoint_verification.json"
)


def parse_layer_values(text, n_layers):
    values = np.fromstring(text, sep=",", dtype=np.double)
    if values.size != n_layers:
        raise ValueError(f"Expected {n_layers} layer values, got {values.size}.")
    return values


def parse_h_values(text):
    parts = [p.strip() for p in text.split(",") if p.strip()]
    if not parts:
        raise ValueError("At least one FD step must be provided.")
    return [float(p) for p in parts]


def component_rel_error(adj, fd):
    denom = np.maximum.reduce([np.abs(adj), np.abs(fd), np.full_like(adj, 1.0e-14)])
    return np.abs(adj - fd) / denom


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


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reduced-stiffness/reduced-horizon fixed-point adjoint verification "
            "(FD sweep + directional derivative check)."
        )
    )
    parser.add_argument("--nx", type=int, default=60)
    parser.add_argument("--ny", type=int, default=12)
    parser.add_argument("--width", type=float, default=15.0)
    parser.add_argument("--height", type=float, default=3.0)
    parser.add_argument("--nsteps", type=int, default=300)
    parser.add_argument("--nsub-steps", type=int, default=1)
    parser.add_argument("--dt", type=float, default=2.5e-3)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--n-sensors", type=int, default=5)
    parser.add_argument("--sensor-distribution-width", type=float, default=13.5)
    parser.add_argument("--impact-window-width", type=float, default=1.35)
    parser.add_argument("--impact-velocity", type=float, default=0.45)

    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--youngs-modulus", type=float, default=3.0)
    parser.add_argument("--poissons-ratio", type=float, default=0.28)
    parser.add_argument("--shear-modulus-maxwell", type=float, default=1.0)
    parser.add_argument("--mass-damping-factor", type=float, default=0.0)

    parser.add_argument("--true-layers", type=str, default="6.0,3.5")
    parser.add_argument("--init-layers", type=str, default="4.8,4.2")
    parser.add_argument("--true-tau", type=float, default=1.0)
    parser.add_argument("--init-tau", type=float, default=1.3)

    parser.add_argument("--overflow-limit", type=int, default=200)
    parser.add_argument("--mat-overflow-limit", type=int, default=200)

    parser.add_argument(
        "--h-values",
        type=str,
        default="1e-2",
        help="Comma-separated FD step sizes.",
    )
    parser.add_argument(
        "--directional-h",
        type=float,
        default=1.0e-2,
        help="Step size for directional derivative FD.",
    )
    parser.add_argument(
        "--min-fd-step",
        type=float,
        default=1.0e-5,
        help="Reject FD steps smaller than this threshold.",
    )
    parser.add_argument("--output-json", type=str, default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()

    h_values = [h for h in parse_h_values(args.h_values) if h >= args.min_fd_step]
    if not h_values:
        raise ValueError(
            "No admissible FD steps remain after min-fd-step filtering. "
            "Provide larger steps."
        )
    if args.directional_h < args.min_fd_step:
        raise ValueError(
            f"directional-h={args.directional_h} is below min-fd-step={args.min_fd_step}."
        )

    true_layers = parse_layer_values(args.true_layers, args.n_layers)
    init_layers = parse_layer_values(args.init_layers, args.n_layers)

    runtime_args = SimpleNamespace(
        nx=args.nx,
        ny=args.ny,
        width=args.width,
        height=args.height,
        nsteps=args.nsteps,
        nsub_steps=args.nsub_steps,
        dt=args.dt,
        n_layers=args.n_layers,
        n_sensors=args.n_sensors,
        sensor_distribution_width=args.sensor_distribution_width,
        impact_window_width=args.impact_window_width,
        impact_velocity=args.impact_velocity,
        integrator_type="fixed_visco",
        density=args.density,
        youngs_modulus=args.youngs_modulus,
        poissons_ratio=args.poissons_ratio,
        shear_modulus_maxwell=args.shear_modulus_maxwell,
        mass_damping_factor=args.mass_damping_factor,
        overflow_limit=args.overflow_limit,
        mat_overflow_limit=args.mat_overflow_limit,
        adjoint_debug_dump=False,
        adjoint_debug_threshold=0.0,
        adjoint_debug_max_rows=1000,
        adjoint_debug_stride=1,
        dual_overflow_warn=False,
        dual_overflow_warn_fraction=0.95,
        dual_overflow_warn_limit=20,
    )

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
        args.true_tau,
        runtime_args,
        observed_history=None,
        compute_gradients=False,
    )["sensor_history"]

    run = inv.run_forward_or_adjoint(
        problem,
        init_layers,
        args.init_tau,
        runtime_args,
        observed_history=observed,
        compute_gradients=True,
    )
    g_adj = np.concatenate([run["grad_layers"], np.array([run["grad_tau"]], dtype=np.double)])

    fd_rows = []
    for h in h_values:
        fd_layers, fd_tau = inv.finite_difference_gradients(
            problem,
            init_layers,
            args.init_tau,
            runtime_args,
            observed,
            h,
        )
        g_fd = np.concatenate([fd_layers, np.array([fd_tau], dtype=np.double)])
        rel = component_rel_error(g_adj, g_fd)
        fd_rows.append(
            {
                "h": float(h),
                "fd_gradient": g_fd.tolist(),
                "rel_error": rel.tolist(),
                "max_rel_error": float(np.max(rel)),
                "max_abs_error": float(np.max(np.abs(g_adj - g_fd))),
            }
        )

    best_fd = min(fd_rows, key=lambda r: (r["max_rel_error"], r["max_abs_error"]))

    direction = np.array([0.8427009716003844, -0.48154341234307685, 0.24077170617153842], dtype=np.double)
    dir_adj = float(np.dot(g_adj, direction))
    dir_fd = directional_fd(
        problem,
        init_layers,
        args.init_tau,
        runtime_args,
        observed,
        direction,
        args.directional_h,
    )
    dir_rel = abs(dir_adj - dir_fd) / max(abs(dir_adj), abs(dir_fd), 1.0e-14)

    payload = {
        "config": {
            "mesh": [args.nx, args.ny],
            "domain": [args.width, args.height],
            "time": {"nsteps": args.nsteps, "dt": args.dt, "nsub_steps": args.nsub_steps},
            "n_layers": args.n_layers,
            "n_sensors": args.n_sensors,
            "impact_velocity": args.impact_velocity,
            "true_layers": true_layers.tolist(),
            "init_layers": init_layers.tolist(),
            "true_tau": float(args.true_tau),
            "init_tau": float(args.init_tau),
            "integrator": "fixed_visco",
            "overflow_limit": int(args.overflow_limit),
            "mat_overflow_limit": int(args.mat_overflow_limit),
            "h_values": [float(h) for h in h_values],
            "directional_h": float(args.directional_h),
            "min_fd_step": float(args.min_fd_step),
        },
        "loss": float(run["loss"]),
        "adjoint_gradient": g_adj.tolist(),
        "fd_sweep": fd_rows,
        "best_fd_row": best_fd,
        "directional_check": {
            "direction": direction.tolist(),
            "adjoint_directional": dir_adj,
            "fd_directional": dir_fd,
            "rel_error": float(dir_rel),
        },
    }

    output_path = Path(args.output_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    print(f"Saved: {output_path}")
    print(
        "Fixed-point control summary: "
        f"best h={best_fd['h']:.3e}, "
        f"best max_rel={best_fd['max_rel_error']:.6f}, "
        f"directional rel={payload['directional_check']['rel_error']:.6f}"
    )


if __name__ == "__main__":
    main()
