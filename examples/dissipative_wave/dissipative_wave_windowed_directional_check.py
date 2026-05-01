import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np

import dissipative_wave_inverse as inv
import REMAT


THIS_DIR = Path(__file__).resolve().parent
DEFAULT_OUTPUT = (
    THIS_DIR
    / "inverse_case_outputs"
    / "case3_adjoint_diag"
    / "windowed_directional_check.json"
)


def parse_layer_values(text, n_layers):
    values = np.fromstring(text, sep=",", dtype=np.double)
    if values.size != n_layers:
        raise ValueError(f"Expected {n_layers} layer values, got {values.size}.")
    return values


def make_runtime_args(args):
    return SimpleNamespace(
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


def build_windows(nsteps, window_frac):
    win_len = max(1, int(round(window_frac * nsteps)))
    win_len = min(win_len, nsteps)

    early = (0, win_len)
    mid_start = max(0, (nsteps // 2) - (win_len // 2))
    mid_end = min(nsteps, mid_start + win_len)
    mid_start = max(0, mid_end - win_len)
    mid = (mid_start, mid_end)
    late = (nsteps - win_len, nsteps)

    return [
        {"label": "early", "k_start": int(early[0]), "k_end": int(early[1])},
        {"label": "mid", "k_start": int(mid[0]), "k_end": int(mid[1])},
        {"label": "late", "k_start": int(late[0]), "k_end": int(late[1])},
        {"label": "full", "k_start": 0, "k_end": int(nsteps)},
    ]


def _sum_by_layer(values, elem_layer_ids, n_layers):
    out = np.zeros(n_layers, dtype=np.double)
    for lid in range(n_layers):
        out[lid] = float(np.sum(values[elem_layer_ids == lid]))
    return out


def run_forward_sensor(problem, layer_coeffs, tau, runtime_args):
    inv.configure_run(problem, layer_coeffs, tau, runtime_args)
    sensor_nodes = problem["sensor_nodes"]
    sensor_history = np.zeros((runtime_args.nsteps, sensor_nodes.size), dtype=np.double)
    for k in range(runtime_args.nsteps):
        REMAT.API.update_state(runtime_args.dt, runtime_args.nsub_steps, REMAT.PASS_FORWARD)
        sensor_history[k, :] = REMAT.get_field(b"node", "velocity_Y")[sensor_nodes]
    return sensor_history


def window_loss(sensor_history, observed_history, k_start, k_end):
    residual = sensor_history[k_start:k_end, :] - observed_history[k_start:k_end, :]
    return 0.5 * float(np.sum(residual * residual))


def run_windowed_adjoint_gradient(problem, layer_coeffs, tau, runtime_args, observed_history, k_start, k_end):
    inv.configure_run(problem, layer_coeffs, tau, runtime_args)
    sensor_nodes = problem["sensor_nodes"]
    nsensors = sensor_nodes.size
    sensor_history = np.zeros((runtime_args.nsteps, nsensors), dtype=np.double)

    for k in range(runtime_args.nsteps):
        REMAT.API.update_state(runtime_args.dt, runtime_args.nsub_steps, REMAT.PASS_FORWARD)
        sensor_history[k, :] = REMAT.get_field(b"node", "velocity_Y")[sensor_nodes]

    REMAT.clear_adjoint_state()
    for rev in range(runtime_args.nsteps):
        k = runtime_args.nsteps - 1 - rev
        if k_start <= k < k_end:
            residual = sensor_history[k, :] - observed_history[k, :]
            seed_xy = np.zeros((nsensors, 2), dtype=np.double)
            seed_xy[:, 1] = residual
            REMAT.add_nodal_velocity_adjoint_seed(sensor_nodes, seed_xy)
        REMAT.API.update_state(runtime_args.dt, runtime_args.nsub_steps, REMAT.PASS_BACKWARD_ADJOINT)

    grad_tau = float(REMAT.get_field(b"global", "dL_dparam_relaxation_time")[0])
    elem_grad_scaling = REMAT.get_field(b"element", "dparam_stiffness_scaling_factor")
    grad_layers = _sum_by_layer(elem_grad_scaling, problem["elem_layer_ids"], runtime_args.n_layers)
    grad = np.concatenate([grad_layers, np.array([grad_tau], dtype=np.double)])

    return {
        "gradient": grad,
        "loss": window_loss(sensor_history, observed_history, k_start, k_end),
    }


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Windowed directional-derivative localization check for the Case-3-style fixed-point inverse setup."
        )
    )
    parser.add_argument("--width", type=float, default=15.0)
    parser.add_argument("--height", type=float, default=3.0)
    parser.add_argument("--sensor-distribution-width", type=float, default=13.5)
    parser.add_argument("--nx", type=int, default=150)
    parser.add_argument("--ny", type=int, default=30)
    parser.add_argument("--nsteps", type=int, default=1000)
    parser.add_argument("--nsub-steps", type=int, default=1)
    parser.add_argument("--dt", type=float, default=4.0e-3)
    parser.add_argument("--n-layers", type=int, default=2)
    parser.add_argument("--n-sensors", type=int, default=5)
    parser.add_argument("--impact-velocity", type=float, default=1.0)
    parser.add_argument("--impact-window-width", type=float, default=1.35)
    parser.add_argument("--true-layers", type=str, default="15,7.7")
    parser.add_argument("--init-layers", type=str, default="10.0,4.0")
    parser.add_argument("--true-tau", type=float, default=0.08)
    parser.add_argument("--init-tau", type=float, default=0.10)
    parser.add_argument("--overflow-limit", type=int, default=10)
    parser.add_argument("--mat-overflow-limit", type=int, default=10)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--youngs-modulus", type=float, default=5.0)
    parser.add_argument("--poissons-ratio", type=float, default=0.28)
    parser.add_argument("--shear-modulus-maxwell", type=float, default=2.0)
    parser.add_argument("--mass-damping-factor", type=float, default=0.0)
    parser.add_argument("--h", type=float, default=1.0e-2)
    parser.add_argument(
        "--window-frac",
        type=float,
        default=0.2,
        help="Fraction of total timesteps for each short window (early/mid/late).",
    )
    parser.add_argument("--output-json", type=str, default=str(DEFAULT_OUTPUT))
    args = parser.parse_args()

    if args.h <= 0.0:
        raise ValueError("h must be positive.")
    if not (0.0 < args.window_frac <= 1.0):
        raise ValueError("window-frac must lie in (0,1].")

    true_layers = parse_layer_values(args.true_layers, args.n_layers)
    init_layers = parse_layer_values(args.init_layers, args.n_layers)
    runtime_args = make_runtime_args(args)

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

    windows = build_windows(runtime_args.nsteps, args.window_frac)
    direction = np.array([0.8427009716003844, -0.48154341234307685, 0.24077170617153842], dtype=np.double)
    direction /= np.linalg.norm(direction)

    observed_history = run_forward_sensor(problem, true_layers, args.true_tau, runtime_args)

    p0 = np.concatenate([init_layers, np.array([args.init_tau], dtype=np.double)])
    p_plus = p0 + args.h * direction
    p_minus = p0 - args.h * direction
    sensor_plus = run_forward_sensor(problem, p_plus[: runtime_args.n_layers], float(p_plus[runtime_args.n_layers]), runtime_args)
    sensor_minus = run_forward_sensor(problem, p_minus[: runtime_args.n_layers], float(p_minus[runtime_args.n_layers]), runtime_args)

    rows = []
    for w in windows:
        k_start = int(w["k_start"])
        k_end = int(w["k_end"])
        adj = run_windowed_adjoint_gradient(
            problem,
            init_layers,
            args.init_tau,
            runtime_args,
            observed_history,
            k_start,
            k_end,
        )
        dir_adj = float(np.dot(adj["gradient"], direction))
        loss_plus = window_loss(sensor_plus, observed_history, k_start, k_end)
        loss_minus = window_loss(sensor_minus, observed_history, k_start, k_end)
        dir_fd = float((loss_plus - loss_minus) / (2.0 * args.h))
        rel_err = abs(dir_adj - dir_fd) / max(abs(dir_adj), abs(dir_fd), 1.0e-14)

        rows.append(
            {
                "label": w["label"],
                "k_start": k_start,
                "k_end": k_end,
                "t_start": float(k_start * runtime_args.dt),
                "t_end": float(k_end * runtime_args.dt),
                "window_loss": float(adj["loss"]),
                "adjoint_directional": dir_adj,
                "fd_directional": dir_fd,
                "rel_error": float(rel_err),
            }
        )

    payload = {
        "config": {
            "mesh": [runtime_args.nx, runtime_args.ny],
            "domain": [runtime_args.width, runtime_args.height],
            "time": {
                "nsteps": runtime_args.nsteps,
                "dt": runtime_args.dt,
                "nsub_steps": runtime_args.nsub_steps,
            },
            "true_layers": true_layers.tolist(),
            "init_layers": init_layers.tolist(),
            "true_tau": float(args.true_tau),
            "init_tau": float(args.init_tau),
            "overflow_limit": int(runtime_args.overflow_limit),
            "mat_overflow_limit": int(runtime_args.mat_overflow_limit),
            "integrator": runtime_args.integrator_type,
            "h": float(args.h),
            "window_frac": float(args.window_frac),
            "direction": direction.tolist(),
        },
        "windows": rows,
    }

    out_path = Path(args.output_json)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    csv_path = out_path.with_suffix(".csv")
    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(
            [
                "label",
                "k_start",
                "k_end",
                "t_start",
                "t_end",
                "window_loss",
                "adjoint_directional",
                "fd_directional",
                "rel_error",
            ]
        )
        for row in rows:
            writer.writerow(
                [
                    row["label"],
                    row["k_start"],
                    row["k_end"],
                    f"{row['t_start']:.12e}",
                    f"{row['t_end']:.12e}",
                    f"{row['window_loss']:.12e}",
                    f"{row['adjoint_directional']:.12e}",
                    f"{row['fd_directional']:.12e}",
                    f"{row['rel_error']:.12e}",
                ]
            )

    print(f"Wrote: {out_path}")
    print(f"Wrote: {csv_path}")
    print("Windowed directional-derivative localization (h=1e-2):")
    for row in rows:
        print(
            f"  {row['label']:>5s}  "
            f"k=[{row['k_start']},{row['k_end']})  "
            f"adj={row['adjoint_directional']:.6e}  "
            f"fd={row['fd_directional']:.6e}  "
            f"rel={row['rel_error']:.3e}"
        )


if __name__ == "__main__":
    main()
