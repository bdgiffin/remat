import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
sys.path.append(str(REPO_ROOT / "install" / "package"))

import REMAT


_ACTIVE_ELEM_COEFFS = np.ones(1, dtype=np.double)
_ELEM_DX = 1.0
_ELEM_DY = 1.0
_ELEM_NX = 1
_ELEM_NY = 1
_ELEM_WIDTH = 1.0
_ELEM_HEIGHT = 1.0


def elementwise_stiffness_scaling(x, y):
    x_clamped = min(max(float(x), 0.0), _ELEM_WIDTH - 1.0e-12)
    y_clamped = min(max(float(y), 0.0), _ELEM_HEIGHT - 1.0e-12)
    ix = int(np.clip(np.floor(x_clamped / _ELEM_DX), 0, _ELEM_NX - 1))
    iy = int(np.clip(np.floor(y_clamped / _ELEM_DY), 0, _ELEM_NY - 1))
    eid = iy * _ELEM_NX + ix
    return float(_ACTIVE_ELEM_COEFFS[eid])


def make_structured_quad_problem(
    nx,
    ny,
    width,
    height,
    impact_velocity,
    impact_window_width,
    nsensors,
    sensor_distribution_width=None,
):
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
    top_nodes = np.where(np.abs(node_y - height) < eps)[0]

    fixity[bottom_nodes, :] = True

    xmid = 0.5 * width
    impact_span = min(float(impact_window_width), width)
    impact_half_width = 0.5 * impact_span
    impact_x_lo = xmid - impact_half_width
    impact_x_hi = xmid + impact_half_width
    impact_nodes = top_nodes[
        (node_x[top_nodes] >= impact_x_lo - eps) & (node_x[top_nodes] <= impact_x_hi + eps)
    ]
    if impact_nodes.size == 0:
        nearest = int(np.argmin(np.abs(node_x[top_nodes] - xmid)))
        impact_nodes = top_nodes[np.array([nearest], dtype=np.int64)]
    velocities[impact_nodes, 1] = -abs(float(impact_velocity))

    top_interior = top_nodes[(node_x[top_nodes] > eps) & (node_x[top_nodes] < width - eps)]
    if sensor_distribution_width is None:
        sensor_distribution_width = width
    sensor_span = min(float(sensor_distribution_width), width)
    sensor_x_lo = xmid - 0.5 * sensor_span
    sensor_x_hi = xmid + 0.5 * sensor_span
    sensor_candidates = top_interior[
        (node_x[top_interior] >= sensor_x_lo - eps) & (node_x[top_interior] <= sensor_x_hi + eps)
    ]
    if nsensors > sensor_candidates.size:
        nsensors = sensor_candidates.size
    sensor_pick = np.unique(np.round(np.linspace(0, sensor_candidates.size - 1, nsensors)).astype(int))
    sensor_nodes = sensor_candidates[sensor_pick]

    return {
        "coordinates": coordinates,
        "velocities": velocities,
        "fixity": fixity,
        "connectivity": connectivity,
        "sensor_nodes": sensor_nodes.astype(np.int32),
        "width": float(width),
        "height": float(height),
        "nx": int(nx),
        "ny": int(ny),
    }


def configure_run(problem, elem_coeffs, tau, args):
    global _ACTIVE_ELEM_COEFFS, _ELEM_DX, _ELEM_DY, _ELEM_NX, _ELEM_NY, _ELEM_WIDTH, _ELEM_HEIGHT
    _ACTIVE_ELEM_COEFFS = np.asarray(elem_coeffs, dtype=np.double).copy()
    _ELEM_NX = int(problem["nx"])
    _ELEM_NY = int(problem["ny"])
    _ELEM_WIDTH = float(problem["width"])
    _ELEM_HEIGHT = float(problem["height"])
    _ELEM_DX = _ELEM_WIDTH / float(_ELEM_NX)
    _ELEM_DY = _ELEM_HEIGHT / float(_ELEM_NY)

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
    REMAT.API.define_parameter(b"relaxation_time", float(tau))
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
    REMAT.define_variable_properties(elementwise_stiffness_scaling)
    REMAT.API.initialize()


def run_forward_or_adjoint(problem, elem_coeffs, tau, args, observed_history=None, compute_gradients=False):
    configure_run(problem, elem_coeffs, tau, args)
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

    elem_grad = np.zeros(problem["connectivity"].shape[0], dtype=np.double)
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
        elem_grad = np.asarray(
            REMAT.get_field(b"element", "dparam_stiffness_scaling_factor"), dtype=np.double
        ).reshape(-1)

    return {
        "loss": float(loss),
        "sensor_history": sensor_history,
        "elem_grad": elem_grad,
    }


def fd_element_gradient(problem, base_coeffs, tau, args, observed_history, h):
    nelem = base_coeffs.size
    out = np.zeros(nelem, dtype=np.double)
    for e in range(nelem):
        plus = base_coeffs.copy()
        minus = base_coeffs.copy()
        plus[e] += h
        minus[e] -= h
        lp = run_forward_or_adjoint(problem, plus, tau, args, observed_history, compute_gradients=False)["loss"]
        lm = run_forward_or_adjoint(problem, minus, tau, args, observed_history, compute_gradients=False)["loss"]
        out[e] = (lp - lm) / (2.0 * h)
    return out


def choose_true_coeffs(base_coeffs, nx, ny):
    true_coeffs = base_coeffs.copy()
    ix_center = nx // 2
    iy_top = ny - 1
    iy_bottom = 0
    center_top = iy_top * nx + ix_center
    center_bottom = iy_bottom * nx + ix_center
    true_coeffs[center_top] = 1.35
    true_coeffs[center_bottom] = 0.75
    if ix_center - 1 >= 0:
        true_coeffs[iy_top * nx + (ix_center - 1)] = 1.20
    if ix_center + 1 < nx:
        true_coeffs[iy_bottom * nx + (ix_center + 1)] = 0.85
    return true_coeffs


def summarize(adj, fd):
    eps = 1.0e-14
    rel = np.abs(adj - fd) / np.maximum.reduce([np.abs(adj), np.abs(fd), np.full_like(adj, eps)])
    sign_ok = np.sign(adj) == np.sign(fd)
    dot = float(np.dot(adj, fd))
    na = float(np.linalg.norm(adj))
    nf = float(np.linalg.norm(fd))
    cos = dot / max(na * nf, eps)
    nd = float(np.linalg.norm(adj - fd))
    metrics = {
        "max_rel_error": float(np.max(rel)),
        "mean_rel_error": float(np.mean(rel)),
        "median_rel_error": float(np.median(rel)),
        "sign_agreement_fraction": float(np.mean(sign_ok.astype(np.double))),
        "cosine_similarity": float(cos),
        "adjoint_norm_l2": na,
        "fd_norm_l2": nf,
        "diff_norm_l2": nd,
        "relative_l2_mismatch": nd / max(na, nf, eps),
        "dot_product": dot,
    }

    active_thresholds = [1.0e-2, 1.0e-3, 1.0e-4]
    active_stats = {}
    amp = np.maximum(np.abs(adj), np.abs(fd))
    for thr in active_thresholds:
        mask = amp >= thr
        if not np.any(mask):
            continue
        adj_m = adj[mask]
        fd_m = fd[mask]
        nd_m = float(np.linalg.norm(adj_m - fd_m))
        na_m = float(np.linalg.norm(adj_m))
        nf_m = float(np.linalg.norm(fd_m))
        sign_m = np.sign(adj_m) == np.sign(fd_m)
        active_stats[f"thr_{thr:.0e}"] = {
            "count": int(np.count_nonzero(mask)),
            "relative_l2_mismatch": nd_m / max(na_m, nf_m, eps),
            "sign_agreement_fraction": float(np.mean(sign_m.astype(np.double))),
        }

    metrics["active_stats"] = active_stats
    return metrics, rel, sign_ok


def main():
    parser = argparse.ArgumentParser(
        description="Element-level coupling audit: compare element adjoint sensitivities to element FD sensitivities."
    )
    parser.add_argument("--nx", type=int, default=6)
    parser.add_argument("--ny", type=int, default=2)
    parser.add_argument("--width", type=float, default=15.0)
    parser.add_argument("--height", type=float, default=3.0)
    parser.add_argument("--nsteps", type=int, default=220)
    parser.add_argument("--nsub-steps", type=int, default=1)
    parser.add_argument("--dt", type=float, default=2.5e-3)
    parser.add_argument("--impact-velocity", type=float, default=0.8)
    parser.add_argument("--impact-window-width", type=float, default=4.0)
    parser.add_argument("--n-sensors", type=int, default=3)
    parser.add_argument("--sensor-distribution-width", type=float, default=10.0)
    parser.add_argument("--tau", type=float, default=1.0)
    parser.add_argument("--h", type=float, default=1.0e-2)
    parser.add_argument("--integrator-type", type=str, default="fixed_visco")
    parser.add_argument("--overflow-limit", type=int, default=200)
    parser.add_argument("--mat-overflow-limit", type=int, default=200)
    parser.add_argument("--mass-damping-factor", type=float, default=0.0)
    parser.add_argument("--density", type=float, default=1.0)
    parser.add_argument("--youngs-modulus", type=float, default=4.0)
    parser.add_argument("--poissons-ratio", type=float, default=0.25)
    parser.add_argument("--shear-modulus-maxwell", type=float, default=2.0)
    parser.add_argument(
        "--output-dir",
        type=str,
        default="examples/dissipative_wave/inverse_case_outputs/case3_adjoint_diag",
    )
    args = parser.parse_args()

    problem = make_structured_quad_problem(
        nx=args.nx,
        ny=args.ny,
        width=args.width,
        height=args.height,
        impact_velocity=args.impact_velocity,
        impact_window_width=args.impact_window_width,
        nsensors=args.n_sensors,
        sensor_distribution_width=args.sensor_distribution_width,
    )
    nelem = args.nx * args.ny
    base_coeffs = np.ones(nelem, dtype=np.double)
    true_coeffs = choose_true_coeffs(base_coeffs, args.nx, args.ny)

    observed = run_forward_or_adjoint(problem, true_coeffs, args.tau, args, observed_history=None, compute_gradients=False)
    observed_history = observed["sensor_history"]

    audit = run_forward_or_adjoint(problem, base_coeffs, args.tau, args, observed_history=observed_history, compute_gradients=True)
    g_adj = audit["elem_grad"]
    g_fd = fd_element_gradient(problem, base_coeffs, args.tau, args, observed_history, args.h)

    metrics, rel, sign_ok = summarize(g_adj, g_fd)

    top_idx = np.argsort(-np.abs(g_adj - g_fd))[: min(8, nelem)]
    top_rows = []
    for idx in top_idx:
        top_rows.append(
            {
                "element_id": int(idx),
                "adjoint": float(g_adj[idx]),
                "fd": float(g_fd[idx]),
                "abs_diff": float(abs(g_adj[idx] - g_fd[idx])),
                "rel_error": float(rel[idx]),
                "sign_match": bool(sign_ok[idx]),
            }
        )

    payload = {
        "config": {
            "nx": int(args.nx),
            "ny": int(args.ny),
            "width": float(args.width),
            "height": float(args.height),
            "nsteps": int(args.nsteps),
            "dt": float(args.dt),
            "tau": float(args.tau),
            "h": float(args.h),
            "integrator": args.integrator_type,
            "overflow_limit": int(args.overflow_limit),
            "mat_overflow_limit": int(args.mat_overflow_limit),
        },
        "metrics": metrics,
        "top_abs_difference_rows": top_rows,
    }

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    json_path = out_dir / "element_coupling_audit.json"
    csv_path = out_dir / "element_coupling_audit_details.csv"

    with json_path.open("w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["element_id", "adjoint", "fd", "abs_diff", "rel_error", "sign_match"])
        for i in range(nelem):
            writer.writerow(
                [
                    int(i),
                    f"{float(g_adj[i]):.12e}",
                    f"{float(g_fd[i]):.12e}",
                    f"{float(abs(g_adj[i] - g_fd[i])):.12e}",
                    f"{float(rel[i]):.12e}",
                    int(sign_ok[i]),
                ]
            )

    print(f"Wrote: {json_path}")
    print(f"Wrote: {csv_path}")
    print("Element coupling audit metrics:")
    for k, v in metrics.items():
        if isinstance(v, dict):
            print(f"  {k}:")
            for kk, vv in v.items():
                print(f"    {kk}: {vv}")
        else:
            print(f"  {k}: {v:.12e}")


if __name__ == "__main__":
    main()
