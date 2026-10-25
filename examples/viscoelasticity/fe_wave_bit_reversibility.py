"""Forward--backward reconstruction of a viscoelastic finite-element wave problem.

Uses the layered-medium problem of examples/dissipative_wave (true material
fields, one surface impact) and runs it forward and then backward with
  fixed_visco : fixed-point leapfrog + bit-reversible Maxwell update (Algorithm 3)
  float_visco : the same model in double precision with the explicit inverse.

For every time step the backward state is compared with the forward state:
  - a hash of the full persistent state (nodal displacement, velocity and their
    ancillary pairs, and all 29 state values at every integration point),
  - entrywise comparison of nodal u, v and integration-point viscous strain and
    its ancillary variable.
The element-average output slot is excluded because it is a diagnostic copy that
is recomputed at each update and is not part of the persistent state.

The library is compiled here with floating-point contraction disabled
(-ffp-contract=off). With contraction enabled, the compiler may fuse the stress
evaluation into multiply--add instructions differently in the forward and
backward code paths, so the stress (which feeds the nodal forces) can differ by
one unit in the last place between the two passes. A check with the default
installed build is included to document this.

Usage
  python fe_wave_bit_reversibility.py            # runs both modes, writes figure/json
  python fe_wave_bit_reversibility.py --plot-only

Outputs
  fe_wave_bit_reversibility.json
  fe_wave_bit_reversibility.{pdf,png}
"""

import argparse
import hashlib
import json
import math
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
REPO_ROOT = THIS_DIR.parent.parent
WAVE_DIR = THIS_DIR.parent / "dissipative_wave"
sys.path.insert(0, str(WAVE_DIR))
if "--package" in sys.argv:
    # Select the REMAT build before dissipative_wave_inverse imports REMAT
    sys.path.insert(0, sys.argv[sys.argv.index("--package") + 1])

import dissipative_wave_inverse as inverse  # noqa: E402  (also imports REMAT)

STRICT_FLAGS = ["-std=c++17", "-fPIC", "-O3", "-march=native", "-ffp-contract=off"]

REMAT = inverse.REMAT

N_STEPS = 2500            # 10 s with dt = 4e-3
IMPACT_CENTER = 6.0
INT32_MAX = 2**31 - 1
M_MAX = 2**30             # largest admissible ancillary mantissa, half of the int32 range
NUM_GP = 4                # 2x2 Gauss quadrature per element
NUM_MAT_STATE = 29        # ViscoElasticity::num_state_vars()
VISCOUS = slice(9, 12)    # viscous_strain_{xx,yy,xy}
ANCILLARY = slice(12, 15) # viscous_strain_{xx,yy,xy}_dual
FIXED_H = 1.0e-6          # lattice spacing of Fixed_E (Fixed<10,-6>)
NODE_COLS = [0, 1, 3, 4, 9, 10, 12, 13]  # u_x, u_y, v_x, v_y and their ancillary pairs

OUT_JSON = THIS_DIR / "fe_wave_bit_reversibility.json"


def flush_interval(dt, tau_min):
    """Largest flush interval allowed by eq:flush-time for the shortest relaxation time.

    The ancillary mantissa obeys m_{n+1} <= (q/p) m_n + (q - 1), so from m_0 = 0 it stays
    below M_MAX for N < ln(1 + (q/p - 1) M_MAX / (q - 1)) / ln(q/p).
    """
    q = int(math.sqrt(INT32_MAX))
    p = int(math.exp(-dt / tau_min) * q)
    ratio = q / p
    bound = math.log(1.0 + (ratio - 1.0) * M_MAX / (q - 1)) / math.log(ratio)
    return math.ceil(bound) - 1


def setup(mode, mat_limit, v_limit):
    problem = inverse.make_experiment(IMPACT_CENTER)
    labels = inverse.make_labels(problem)
    stiffness = inverse.region_stiffness_values()[labels]
    tau = inverse.region_tau_values()[labels]
    inverse.INTEGRATOR_TYPE = mode
    inverse.OVERFLOW_LIMIT = float(v_limit)
    inverse.MAT_OVERFLOW_LIMIT = float(mat_limit)
    with inverse.suppress_solver_output():
        inverse.configure_run(problem, stiffness, tau)
    return problem, tau


def capture(n_nodes, n_node_fields, n_raw):
    node = np.zeros((n_nodes, n_node_fields), dtype=np.double)
    REMAT.API.get_fields(b"node", node)
    raw = np.zeros(n_raw, dtype=np.double)
    REMAT.API.get_raw_element_state(raw)
    gp = raw.reshape(-1, NUM_GP + 1, NUM_MAT_STATE)[:, 1:, :]   # drop element-average slot
    node_key = np.ascontiguousarray(node[:, NODE_COLS])
    digest = hashlib.blake2b(node_key.tobytes() + np.ascontiguousarray(gp).tobytes(), digest_size=16).hexdigest()
    return node_key, np.ascontiguousarray(gp[:, :, VISCOUS]), np.ascontiguousarray(gp[:, :, ANCILLARY]), digest


def worker_verify(mode, n_steps, mat_limit, v_limit, out_npz):
    problem, tau = setup(mode, mat_limit, v_limit)
    dt = inverse.DT
    n_nodes = REMAT.API.get_num_entities(b"node")
    n_node_fields = REMAT.API.get_num_fields(b"node")
    n_raw = REMAT.API.get_num_raw_element_state()

    fwd_node, fwd_q, fwd_anc, fwd_hash = [], [], [], []
    record_bytes = []
    node_key, q, anc, h = capture(n_nodes, n_node_fields, n_raw)
    fwd_node.append(node_key); fwd_q.append(q); fwd_anc.append(anc); fwd_hash.append(h)
    max_anc_pre_flush = []
    with inverse.suppress_solver_output():
        for _ in range(n_steps):
            REMAT.API.update_state(dt, 1, REMAT.PASS_FORWARD)
            node_key, q, anc, h = capture(n_nodes, n_node_fields, n_raw)
            fwd_node.append(node_key); fwd_q.append(q); fwd_anc.append(anc); fwd_hash.append(h)
            record_bytes.append(REMAT.API.get_overflow_record_bytes())
    time_history_bytes = REMAT.API.get_time_history_bytes()
    peak_record_bytes = REMAT.API.get_overflow_record_bytes()

    hash_mismatch = np.zeros(n_steps + 1, dtype=np.int64)
    entry_mismatch = np.zeros(n_steps + 1, dtype=np.int64)
    err_v = np.zeros(n_steps + 1)
    err_q = np.zeros(n_steps + 1)
    v_scale = max(float(np.max(np.abs(a[:, 2:4]))) for a in fwd_node)
    q_scale = max(float(np.max(np.abs(a))) for a in fwd_q)

    def compare(n, node_key, q, anc, h):
        hash_mismatch[n] = int(h != fwd_hash[n])
        entry_mismatch[n] = (
            int(np.count_nonzero(node_key != fwd_node[n]))
            + int(np.count_nonzero(q != fwd_q[n]))
            + int(np.count_nonzero(anc != fwd_anc[n]))
        )
        err_v[n] = float(np.max(np.abs(node_key[:, 2:4] - fwd_node[n][:, 2:4]))) / v_scale
        err_q[n] = float(np.max(np.abs(q - fwd_q[n]))) / q_scale

    compare(n_steps, *capture(n_nodes, n_node_fields, n_raw))
    with inverse.suppress_solver_output():
        for n in range(n_steps - 1, -1, -1):
            REMAT.API.update_state(dt, 1, REMAT.PASS_BACKWARD)
            compare(n, *capture(n_nodes, n_node_fields, n_raw))

    anc_arr = np.stack(fwd_anc)                    # (N+1, Ne, 4, 3)
    if mode.startswith("fixed"):
        anc_mantissa = np.max(np.abs(anc_arr), axis=(1, 2, 3))   # stored as mantissas
    else:
        anc_mantissa = np.zeros(n_steps + 1)
    v_anc_max = float(np.max(np.abs(np.stack(fwd_node)[:, :, 4:8])))

    np.savez(
        out_npz,
        time=np.arange(n_steps + 1) * dt,
        hash_mismatch=hash_mismatch,
        entry_mismatch=entry_mismatch,
        err_v=err_v,
        err_q=err_q,
        anc_mantissa=anc_mantissa,
        record_bytes=np.asarray([0] + record_bytes),
        tau_elem=tau,
        elem_centers=problem["elem_centers"],
        impact_center=IMPACT_CENTER,
    )
    meta = {
        "mode": mode,
        "n_steps": n_steps,
        "dt": dt,
        "n_nodes": n_nodes,
        "n_elements": int(n_raw // ((NUM_GP + 1) * NUM_MAT_STATE)),
        "n_integration_points": int(n_raw // ((NUM_GP + 1) * NUM_MAT_STATE)) * NUM_GP,
        "mat_flush_interval": mat_limit,
        "velocity_flush_interval": v_limit,
        "steps_with_hash_mismatch": int(np.count_nonzero(hash_mismatch)),
        "total_entry_mismatches": int(np.sum(entry_mismatch)),
        "entries_compared_per_step": int(fwd_node[0].size + fwd_q[0].size + fwd_anc[0].size),
        "max_rel_error_velocity": float(np.max(err_v)),
        "max_rel_error_viscous_strain": float(np.max(err_q)),
        "max_ancillary_mantissa": float(np.max(anc_mantissa)),
        "max_abs_velocity_ancillary": v_anc_max,
        "peak_overflow_record_bytes": int(peak_record_bytes),
        "time_history_bytes": int(time_history_bytes),
    }
    print(json.dumps(meta))


def timed_step(dt, phase):
    """Advance one step and return the time spent in the solver call."""
    t0 = time.perf_counter()
    REMAT.API.update_state(dt, 1, phase)
    return time.perf_counter() - t0


def worker_timing(mode, n_steps, mat_limit, v_limit, repeats):
    """Solver time of the forward and backward passes (bit-reversible strategy)."""
    t_fwd, t_bwd = [], []
    for _ in range(repeats):
        setup(mode, mat_limit, v_limit)
        dt = inverse.DT
        with inverse.suppress_solver_output():
            t_fwd.append(sum(timed_step(dt, REMAT.PASS_FORWARD) for _ in range(n_steps)))
            t_bwd.append(sum(timed_step(dt, REMAT.PASS_BACKWARD) for _ in range(n_steps)))
    print(json.dumps({"mode": mode, "t_forward": float(np.median(t_fwd)), "t_backward": float(np.median(t_bwd))}))


# ----------------------------------------------------------------- checkpointing
# A checkpoint holds the restart state of the double-precision calculation: time,
# step and overflow counters, nodal displacement, velocity, and acceleration, and
# the viscous strain at every integration point. The remaining entries of the
# solver snapshot are material constants, which are taken from the initial
# snapshot, or quantities recomputed at the next update.

def snapshot_layout():
    n_dofs = 2 * REMAT.API.get_num_entities(b"node")
    n_raw = REMAT.API.get_num_raw_element_state()
    dof0 = 3
    state0 = dof0 + 7 * n_dofs
    gp = np.arange(n_raw).reshape(-1, NUM_GP + 1, NUM_MAT_STATE)[:, 1:, VISCOUS].ravel()
    dof = dof0 + 7 * np.arange(n_dofs)
    return {
        "size": REMAT.API.get_snapshot_size(),
        "head": np.arange(3),
        "u": dof + 0, "v": dof + 2, "a": dof + 4,
        "viscous": state0 + gp,
    }


def save_snapshot(layout):
    snap = np.zeros(layout["size"])
    REMAT.API.save_snapshot(snap)
    return snap


def restart_state(snap, layout):
    return np.concatenate([snap[layout[k]] for k in ("head", "u", "v", "a", "viscous")])


def level_state(snap, layout):
    """State recovered at each time level: nodal u and v and integration-point viscous strain."""
    return np.concatenate([snap[layout[k]] for k in ("u", "v", "viscous")])


def load_restart(template, restart, layout):
    snap = template.copy()
    k = 0
    for key in ("head", "u", "v", "a", "viscous"):
        n = layout[key].size
        snap[layout[key]] = restart[k:k + n]
        k += n
    REMAT.API.load_snapshot(snap)


def checkpoint_levels(n_steps):
    k = math.ceil(math.sqrt(n_steps))
    starts = list(range(0, n_steps, k))
    return k, [(n0, min(n0 + k, n_steps)) for n0 in starts]


def worker_checkpoint(n_steps, mat_limit, v_limit, repeats):
    """Uniform sqrt(N) checkpointing of the double-precision calculation.

    The verification pass stores every level, restarts each segment from its
    checkpoint restart state alone, and compares every recomputed level with the
    stored forward level. The timed passes record the solver time of the forward
    pass and of the recomputation; copies of states to and from storage are not
    timed for any strategy.
    """
    dt = inverse.DT
    k, segments = checkpoint_levels(n_steps)

    # Verification (floating point needs no ancillary flushes)
    setup("float_visco", v_limit, v_limit)
    layout = snapshot_layout()
    with inverse.suppress_solver_output():
        template = save_snapshot(layout)
        reference = [restart_state(template, layout)]
        for _ in range(n_steps):
            REMAT.API.update_state(dt, 1, REMAT.PASS_FORWARD)
            reference.append(restart_state(save_snapshot(layout), layout))
        mismatches, recomputed = 0, 0
        for n0, n1 in reversed(segments):
            load_restart(template, reference[n0], layout)
            mismatches += int(np.count_nonzero(restart_state(save_snapshot(layout), layout) != reference[n0]))
            for n in range(n0 + 1, n1):
                REMAT.API.update_state(dt, 1, REMAT.PASS_FORWARD)
                recomputed += 1
                mismatches += int(np.count_nonzero(restart_state(save_snapshot(layout), layout) != reference[n]))
    restart_values = int(reference[0].size)
    level_values = int(level_state(template, layout).size)
    del reference

    # Timing
    t_fwd, t_rec, t_restore = [], [], []
    for _ in range(repeats):
        setup("float_visco", v_limit, v_limit)
        with inverse.suppress_solver_output():
            template = save_snapshot(layout)
            checkpoints = {}
            tf = 0.0
            for n in range(n_steps):
                if n % k == 0:
                    checkpoints[n] = restart_state(save_snapshot(layout), layout)
                tf += timed_step(dt, REMAT.PASS_FORWARD)
            tr, tl = 0.0, 0.0
            for n0, n1 in reversed(segments):
                t0 = time.perf_counter()
                load_restart(template, checkpoints[n0], layout)
                tl += time.perf_counter() - t0
                for _ in range(n0 + 1, n1):
                    tr += timed_step(dt, REMAT.PASS_FORWARD)
        t_fwd.append(tf); t_rec.append(tr); t_restore.append(tl)

    print(json.dumps({
        "segment_length": k,
        "n_checkpoints": len(segments),
        "recomputed_steps": recomputed,
        "restart_mismatches": mismatches,
        "restart_values": restart_values,
        "level_values": level_values,
        "t_forward": float(np.median(t_fwd)),
        "t_recompute": float(np.median(t_rec)),
        "t_restore": float(np.median(t_restore)),
    }))


def build_strict_package(dest):
    """Compile REMAT with floating-point contraction disabled into dest."""
    dest.mkdir(parents=True, exist_ok=True)
    shutil.copy(REPO_ROOT / "src" / "REMAT.py", dest / "REMAT.py")
    subprocess.run(
        ["g++", "-DREMAT_EXPORTS", *STRICT_FLAGS, "-shared", "-o", str(dest / "libREMAT.so"),
         str(REPO_ROOT / "src" / "REMAT_API.cpp")],
        check=True,
    )
    return dest


def spawn(args, package=None):
    extra = ["--package", str(package)] if package else []
    out = subprocess.run(
        ["taskset", "-c", "2", sys.executable, str(Path(__file__).resolve())] + args + extra,
        check=True, capture_output=True, text=True, cwd=THIS_DIR,
    ).stdout.strip().splitlines()
    return json.loads(out[-1])


def storage_summary(meta, timing_fixed, ckpt):
    """Bytes and solver times for backward access to every forward state, for each strategy."""
    n = meta["n_steps"]
    n_ip = meta["n_integration_points"]
    level_bytes = 8 * ckpt["level_values"]        # nodal u and v, integration-point viscous strain
    restart_bytes = 8 * ckpt["restart_values"]    # level state plus nodal acceleration and counters
    full_time = ckpt["t_forward"]
    ckpt_time = ckpt["t_forward"] + ckpt["t_recompute"] + ckpt["t_restore"]
    rev_bytes = meta["peak_overflow_record_bytes"] + 8 * 3 * n_ip + meta["time_history_bytes"]
    rev_time = timing_fixed["t_forward"] + timing_fixed["t_backward"]
    return {
        "level_bytes": level_bytes,
        "restart_bytes": restart_bytes,
        "full_storage_bytes": level_bytes * (n + 1),
        "checkpoint_bytes": restart_bytes * ckpt["n_checkpoints"] + level_bytes * ckpt["segment_length"],
        "checkpoint_recomputed_steps": ckpt["recomputed_steps"],
        "reversible_record_bytes": meta["peak_overflow_record_bytes"],
        "reversible_ancillary_bytes": 8 * 3 * n_ip,
        "reversible_time_history_bytes": meta["time_history_bytes"],
        "reversible_bytes": rev_bytes,
        "full_storage_time": full_time,
        "checkpoint_forward_time": ckpt["t_forward"],
        "checkpoint_backward_time": ckpt["t_recompute"] + ckpt["t_restore"],
        "checkpoint_time": ckpt_time,
        "reversible_forward_time": timing_fixed["t_forward"],
        "reversible_backward_time": timing_fixed["t_backward"],
        "reversible_time": rev_time,
        "checkpoint_time_ratio": ckpt_time / full_time,
        "reversible_time_ratio": rev_time / full_time,
        "reversible_over_checkpoint_time": rev_time / ckpt_time,
    }


def plot(results):
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
            "mathtext.fontset": "cm",
        }
    )
    fx = np.load(THIS_DIR / "fe_wave_fixed_visco.npz")
    fl = np.load(THIS_DIR / "fe_wave_float_visco.npz")
    meta = results["fixed_visco"]
    t = fx["time"]
    t_end = meta["n_steps"] * meta["dt"]

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    fig = plt.figure(figsize=(7.6, 5.0))
    ax0 = fig.add_subplot(2, 1, 1)
    nx, nz = inverse.NX, inverse.NZ
    tau_grid = np.asarray(fx["tau_elem"]).reshape(nz, nx)
    im = ax0.imshow(tau_grid, origin="lower", extent=[0, inverse.WIDTH, 0, inverse.DEPTH],
                    cmap="magma", aspect="equal")
    half = 0.5 * inverse.IMPACT_WINDOW_WIDTH
    ax0.plot([IMPACT_CENTER - half, IMPACT_CENTER + half], [inverse.DEPTH, inverse.DEPTH],
             color="#2b738eff", linewidth=4.0, solid_capstyle="butt", clip_on=False)
    ax0.text(IMPACT_CENTER + half + 0.15, inverse.DEPTH + 0.12, "impact", fontsize="medium",
             color="#2b738eff", va="bottom")
    ax0.set_xlabel("$x$", fontsize="large")
    ax0.set_ylabel("$y$", fontsize="large")
    cax = make_axes_locatable(ax0).append_axes("right", size="1.5%", pad=0.08)
    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r"$\tau$", fontsize="large")
    ax0.set_title("(a) relaxation-time field, bottom edge fixed", fontsize="medium", pad=16)

    ax1 = fig.add_subplot(2, 2, 3)
    ax1.semilogy(t, np.maximum(fx["anc_mantissa"], 1.0), color="#2b738eff", linewidth=1.6)
    ax1.axhline(INT32_MAX, color="k", linestyle="--", linewidth=1.0)
    ax1.axhline(M_MAX, color="k", linestyle=":", linewidth=1.0)
    ax1.text(0.02 * t_end, INT32_MAX * 1.6, "int32 limit", fontsize="small", va="bottom")
    ax1.text(0.02 * t_end, M_MAX * 0.6, r"sizing level $m_{\max}$", fontsize="small", va="top")
    ax1.set_ylim(top=INT32_MAX * 12)
    ax1.set_xlim(0.0, t_end)
    ax1.set_xlabel("time", fontsize="large")
    ax1.set_ylabel("max ancillary mantissa", fontsize="large")
    ax1.set_title("(b) ancillary growth with flushes", fontsize="medium")

    ax2 = fig.add_subplot(2, 2, 4)
    ax2.semilogy(t, np.maximum(fl["err_v"], 1e-18), color="#f9826bff", linewidth=1.6, label="velocity")
    ax2.semilogy(t, np.maximum(fl["err_q"], 1e-18), color="#f9826bff", linewidth=1.6, linestyle="--",
                 label="viscous strain")
    ax2.set_xlim(0.0, t_end)
    ax2.set_xlabel("time", fontsize="large")
    ax2.set_ylabel("relative backward error", fontsize="large")
    ax2.set_title("(c) floating point (fixed point: 0)", fontsize="medium")
    ax2.legend(fontsize="small", frameon=False)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(THIS_DIR / f"fe_wave_bit_reversibility.{ext}", dpi=200)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--worker", choices=["verify", "timing", "checkpoint"])
    parser.add_argument("--mode", default="fixed_visco")
    parser.add_argument("--steps", type=int, default=N_STEPS)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--plot-only", action="store_true")
    parser.add_argument("--package", default=None)
    parser.add_argument("--tag", default="")
    args = parser.parse_args()

    tau_min = float(np.min(inverse.region_tau_values()))
    mat_limit = flush_interval(inverse.DT, tau_min)
    v_limit = args.steps + 1   # no mass damping: the velocity ancillary stays zero (checked)

    if args.worker == "verify":
        worker_verify(args.mode, args.steps, mat_limit, v_limit,
                      THIS_DIR / f"fe_wave_{args.mode}{args.tag}.npz")
        return
    if args.worker == "timing":
        worker_timing(args.mode, args.steps, mat_limit, v_limit, args.repeats)
        return
    if args.worker == "checkpoint":
        worker_checkpoint(args.steps, mat_limit, v_limit, args.repeats)
        return

    if args.plot_only:
        results = json.loads(OUT_JSON.read_text())
    else:
        results = {"compile_flags": " ".join(STRICT_FLAGS)}
        with tempfile.TemporaryDirectory() as tmpdir:
            pkg = build_strict_package(Path(tmpdir) / "package")
            for mode in ["fixed_visco", "float_visco"]:
                results[mode] = spawn(["--worker", "verify", "--mode", mode, "--steps", str(args.steps)], pkg)
            timing = {"fixed_visco": spawn(["--worker", "timing", "--mode", "fixed_visco", "--steps",
                                            str(args.steps), "--repeats", str(args.repeats)], pkg)}
            checkpointing = spawn(["--worker", "checkpoint", "--steps", str(args.steps),
                                   "--repeats", str(args.repeats)], pkg)
        # Same check with the default installed build (floating-point contraction enabled)
        results["fixed_visco_default_build"] = spawn(
            ["--worker", "verify", "--mode", "fixed_visco", "--steps", str(args.steps), "--tag", "_default_build"])
        results["timing"] = timing
        results["checkpointing"] = checkpointing
        results["storage"] = storage_summary(results["fixed_visco"], timing["fixed_visco"], checkpointing)
        OUT_JSON.write_text(json.dumps(results, indent=2))
    print(json.dumps(results, indent=2))
    plot(results)


if __name__ == "__main__":
    main()
