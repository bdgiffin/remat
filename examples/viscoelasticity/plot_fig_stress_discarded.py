"""Three-column fig:stress-backward-comparison of Paper 1 (step and constant-rate histories).

Columns: double-precision floating point, the same fixed-point lattice with the
Euclidean remainders discarded, and the bit-reversible fixed-point update.

The floating-point and bit-reversible columns are REMAT runs
(plot_fig_stress.run_truss_relaxation). REMAT has no mode that discards the
remainders, so that column comes from an exact integer replica of the REMAT
fixed-point update (Rational: q = floor(sqrt(2^31 - 1)), p = trunc(A q); squeezeE;
Fixed<10,-6>), driven by the strains REMAT computes. The ancillary is set to zero
before every squeeze map in both directions, so each step keeps only the Euclidean
quotient. Before it is used, the replica is checked against REMAT: with the
ancillary kept and flushed every N_max steps, it must reproduce the REMAT
bit-reversible mantissas at every time level of both passes.

The primary mantissa of the discarded update is an unbounded Python integer; its
backward curve is drawn only while it fits the 32-bit mantissa of the
implementation.

The panels are 4.0 x 3.5 in scaled by 0.31/0.45, the ratio of the panel widths in
the three-column figure and in the two-column stress figures, so their text has the
same size in the paper.

Usage
  python plot_fig_stress_discarded.py
      writes stress_backward_{float,discarded,fixed}__{step,constant_rate}.{pdf,png}
  python plot_fig_stress_discarded.py --duration 15 6 --outdir fig5_preview
      previews, with the duration appended to each file name
"""

import argparse
import contextlib
import math
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(THIS_DIR))
os.chdir(THIS_DIR)

import plot_fig_stress as ps
from plot_fig_flt_fxd_dt_convergence import fixed_flush_interval

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)

TAU = 0.3
DT = 1.0e-3
E = 1.0
Q = math.isqrt(2**31 - 1)        # Rational(Real): denominator = sqrt(INT32_MAX) cast to int32
DECIMAL = math.pow(10.0, -6.0)   # Fixed<10,-6>::decimal
INV_DECIMAL = math.pow(10.0, 6.0)
INT32_MAX = 2**31 - 1
HISTORIES = [("right_node_step", 0.1, "step", (-0.02, 0.12)),
             ("right_node_constant_rate", 0.1, "constant_rate", (0.0, 0.04))]
FORWARD_COLOR = "#2b738eff"
PANEL_SIZE = (4.0 * 0.31 / 0.45, 3.5 * 0.31 / 0.45)
BACKWARD_COLOR = "#f9826bff"


@contextlib.contextmanager
def quiet():
    sys.stdout.flush()
    fd = os.dup(1)
    with open(os.devnull, "w") as devnull:
        os.dup2(devnull.fileno(), 1)
        try:
            yield
        finally:
            sys.stdout.flush()
            os.dup2(fd, 1)
            os.close(fd)


def squeeze(first, second, p, q):
    """squeezeE of src/arithmetic.h; Python divmod is Euclidean for p, q > 0."""
    quotient1, remainder1 = divmod(second, p)
    first = first * p + remainder1
    quotient2, remainder2 = divmod(first, q)
    return quotient2, quotient1 * q + remainder2


def b_mantissa(strain, one_minus_pq):
    """Fixed<10,-6>(strain*(1 - A_rat)): truncation toward zero of the scaled product."""
    return int((strain * one_minus_pq) * INV_DECIMAL)


def replica(fwd_strain, bwd_strain, flush, keep_remainders):
    """Primary mantissas of the replica at levels 0..N, forward and backward.

    fwd_strain[n-1] is the strain REMAT used at forward step n (level n);
    bwd_strain[j] is the strain REMAT used at backward step j (level N-1-j).
    """
    A = math.exp(-DT / TAU)
    p = int(A * Q)
    one_minus_pq = 1 - p / Q
    N = len(fwd_strain)
    # REMAT's initialize() calls the material update once with dt = 0, which leaves the
    # state unchanged but advances the flush counter to 1
    x, xs, prev, counter, records = 0, 0, 0.0, 1, []
    fwd = [0]
    for n in range(N):
        if not keep_remainders:
            xs = 0
        x, xs = squeeze(x, xs, p, Q)
        prev = fwd_strain[n]
        x += b_mantissa(prev, one_minus_pq)
        fwd.append(x)
        counter += 1
        if keep_remainders and counter == flush:
            records.append(xs)
            xs, counter = 0, 0
    bwd = [None] * (N + 1)
    bwd[N] = x
    for j in range(N):
        if keep_remainders and counter == 0 and records:
            counter, xs = flush, records.pop()
        if not keep_remainders:
            xs = 0
        x -= b_mantissa(prev, one_minus_pq)
        prev = bwd_strain[j]
        x, xs = squeeze(x, xs, Q, p)
        counter -= 1
        bwd[N - 1 - j] = x
    return np.array(fwd, dtype=object), np.array(bwd, dtype=object)


def remat_run(mode, bc, eps0, n_steps, flush):
    with quiet():
        r = ps.run_truss_relaxation(
            dt=DT, Nsteps=n_steps, Nsub_steps=1, epsilon0=eps0, bc_name=bc, relaxation_time=TAU,
            overflow_limit=flush, include_backward=True,
            record_states=("axial_stress", "axial_strain", "viscous_strain"),
            set_integrator_type=mode,
        )
    h = r["state_history"]
    return r["forward_time"], {k: (h[k]["forward"], h[k]["backward"]) for k in h}


def draw(ax, t, fwd, bwd, ylim, duration, legend_loc="upper right"):
    """Same layout as plot_fig_stress.py: forward solid, backward dashed."""
    ax.plot(t, fwd, linewidth=1.6, label="forward", color=FORWARD_COLOR)
    ax.plot(t[:-1], bwd, linewidth=1.6, dashes=(4, 3.5), label="backward", color=BACKWARD_COLOR)
    ax.set_xlabel("time (s)", fontsize="large")
    ax.set_ylabel(r"axial stress ($\sigma_{xx}$)", fontsize="large")
    legend = ax.legend(loc=legend_loc, fontsize="small")
    legend.get_texts()[0].set_color(FORWARD_COLOR)
    legend.get_texts()[1].set_color(BACKWARD_COLOR)
    legend.get_lines()[1].set_linestyle("dashed")
    ax.set_xlim(0, duration)
    ax.set_ylim(*ylim)


def departure_time(t, fwd, bwd, tol=1.0e-3):
    """Latest time at which |backward - forward| exceeds tol (backward runs from T to 0)."""
    diff = np.abs(np.asarray(bwd, dtype=float) - np.asarray(fwd[:-1], dtype=float))
    bad = np.where(~(diff <= tol))[0]
    return float(t[bad[-1]]) if bad.size else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--duration", type=float, nargs="+", default=[15.0])
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()
    outdir = THIS_DIR / args.outdir
    outdir.mkdir(exist_ok=True)
    flush = int(fixed_flush_interval(DT, TAU))

    for duration in args.duration:
        n_steps = int(round(duration / DT))
        tag = f"__T{duration:g}" if len(args.duration) > 1 else ""
        for bc, eps0, name, ylim in HISTORIES:
            t, fx = remat_run(b"fixed_truss_visco", bc, eps0, n_steps, flush)
            _, fl = remat_run(b"float_truss_visco", bc, eps0, n_steps, flush)
            fwd_strain, bwd_strain = fx["axial_strain"]
            # strain at levels 1..N (forward) and N-1..0 (backward)
            fwd_eps = np.concatenate([[0.0], fwd_strain])
            bwd_eps = np.empty(n_steps + 1)
            bwd_eps[n_steps] = fwd_strain[-1]
            bwd_eps[:n_steps] = bwd_strain[::-1]

            # replica check: with the ancillary kept it must equal REMAT bit for bit
            kf, kb = replica(fwd_strain, bwd_strain, flush, keep_remainders=True)
            remat_fwd = np.rint(np.asarray(fx["viscous_strain"][0]) / DECIMAL).astype(np.int64)
            remat_bwd = np.rint(np.asarray(fx["viscous_strain"][1][::-1]) / DECIMAL).astype(np.int64)
            assert np.array_equal(kf[1:].astype(np.int64), remat_fwd), "replica forward differs from REMAT"
            assert np.array_equal(kb[:-1].astype(np.int64), remat_bwd), "replica backward differs from REMAT"

            # remainders discarded: same lattice, quotient only, in both passes
            df, db = replica(fwd_strain, bwd_strain, flush, keep_remainders=False)
            in_range = np.array([abs(int(v)) <= INT32_MAX for v in db])
            first_over = np.where(~in_range)[0]
            d_fwd = E * (fwd_eps - np.array([float(v) * DECIMAL for v in df]))
            d_bwd = E * (bwd_eps - np.array([float(v) * DECIMAL for v in db]))
            if first_over.size:
                d_bwd[: first_over[-1] + 1] = np.nan  # backward runs from level N down to 0

            columns = [
                ("float", fl["axial_stress"][0], fl["axial_stress"][1][::-1][1:]),
                ("discarded", d_fwd[1:], d_bwd[1:-1]),
                ("fixed", fx["axial_stress"][0], fx["axial_stress"][1][::-1][1:]),
            ]
            for label, fwd, bwd in columns:
                fig, ax = plt.subplots(figsize=PANEL_SIZE)
                # the discarded backward curve leaves near the right edge, so its legend moves
                draw(ax, t, fwd, bwd, ylim, duration, "center" if label == "discarded" else "upper right")
                fig.tight_layout()
                stem = outdir / f"stress_backward_{label}__{name}{tag}"
                fig.savefig(f"{stem}.pdf", dpi=200)
                fig.savefig(f"{stem}.png", dpi=200)
                plt.close(fig)
                dep = departure_time(t, fwd, bwd)
                print(f"T={duration:g} {name:13s} {label:9s} backward departs (|diff|>1e-3) at t = "
                      f"{'never' if dep is None else f'{dep:.2f} s'}; "
                      f"|diff| at t=dt: {abs(float(bwd[0]) - float(fwd[0])) if np.isfinite(float(bwd[0])) else float('nan'):.2e}")
            over_t = t[first_over[-1] - 1] if first_over.size and first_over[-1] >= 1 else None
            print(f"T={duration:g} {name:13s} replica matches REMAT bit for bit; discarded 32-bit overflow below t = "
                  f"{'none' if over_t is None else f'{over_t:.2f} s'}")


if __name__ == "__main__":
    main()
