"""Finite-precision contributions to the time-step study of Section 8.

The fixed-point Maxwell update has two truncation errors that grow as dt/tau
decreases:
  - the inhomogeneous term b is truncated toward zero on the lattice of spacing h,
    which lowers the viscous strain and raises the stress;
  - the numerator p = floor(A q) makes p/q slightly smaller than A, which shortens
    the effective relaxation time and lowers the stress.
To separate them, REMAT is recompiled with a 100 times smaller lattice spacing
and/or a larger common denominator q, and the forward constant-rate study of
plot_fig_flt_fxd_dt_convergence.py is repeated for each build.

Output: dt_convergence_finite_precision.json
"""

import contextlib
import json
import multiprocessing as mp
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

THIS_DIR = Path(__file__).resolve().parent
SRC_DIR = THIS_DIR.parent.parent / "src"

DEFAULT_Q = "std::sqrt(std::numeric_limits<Integer>::max())"
LARGE_Q = "(Integer(1)<<24)"
VARIANTS = {
    "baseline": {"q": DEFAULT_Q, "exponent": -6},
    "larger_q": {"q": LARGE_Q, "exponent": -6},
    "smaller_h": {"q": DEFAULT_Q, "exponent": -8},
    "both": {"q": LARGE_Q, "exponent": -8},
}
COMPILE_FLAGS = ["-std=c++17", "-fPIC", "-O3", "-march=native"]


def build_variant(dest, q_expr, exponent):
    """Copy the REMAT sources, set q and the decimal exponent of h, and compile."""
    src = dest / "src"
    pkg = dest / "package"
    src.mkdir(parents=True)
    pkg.mkdir()
    for f in list(SRC_DIR.glob("*.h")) + list(SRC_DIR.glob("*.cpp")):
        shutil.copy(f, src / f.name)
    shutil.copy(SRC_DIR / "REMAT.py", pkg / "REMAT.py")
    rational = (src / "Rational.h").read_text()
    assert rational.count(f"denominator = {DEFAULT_Q};") == 1
    (src / "Rational.h").write_text(rational.replace(f"denominator = {DEFAULT_Q};", f"denominator = {q_expr};"))
    fixed = (src / "Fixed.h").read_text()
    (src / "Fixed.h").write_text(fixed.replace("Fixed<10,-6>", f"Fixed<10,{exponent}>"))
    subprocess.run(["g++", "-DREMAT_EXPORTS", *COMPILE_FLAGS, "-shared", "-o", str(pkg / "libREMAT.so"),
                    str(src / "REMAT_API.cpp")], check=True)
    return pkg


def init_worker(pkg):
    sys.path.insert(0, str(pkg))
    sys.path.insert(1, str(THIS_DIR))
    os.chdir(THIS_DIR)


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


def run_case(args):
    variant, dt, mode = args
    import plot_fig_flt_fxd_dt_convergence as conv
    sc = conv.SCENARIO
    n_steps = max(1, int(np.round(conv.DURATION / dt)))
    fixed = mode == "fixed"
    with quiet():
        r = conv.run_truss_relaxation(
            dt=dt, Nsteps=n_steps, Nsub_steps=1, epsilon0=sc["epsilon0"], bc_name=sc["bc_name"],
            relaxation_time=sc["relaxation_time"],
            overflow_limit=conv.fixed_flush_interval(dt, sc["relaxation_time"]) if fixed else conv.FLOAT_OVERFLOW_LIMIT,
            record_states=(conv.STATE_TO_PLOT,),
            set_integrator_type=b"fixed_truss_visco" if fixed else b"float_truss_visco",
            include_backward=False,
        )
    stress = r["state_history"][conv.STATE_TO_PLOT]["forward"]
    analytical = conv.compute_analytical_stress(r["forward_time"], {**sc, "dt": dt})
    return {
        "variant": variant,
        "mode": mode,
        "dt": dt,
        "dt_over_tau": dt / sc["relaxation_time"],
        "relative_l2_error": conv.compute_relative_l2_error(analytical, stress),
        "signed_mean_error": float(np.mean(stress - analytical) / np.mean(analytical)),
    }


def main():
    sys.path.insert(0, str(THIS_DIR))
    import plot_fig_flt_fxd_dt_convergence as conv
    dts = [float(dt) for dt in conv.DT_VALUES]
    results = []
    ctx = mp.get_context("spawn")
    with tempfile.TemporaryDirectory() as tmpdir:
        for variant, cfg in VARIANTS.items():
            pkg = build_variant(Path(tmpdir) / variant, cfg["q"], cfg["exponent"])
            jobs = [(variant, dt, "fixed") for dt in dts]
            if variant == "baseline":
                jobs += [(variant, dt, "float") for dt in dts]
            with ProcessPoolExecutor(max_workers=min(len(jobs), os.cpu_count()), mp_context=ctx,
                                     initializer=init_worker, initargs=(pkg,)) as ex:
                results += list(ex.map(run_case, jobs))
    (THIS_DIR / "dt_convergence_finite_precision.json").write_text(json.dumps(results, indent=2))

    print(f"{'dt/tau':>10} " + " ".join(f"{v:>22}" for v in VARIANTS) + f" {'double':>22}")
    for dt in dts:
        row = {r["variant"]: r for r in results if r["dt"] == dt and r["mode"] == "fixed"}
        fl = next(r for r in results if r["dt"] == dt and r["mode"] == "float")
        cells = [f"{row[v]['relative_l2_error']:.3e} ({row[v]['signed_mean_error']:+.2e})" for v in VARIANTS]
        print(f"{dt / conv.SCENARIO['relaxation_time']:10.3e} " + " ".join(f"{c:>22}" for c in cells)
              + f" {fl['relative_l2_error']:.3e} ({fl['signed_mean_error']:+.2e})")


if __name__ == "__main__":
    main()
