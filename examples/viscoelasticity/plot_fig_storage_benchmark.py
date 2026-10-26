"""Storage and runtime of backward access to Maxwell history variables.

Compiles and runs storage_benchmark.cpp, which compares full state storage,
uniform sqrt(N) checkpointing, and the bit-reversible fixed-point update
(Algorithm 3) with ancillary flushes, all in the same C++ binary.

Outputs
  storage_benchmark.json            raw results (table case and N sweep)
  storage_benchmark_table.tex       LaTeX table body for the paper
  storage_vs_steps__{bytes,time}.{pdf,png}   storage (a) and wall time (b) versus number of steps

With --plot-only the figure is redrawn from storage_benchmark.json without rerunning the benchmark.
"""

import argparse
import json
import math
import subprocess
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "pdf.fonttype": 42,
    }
)

THIS_DIR = Path(__file__).resolve().parent
SRC_DIR = THIS_DIR.parent.parent / "src"

DT = 1.0e-3
TAU = 0.3
INT32_MAX = 2**31 - 1
M_MAX = 2**30  # largest admissible ancillary mantissa, half of the int32 range

TABLE_NP = 10368  # 864 elements x 4 integration points x 3 deviatoric components
TABLE_N = 15000
SWEEP_NP = 1024
SWEEP_N = [1000, 3000, 10000, 30000, 100000]
REPEATS = 5

COLORS = {"store_all": "#2b738eff", "checkpoint": "#f9826bff", "reversible": "k"}
LABELS = {
    "store_all": "full storage",
    "checkpoint": r"$\sqrt{N}$ checkpointing",
    "reversible": "bit-reversible",
}


def flush_interval(dt, tau):
    """Largest flush interval allowed by eq:flush-time.

    The ancillary mantissa obeys m_{n+1} <= (q/p) m_n + (q - 1), so from m_0 = 0 it stays
    below M_MAX for N < ln(1 + (q/p - 1) M_MAX / (q - 1)) / ln(q/p).
    """
    q = int(math.sqrt(INT32_MAX))
    p = int(math.exp(-dt / tau) * q)
    ratio = q / p
    bound = math.log(1.0 + (ratio - 1.0) * M_MAX / (q - 1)) / math.log(ratio)
    return math.ceil(bound) - 1, p, q


def build(tmpdir):
    exe = Path(tmpdir) / "storage_benchmark"
    subprocess.run(
        ["g++", "-O3", "-march=native", "-ffp-contract=off", "-std=c++17", f"-I{SRC_DIR}",
         str(THIS_DIR / "storage_benchmark.cpp"), "-o", str(exe)],
        check=True,
    )
    return exe


def run(exe, n_points, n_steps, flush, repeats):
    out = subprocess.run(
        ["taskset", "-c", "2", str(exe), str(n_points), str(n_steps), str(flush), str(repeats)],
        check=True, capture_output=True, text=True,
    ).stdout.splitlines()
    verify = out[0]
    header = out[1].split(",")
    rows = [dict(zip(header, line.split(","))) for line in out[2:]]
    result = {"verify": verify, "methods": {}}
    for name in ["store_all", "checkpoint", "reversible", "explicit_inverse"]:
        sel = [r for r in rows if r["method"] == name]
        tf = np.median([float(r["t_forward"]) for r in sel])
        tb = np.median([float(r["t_backward"]) for r in sel])
        result["methods"][name] = {
            "t_forward": tf,
            "t_backward": tb,
            "t_total": tf + tb,
            "bytes": int(sel[0]["bytes"]),
            "extra_forward_steps": int(sel[0]["extra_forward_steps"]),
            "checksum": sel[0]["checksum"],
            "max_error": float(sel[0]["max_error"]),
            "flushes": int(sel[0]["flushes"]),
            "max_ancillary": int(sel[0]["max_ancillary"]),
        }
    return result


def fmt_bytes(b):
    if b >= 1e9:
        return f"{b / 1e9:.2f}~GB"
    if b >= 1e6:
        return f"{b / 1e6:.2f}~MB"
    return f"{b / 1e3:.1f}~kB"


def write_table(res, path):
    m = res["methods"]
    ref = m["store_all"]
    lines = []
    for name, label in [("store_all", "Full state storage"),
                        ("checkpoint", r"Uniform $\sqrt{N}$ checkpointing"),
                        ("reversible", "Bit-reversible (\\autoref{alg:bit_reversible_maxwell})")]:
        r = m[name]
        lines.append(
            f"{label} & {fmt_bytes(r['bytes'])} & {r['bytes'] / ref['bytes']:.2e} & "
            f"{r['extra_forward_steps']} & {r['t_forward']:.2f} & {r['t_backward']:.2f} & "
            f"{r['t_total'] / ref['t_total']:.2f} \\\\"
        )
    path.write_text("\n".join(lines) + "\n")


def plot(sweep):
    """Storage (a) and total wall time (b) versus number of steps, one figure file per panel.

    Both panels are printed side by side at 0.45\\linewidth (6.5 in) in paper_1.tex; an 8.6 pt
    printed base size keeps superscripts at 6 pt or more (Elsevier artwork sizing).
    """
    panel_size = (3.8, 3.2)
    plt.rcParams["font.size"] = 8.6 * panel_size[0] / (0.45 * 6.5)
    n = np.asarray(SWEEP_N, dtype=float)
    panels = [
        ("bytes", lambda r: r["bytes"] / SWEEP_NP, "stored bytes per\nhistory variable"),
        ("time", lambda r: r["t_total"] / SWEEP_NP * 1e6, r"forward + backward time ($\mu$s)"),
    ]
    for tag, value, ylabel in panels:
        fig, ax = plt.subplots(figsize=panel_size)
        for name in ["store_all", "checkpoint", "reversible"]:
            y = [value(sweep[k]["methods"][name]) for k in SWEEP_N]
            ax.loglog(n, y, color=COLORS[name], linewidth=1.6, marker="o", markersize=3.5, label=LABELS[name])
        ax.set_xlabel("number of time steps $N$", fontsize="large")
        ax.set_ylabel(ylabel, fontsize="large")
        if tag == "bytes":
            ax.set_ylim(top=1.0e8)  # room for the legend above the full-storage line
            ax.legend(loc="upper left", fontsize="small", frameon=False)
        fig.tight_layout()
        for ext in ("pdf", "png"):
            fig.savefig(THIS_DIR / f"storage_vs_steps__{tag}.{ext}", dpi=200)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plot-only", action="store_true",
                        help="redraw the figure from storage_benchmark.json without rerunning the benchmark")
    args = parser.parse_args()
    if args.plot_only:
        out = json.loads((THIS_DIR / "storage_benchmark.json").read_text())
        plot({int(n): r for n, r in out["sweep"]["results"].items()})
        return

    flush, p, q = flush_interval(DT, TAU)
    with tempfile.TemporaryDirectory() as tmpdir:
        exe = build(tmpdir)
        table = run(exe, TABLE_NP, TABLE_N, flush, REPEATS)
        sweep = {n: run(exe, SWEEP_NP, n, flush, REPEATS) for n in SWEEP_N}

    out = {
        "dt": DT, "tau": TAU, "p": p, "q": q, "flush_interval": flush, "m_max": M_MAX,
        "table": {"n_points": TABLE_NP, "n_steps": TABLE_N, **table},
        "sweep": {"n_points": SWEEP_NP, "results": {str(n): r for n, r in sweep.items()}},
    }
    (THIS_DIR / "storage_benchmark.json").write_text(json.dumps(out, indent=2))
    write_table(table, THIS_DIR / "storage_benchmark_table.tex")
    print(table["verify"])
    for name, r in table["methods"].items():
        print(f"{name:17s} bytes={r['bytes']:>12d} fwd={r['t_forward']:.3f}s bwd={r['t_backward']:.3f}s "
              f"extra={r['extra_forward_steps']} flushes={r['flushes']} max_err={r['max_error']:.3e}")
    plot(sweep)


if __name__ == "__main__":
    main()
