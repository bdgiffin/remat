"""
Finite-difference Jacobian study for integer squeeze multiplication by p
only (no division by q). This mirrors the squeeze step in `squeeze_jacobian.py`
but omits the inverse step. We retain the dual-variable format (first, second)
and compare discrete forward FD Jacobians to the smoothed continuous assumption.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
    }
)


# --------------------------------------------------------------------------- #
# Integer arithmetic (Euclidean divmod) as in arithmetic.h

def divmod_truncated(dividend: int, divisor: int) -> Tuple[int, int]:
    abs_q = abs(dividend) // abs(divisor)
    quotient = -abs_q if (dividend < 0) ^ (divisor < 0) else abs_q
    remainder = dividend - quotient * divisor
    return quotient, remainder


def divmod_euclidean(dividend: int, divisor: int) -> Tuple[int, int]:
    quotient, remainder = divmod_truncated(dividend, divisor)
    if remainder < 0:
        if divisor > 0:
            quotient -= 1
            remainder += divisor
        else:
            quotient += 1
            remainder -= divisor
    return quotient, remainder


# --------------------------------------------------------------------------- #
# Mapping: multiplication by p (squeeze step only)

def squeeze_step(first: int, second: int, p: int) -> Tuple[int, int]:
    """
    S_p(first, second): Euclidean divmod on second, then combine with first.
      (q1, r1) = divmodE(second, p)
      new_first  = first * p + r1
      new_second = q1
    """
    q1, r1 = divmod_euclidean(second, p)
    new_first = first * p + r1
    new_second = q1
    return new_first, new_second


# --------------------------------------------------------------------------- #
# Jacobians

def forward_fd_jacobian(
    mapping: Callable[[int, int, int], Tuple[int, int]],
    first: int,
    second: int,
    p: int,
    step: int = 1,
) -> np.ndarray:
    base_first, base_second = mapping(first, second, p)
    f_dx0, s_dx0 = mapping(first + step, second, p)
    f_dy0, s_dy0 = mapping(first, second + step, p)

    j00 = (f_dx0 - base_first) / step
    j10 = (s_dx0 - base_second) / step
    j01 = (f_dy0 - base_first) / step
    j11 = (s_dy0 - base_second) / step
    return np.array([[j00, j01], [j10, j11]], dtype=float)


def smoothed_jacobian(p: int) -> np.ndarray:
    """
    Continuous proxy: S_p(x, x*) ≈ (p x + x*/p, x*/p)
    """
    return np.array([[p, 1.0 / p], [0.0, 1.0 / p]], dtype=float)


# --------------------------------------------------------------------------- #
# Sweep

def run_sweep(
    p_values: Iterable[int],
    mantissa_values: Iterable[int],
    step: int,
) -> pd.DataFrame:
    records: List[Dict[str, float]] = []
    for p in p_values:
        if p <= 0:
            continue
        for first in mantissa_values:
            for second in mantissa_values:
                out_first, out_second = squeeze_step(first, second, p)
                j_fd = forward_fd_jacobian(squeeze_step, first, second, p, step)
                j_smooth = smoothed_jacobian(p)
                diff = j_fd - j_smooth
                abs_j01 = abs(j_fd[0, 1])
                abs_j10 = abs(j_fd[1, 0])
                records.append(
                    {
                        "p": p,
                        "x": first,
                        "x_star": second,
                        "y": out_first,
                        "y_star": out_second,
                        "J00": j_fd[0, 0],
                        "J01": j_fd[0, 1],
                        "J10": j_fd[1, 0],
                        "J11": j_fd[1, 1],
                        "J00_smooth": j_smooth[0, 0],
                        "J01_smooth": j_smooth[0, 1],
                        "J10_smooth": j_smooth[1, 0],
                        "J11_smooth": j_smooth[1, 1],
                        "J00_delta": diff[0, 0],
                        "J01_delta": diff[0, 1],
                        "J10_delta": diff[1, 0],
                        "J11_delta": diff[1, 1],
                        "abs_J01": abs_j01,
                        "abs_J10": abs_j10,
                    }
                )
    return pd.DataFrame.from_records(records)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    return (
        df.groupby("p")
        .agg(
            max_abs_J01=("abs_J01", "max"),
            frac_nonzero_J01=("abs_J01", lambda s: float((s != 0).mean())),
            max_diag_err=("J00_delta", lambda s: max(abs(s.max()), abs(s.min()))),
        )
        .reset_index()
    )


# --------------------------------------------------------------------------- #
# Plots

def heatmap_j01_by_p(summary: pd.DataFrame, out_path: Path) -> None:
    data = summary.set_index("p")["max_abs_J01"]
    fig, ax = plt.subplots(figsize=(8, 1.6), constrained_layout=True)
    im = ax.imshow([data.values], origin="lower", cmap="magma")
    ax.set_xticks(range(len(data)))
    ax.set_xticklabels(data.index)
    ax.set_yticks([0])
    ax.set_yticklabels(["max |J01|"])
    for j, val in enumerate(data.values):
        ax.text(j, 0, f"{val:.1f}", ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(im, ax=ax, shrink=0.8, label="max |J01| over states")
    ax.set_xlabel("p")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def heatmap_j01_by_state(df: pd.DataFrame, out_path: Path) -> None:
    mantissas = sorted(df["x"].unique())
    state_agg = (
        df.groupby(["x", "x_star"])["abs_J01"]
        .max()
        .unstack(fill_value=0)
        .reindex(index=mantissas, columns=mantissas, fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
    im = ax.imshow(
        state_agg.values,
        origin="lower",
        extent=[mantissas[0] - 0.5, mantissas[-1] + 0.5, mantissas[0] - 0.5, mantissas[-1] + 0.5],
        cmap="viridis",
        aspect="auto",
    )
    ax.set_xlabel("x* (second mantissa)")
    ax.set_ylabel("x (first mantissa)")
    ax.set_title("Max |J01| over p")
    fig.colorbar(im, ax=ax, shrink=0.8, label="max |J01|")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# CLI

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="FD Jacobians for integer squeeze multiplication by p."
    )
    parser.add_argument("--min-mantissa", type=int, default=-30, help="minimum value for x and x*")
    parser.add_argument("--max-mantissa", type=int, default=30, help="maximum value for x and x*")
    parser.add_argument(
        "--p",
        type=int,
        nargs="+",
        default=np.ones(29,dtype=int)*4000+np.array([2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]),
        # default=[2, 3, 5, 7],
        help="squeeze parameter p values (non-negative)",
    )
    parser.add_argument("--step", type=int, default=1, help="finite-difference step Δ")
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="directory to write CSVs and figures",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    mantissa_values = list(range(args.min_mantissa, args.max_mantissa + 1))

    df = run_sweep(args.p, mantissa_values, args.step)
    summary = summarize(df)

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    results_path = out_dir / "squeeze_mul_p_results.csv"
    summary_path = out_dir / "squeeze_mul_p_summary.csv"
    plot_p_path = out_dir / "squeeze_mul_p_by_p.svg"
    plot_state_path = out_dir / "squeeze_mul_p_by_state.svg"

    df.to_csv(results_path, index=False)
    summary.to_csv(summary_path, index=False)
    heatmap_j01_by_p(summary, plot_p_path)
    heatmap_j01_by_state(df, plot_state_path)

    print(f"[write] {results_path}")
    print(f"[write] {summary_path}")
    print(f"[write] {plot_p_path}")
    print(f"[write] {plot_state_path}")


if __name__ == "__main__":
    main()
