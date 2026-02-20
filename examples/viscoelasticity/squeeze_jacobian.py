"""
Finite-difference Jacobian study for the integer squeeze mapping defined in
`src/arithmetic.h` (function `squeezeE`).  
Mirror that implementation exactly, then compute forward FD Jacobians (Δ configurable, default 1) for two
orderings:

  - Sqinv_Sp : S_q^{-1} ∘ S_p (matches current squeezeE implementation in C++)
  - Sp_Sqinv : S_p ∘ S_q^{-1}

We sweep over integer states (x, x*) and squeeze parameters (p, q), log the
full FD Jacobians, compare side by side for the two orderings, and visualize the off-diagonal sensitivity J01 as a

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
# Integer arithmetic helpers (mirror arithmetic.h)

def divmod_truncated(dividend: int, divisor: int) -> Tuple[int, int]:
    """
    Truncated division (quotient toward zero) with remainder carrying the
    dividend's sign. Matches C++ integer / and %.
    """
    if divisor == 0:
        raise ZeroDivisionError("divisor must be non-zero")
    abs_q = abs(dividend) // abs(divisor)
    quotient = -abs_q if (dividend < 0) ^ (divisor < 0) else abs_q
    remainder = dividend - quotient * divisor
    return quotient, remainder


def divmod_euclidean(dividend: int, divisor: int) -> Tuple[int, int]:
    """
    Euclidean division: remainder shares the divisor's sign, adjusted from the
    truncated result to ensure 0 <= |remainder| < |divisor|.
    """
    quotient, remainder = divmod_truncated(dividend, divisor)
    if remainder < 0:
        if divisor > 0:
            quotient -= 1
            remainder += divisor
        else:
            quotient += 1
            remainder -= divisor
    return quotient, remainder


def squeeze_step(first: int, second: int, p: int) -> Tuple[int, int]:
    """
    S_p(first, second): squeeze mapping using parameter p.
    """
    q1, r1 = divmod_euclidean(second, p)
    new_first = first * p + r1
    new_second = q1
    return new_first, new_second


def inverse_squeeze_step(first: int, second: int, q: int) -> Tuple[int, int]:
    """
    S_q^{-1}(first, second): inverse squeeze mapping using parameter q.
    """
    q2, r2 = divmod_euclidean(first, q)
    new_second = second * q + r2
    new_first = q2
    return new_first, new_second


def compose_sqinv_sp(first: int, second: int, p: int, q: int) -> Tuple[int, int]:
    """
    Ordering implemented by squeezeE: S_q^{-1} ∘ S_p.
    """
    a, b = squeeze_step(first, second, p)
    return inverse_squeeze_step(a, b, q)


def compose_sp_sqinv(first: int, second: int, p: int, q: int) -> Tuple[int, int]:
    """
    Alternate ordering: S_p ∘ S_q^{-1}.
    """
    a, b = inverse_squeeze_step(first, second, q)
    return squeeze_step(a, b, p)


# --------------------------------------------------------------------------- #
# Jacobian calculations

def forward_fd_jacobian(
    mapping: Callable[[int, int, int, int], Tuple[int, int]],
    first: int,
    second: int,
    p: int,
    q: int,
    step: float = 1.0,
) -> np.ndarray:
    """
    Forward finite-difference Jacobian with Δ=step (default 1).
    """
    base_first, base_second = mapping(first, second, p, q)
    f_dx0, s_dx0 = mapping(first + step, second, p, q)
    f_dy0, s_dy0 = mapping(first, second + step, p, q)

    j00 = (f_dx0 - base_first) / step
    j10 = (s_dx0 - base_second) / step
    j01 = (f_dy0 - base_first) / step
    j11 = (s_dy0 - base_second) / step

    return np.array([[j00, j01], [j10, j11]], dtype=float)


def smoothed_jacobian(p: int, q: int) -> np.ndarray:
    """
    Jacobian if div/mod were smooth (continuous) operations.
    """
    return np.array([[p / q, 0.0], [0.0, q / p]], dtype=float)


# --------------------------------------------------------------------------- #
# Sweeps and aggregation

ORDERINGS: Dict[str, Callable[[int, int, int, int], Tuple[int, int]]] = {
    "Sqinv_Sp": compose_sqinv_sp,
    "Sp_Sqinv": compose_sp_sqinv,
}

# Pretty mathtext labels for plot titles (order preserved via ORDERINGS keys)
ORDERING_TITLES = {
    "Sqinv_Sp": r"($S_q^{-1} \circ S_p)(\mathbf{x})$",
    "Sp_Sqinv": r"($S_p \circ S_q^{-1})(\mathbf{x})$",
}


def run_sweep(
    p_values: Iterable[int],
    q_values: Iterable[int],
    mantissa_values: Iterable[int],
    step: int,
) -> pd.DataFrame:
    records: List[Dict[str, float]] = []
    for p in p_values:
        for q in q_values:
            if p <= 0 or q <= 0:
                continue  # squeezeE assumes non-negative parameters
            for first in mantissa_values:
                for second in mantissa_values:
                    for ordering, mapping in ORDERINGS.items():
                        out_first, out_second = mapping(first, second, p, q)
                        j_fd = forward_fd_jacobian(mapping, first, second, p, q, step)
                        j_smooth = smoothed_jacobian(p, q)
                        diff = j_fd - j_smooth
                        abs_j01 = abs(j_fd[0, 1])
                        abs_j10 = abs(j_fd[1, 0])
                        offdiag = max(abs_j01, abs_j10)
                        diag_error = max(abs(diff[0, 0]), abs(diff[1, 1]))
                        records.append(
                            {
                                "p": p,
                                "q": q,
                                "x": first,
                                "x_star": second,
                                "ordering": ordering,
                                "y": out_first,
                                "y_star": out_second,
                                "J00": j_fd[0, 0],
                                "J01": j_fd[0, 1],
                                "J10": j_fd[1, 0],
                                "J11": j_fd[1, 1],
                                "J00_smooth": j_smooth[0, 0],
                                "J11_smooth": j_smooth[1, 1],
                                "J00_delta": diff[0, 0],
                                "J11_delta": diff[1, 1],
                                "abs_J01": abs_j01,
                                "abs_J10": abs_j10,
                                "offdiag_max": offdiag,
                                "diag_error_max": diag_error,
                            }
                        )
    return pd.DataFrame.from_records(records)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    grouped = (
        df.groupby(["ordering", "p", "q"])
        .agg(
            offdiag01_max=("abs_J01", "max"),
            offdiag01_count=("abs_J01", lambda s: int((s != 0).sum())),
            offdiag01_fraction=("abs_J01", lambda s: float((s != 0).mean())),
            diag_error_max=("diag_error_max", "max"),
        )
        .reset_index()
    )
    return grouped


# --------------------------------------------------------------------------- #
# Plotting

def heatmap_offdiag_by_pq(summary: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True, constrained_layout=True)
    for ax, ordering in zip(axes, ORDERINGS.keys()):
        data = (
            summary[summary["ordering"] == ordering]
            .pivot(index="p", columns="q", values="offdiag01_max")
            .sort_index()
        )
        im = ax.imshow(data.values, origin="lower", cmap="magma")
        # show every 4th tick on p and q to declutter; always include the last tick
        x_positions = list(range(0, len(data.columns), 4))
        if (len(data.columns) - 1) not in x_positions:
            x_positions.append(len(data.columns) - 1)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([data.columns[i] for i in x_positions])

        y_positions = list(range(0, len(data.index), 4))
        if (len(data.index) - 1) not in y_positions:
            y_positions.append(len(data.index) - 1)
        ax.set_yticks(y_positions)
        ax.set_yticklabels([data.index[i] for i in y_positions])
        ax.set_xlabel("q", fontsize="large")
        ax.set_ylabel("p", fontsize="large")
        ax.set_title(ORDERING_TITLES.get(ordering, ordering), fontsize="large")
        # for (i, j), val in np.ndenumerate(data.values):
        #     ax.text(j, i, f"{val:.1f}", ha="center", va="center", color="white", fontsize=8)
    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8, label=r"Cross-sensitivity $|J_{12}^{\mathrm{FD}}|$")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def hist_offdiag_by_pq(summary: pd.DataFrame, out_path: Path) -> None:
    """
    Histogram of max |J01| values across all (p,q) cells for each ordering.
    Interpreted as an empirical probability distribution over heatmap values.
    """
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True, constrained_layout=True)
    for ax, ordering in zip(axes, ORDERINGS.keys()):
        vals = summary.loc[summary["ordering"] == ordering, "offdiag01_max"].values
        if vals.size == 0:
            continue
        upper = max(vals)
        bins = np.arange(-0.5, upper + 1.5, 1.0)  # integer bins centered on counts
        ax.hist(vals, bins=bins, density=True, alpha=0.85, color="tab:blue", edgecolor="black")
        ax.set_xlabel("max |J01| per (p,q) cell", fontsize="large")
        ax.set_ylabel("Probability density", fontsize="large")
        # ax.set_title(f"Distribution of max |J01| ({ordering})", fontsize="large")
        ax.set_title(ORDERING_TITLES.get(ordering, ordering), fontsize="large")

    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def heatmap_offdiag_by_state(df: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(
        1, 2, figsize=(10, 4), sharey=True, sharex=True, constrained_layout=True
    )
    mantissas = sorted(df["x"].unique())
    for ax, ordering in zip(axes, ORDERINGS.keys()):
        state_agg = (
            df[df["ordering"] == ordering]
            .groupby(["x", "x_star"])["abs_J01"]
            .max()
            .unstack(fill_value=0)
            .reindex(index=mantissas, columns=mantissas, fill_value=0)
        )
        im = ax.imshow(
            state_agg.values,
            origin="lower",
            extent=[mantissas[0] - 0.5, mantissas[-1] + 0.5, mantissas[0] - 0.5, mantissas[-1] + 0.5],
            cmap="viridis",
            aspect="auto",
        )
        ax.set_xlabel("x*", fontsize="large")
        ax.set_ylabel("x", fontsize="large")
        # ax.set_title(f"Max |J01| over p,q ({ordering})", fontsize="large")
        ax.set_title(ORDERING_TITLES.get(ordering, ordering), fontsize="large")

    fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8, label=r"Cross-sensitivity $|J_{12}^{\mathrm{FD}}|$")
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def violin_offdiag_by_delta(
    p_values: Iterable[int],
    q_values: Iterable[int],
    mantissa_values: Iterable[int],
    deltas: Iterable[float],
    out_path: Path,
) -> None:
    """
    For each finite-difference delta, sweep (p,q,x,x*) and collect max |J01|
    per (p,q) cell, then visualize distributions across cells with violins.
    """
    ordering_colors = {"Sqinv_Sp": "tab:purple", "Sp_Sqinv": "tab:orange"}
    delta_list = list(deltas)
    all_data = {o: [] for o in ORDERINGS.keys()}

    for delta in delta_list:
        df_delta = run_sweep(p_values, q_values, mantissa_values, step=delta)
        summary_delta = summarize(df_delta)
        for ordering in ORDERINGS.keys():
            vals = summary_delta.loc[summary_delta["ordering"] == ordering, "offdiag01_max"].values
            all_data[ordering].append(vals)

    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
    positions = np.arange(len(delta_list))
    width = 0.35

    parts_left = ax.violinplot(
        all_data["Sqinv_Sp"],
        positions=positions - width / 2,
        widths=0.28,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for pc in parts_left["bodies"]:
        pc.set_facecolor(ordering_colors["Sqinv_Sp"])
        pc.set_alpha(0.6)
    parts_left["cmedians"].set_color("black")

    parts_right = ax.violinplot(
        all_data["Sp_Sqinv"],
        positions=positions + width / 2,
        widths=0.28,
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )
    for pc in parts_right["bodies"]:
        pc.set_facecolor(ordering_colors["Sp_Sqinv"])
        pc.set_alpha(0.6)
    parts_right["cmedians"].set_color("black")

    ax.set_xticks(positions)
    ax.set_xticklabels([f"{d:g}" for d in delta_list])
    ax.set_xlabel(r"Finite-difference $\Delta$", fontsize="large")
    ax.set_ylabel(r"Cross-sensitivity $|J_{12}^{\mathrm{FD}}|$", fontsize="large")
    # ax.set_title(r"Sensitivity vs $\Delta$ (violins = PDF over $(p,q)$)", fontsize="large")
    ax.legend(
        handles=[
            plt.Line2D(
                [0], [0],
                color=ordering_colors["Sqinv_Sp"],
                lw=6,
                alpha=0.6,
                label=ORDERING_TITLES["Sqinv_Sp"],
            ),
            plt.Line2D(
                [0], [0],
                color=ordering_colors["Sp_Sqinv"],
                lw=6,
                alpha=0.6,
                label=ORDERING_TITLES["Sp_Sqinv"],
            ),
        ],
        loc="upper right",
        fontsize="large",
    )
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# CLI

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Forward FD Jacobians for integer squeeze mapping (squeezeE)."
    )
    parser.add_argument("--min-mantissa", type=int, default=-20, help="minimum value for x and x*")
    parser.add_argument("--max-mantissa", type=int, default=20, help="maximum value for x and x*")
    parser.add_argument(
        "--random-mantissa",
        action="store_true",
        help="sample random x/x* values instead of using a contiguous range",
    )
    parser.add_argument(
        "--random-count",
        type=int,
        default=40,
        help="number of random x (and x*) values to sample when --random-mantissa is set",
    )
    parser.add_argument(
        "--random-low",
        type=int,
        default=-2000,
        help="minimum random x/x* value when --random-mantissa is set",
    )
    parser.add_argument(
        "--random-high",
        type=int,
        default=2000,
        help="maximum random x/x* value when --random-mantissa is set",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=None,
        help="optional RNG seed for reproducible random x/x* samples",
    )
    parser.add_argument(
        "--p",
        type=int,
        nargs="+",
        default=np.array([2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]),
        # default=[2, 3, 5, 7],
        help="squeeze parameter p values (non-negative)",
    )
    parser.add_argument(
        "--q",
        type=int,
        nargs="+",
        # default=[2, 3, 5, 7],
        default=np.array([2, 3, 4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30]),
        help="squeeze parameter q values (non-negative)",
    )
    parser.add_argument(
        "--fd-delta",
        "--step",
        dest="fd_delta",
        type=float,
        default=1.0,
        help="finite-difference perturbation Δ for x and x* (can be >1)",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="directory to write CSVs and figures",
    )
    parser.add_argument(
        "--fd-delta-violin",
        action="store_true",
        help="if set, sweep Δ=1..10 and add a violin plot of max |J01| distributions over (p,q).",
    )
    return parser.parse_args()


def build_mantissa_values(args: argparse.Namespace) -> List[int]:
    if not args.random_mantissa:
        return list(range(args.min_mantissa, args.max_mantissa + 1))

    if args.random_count <= 0:
        raise ValueError("--random-count must be positive")

    low = min(args.random_low, args.random_high)
    high = max(args.random_low, args.random_high)
    total_available = high - low + 1
    if args.random_count > total_available:
        raise ValueError(
            f"--random-count={args.random_count} exceeds available integers in "
            f"[{low}, {high}] ({total_available})"
        )

    rng = np.random.default_rng(args.random_seed)
    sampled = rng.choice(np.arange(low, high + 1), size=args.random_count, replace=False)
    return sorted(int(v) for v in sampled)


def main() -> None:
    args = parse_args()
    mantissa_values = build_mantissa_values(args)

    df = run_sweep(args.p, args.q, mantissa_values, args.fd_delta)
    summary = summarize(df)

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    results_path = out_dir / "squeeze_jacobian_results.csv"
    summary_path = out_dir / "squeeze_jacobian_summary.csv"
    plot_pq_path = out_dir / "squeeze_offdiag_by_pq.svg"
    plot_pq_hist_path = out_dir / "squeeze_offdiag_hist.svg"
    plot_state_path = out_dir / "squeeze_offdiag_by_state.svg"
    plot_violin_path = out_dir / "squeeze_offdiag_violin_by_delta.svg"

    # df.to_csv(results_path, index=False)
    # summary.to_csv(summary_path, index=False)

    heatmap_offdiag_by_pq(summary, plot_pq_path)
    hist_offdiag_by_pq(summary, plot_pq_hist_path)
    heatmap_offdiag_by_state(df, plot_state_path)
    if args.fd_delta_violin:
        violin_offdiag_by_delta(
            args.p,
            args.q,
            mantissa_values,
            deltas=range(1, 11),
            out_path=plot_violin_path,
        )

    print(f"[write] {results_path}")
    print(f"[write] {summary_path}")
    print(f"[write] {plot_pq_path}")
    print(f"[write] {plot_pq_hist_path}")
    print(f"[write] {plot_state_path}")
    if args.fd_delta_violin:
        print(f"[write] {plot_violin_path}")


if __name__ == "__main__":
    main()
