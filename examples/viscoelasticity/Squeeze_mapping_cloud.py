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
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

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
# Cloud-based affine study

def sample_cloud_points(
    center_x: int,
    center_x_star: int,
    half_width_x: int,
    half_width_x_star: int,
    sample_mode: str,
    grid_step: int,
    random_count: int,
    rng: Optional[np.random.Generator] = None,
) -> np.ndarray:
    if half_width_x <= 0 or half_width_x_star <= 0:
        raise ValueError("half-widths must be positive")
    if sample_mode == "grid":
        if grid_step <= 0:
            raise ValueError("--cloud-grid-step must be positive")
        x_vals = np.arange(center_x - half_width_x, center_x + half_width_x + 1, grid_step)
        x_star_vals = np.arange(
            center_x_star - half_width_x_star,
            center_x_star + half_width_x_star + 1,
            grid_step,
        )
        xx, xx_star = np.meshgrid(x_vals, x_star_vals, indexing="ij")
        return np.column_stack([xx.ravel(), xx_star.ravel()]).astype(float)

    if sample_mode == "random":
        if random_count <= 0:
            raise ValueError("--cloud-random-count must be positive")
        if rng is None:
            rng = np.random.default_rng()
        x_samples = rng.integers(center_x - half_width_x, center_x + half_width_x + 1, size=random_count)
        x_star_samples = rng.integers(
            center_x_star - half_width_x_star,
            center_x_star + half_width_x_star + 1,
            size=random_count,
        )
        return np.column_stack([x_samples, x_star_samples]).astype(float)

    raise ValueError(f"unsupported cloud sample mode: {sample_mode}")


def apply_mapping_to_cloud(
    points: np.ndarray,
    mapping: Callable[[int, int, int, int], Tuple[int, int]],
    p: int,
    q: int,
) -> np.ndarray:
    mapped = np.empty_like(points, dtype=float)
    for idx, (first, second) in enumerate(points):
        y, y_star = mapping(int(round(first)), int(round(second)), p, q)
        mapped[idx, 0] = y
        mapped[idx, 1] = y_star
    return mapped


def propagate_cloud(
    initial_points: np.ndarray,
    mapping: Callable[[int, int, int, int], Tuple[int, int]],
    p: int,
    q: int,
    max_iterations: int,
) -> List[np.ndarray]:
    states = [initial_points.astype(float)]
    current = initial_points.astype(float)
    for _ in range(max_iterations):
        current = apply_mapping_to_cloud(current, mapping, p, q)
        states.append(current)
    return states


def iterated_mapping(
    mapping: Callable[[int, int, int, int], Tuple[int, int]],
    iterations: int,
) -> Callable[[int, int, int, int], Tuple[int, int]]:
    def mapping_k(first: int, second: int, p: int, q: int) -> Tuple[int, int]:
        cur_first, cur_second = first, second
        for _ in range(iterations):
            cur_first, cur_second = mapping(cur_first, cur_second, p, q)
        return cur_first, cur_second

    return mapping_k


def fit_affine_map(source_points: np.ndarray, target_points: np.ndarray) -> Dict[str, np.ndarray | float]:
    if source_points.shape != target_points.shape:
        raise ValueError("source and target point clouds must have the same shape")
    if source_points.shape[0] < 3:
        raise ValueError("at least 3 points are required for affine fitting")

    source_mean = source_points.mean(axis=0)
    target_mean = target_points.mean(axis=0)
    source_centered = source_points - source_mean
    target_centered = target_points - target_mean

    a_t, _, _, _ = np.linalg.lstsq(source_centered, target_centered, rcond=None)
    a = a_t.T
    b = target_mean - a @ source_mean

    predicted = source_points @ a.T + b
    residual_norm = np.linalg.norm(target_points - predicted)
    scale_norm = np.linalg.norm(target_centered)
    non_affinity = float(residual_norm / scale_norm) if scale_norm > 0 else 0.0

    singular_values = np.linalg.svd(a, compute_uv=False)
    sigma_max = float(singular_values[0])
    sigma_min = float(singular_values[-1])
    condition = float(sigma_max / sigma_min) if sigma_min > 0 else float("inf")

    return {
        "A": a,
        "b": b,
        "non_affinity": non_affinity,
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "condition": condition,
    }


def propagate_box_corners_float_x_scale(
    initial_corners: np.ndarray,
    p: int,
    q: int,
    max_iterations: int,
) -> List[np.ndarray]:
    if q == 0:
        raise ValueError("q must be non-zero for floating-point corner propagation")
    scale = float(p) / float(q)
    states = [initial_corners.astype(float)]
    current = initial_corners.astype(float)
    for _ in range(max_iterations):
        next_state = np.empty_like(current, dtype=float)
        next_state[:, 0] = current[:, 0] * scale
        next_state[:, 1] = current[:, 1]
        states.append(next_state)
        current = next_state
    return states


def box_corners_progression_frame(
    progression: Dict[str, Dict[str, object]],
    box_tracking_mode: str,
) -> pd.DataFrame:
    records: List[Dict[str, float | int | str]] = []
    for ordering in ORDERINGS.keys():
        data = progression[ordering]
        states: List[np.ndarray] = data["states"]  # type: ignore[assignment]
        base_corners: np.ndarray = data["box_corners"]  # type: ignore[assignment]

        corners_by_iteration: List[np.ndarray] = [base_corners]
        if box_tracking_mode == "affine":
            fit_by_iteration: Dict[int, Dict[str, np.ndarray | float]] = data["fit_by_iteration"]  # type: ignore[assignment]
            for iteration in range(1, len(states)):
                fit = fit_by_iteration[iteration]
                a = fit["A"]  # type: ignore[index]
                b = fit["b"]  # type: ignore[index]
                corners_by_iteration.append(base_corners @ a.T + b)
        elif box_tracking_mode == "float-x-scale":
            float_states: List[np.ndarray] = data["float_box_corners_by_iteration"]  # type: ignore[assignment]
            corners_by_iteration = float_states
        else:
            raise ValueError(f"unsupported cloud box tracking mode: {box_tracking_mode}")

        for iteration, corners in enumerate(corners_by_iteration):
            for corner_idx, (x_val, x_star_val) in enumerate(corners):
                records.append(
                    {
                        "ordering": ordering,
                        "box_tracking_mode": box_tracking_mode,
                        "iteration": int(iteration),
                        "corner_index": int(corner_idx),
                        "x": float(x_val),
                        "x_star": float(x_star_val),
                    }
                )
    return pd.DataFrame.from_records(records)


def run_cloud_affine_study(
    p: int,
    q: int,
    center_x: int,
    center_x_star: int,
    box_half_widths: Iterable[int],
    max_iterations: int,
    sample_mode: str,
    grid_step: int,
    random_count: int,
    random_seed: Optional[int],
    fd_delta: float,
    box_tracking_mode: str = "affine",
    reference_half_width: Optional[int] = None,
) -> Tuple[pd.DataFrame, Dict[str, Dict[str, object]], int]:
    if p <= 0 or q <= 0:
        raise ValueError("cloud study requires positive --cloud-p and --cloud-q")
    if max_iterations <= 0:
        raise ValueError("--cloud-max-iterations must be positive")
    if box_tracking_mode not in {"affine", "float-x-scale"}:
        raise ValueError("cloud box tracking mode must be 'affine' or 'float-x-scale'")

    cleaned_half_widths: List[int] = []
    seen = set()
    for width in box_half_widths:
        width_int = int(width)
        if width_int <= 0:
            raise ValueError("all --cloud-box-half-widths must be positive")
        if width_int not in seen:
            cleaned_half_widths.append(width_int)
            seen.add(width_int)
    if not cleaned_half_widths:
        raise ValueError("at least one --cloud-box-half-width is required")

    ref_half_width = cleaned_half_widths[0] if reference_half_width is None else int(reference_half_width)
    if ref_half_width not in cleaned_half_widths:
        cleaned_half_widths.append(ref_half_width)

    rng = np.random.default_rng(random_seed)
    initial_clouds: Dict[int, np.ndarray] = {}
    for width in cleaned_half_widths:
        initial_clouds[width] = sample_cloud_points(
            center_x=center_x,
            center_x_star=center_x_star,
            half_width_x=width,
            half_width_x_star=width,
            sample_mode=sample_mode,
            grid_step=grid_step,
            random_count=random_count,
            rng=rng,
        )

    records: List[Dict[str, float]] = []
    progression: Dict[str, Dict[str, object]] = {}
    for ordering, mapping in ORDERINGS.items():
        fd_jacobians: Dict[int, np.ndarray] = {}
        for iteration in range(1, max_iterations + 1):
            fd_map = iterated_mapping(mapping, iteration)
            fd_jacobians[iteration] = forward_fd_jacobian(
                fd_map,
                center_x,
                center_x_star,
                p,
                q,
                step=fd_delta,
            )

        for width in cleaned_half_widths:
            source = initial_clouds[width]
            states = propagate_cloud(source, mapping, p, q, max_iterations)
            fit_by_iteration: Dict[int, Dict[str, np.ndarray | float]] = {}

            for iteration in range(1, max_iterations + 1):
                fit = fit_affine_map(source, states[iteration])
                fit_by_iteration[iteration] = fit
                fd_j = fd_jacobians[iteration]
                a = fit["A"]
                b = fit["b"]
                records.append(
                    {
                        "ordering": ordering,
                        "p": p,
                        "q": q,
                        "center_x": center_x,
                        "center_x_star": center_x_star,
                        "sample_mode": sample_mode,
                        "box_half_width": width,
                        "n_points": int(source.shape[0]),
                        "iteration": iteration,
                        "A00": float(a[0, 0]),
                        "A01": float(a[0, 1]),
                        "A10": float(a[1, 0]),
                        "A11": float(a[1, 1]),
                        "b0": float(b[0]),
                        "b1": float(b[1]),
                        "non_affinity": float(fit["non_affinity"]),
                        "sigma_max": float(fit["sigma_max"]),
                        "sigma_min": float(fit["sigma_min"]),
                        "condition": float(fit["condition"]),
                        "fd_J00": float(fd_j[0, 0]),
                        "fd_J01": float(fd_j[0, 1]),
                        "fd_J10": float(fd_j[1, 0]),
                        "fd_J11": float(fd_j[1, 1]),
                        "A01_minus_fd_J01": float(a[0, 1] - fd_j[0, 1]),
                    }
                )

            if width == ref_half_width:
                box_corners = np.array(
                    [
                        [center_x - width, center_x_star - width],
                        [center_x + width, center_x_star - width],
                        [center_x + width, center_x_star + width],
                        [center_x - width, center_x_star + width],
                    ],
                    dtype=float,
                )
                progression[ordering] = {
                    "states": states,
                    "fit_by_iteration": fit_by_iteration,
                    "box_corners": box_corners,
                }
                if box_tracking_mode == "float-x-scale":
                    progression[ordering]["float_box_corners_by_iteration"] = (
                        propagate_box_corners_float_x_scale(
                            initial_corners=box_corners,
                            p=p,
                            q=q,
                            max_iterations=max_iterations,
                        )
                    )

    metrics = pd.DataFrame.from_records(records)
    return metrics, progression, ref_half_width


def plot_cloud_deformation_progression(
    progression: Dict[str, Dict[str, object]],
    reference_half_width: int,
    out_path: Path,
    box_tracking_mode: str = "affine",
) -> None:
    n_rows = len(ORDERINGS)
    if n_rows == 0:
        return
    if box_tracking_mode not in {"affine", "float-x-scale"}:
        raise ValueError("cloud box tracking mode must be 'affine' or 'float-x-scale'")
    sample_ordering = next(iter(ORDERINGS.keys()))
    n_cols = len(progression[sample_ordering]["states"])
    fig, axes = plt.subplots(
        n_rows,
        n_cols,
        figsize=(3.0 * n_cols, 2.8 * n_rows),
        squeeze=False,
        constrained_layout=True,
    )
    corner_gradient = LinearSegmentedColormap.from_list(
        "corner_gradient",
        ["#3b0f70", "#2c7fb8", "#41b6c4", "#f1e51d"],
        N=256,
    )
    global_x_min = float("inf")
    global_x_max = float("-inf")
    global_y_min = float("inf")
    global_y_max = float("-inf")

    for row, ordering in enumerate(ORDERINGS.keys()):
        states: List[np.ndarray] = progression[ordering]["states"]  # type: ignore[assignment]
        fit_by_iteration: Dict[int, Dict[str, np.ndarray | float]] = progression[ordering][
            "fit_by_iteration"
        ]  # type: ignore[assignment]
        box_corners: np.ndarray = progression[ordering]["box_corners"]  # type: ignore[assignment]
        float_box_states: Optional[List[np.ndarray]] = None
        if box_tracking_mode == "float-x-scale":
            float_box_states = progression[ordering]["float_box_corners_by_iteration"]  # type: ignore[assignment]
        source_points = states[0]
        corner_scalar = source_points[:, 0] + source_points[:, 1]
        scalar_min = float(corner_scalar.min())
        scalar_max = float(corner_scalar.max())
        if scalar_max > scalar_min:
            point_colors = (corner_scalar - scalar_min) / (scalar_max - scalar_min)
        else:
            point_colors = np.zeros_like(corner_scalar)

        for iteration, points in enumerate(states):
            ax = axes[row, iteration]
            ax.scatter(
                points[:, 0],
                points[:, 1],
                s=9,
                alpha=0.8,
                c=point_colors,
                cmap=corner_gradient,
                vmin=0.0,
                vmax=1.0,
                linewidths=0,
            )

            if iteration == 0:
                mapped_corners = box_corners
                style = {"linestyle": "--", "linewidth": 1.2}
            else:
                if box_tracking_mode == "affine":
                    fit = fit_by_iteration[iteration]
                    a = fit["A"]  # type: ignore[index]
                    b = fit["b"]  # type: ignore[index]
                    mapped_corners = box_corners @ a.T + b
                else:
                    assert float_box_states is not None
                    mapped_corners = float_box_states[iteration]
                style = {"linestyle": "-", "linewidth": 1.2}

            polygon = np.vstack([mapped_corners, mapped_corners[0]])
            ax.plot(polygon[:, 0], polygon[:, 1], color="black", **style)

            global_x_min = min(global_x_min, float(points[:, 0].min()), float(polygon[:, 0].min()))
            global_x_max = max(global_x_max, float(points[:, 0].max()), float(polygon[:, 0].max()))
            global_y_min = min(global_y_min, float(points[:, 1].min()), float(polygon[:, 1].min()))
            global_y_max = max(global_y_max, float(points[:, 1].max()), float(polygon[:, 1].max()))

            if row == 0:
                ax.set_title(f"k={iteration}", fontsize=10)
            if iteration == 0:
                ax.set_ylabel(
                    f"{ORDERING_TITLES.get(ordering, ordering)}\nancillary variable x*",
                    fontsize="large",
                )
            if row == n_rows - 1:
                ax.set_xlabel("mapped x", fontsize="large")
            ax.grid(alpha=0.2, linestyle=":")

    if np.isfinite(global_x_min) and np.isfinite(global_x_max) and np.isfinite(global_y_min) and np.isfinite(global_y_max):
        global_min = min(global_x_min, global_y_min)
        global_max = max(global_x_max, global_y_max)
        span = max(global_max - global_min, 1.0)
        pad = 0.03 * span
        lim_min = global_min - pad
        lim_max = global_max + pad
        for row in range(n_rows):
            for col in range(n_cols):
                ax = axes[row, col]
                ax.set_xlim(lim_min, lim_max)
                ax.set_ylim(lim_min, lim_max)
                ax.set_aspect("equal", adjustable="box")

    mode_label = "fitted affine box" if box_tracking_mode == "affine" else "float box x <- x*p/q, x* fixed"
    # fig.suptitle(f"Cloud deformation and {mode_label} (half-width={reference_half_width})", fontsize=11)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_a01_vs_box(metrics: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True, constrained_layout=True)
    for ax, ordering in zip(axes, ORDERINGS.keys()):
        subset = metrics[metrics["ordering"] == ordering]
        iterations = sorted(subset["iteration"].unique())
        if not iterations:
            continue
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(iterations)))
        for color, iteration in zip(colors, iterations):
            sub_it = subset[subset["iteration"] == iteration].sort_values("box_half_width")
            ax.plot(
                sub_it["box_half_width"],
                sub_it["A01"],
                marker="o",
                color=color,
                label=f"fit k={iteration}",
            )
            ax.scatter(
                sub_it["box_half_width"],
                sub_it["fd_J01"],
                marker="x",
                color=color,
                s=36,
                label=f"FD k={iteration}",
            )
        ax.set_xlabel("box half-width", fontsize="large")
        ax.set_title(ORDERING_TITLES.get(ordering, ordering), fontsize="large")
        ax.grid(alpha=0.25, linestyle=":")

    axes[0].set_ylabel(r"$A_{01}$ (cloud fit) and FD $J_{01}$", fontsize="large")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=8)
    fig.savefig(out_path, dpi=300)
    plt.close(fig)


def plot_non_affinity_vs_box(metrics: pd.DataFrame, out_path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharey=True, constrained_layout=True)
    for ax, ordering in zip(axes, ORDERINGS.keys()):
        subset = metrics[metrics["ordering"] == ordering]
        iterations = sorted(subset["iteration"].unique())
        if not iterations:
            continue
        colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(iterations)))
        for color, iteration in zip(colors, iterations):
            sub_it = subset[subset["iteration"] == iteration].sort_values("box_half_width")
            ax.plot(
                sub_it["box_half_width"],
                sub_it["non_affinity"],
                marker="o",
                color=color,
                label=f"k={iteration}",
            )
        ax.set_xlabel("box half-width", fontsize="large")
        ax.set_title(ORDERING_TITLES.get(ordering, ordering), fontsize="large")
        ax.grid(alpha=0.25, linestyle=":")

    axes[0].set_ylabel(r"Non-affinity $\|Y-\hat{Y}\|/\|Y-\bar{Y}\|$", fontsize="large")
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=4, fontsize=8)
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
    parser.add_argument(
        "--cloud-study",
        action="store_true",
        help="run cloud-based affine sensitivity study instead of the full FD sweep",
    )
    parser.add_argument(
        "--cloud-p",
        type=int,
        default=8,
        help="p parameter for cloud study",
    )
    parser.add_argument(
        "--cloud-q",
        type=int,
        default=7,
        help="q parameter for cloud study",
    )
    parser.add_argument(
        "--cloud-center-x",
        type=int,
        default=0,
        help="center x for cloud study",
    )
    parser.add_argument(
        "--cloud-center-x-star",
        type=int,
        default=0,
        help="center x* for cloud study",
    )
    parser.add_argument(
        "--cloud-box-half-widths",
        type=int,
        nargs="+",
        default=[2, 4, 8, 16],
        help="half-widths of square cloud boxes around the center",
    )
    parser.add_argument(
        "--cloud-max-iterations",
        type=int,
        default=2,
        help="max number of repeated squeeze mappings for cloud study",
    )
    parser.add_argument(
        "--cloud-sample-mode",
        choices=["grid", "random"],
        default="grid",
        help="cloud point sampling mode",
    )
    parser.add_argument(
        "--cloud-grid-step",
        type=int,
        default=1,
        help="integer step for grid cloud sampling",
    )
    parser.add_argument(
        "--cloud-random-count",
        type=int,
        default=400,
        help="number of sampled points when --cloud-sample-mode random",
    )
    parser.add_argument(
        "--cloud-random-seed",
        type=int,
        default=None,
        help="optional RNG seed for cloud random sampling",
    )
    parser.add_argument(
        "--cloud-reference-half-width",
        type=int,
        default=None,
        help="half-width used for deformation progression plot (defaults to first value in --cloud-box-half-widths)",
    )
    parser.add_argument(
        "--cloud-box-tracking-mode",
        choices=["affine", "float-x-scale"],
        default="float-x-scale",
        help="how to propagate the plotted box: affine fit or floating-point x <- x*p/q with x* unchanged",
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
    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.cloud_study:
        cloud_metrics, cloud_progression, ref_half_width = run_cloud_affine_study(
            p=args.cloud_p,
            q=args.cloud_q,
            center_x=args.cloud_center_x,
            center_x_star=args.cloud_center_x_star,
            box_half_widths=args.cloud_box_half_widths,
            max_iterations=args.cloud_max_iterations,
            sample_mode=args.cloud_sample_mode,
            grid_step=args.cloud_grid_step,
            random_count=args.cloud_random_count,
            random_seed=args.cloud_random_seed,
            fd_delta=args.fd_delta,
            box_tracking_mode=args.cloud_box_tracking_mode,
            reference_half_width=args.cloud_reference_half_width,
        )

        cloud_metrics_path = out_dir / "squeeze_cloud_affine_metrics.csv"
        cloud_progress_path = out_dir / "squeeze_cloud_deformation_progression.svg"
        cloud_a01_path = out_dir / "squeeze_cloud_A01_vs_box.svg"
        cloud_non_affinity_path = out_dir / "squeeze_cloud_non_affinity_vs_box.svg"
        cloud_box_corners_path = out_dir / "squeeze_cloud_box_corners.csv"
        cloud_box_corners = box_corners_progression_frame(
            cloud_progression,
            box_tracking_mode=args.cloud_box_tracking_mode,
        )

        cloud_metrics.to_csv(cloud_metrics_path, index=False)
        cloud_box_corners.to_csv(cloud_box_corners_path, index=False)
        plot_cloud_deformation_progression(
            cloud_progression,
            ref_half_width,
            cloud_progress_path,
            box_tracking_mode=args.cloud_box_tracking_mode,
        )
        plot_a01_vs_box(cloud_metrics, cloud_a01_path)
        plot_non_affinity_vs_box(cloud_metrics, cloud_non_affinity_path)

        print(f"[write] {cloud_metrics_path}")
        print(f"[write] {cloud_box_corners_path}")
        print(f"[write] {cloud_progress_path}")
        print(f"[write] {cloud_a01_path}")
        print(f"[write] {cloud_non_affinity_path}")
        return

    mantissa_values = build_mantissa_values(args)
    df = run_sweep(args.p, args.q, mantissa_values, args.fd_delta)
    summary = summarize(df)

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
