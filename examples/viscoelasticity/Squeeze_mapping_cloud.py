"""
Generate the Paper 1 squeeze-map deformation progression figure.

It places each initial fixed-point lattice directly before its image under the
fixed-point realization of S_{p/q}. The integer squeeze map is applied to the
mantissas, and the returned mantissas are rescaled by the lattice spacing h.
Tile color is assigned from the initial z_1 coordinate, so the mapped cloud
records how columns of the original lattice are carried by the fixed-point map.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]

P = 7
Q = 8
REFINEMENTS = (16, 32)
OUTPUT = SCRIPT_DIR / "squeeze_cloud_deformation_progression.pdf"

REFERENCE_COLOR = "#2b738eff"
LOW_COLOR = "#2b738eff"
HIGH_COLOR = "#f9826bff"
POINT_CMAP = LinearSegmentedColormap.from_list(
    "squeeze_progression", [LOW_COLOR, "#b9aaa4", HIGH_COLOR]
)
POINT_NORM = Normalize(vmin=0.0, vmax=1.0)
POINT_SIZE = {16: 7, 32: 3.05, 64: 0.5, 128: 0.30, 256: 0.16}
COMMON_X_LIMITS = (-0.03, 1.03)
COMMON_Y_LIMITS = (-0.03, 1.50)

CORE_DRIVER_SOURCE = SCRIPT_DIR / "squeeze_core_driver.cpp"
CORE_DRIVER_EXE = REPO_ROOT / "build" / "squeeze_core_driver"

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
    }
)


def build_core_driver() -> Path:
    deps = [
        CORE_DRIVER_SOURCE,
        REPO_ROOT / "src/arithmetic.h",
        REPO_ROOT / "src/types.h",
    ]
    if CORE_DRIVER_EXE.exists() and CORE_DRIVER_EXE.stat().st_mtime >= max(
        dep.stat().st_mtime for dep in deps
    ):
        return CORE_DRIVER_EXE

    CORE_DRIVER_EXE.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        [
            os.environ.get("CXX", "g++"),
            "-std=c++17",
            f"-I{REPO_ROOT / 'src'}",
            str(CORE_DRIVER_SOURCE),
            "-o",
            str(CORE_DRIVER_EXE),
        ],
        check=True,
    )
    return CORE_DRIVER_EXE


def core_squeeze(points: np.ndarray) -> np.ndarray:
    stdin = "\n".join(f"{int(x)} {int(x_star)}" for x, x_star in points) + "\n"
    result = subprocess.run(
        [str(build_core_driver()), str(P), str(Q)],
        input=stdin,
        text=True,
        capture_output=True,
        check=True,
    )
    values = np.fromstring(result.stdout, sep=" ", dtype=np.int64)
    if values.size != 2 * len(points):
        raise RuntimeError("unexpected number of values returned by squeeze core helper")
    return values.reshape((-1, 2))


def lattice_image(refinement: int) -> tuple[np.ndarray, np.ndarray]:
    axis = np.arange(refinement + 1, dtype=np.int64)
    x, x_star = np.meshgrid(axis, axis)
    points = np.column_stack((x.ravel(), x_star.ravel()))
    colors = points[:, 0].astype(float) / refinement
    mapped = core_squeeze(points).astype(float) / refinement
    return mapped, colors


def initial_lattice(refinement: int) -> tuple[np.ndarray, np.ndarray]:
    axis = np.arange(refinement + 1, dtype=np.int64)
    x, x_star = np.meshgrid(axis, axis)
    points = np.column_stack((x.ravel(), x_star.ravel())).astype(float) / refinement
    colors = x.ravel().astype(float) / refinement
    return points, colors


def continuous_reference_box() -> np.ndarray:
    x_scale = P / Q
    x_star_scale = Q / P
    return np.array(
        [
            [0, 0],
            [x_scale, 0],
            [x_scale, x_star_scale],
            [0, x_star_scale],
            [0, 0],
        ],
        dtype=float,
    )


def draw_gradient_key(ax) -> None:
    gradient = np.linspace(0.0, 1.0, 256).reshape(1, -1)
    ax.imshow(gradient, cmap=POINT_CMAP, norm=POINT_NORM, aspect="auto")
    ax.set_xticks([0, 255])
    ax.set_xticklabels([r"$0$", r"$1$"])
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize="small", length=0, pad=1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlabel(r"initial column $z_1$", fontsize="medium", labelpad=1)


def draw_points(ax, points: np.ndarray, colors: np.ndarray, refinement: int) -> None:
    ax.scatter(
        points[:, 0],
        points[:, 1],
        s=POINT_SIZE.get(refinement, 1.0),
        c=colors,
        cmap=POINT_CMAP,
        norm=POINT_NORM,
        marker="s",
        alpha=0.82,
        linewidths=0,
        rasterized=True,
    )


def draw_panel_label(ax, text: str, *, left: float = 0.035) -> None:
    ax.text(
        left,
        0.965,
        text,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize="medium",
        bbox={
            "boxstyle": "round,pad=0.2",
            "facecolor": "white",
            "edgecolor": "0.2",
            "linewidth": 0.7,
            "alpha": 0.92,
        },
    )


def style_axes(
    ax,
    *,
    x_limits: tuple[float, float],
    y_limits: tuple[float, float],
    xticks: list[float],
    yticks: list[float],
    xticklabels: list[str],
    yticklabels: list[str],
    show_ylabel: bool,
) -> None:
    ax.set_xlim(*x_limits)
    ax.set_ylim(*y_limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"$z_1$", fontsize="large", labelpad=0)
    ax.set_ylabel(r"$z_2$" if show_ylabel else "", fontsize="large", labelpad=1)
    if show_ylabel:
        square_center = (0.5 - y_limits[0]) / (y_limits[1] - y_limits[0])
        ax.yaxis.label.set_y(square_center)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.set_xticklabels(xticklabels)
    ax.set_yticklabels(yticklabels)
    ax.tick_params(axis="both", labelsize="small", length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)


def draw_initial(ax, refinement: int, *, show_ylabel: bool) -> None:
    points, colors = initial_lattice(refinement)
    exponent = refinement.bit_length() - 1
    draw_points(ax, points, colors, refinement)
    draw_panel_label(ax, rf"$h=2^{{-{exponent}}}$")
    style_axes(
        ax,
        x_limits=COMMON_X_LIMITS,
        y_limits=COMMON_Y_LIMITS,
        xticks=[0.0, 1.0],
        yticks=[0.0, 1.0],
        xticklabels=[r"$0$", r"$1$"],
        yticklabels=[r"$0$", r"$1$"],
        show_ylabel=show_ylabel,
    )


def draw_mapped(ax, refinement: int, *, show_ylabel: bool) -> None:
    mapped, colors = lattice_image(refinement)
    reference_box = continuous_reference_box()
    x_scale = P / Q
    x_star_scale = Q / P

    draw_points(ax, mapped, colors, refinement)
    ax.plot(
        reference_box[:, 0],
        reference_box[:, 1],
        "--",
        color=REFERENCE_COLOR,
        linewidth=1.6,
        dashes=(6, 5),
    )
    # draw_panel_label(ax, r"$S_{p/q2}$", left=0.04)
    style_axes(
        ax,
        x_limits=COMMON_X_LIMITS,
        y_limits=COMMON_Y_LIMITS,
        xticks=[0.0, x_scale],
        yticks=[0.0, x_star_scale],
        xticklabels=[r"$0$", r"$p/q$"],
        yticklabels=[r"$0$", r"$q/p$"],
        show_ylabel=show_ylabel,
    )


def add_scenario_divider(fig, left_ax, right_ax) -> None:
    left_box = left_ax.get_position()
    right_box = right_ax.get_position()
    x_position = 0.5 * (left_box.x1 + right_box.x0)
    fig.add_artist(
        Line2D(
            [x_position, x_position],
            [left_box.y0-.1, left_box.y1],
            transform=fig.transFigure,
            color="0.82",
            linewidth=0.8,
        )
    )


def add_transform_arrow(fig, left_ax, right_ax) -> None:
    left_box = left_ax.get_position()
    right_box = right_ax.get_position()
    y_position = max(left_box.y1, right_box.y1) - 0.35
    x_start = left_box.x1 + 0.012
    x_end = right_box.x0 - 0.012
    arrow = FancyArrowPatch(
        (x_start, y_position),
        (x_end, y_position),
        transform=fig.transFigure,
        arrowstyle="->",
        mutation_scale=9,
        linewidth=0.9,
        color="0.35",
    )
    fig.add_artist(arrow)
    fig.text(
        0.5 * (x_start + x_end),
        y_position +.02,
        r"$\widehat{S}_{p/q,h}$",
        ha="center",
        va="bottom",
        fontsize="small",
    )


def main() -> None:
    if P <= 0 or Q <= 0:
        raise ValueError("P and Q must be positive for squeezeE")

    fig = plt.figure(figsize=(6.5, 3.0), constrained_layout=False)
    grid = fig.add_gridspec(
        2,
        4,
        height_ratios=[1.0, 0.055],
        width_ratios=[1.0, 0.88, 1.0, 0.88],
        left=0.07,
        right=0.98,
        bottom=0.14,
        top=0.86,
        hspace=0.85,
        wspace=0.18,
    )
    axes = [fig.add_subplot(grid[0, index]) for index in range(4)]
    for pair_index, refinement in enumerate(REFINEMENTS):
        first_col = 2 * pair_index
        draw_initial(axes[first_col], refinement, show_ylabel=first_col == 0)
        draw_mapped(axes[first_col + 1], refinement, show_ylabel=False)
    add_scenario_divider(fig, axes[1], axes[2])
    add_transform_arrow(fig, axes[0], axes[1])
    add_transform_arrow(fig, axes[2], axes[3])

    # key_ax = fig.add_subplot(grid[1, :])
    # draw_gradient_key(key_ax)
    fig.legend(
        [
            Line2D(
                [0],
                [0],
                linestyle=(0, (6, 5)),
                color=REFERENCE_COLOR,
                linewidth=1.6,
            )
        ],
        [r"continuous reference $S^{\mathbb{R}}_{p/q}$"],
        loc="lower center",
        frameon=False,
        fontsize="small",
        bbox_to_anchor=(0.5, 0.15),
        borderaxespad=0.0,
    )
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=200, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
