"""
Generate the Paper 1 squeeze-map lattice figure.

The finite points are produced by `squeezeE` in `src/arithmetic.h` through a
small C++ helper. The dashed boundary is the real-valued squeeze reference.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch


SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]

P = 7
Q = 8
REFINEMENTS = (32, 256)
OUTPUT = SCRIPT_DIR / "squeeze_cloud_deformation_progression.pdf"

LATTICE_COLOR = "#f9826bff"
REFERENCE_COLOR = "#2b738eff"
POINT_SIZE = {16: 10, 32: 5, 64: 0.5, 128: 0.25, 256: 0.5}

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
    deps = [CORE_DRIVER_SOURCE, REPO_ROOT / "src/arithmetic.h", REPO_ROOT / "src/types.h"]
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


def lattice_pair(refinement: int) -> tuple[np.ndarray, np.ndarray]:
    axis = np.arange(refinement + 1, dtype=np.int64)
    x, x_star = np.meshgrid(axis, axis)
    points = np.column_stack((x.ravel(), x_star.ravel()))
    return points.astype(float) / refinement, core_squeeze(points).astype(float) / refinement


def set_square_axes(ax, xlim, ylim) -> None:
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")


def main() -> None:
    if P <= 0 or Q <= 0:
        raise ValueError("P and Q must be positive for squeezeE")
    if len(REFINEMENTS) != 2:
        raise ValueError("this layout expects exactly two refinement levels")

    x_scale = P / Q
    x_star_scale = Q / P
    initial_box = np.array([[0, 0], [1, 0], [1, 1], [0, 1], [0, 0]], dtype=float)
    reference_box = np.array(
        [[0, 0], [x_scale, 0], [x_scale, x_star_scale], [0, x_star_scale], [0, 0]],
        dtype=float,
    )
    lattices = {refinement: lattice_pair(refinement) for refinement in REFINEMENTS}

    max_x = max(1.0, reference_box[:, 0].max(), *(mapped[:, 0].max() for _, mapped in lattices.values()))
    max_y = max(1.0, reference_box[:, 1].max(), *(mapped[:, 1].max() for _, mapped in lattices.values()))
    margin = 0.06 * max(max_x, max_y)
    common_xlim = (-margin, max(max_x, max_y) + margin)
    common_ylim = common_xlim

    fig = plt.figure(figsize=(8.1, 7.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.48, wspace=0.72)
    initial_axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[1, 0])]
    mapped_axes = [fig.add_subplot(grid[0, 1]), fig.add_subplot(grid[1, 1])]

    for initial_ax, mapped_ax, refinement in zip(initial_axes, mapped_axes, REFINEMENTS):
        initial, mapped = lattices[refinement]
        size = POINT_SIZE.get(refinement, 5)

        initial_ax.scatter(
            initial[:, 0],
            initial[:, 1],
            s=size,
            color=LATTICE_COLOR,
            alpha=0.72,
            linewidths=0,
            rasterized=True,
        )
        initial_ax.plot(initial_box[:, 0], initial_box[:, 1], color="0.2", linewidth=1.3)
        initial_ax.text(
            0.05,
            common_ylim[1] - 0.04 * (common_ylim[1] - common_ylim[0]),
            rf"spacing $=1/{refinement}$",
            ha="left",
            va="top",
            fontsize="small",
            bbox={
                "boxstyle": "round,pad=0.2",
                "facecolor": "white",
                "edgecolor": "0.2",
                "linewidth": 0.7,
                "alpha": 0.9,
            },
        )
        initial_ax.set_xlabel(r"$x$", fontsize="large")
        initial_ax.set_ylabel(r"$x^*$", fontsize="large")
        initial_ax.set_xticks([0.0, 1.0])
        initial_ax.set_yticks([0.0, 1.0])
        initial_ax.set_xticklabels([r"$0$", r"$1$"])
        initial_ax.set_yticklabels([r"$0$", r"$1$"])
        set_square_axes(initial_ax, common_xlim, common_ylim)

        mapped_ax.scatter(
            mapped[:, 0],
            mapped[:, 1],
            s=size,
            color=LATTICE_COLOR,
            alpha=0.72,
            linewidths=0,
            rasterized=True,
        )
        mapped_ax.plot(
            reference_box[:, 0],
            reference_box[:, 1],
            "--",
            color=REFERENCE_COLOR,
            linewidth=1.6,
            dashes=(6, 5),
        )
        mapped_ax.set_xlabel(r"$x$", fontsize="large")
        mapped_ax.set_ylabel(r"$x^*$", fontsize="large")
        mapped_ax.yaxis.set_label_position("right")
        mapped_ax.set_xticks([0.0, x_scale])
        mapped_ax.set_yticks([0.0, x_star_scale])
        mapped_ax.set_xticklabels([r"$0$", r"$1\times \frac{p}{q}$"])
        mapped_ax.set_yticklabels([r"$0$", r"$1\div(\frac{p}{q})$"])
        set_square_axes(mapped_ax, common_xlim, common_ylim)

    for initial_ax, mapped_ax in zip(initial_axes, mapped_axes):
        fig.add_artist(
            ConnectionPatch(
                xyA=(1.02, 0.42),
                coordsA=initial_ax.transAxes,
                xyB=(-0.05, 0.42),
                coordsB=mapped_ax.transAxes,
                arrowstyle="->",
                mutation_scale=13,
                linewidth=1.2,
                color="0.25",
                shrinkA=2,
                shrinkB=2,
            )
        )

    arrow_label = r"$S_{p/q}(x,x^*)$"
    fig.text(0.52, 0.775, arrow_label, ha="center", va="center", fontsize=12)
    fig.text(0.52, 0.285, arrow_label, ha="center", va="center", fontsize=12)
    fig.legend(
        [
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="None",
                color=LATTICE_COLOR,
                alpha=0.72,
                markersize=max(4.0, POINT_SIZE.get(REFINEMENTS[0], 5)),
                markeredgewidth=0,
            ),
            Line2D(
                [0],
                [0],
                linestyle=(0, (6, 5)),
                color=REFERENCE_COLOR,
                linewidth=1.6,
            ),
        ],
        ["finite integer-lattice map", "real-valued squeeze reference"],
        loc="lower center",
        ncol=2,
        frameon=False,
        fontsize="medium",
        bbox_to_anchor=(0.5, 0.02),
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.13, top=0.95)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT, dpi=200)
    plt.close(fig)
    print(f"Wrote {OUTPUT}")


if __name__ == "__main__":
    main()
