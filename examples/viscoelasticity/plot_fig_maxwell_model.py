"""Schematic of the one-dimensional Maxwell element (paper Fig. 1, fig:maxwell_model).

An elastic spring (modulus E) and a viscous dashpot (viscosity eta) in series, loaded by
the stress sigma at both ends. Dimension lines mark the total strain epsilon and its split
into the elastic strain epsilon^e of the spring and the viscous strain epsilon^v of the
dashpot.

Output
  maxwell_model.{pdf,png}
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

THIS_DIR = Path(__file__).resolve().parent

FIG_SIZE = (6.0, 2.35)  # aspect of the drawing limits set in main()
# Printed at 0.5\textwidth (6.5 in) in paper_1.tex; an 8.6 pt printed base size keeps
# superscripts at 6 pt or more (Elsevier artwork sizing).
FONT_SIZE = 8.6 * FIG_SIZE[0] / (0.5 * 6.5)

plt.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm",
        "font.size": FONT_SIZE,
        "pdf.fonttype": 42,
    }
)

ELEMENT_COLOR = "k"
ANNOTATION_COLOR = "0.45"
ELEMENT_WIDTH = 2.4    # about 1.3 pt at printed size
ANNOTATION_WIDTH = 1.0

X_LEFT, X_MID, X_RIGHT = 0.0, 4.0, 8.0   # nodes: left end, spring--dashpot junction, right end
Y_DIM_PARTS = 1.25                       # dimension line for epsilon^e and epsilon^v
Y_DIM_TOTAL = 2.15                       # dimension line for epsilon


def draw_spring(ax, x0, x1, teeth=4, lead=0.55, height=0.5):
    zig_x = np.linspace(x0 + lead, x1 - lead, 2 * teeth + 1)
    zig_y = np.zeros_like(zig_x)
    zig_y[1:-1] = height * (-1.0) ** np.arange(2 * teeth - 1)  # +h, -h, ..., +h between the leads
    x = np.concatenate([[x0], zig_x, [x1]])
    y = np.concatenate([[0.0], zig_y, [0.0]])
    ax.plot(x, y, color=ELEMENT_COLOR, linewidth=ELEMENT_WIDTH, solid_joinstyle="miter")


def draw_dashpot(ax, x0, x1, half_height=0.55, cylinder=(1.25, 2.85), piston=2.1, piston_half=0.4):
    c0, c1 = x0 + cylinder[0], x0 + cylinder[1]
    p = x0 + piston
    kw = dict(color=ELEMENT_COLOR, linewidth=ELEMENT_WIDTH, solid_capstyle="butt")
    ax.plot([x0, c0], [0.0, 0.0], **kw)
    ax.plot([c1, c0, c0, c1], [half_height, half_height, -half_height, -half_height], **kw)
    ax.plot([p, p], [-piston_half, piston_half], **kw)
    ax.plot([p, x1], [0.0, 0.0], **kw)


def draw_dimension(ax, x0, x1, y, label, tick=0.18):
    kw = dict(color=ANNOTATION_COLOR, linewidth=ANNOTATION_WIDTH)
    ax.plot([x0, x1], [y, y], **kw)
    for x in (x0, x1):
        ax.plot([x, x], [y - tick, y + tick], **kw)
    ax.text(0.5 * (x0 + x1), y + 0.12, label, color=ANNOTATION_COLOR, fontsize="x-large",
            ha="center", va="bottom")


def draw_load(ax, x, direction, label, length=0.55, gap=0.22):
    start = x + direction * gap
    end = start + direction * length
    ax.annotate("", xy=(end, 0.0), xytext=(start, 0.0),
                arrowprops=dict(arrowstyle="-|>", color=ANNOTATION_COLOR, linewidth=2.0,
                                mutation_scale=18, shrinkA=0, shrinkB=0))
    ax.text(end + direction * 0.15, 0.0, label, color=ANNOTATION_COLOR, fontsize="x-large",
            ha="left" if direction > 0 else "right", va="center")


def main():
    fig, ax = plt.subplots(figsize=FIG_SIZE)

    draw_spring(ax, X_LEFT, X_MID)
    draw_dashpot(ax, X_MID, X_RIGHT)
    ax.plot([X_LEFT, X_MID, X_RIGHT], [0.0, 0.0, 0.0], linestyle="none", marker="o", markersize=7,
            markerfacecolor="white", markeredgecolor=ELEMENT_COLOR, markeredgewidth=2.0, zorder=3)

    for x in (X_LEFT, X_MID, X_RIGHT):
        top = Y_DIM_TOTAL if x != X_MID else Y_DIM_PARTS
        ax.plot([x, x], [0.35, top + 0.18], color=ANNOTATION_COLOR, linewidth=0.6)
    draw_dimension(ax, X_LEFT, X_MID, Y_DIM_PARTS, r"$\varepsilon^e$")
    draw_dimension(ax, X_MID, X_RIGHT, Y_DIM_PARTS, r"$\varepsilon^v$")
    draw_dimension(ax, X_LEFT, X_RIGHT, Y_DIM_TOTAL, r"$\varepsilon$")

    draw_load(ax, X_LEFT, -1.0, r"$\sigma$")
    draw_load(ax, X_RIGHT, 1.0, r"$\sigma$")

    ax.text(0.5 * (X_LEFT + X_MID), -0.75, r"$E$", fontsize="x-large", ha="center", va="top")
    ax.text(X_MID + 2.1, -0.75, r"$\eta$", fontsize="x-large", ha="center", va="top")

    ax.set_aspect("equal")
    ax.set_xlim(X_LEFT - 1.45, X_RIGHT + 1.45)
    ax.set_ylim(-1.35, Y_DIM_TOTAL + 0.75)
    ax.axis("off")
    fig.subplots_adjust(left=0.0, right=1.0, bottom=0.0, top=1.0)
    for ext in ("pdf", "png"):
        fig.savefig(THIS_DIR / f"maxwell_model.{ext}", dpi=200)


if __name__ == "__main__":
    main()
