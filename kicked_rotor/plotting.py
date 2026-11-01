"""Plotting helpers shared by the experiment scripts.

Typical use in a script:

    plotting.use_backend(show)   # before any figure is made
    plotting.apply_style()
    fig, ax = plt.subplots()
    ...
    plotting.save_figure(fig, figures_dir, "01_phase_portraits")
"""

import pathlib
from fractions import Fraction

import matplotlib
import matplotlib.pyplot as plt
import numpy as np


def use_backend(show):
    """Draw off screen unless the figures will be shown.

    With the non-interactive Agg backend a script never opens a window and
    never waits for one to be closed.
    """
    if not show:
        matplotlib.use("Agg")


def apply_style():
    """Set one consistent look for every figure in the project."""
    plt.rcParams.update(
        {
            "font.family": "STIXGeneral",
            "mathtext.fontset": "stix",
            "font.size": 10,
            "axes.titlesize": 10,
            "axes.labelsize": 10,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "legend.fontsize": 9,
            "axes.linewidth": 0.8,
            "xtick.direction": "in",
            "ytick.direction": "in",
            "xtick.top": True,
            "ytick.right": True,
            "savefig.dpi": 200,
        }
    )


def pi_label(numerator, denominator):
    """TeX label for the angle (numerator / denominator) * pi, such as 3pi/2."""
    fraction = Fraction(numerator, denominator)
    if fraction == 0:
        return "$0$"
    sign = "-" if fraction < 0 else ""
    top = abs(fraction.numerator)
    top_text = r"\pi" if top == 1 else rf"{top}\pi"
    if fraction.denominator == 1:
        return f"${sign}{top_text}$"
    return f"${sign}{top_text}/{fraction.denominator}$"


def set_pi_ticks(axis, numerators, denominator):
    """Put ticks at numerator * pi / denominator on ``axis`` (ax.xaxis or ax.yaxis)."""
    numerators = list(numerators)
    axis.set_ticks([n * np.pi / denominator for n in numerators])
    axis.set_ticklabels([pi_label(n, denominator) for n in numerators])


def format_torus_axes(ax):
    """Axes for one cell of the phase-space torus: theta in [0, 2pi), p in [-pi, pi)."""
    ax.set_xlim(0.0, 2.0 * np.pi)
    ax.set_ylim(-np.pi, np.pi)
    ax.set_aspect("equal")
    set_pi_ticks(ax.xaxis, range(0, 5), 2)
    set_pi_ticks(ax.yaxis, range(-2, 3), 2)


def draw_phase_portrait(
    ax,
    theta,
    p,
    chaotic,
    color_values,
    cmap="viridis",
    color_range=(0.0, 0.85),
    chaotic_color="0.6",
    marker_size=0.6,
):
    """Draw orbits on the torus: chaotic ones in gray, regular ones in color.

    Parameters
    ----------
    theta, p : arrays of shape (n_points, n_orbits)
        Orbit points, with p already folded into [-pi, pi).
    chaotic : bool array of length n_orbits
        Which orbits to draw in ``chaotic_color``. They go underneath.
    color_values : array of length n_orbits with values in [0, 1]
        Picks each regular orbit's color from ``cmap``, rescaled into
        ``color_range`` (the default skips the palest end of viridis, which is
        hard to see on white).

    The dots are rasterized, so a PDF stays small while the axes and labels
    stay as vector graphics.
    """
    colormap = matplotlib.colormaps[cmap]
    chaotic = np.asarray(chaotic, dtype=bool)
    dots = dict(linestyle="none", marker=".", markersize=marker_size, markeredgewidth=0, rasterized=True)
    if chaotic.any():
        ax.plot(theta[:, chaotic].ravel(), p[:, chaotic].ravel(), color=chaotic_color, **dots)
    low, high = color_range
    for j in np.flatnonzero(~chaotic):
        ax.plot(theta[:, j], p[:, j], color=colormap(low + (high - low) * color_values[j]), **dots)


def save_figure(fig, directory, stem):
    """Save ``fig`` as ``<directory>/<stem>.png`` and ``.pdf``; return both paths.

    The PDF gets no creation date, so rerunning a script with unchanged
    results rewrites identical files and git sees no change.
    """
    directory = pathlib.Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    png, pdf = directory / f"{stem}.png", directory / f"{stem}.pdf"
    fig.savefig(png, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight", metadata={"CreationDate": None})
    return [png, pdf]
