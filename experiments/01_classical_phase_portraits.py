"""Experiment 01: phase portraits of the classical kicked rotor.

Draws the standard map

    p_{n+1} = p_n + K sin(theta_n),   theta_{n+1} = theta_n + p_{n+1}   (mod 2 pi)

on one cell of the phase-space torus for four kick strengths:

    K = 0.5    mostly regular: invariant circles (KAM tori) and island chains
    K = 0.97   just below K_c = 0.9716, where the last invariant circle that
               runs all the way around in theta breaks up
    K = 2      a large chaotic sea around the main island
    K = 5      almost all chaotic (the default K of the project)

Each orbit is labelled chaotic or regular by its finite-time Lyapunov
exponent. The script also measures the chaotic area fraction and the Lyapunov
exponent of the chaotic orbits on a grid of starts covering the torus, and
how far the momentum of the portrait orbits gets from its starting value.

Run from the repo root:
    python experiments/01_classical_phase_portraits.py           # save figures
    python experiments/01_classical_phase_portraits.py --show    # also open a window
"""

import argparse
import pathlib
import sys
import time

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

# "python experiments/<script>.py" puts experiments/ on the import path, not
# the repo root, so add the root to make the kicked_rotor package importable.
REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from kicked_rotor import classical, plotting  # noqa: E402

FIGURE_STEM = "01_phase_portraits"
K_VALUES = (0.5, 0.97, 2.0, 5.0)
K_CRITICAL = 0.9716  # the last invariant circle spanning all theta breaks up here

# Portraits: N_PER_LINE starts on each of the four symmetry lines of the map.
# Each orbit is drawn over the same kicks the chaos test measures (kicks
# N_TRANSIENT to N_TRANSIENT + N_WINDOW), so its color matches what is drawn.
N_PER_LINE = 30

# Chaos test. The Lyapunov exponent is measured over N_WINDOW kicks that follow
# N_TRANSIENT uncounted kicks. On a regular orbit a small displacement grows
# only linearly in time, so its log grows by about ln 2 from kick n to kick 2n.
# The displacement also stretches and shrinks as it goes around an island or
# an invariant circle, which adds up to about +-6 to the log. (Only this wobble
# can make lambda negative, so the most negative lambda on the grid measures
# it; the summary prints that value, -0.0030 = -6 / 2000.) So regular orbits
# give |lambda| up to about 0.003. An orbit counts as chaotic if
# lambda > 0.01, that is, if displacements grow by more than e^20 over the
# window. Sticky chaotic orbits that spend the whole window near an island can
# fall below the threshold and count as regular.
N_TRANSIENT = 2000
N_WINDOW = 2000
CHAOS_THRESHOLD = 0.01

# Area statistics use an N_GRID x N_GRID midpoint grid of starts on the torus.
# With 6400 starts, the binomial error of an area fraction f is
# sqrt(f (1 - f) / 6400), at most 0.6%.
N_GRID = 80


def grid_starts(n_grid):
    """Midpoint grid of n_grid * n_grid starts on theta in [0, 2pi), p in [-pi, pi)."""
    centers = (np.arange(n_grid) + 0.5) / n_grid
    theta0, p0 = np.meshgrid(classical.TWO_PI * centers, -np.pi + classical.TWO_PI * centers)
    return theta0.ravel(), p0.ravel()


def portrait_data(K, n_per_line, n_transient, n_window, threshold):
    """Orbits for one panel, each labelled chaotic or regular."""
    theta0, p0 = classical.symmetry_line_starts(n_per_line)
    theta, p = classical.trajectory(theta0, p0, K, n_transient + n_window)
    lam = classical.lyapunov_exponent(theta0, p0, K, n_window, n_transient=n_transient)
    drawn = slice(n_transient, None)  # the kicks the chaos test measured
    return {
        "K": K,
        "theta": theta[drawn],
        "p": classical.wrap_centered(p[drawn]),  # folded for drawing only
        "chaotic": lam > threshold,
        "color_values": (p0 + np.pi) / classical.TWO_PI,  # starting momentum, scaled to [0, 1)
        "max_dp": np.abs(p - p0).max(),  # largest |p_n - p_0| over the whole run
    }


def grid_statistics(K, n_grid, n_transient, n_window, threshold):
    """Chaotic area fraction and Lyapunov exponents on the grid of starts."""
    theta0, p0 = grid_starts(n_grid)
    lam = classical.lyapunov_exponent(theta0, p0, K, n_window, n_transient=n_transient)
    chaotic = lam > threshold
    return {
        "K": K,
        "chaotic_fraction": chaotic.mean(),
        "lambda_chaotic": np.median(lam[chaotic]) if chaotic.any() else np.nan,
        "lambda_min": lam.min(),
    }


def make_figure(portraits, threshold):
    """2 x 2 grid of phase portraits."""
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 7.0), sharex=True, sharey=True, layout="constrained")
    for ax, data in zip(axes.flat, portraits):
        plotting.draw_phase_portrait(ax, data["theta"], data["p"], data["chaotic"], data["color_values"])
        plotting.format_torus_axes(ax)
        ax.set_title(f"$K = {data['K']:g}$")
    for ax in axes[1]:
        ax.set_xlabel(r"$\theta$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$p$ (mod $2\pi$)")
    marker = dict(linestyle="none", marker="o", markersize=4)
    handles = [
        Line2D([], [], color="0.6", label=rf"chaotic orbits ($\lambda > {threshold:g}$)", **marker),
        Line2D([], [], color=matplotlib.colormaps["viridis"](0.4), label="regular orbits (colored by starting $p$)", **marker),
    ]
    fig.legend(handles=handles, loc="outside lower center", ncols=2, frameon=False)
    return fig


def print_summary(portraits, stats, paths, elapsed, n_transient, n_window, n_grid, threshold):
    n_total = n_transient + n_window
    n_orbits = portraits[0]["chaotic"].size
    print("Experiment 01: phase portraits of the standard map")
    print("  p' = p + K sin(theta),  theta' = theta + p'  (mod 2 pi)")
    print()
    print(f"lambda = finite-time Lyapunov exponent over kicks {n_transient} to {n_total}; an orbit")
    print(f"is chaotic if lambda > {threshold:g}. Areas come from an {n_grid} x {n_grid} grid of starts")
    print(f"on the torus; the portraits use {n_orbits} starts on the symmetry lines.")
    print()
    print("     K   chaotic area   lambda of chaotic   largest |p_n - p_0|")
    print("                        orbits (median)     in the portraits")
    for portrait, s in zip(portraits, stats):
        print(
            f"  {s['K']:4.2f}   {100 * s['chaotic_fraction']:10.1f}%"
            f"   {s['lambda_chaotic']:17.2f}   {portrait['max_dp']:19.1f}"
        )
    print()
    lam_floor = min(s["lambda_min"] for s in stats)
    print(f"Threshold check: the most negative lambda on the grids is {lam_floor:.4f}.")
    print("Only the wobble of regular orbits gives negative values, so regular orbits")
    print(f"stay about {threshold / abs(lam_floor):.0f} times below the threshold.")
    for s in stats:
        if s["K"] >= 4:
            print(
                f"K = {s['K']:g}: lambda = {s['lambda_chaotic']:.2f}, compared with the large-K "
                f"estimate ln(K/2) = {np.log(s['K'] / 2):.2f}."
            )
    below = [p["K"] for p in portraits if p["max_dp"] < classical.TWO_PI]
    above = [p["K"] for p in portraits if p["max_dp"] >= classical.TWO_PI]
    if below:
        print(f"|p_n - p_0| stays below 2 pi for K = {' and '.join(f'{K:g}' for K in below)}.")
    if above:
        print(f"It grows far past 2 pi for K = {' and '.join(f'{K:g}' for K in above)}.")
    print(f"Below K_c = {K_CRITICAL}, invariant circles that wrap around the torus act as walls in p.")
    print()
    print("Saved " + " and ".join(path.relative_to(REPO_ROOT).as_posix() for path in paths))
    print(f"Run time {elapsed:.1f} s")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Experiment 01: phase portraits of the standard map.")
    parser.add_argument("--show", action="store_true", help="open the figure in a window after saving it")
    args = parser.parse_args(argv)

    plotting.use_backend(args.show)
    plotting.apply_style()
    start = time.perf_counter()

    portraits = [portrait_data(K, N_PER_LINE, N_TRANSIENT, N_WINDOW, CHAOS_THRESHOLD) for K in K_VALUES]
    stats = [grid_statistics(K, N_GRID, N_TRANSIENT, N_WINDOW, CHAOS_THRESHOLD) for K in K_VALUES]
    fig = make_figure(portraits, CHAOS_THRESHOLD)
    paths = plotting.save_figure(fig, REPO_ROOT / "figures", FIGURE_STEM)
    print_summary(
        portraits, stats, paths, time.perf_counter() - start, N_TRANSIENT, N_WINDOW, N_GRID, CHAOS_THRESHOLD
    )

    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
