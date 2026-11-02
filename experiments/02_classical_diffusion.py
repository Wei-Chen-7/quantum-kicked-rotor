"""Experiment 02: classical momentum diffusion and its Bessel-function corrections.

For K well above the critical value, the momentum of the standard map makes a
random walk, <p_n^2> ~ D n. If the kick phases were independent from kick to
kick, D would be the quasilinear value K^2/2 (Chirikov 1979, Eq. 5.46).
Correlations between kicks change it. With C_j = <sin(theta_0) sin(theta_j)>,

    D = K^2 [C_0 + 2 sum_{j >= 1} C_j],   C_0 = 1/2, C_1 = 0,
    C_2 = -J_2(K)/2,   C_3 = (J_3(K)^2 - J_1(K)^2)/2,  ...

and summing such terms gives the corrected formulas (PLAN.md note 3):

    Rechester and White (1980):                D/(K^2/2) = 1 - 2J2 - 2J1^2 + 2J3^2
    Rechester, Rosenbluth and White (1981):    the same + 2 J2^2 (first five terms of their Eq. 26)

The figure shows (a) <p^2> against n for the project's ensemble (p = 0,
theta uniform) at K = 5 and 10, (b) the measured D against K together with
the formulas, and (c) the lag correlations C_j at K = 22 and 22.3. Away from
accelerator modes the formulas miss by more than the noise at a few values
of K (near 9 to 11, and at 22); panel (c) shows that at K = 22 the miss
comes from lag 4 and beyond, which the formulas only partly include.

Run from the repo root:
    python experiments/02_classical_diffusion.py           # save figures
    python experiments/02_classical_diffusion.py --show    # also open a window
"""

import argparse
import pathlib
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.special import jv

# "python experiments/<script>.py" puts experiments/ on the import path, not
# the repo root, so add the root to make the kicked_rotor package importable.
REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from kicked_rotor import classical, plotting  # noqa: E402

FIGURE_STEM = "02_classical_diffusion"
SEED = 2026

# (a) <p^2> against n for the project's initial ensemble
GROWTH_K = (5.0, 10.0)
GROWTH_PARTICLES = 40_000
GROWTH_KICKS = 1000
SLOPE_START = 200  # fit D on kicks 200..1000, after the start-up correlations have died out

# (b) D against K, measured with the torus ensemble from the slope of the
# mean squared displacement between kicks SCAN_N1 and SCAN_N2
SCAN_K = np.arange(4.5, 30.0 + 1e-9, 0.25)
SCAN_PARTICLES = 10_000
SCAN_N1, SCAN_N2 = 100, 400

# (c) lag correlations C_j at K = 22.0 and a neighbor
LAG_K = (22.0, 22.3)
LAG_PARTICLES = 400_000
LAG_MAX = 16
LAG_ORIGINS = 4  # average each C_j over this many starting kicks


def quasilinear_ratio_rw1980(K):
    """Rechester and White (1980), Eq. (20), with sigma -> 0, as D / (K^2/2)."""
    return 1 - 2 * jv(2, K) - 2 * jv(1, K) ** 2 + 2 * jv(3, K) ** 2


def quasilinear_ratio_rrw1981(K):
    """First five terms of Rechester, Rosenbluth and White (1981), Eq. (26), as D / (K^2/2)."""
    return quasilinear_ratio_rw1980(K) + 2 * jv(2, K) ** 2


def accelerator_window(K, margin=0.0):
    """True where K lies in an accelerator-mode window 2 pi l <= K <= sqrt((2 pi l)^2 + 16)."""
    K = np.asarray(K, dtype=float)
    inside = np.zeros(K.shape, dtype=bool)
    for l in range(1, 6):
        low, high = 2 * np.pi * l, np.sqrt((2 * np.pi * l) ** 2 + 16)
        inside |= (K >= low - margin) & (K <= high + margin)
    return inside


def growth_curves(K_values, n_particles, n_kicks):
    """<p_n^2> for the ensemble p0 = 0, theta0 uniform (midpoint grid), and the late-time slope D."""
    results = {}
    for K in K_values:
        theta, p = classical.initial_ensemble(n_particles, sampling="grid")
        p2 = np.empty(n_kicks)
        squares_at = {}
        for n in range(1, n_kicks + 1):
            theta, p = classical.standard_map(theta, p, K)
            p2[n - 1] = np.mean(p**2)
            if n in (SLOPE_START, n_kicks):
                squares_at[n] = p**2
        # Slope between the two kicks, with its statistical error from the
        # spread over particles.
        diff = (squares_at[n_kicks] - squares_at[SLOPE_START]) / (n_kicks - SLOPE_START)
        results[K] = {"p2": p2, "D": diff.mean(), "D_err": diff.std() / np.sqrt(n_particles)}
    return results


def diffusion_scan(K_values, n_particles, n1, n2, seed):
    """D(K) from the torus ensemble, all K at once (one row per K)."""
    K_col = np.asarray(K_values, dtype=float)[:, None]
    theta, p0 = classical.torus_ensemble(K_col.size * n_particles, seed)
    theta, p0 = theta.reshape(K_col.size, n_particles), p0.reshape(K_col.size, n_particles)
    p = p0.copy()
    for n in range(1, n2 + 1):
        theta, p = classical.standard_map(theta, p, K_col)
        if n == n1:
            d1 = (p - p0) ** 2
    diff = ((p - p0) ** 2 - d1) / (n2 - n1)
    return diff.mean(axis=1), diff.std(axis=1) / np.sqrt(n_particles)


def lag_correlations(K, n_particles, max_lag, n_origins, seed):
    """C_j = <sin(theta_t) sin(theta_{t+j})> on the torus ensemble, averaged over n_origins start times t."""
    theta, p = classical.torus_ensemble(n_particles, seed)
    sines = [np.sin(theta)]
    for _ in range(max_lag + n_origins):
        theta, p = classical.standard_map(theta, p, K)
        sines.append(np.sin(theta))
    sines = np.array(sines)
    C = np.array([np.mean([np.mean(sines[t] * sines[t + j]) for t in range(n_origins)]) for j in range(max_lag + 1)])
    # Each product has variance about 1/4, so the mean of n_particles * n_origins of them has this error.
    return C, 0.5 / np.sqrt(n_particles * n_origins)


def make_figure(growth, scan_K, D_ratio, D_ratio_err, lags):
    fig = plt.figure(figsize=(6.5, 6.4), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(1, 1.05))
    ax_a, ax_c = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])
    ax_b = fig.add_subplot(grid[1, :])

    for (K, data), color in zip(growth.items(), ("C0", "C1")):
        n = np.arange(1, data["p2"].size + 1)
        ax_a.plot(n, data["p2"] / (K**2 / 2), color=color, lw=1.2, label=rf"$K = {K:g}$")
        slope = data["D"] / (K**2 / 2)
        fit_n = n[SLOPE_START - 1:]
        offset = data["p2"][-1] / (K**2 / 2) - slope * n[-1]
        ax_a.plot(fit_n, slope * fit_n + offset, "--", color="0.2", lw=0.8)
        above = K == GROWTH_K[0]
        ax_a.annotate(rf"$D/(K^2/2) = {slope:.2f}$", xy=(650, slope * 650 + offset),
                      xytext=(-4, 6) if above else (4, -6), textcoords="offset points",
                      ha="right" if above else "left", va="bottom" if above else "top", fontsize=8)
    ax_a.plot([0, GROWTH_KICKS], [0, GROWTH_KICKS], ":", color="0.4", lw=1, label=r"$K^2 n/2$")
    ax_a.set_xlabel("kicks $n$")
    ax_a.set_ylabel(r"$\langle p^2\rangle / (K^2/2)$")
    ax_a.set_title(r"(a) Growth from $p = 0$")
    ax_a.legend(frameon=False, loc="upper left")
    ax_a.set_xlim(0, GROWTH_KICKS)
    ax_a.set_ylim(0, None)

    fine_K = np.linspace(scan_K[0], scan_K[-1], 2000)
    for l in range(1, 5):
        low, high = 2 * np.pi * l, np.sqrt((2 * np.pi * l) ** 2 + 16)
        if low < scan_K[-1]:
            ax_b.axvspan(low, high, color="0.85", lw=0)
    ax_b.axhline(1.0, color="0.4", lw=0.8, ls=":")
    ax_b.plot(fine_K, quasilinear_ratio_rw1980(fine_K), color="C2", lw=0.9, ls="--", label="Rechester and White (1980)")
    ax_b.plot(fine_K, quasilinear_ratio_rrw1981(fine_K), color="C3", lw=1.2,
              label="Rechester, Rosenbluth and White (1981)")
    window = accelerator_window(scan_K)
    ax_b.errorbar(scan_K[~window], D_ratio[~window], yerr=D_ratio_err[~window], fmt="o", ms=2.5, color="k",
                  elinewidth=0.6, capsize=0, label="simulation")
    top = 2.4
    clipped = window | (D_ratio > top)
    ax_b.plot(scan_K[clipped], np.minimum(D_ratio[clipped], top - 0.05), "^", ms=3.5, mfc="none", color="0.3",
              label="accelerator modes (shaded; off scale if at top)")
    ax_b.set_ylim(0.3, top)
    ax_b.set_xlim(scan_K[0], scan_K[-1])
    ax_b.set_xlabel("$K$")
    ax_b.set_ylabel(r"$D / (K^2/2)$")
    ax_b.set_title(r"(b) Diffusion rate against $K$")
    ax_b.legend(frameon=False, loc="upper right", fontsize=8, ncols=2)

    for (K, (C, err)), color, marker in zip(lags.items(), ("C3", "C0"), ("o", "s")):
        j = np.arange(1, C.size)
        ax_c.plot(j, C[1:], marker=marker, ms=3, lw=0.8, color=color, label=rf"$K = {K:g}$")
    for K, color in zip(lags, ("C3", "C0")):
        ax_c.plot([2, 3], [-jv(2, K) / 2, (jv(3, K) ** 2 - jv(1, K) ** 2) / 2], "x", color=color, ms=6)
    ax_c.axhline(0, color="0.4", lw=0.6)
    ax_c.set_xlabel("lag $j$")
    ax_c.set_ylabel(r"$C_j = \langle \sin\theta_0 \sin\theta_j\rangle$")
    ax_c.set_title(r"(c) Correlations near $K = 22$")
    ax_c.legend(frameon=False, loc="lower right", title=r"$\times$: Bessel formulas", title_fontsize=8)
    return fig


def print_summary(growth, scan_K, D_ratio, D_ratio_err, lags, paths, elapsed):
    print("Experiment 02: classical diffusion and its Bessel-function corrections")
    print()
    print(f"(a) Ensemble p = 0, theta uniform ({GROWTH_PARTICLES} particles), slope of <p^2> over kicks"
          f" {SLOPE_START} to {GROWTH_KICKS}:")
    for K, data in growth.items():
        print(f"  K = {K:g}: D/(K^2/2) = {data['D'] / (K**2 / 2):.3f} +- {data['D_err'] / (K**2 / 2):.3f};"
              f" formula (RRW 1981, five terms) {quasilinear_ratio_rrw1981(K):.3f}")
    clear = ~accelerator_window(scan_K, margin=1.2)
    resid = D_ratio - quasilinear_ratio_rrw1981(scan_K)
    big = clear & (np.abs(resid) > 4 * D_ratio_err)
    print()
    print(f"(b) Torus ensemble ({SCAN_PARTICLES} particles per K), slope over kicks {SCAN_N1} to {SCAN_N2},"
          f" {scan_K.size} values of K:")
    print(f"  typical statistical error of D/(K^2/2): {np.median(D_ratio_err):.3f}")
    print(f"  away from accelerator windows (more than 1.2 from their edges), the five-term formula misses by"
          f" {np.sqrt(np.mean(resid[clear] ** 2)):.3f} rms;")
    print("  larger than 4 errors at K = " + (", ".join(f"{K:g} ({r:+.3f})" for K, r in zip(scan_K[big], resid[big]))
                                              or "none"))
    print("  (misses this large are very unlikely to be noise; panel (c) shows where they come from. Several sit")
    print("  near the 4-error line, so this list can change from run to run and from computer to computer)")
    near = accelerator_window(scan_K, margin=1.2) & ~accelerator_window(scan_K)
    print(f"  just outside the windows D is enhanced: largest D/(K^2/2) there = {D_ratio[near].max():.2f}"
          f" at K = {scan_K[near][np.argmax(D_ratio[near])]:g}")
    print()
    print("(c) Lag correlations C_j (torus ensemble) and the running sum 1 + 4 sum_{i<=j} C_i = D/(K^2/2):")
    for K, (C, err) in lags.items():
        running = 1 + 4 * np.cumsum(C[1:])
        print(f"  K = {K:g}: C2 = {C[2]:+.4f} (formula {-jv(2, K) / 2:+.4f}), C3 = {C[3]:+.4f}"
              f" (formula {(jv(3, K) ** 2 - jv(1, K) ** 2) / 2:+.4f}), C4 = {C[4]:+.4f}, error {err:.4f}")
        print("     running sum at j = 3, 4, 8, 16: " + ", ".join(f"{running[j - 1]:.3f}" for j in (3, 4, 8, 16)))
    print("  The Bessel formulas for C2 and C3 hold. At K = 22.0 the extra diffusion comes from lag 4 and a")
    print("  slowly decaying tail with a period-4 sign pattern. The five-term formula has lags up to 3 exactly")
    print("  but only the leading part of lag 4 (its J2^2 term), so it misses wherever lag 4 and beyond are")
    print("  large (the same pattern appears near K = 9.25).")
    print()
    print("Saved " + " and ".join(path.relative_to(REPO_ROOT).as_posix() for path in paths))
    print(f"Run time {elapsed:.1f} s")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Experiment 02: classical diffusion and Bessel corrections.")
    parser.add_argument("--show", action="store_true", help="open the figure in a window after saving it")
    args = parser.parse_args(argv)

    plotting.use_backend(args.show)
    plotting.apply_style()
    start = time.perf_counter()

    growth = growth_curves(GROWTH_K, GROWTH_PARTICLES, GROWTH_KICKS)
    D, D_err = diffusion_scan(SCAN_K, SCAN_PARTICLES, SCAN_N1, SCAN_N2, SEED)
    D_ratio, D_ratio_err = D / (SCAN_K**2 / 2), D_err / (SCAN_K**2 / 2)
    lags = {K: lag_correlations(K, LAG_PARTICLES, LAG_MAX, LAG_ORIGINS, SEED + i) for i, K in enumerate(LAG_K)}

    fig = make_figure(growth, SCAN_K, D_ratio, D_ratio_err, lags)
    paths = plotting.save_figure(fig, REPO_ROOT / "figures", FIGURE_STEM)
    print_summary(growth, SCAN_K, D_ratio, D_ratio_err, lags, paths, time.perf_counter() - start)

    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
