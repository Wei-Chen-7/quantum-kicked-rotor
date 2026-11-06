"""Experiment 07 (bonus): how the localization length scales with K and hbar.

Experiment 04 measured xi at K = 5 for three values of hbar. Here the
prediction xi ~ D(K)/hbar^2 (PLAN.md note 5) is tested over a range:

(a), (b) K = 5 and hbar from 0.2 to 1 in steps of 0.02 (41 runs): xi should
    scale like 1/hbar^2, so xi hbar^2 / D should stay roughly constant.
(c) K from 4.5 to 11 at hbar = 0.4, 0.5, 0.6, 0.7, 0.8: xi should follow
    the measured classical diffusion rate D(K), including its
    Bessel-function dip near K = 10, not the naive K^2/2.

Single runs scatter. Over hbar the ratio xi hbar^2 / D has an rms spread of
about 0.2, like the spread in the papers, and it wanders in broad bands:
values of hbar 0.01 apart give strongly correlated ratios (PLAN.md note 5),
so nearby hbar values are not independent samples. Panel (c) therefore takes
the median over hbar values 0.1 apart.

The classical D must be well defined, which rules out some K. K stays away
from the accelerator-mode window (6.28 to 7.45 for l = 1) and its sticky
surroundings (at K = 5.5 the diffusion is already anomalous), and K = 9.25
is left out: there the line p = 0 crosses an island that holds about 9% of
the matched ensemble, so D depends on the ensemble (PLAN.md note 3).

Each run lasts n = max(5000, 50 D_f/hbar^2) kicks, averages |psi_m|^2 over
the second half, and fits Izrailev's shape as in experiment 04. Here D_f is
the five-term formula for D (PLAN.md note 3), used instead of the measured D
so that the quantum runs are identical on every computer. The grid starts at
N >= 40 D_f/hbar^2 states and is doubled until the edge probability stays
below 1e-10 (the limit of validation check 4) and the fitted range lies
inside |m| < 3N/8. A run that still fails at N = 2^16 has not localized on
the grid, as happens right next to a quantum resonance; it is reported and
left out.

Run from the repo root (takes a few minutes):
    python experiments/07_bonus_xi_scaling.py           # save figures
    python experiments/07_bonus_xi_scaling.py --show    # also open a window
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

from kicked_rotor import classical, observables, plotting, quantum  # noqa: E402

FIGURE_STEM = "07_bonus_xi_scaling"
K_FIXED = 5.0
HBAR_SCAN = np.round(np.arange(0.20, 1.0 + 1e-9, 0.02), 2)  # panels (a) and (b)
K_VALUES = (4.5, 5.0, 9.75, 10.5, 11.0)  # panel (c)
HBARS_FOR_K = (0.4, 0.5, 0.6, 0.7, 0.8)  # panel (c): 0.1 apart, so roughly independent
MIN_KICKS = 5000
KICKS_PER_XI = 50  # each run lasts at least 50 D_f/hbar^2 kicks and averages over its second half
N_PARTICLES = 40_000  # classical ensemble, as in experiments 03 and 04
D_WINDOW = (200, 2000)  # classical D from the slope of <p^2> over these kicks, as in experiments 03 and 04
START_GRID_FACTOR = 40  # first grid has N >= 40 D_f/hbar^2 states
N_MAX = 2**16
EDGE_LIMIT = 1e-10  # validation check 4's limit on the probability in |m| > 3N/8
TAIL_LEVEL = 1e-14  # fit out to where the average falls to this fraction of its peak (as in experiment 04)
M_MIN = 5  # leave out the narrow peak at m = 0 (as in experiment 04)
PAPERS = (1.04, 0.20)  # Shepelyansky (1986): mean and rms spread of xi / (D/hbar^2)


def classical_D(K_values, n_particles, d_window):
    """D for each K from the matched ensemble (p = 0, theta on a uniform grid), one row per K."""
    K_col = np.asarray(K_values, dtype=float)[:, None]
    theta0, p0 = classical.initial_ensemble(n_particles, sampling="grid")
    theta, p = np.tile(theta0, (K_col.size, 1)), np.tile(p0, (K_col.size, 1))
    first, last = d_window
    for n in range(1, last + 1):
        theta, p = classical.standard_map(theta, p, K_col)
        if n == first:
            p2_first = p**2
    growth = (p**2 - p2_first) / (last - first)
    return growth.mean(axis=1), growth.std(axis=1) / np.sqrt(n_particles)


def averaged_distribution(K, hbar, N, n_kicks, window):
    """Time average of |psi_m|^2 over kicks (first, last], sorted by m, and the largest edge probability."""
    m = quantum.integer_momenta(N)
    outer = np.abs(m) > 3 * N / 8
    first, last = window
    total, edge = np.zeros(N), 0.0
    for n, state in enumerate(quantum.iterate(quantum.momentum_eigenstate(N), K, hbar, n_kicks), start=1):
        prob = np.abs(quantum.to_momentum(state)) ** 2
        edge = max(edge, prob[outer].sum())
        if first < n <= last:
            total += prob
    return np.fft.fftshift(m), np.fft.fftshift(total) / (last - first), edge


def localization_run(K, hbar, D):
    """One run from m = 0: the fitted xi, or None if the state did not localize on grids up to N_MAX.

    The run length and the first grid are set from the five-term formula for D,
    not from the measured D, because the measured D varies slightly from one
    computer to another (README). This way the quantum runs are the same
    everywhere; the measured D only enters the ratio xi hbar^2 / D.
    """
    xi_guess = five_term_formula(K) * (K**2 / 2) / hbar**2
    n_kicks = max(MIN_KICKS, int(round(KICKS_PER_XI * xi_guess)))
    N = max(2**11, 2 ** int(np.ceil(np.log2(START_GRID_FACTOR * xi_guess))))
    while True:
        m, avg, edge = averaged_distribution(K, hbar, N, n_kicks, (n_kicks // 2, n_kicks))
        m_tail = np.abs(m[avg > TAIL_LEVEL * avg.max()]).max()
        if edge < EDGE_LIMIT and m_tail < 3 * N // 8:
            xi, _ = observables.fit_averaged_localization_length(m, avg, M_MIN, m_tail)
            return {"K": K, "hbar": hbar, "xi": xi, "ratio": xi * hbar**2 / D, "N": N, "edge": edge, "kicks": n_kicks}
        if N >= N_MAX:
            return {"K": K, "hbar": hbar, "xi": None, "ratio": None, "N": N, "edge": edge, "kicks": n_kicks}
        N *= 2


def five_term_formula(K):
    """D/(K^2/2) from the first five terms of Rechester, Rosenbluth and White (1981), Eq. 26 (PLAN.md note 3)."""
    return 1 - 2 * jv(2, K) - 2 * jv(1, K) ** 2 + 2 * jv(2, K) ** 2 + 2 * jv(3, K) ** 2


def good(runs):
    return [run for run in runs if run["xi"] is not None]


def make_figure(scan, k_runs, D):
    fig = plt.figure(figsize=(6.5, 4.2), layout="constrained")
    grid = fig.add_gridspec(2, 2, height_ratios=(1.5, 1), width_ratios=(1, 1.05))
    ax_a = fig.add_subplot(grid[0, 0])
    ax_b = fig.add_subplot(grid[1, 0], sharex=ax_a)
    ax_c = fig.add_subplot(grid[:, 1])
    mean_ratio, spread = PAPERS
    ok = good(scan)
    h = np.array([run["hbar"] for run in ok])
    xi = np.array([run["xi"] for run in ok])
    ratio = np.array([run["ratio"] for run in ok])
    failed = [run["hbar"] for run in scan if run["xi"] is None]

    # (a) xi against hbar at K = 5
    h_line = np.geomspace(0.19, 1.05, 100)
    ax_a.fill_between(h_line, (mean_ratio - spread) * D[K_FIXED] / h_line**2,
                      (mean_ratio + spread) * D[K_FIXED] / h_line**2, color="0.88", lw=0,
                      label=rf"papers: $({mean_ratio} \pm {spread})\,D/\hbar^2$")
    ax_a.loglog(h_line, D[K_FIXED] / h_line**2, "k-", lw=0.9, label=r"$\xi = D/\hbar^2$")
    ax_a.loglog(h, xi, "o", color="C0", ms=3, label="single runs")
    ax_a.set_ylabel(r"$\xi$")
    ax_a.set_title(r"(a) $K = 5$: $\xi$ against $\hbar$")
    ax_a.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax_a.tick_params(labelbottom=False)

    # (b) the ratio xi hbar^2 / D against hbar
    ax_b.axhspan(mean_ratio - spread, mean_ratio + spread, color="0.88", lw=0)
    ax_b.semilogx(h, ratio, "o-", color="C0", ms=3, lw=0.6)
    ax_b.axhline(np.median(ratio), color="k", ls="--", lw=0.8, label=f"median {np.median(ratio):.2f}")
    if failed:
        ax_b.plot(failed, np.full(len(failed), 2.05), "x", color="C3", ms=5, label="not localized on the grid")
    ax_b.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax_b.set_ylim(0.6, 2.2)
    ax_b.set_xlim(0.19, 1.05)
    ticks = (0.2, 0.3, 0.4, 0.5, 0.7, 1.0)
    ax_b.set_xticks(ticks, [f"{t:g}" for t in ticks])
    ax_b.set_xticks([], minor=True)
    ax_b.set_xlabel(r"$\hbar$")
    ax_b.set_ylabel(r"$\xi\hbar^2/D$")
    ax_b.set_title(r"(b) the ratio $\xi\hbar^2/D$", fontsize=9)

    # (c) xi hbar^2 / (K^2/2) against K
    K_line = np.linspace(4.2, 11.5, 300)
    window_lo, window_hi = 2 * np.pi, np.sqrt(4 * np.pi**2 + 16)
    ax_c.axvspan(window_lo, window_hi, color="0.88", lw=0)
    ax_c.text(0.5 * (window_lo + window_hi), 0.43, "accelerator\nmodes", ha="center", va="bottom", fontsize=7,
              color="0.35")
    ax_c.plot(K_line, five_term_formula(K_line), "--", color="0.5", lw=0.8, label="five-term formula")
    ax_c.plot(K_VALUES, [D[K] / (K**2 / 2) for K in K_VALUES], "s", mfc="none", mec="k", ms=6, mew=0.9,
              label=r"classical $D/(K^2/2)$")
    for i, K in enumerate(K_VALUES):
        scaled = np.array([run["xi"] * run["hbar"] ** 2 / (K**2 / 2) for run in good(k_runs[K])])
        ax_c.plot(np.full(scaled.size, K), scaled, ".", color="0.45", ms=4, zorder=3,
                  label="quantum, single runs" if i == 0 else None)
        ax_c.plot(K, np.median(scaled), "_", color="C1", ms=16, mew=2.2, zorder=2,
                  label="quantum, median" if i == 0 else None)
    ax_c.set_xlim(4.2, 11.5)
    ax_c.set_ylim(0.4, 2.3)  # room for the legend above the formula's peak (1.85 at K = 6.9)
    ax_c.set_xlabel("$K$")
    ax_c.set_ylabel(r"$\xi\hbar^2/(K^2/2)$  and  $D/(K^2/2)$")
    ax_c.set_title(r"(c) $\hbar = 0.4$ to $0.8$: $\xi$ follows $D(K)$")
    ax_c.legend(frameon=False, fontsize=7.5, loc="upper center", ncols=2, columnspacing=1.2, handlelength=1.6)
    return fig


def print_summary(scan, k_runs, D, D_err, paths, elapsed):
    print("Experiment 07 (bonus): scaling of the localization length with K and hbar")
    print(f"Each run: start in m = 0, n = max({MIN_KICKS}, {KICKS_PER_XI} D_f/hbar^2) kicks (D_f: five-term formula),"
          " Izrailev-shape fit to the")
    print(f"time average over the second half (as in experiment 04); grid doubled until the edge probability < {EDGE_LIMIT:g}.")
    print()
    ok = good(scan)
    ratio = np.array([run["ratio"] for run in ok])
    print(f"(a), (b) K = {K_FIXED:g}, classical D = {D[K_FIXED]:.2f} +- {D_err[K_FIXED]:.2f};"
          f" {len(scan)} runs at hbar = {HBAR_SCAN[0]:g}, {HBAR_SCAN[1]:g}, ..., {HBAR_SCAN[-1]:g}")
    print(f"  xi hbar^2/D: median {np.median(ratio):.2f}, mean {ratio.mean():.2f}, rms spread {ratio.std(ddof=1):.2f},"
          f" range {ratio.min():.2f} to {ratio.max():.2f}")
    print(f"  (Shepelyansky 1986: mean {PAPERS[0]}, rms spread {PAPERS[1]})")
    for lo, hi in ((0.2, 0.3), (0.3, 0.5), (0.5, 0.75), (0.75, 1.01)):
        part = [run["ratio"] for run in ok if lo <= run["hbar"] < hi]
        print(f"  median for {lo:g} <= hbar < {min(hi, 1.0):g}: {np.median(part):.2f} ({len(part)} runs)")
    pairs = [f"{run['hbar']:.2f}: {run['ratio']:.2f}" for run in ok]
    for i in range(0, len(pairs), 7):
        print("    " + "   ".join(pairs[i:i + 7]))
    failed = [run for run in scan if run["xi"] is None]
    print(f"  runs that did not localize on grids up to N = {N_MAX}: "
          + (", ".join(f"hbar = {run['hbar']:g}" for run in failed) or "none"))
    print()
    print(f"(c) K dependence, medians over hbar = {', '.join(f'{h:g}' for h in HBARS_FOR_K)}")
    print("      K   classical D/(K^2/2)   five-term formula   xi hbar^2/(K^2/2) (median, range)   xi hbar^2/D (median)")
    for K in K_VALUES:
        runs = good(k_runs[K])
        scaled = np.array([run["xi"] * run["hbar"] ** 2 / (K**2 / 2) for run in runs])
        ratios = np.array([run["ratio"] for run in runs])
        print(f"  {K:5.2f}   {D[K] / (K**2 / 2):.3f} +- {D_err[K] / (K**2 / 2):.3f}         {five_term_formula(K):.3f}"
              f"            {np.median(scaled):.2f}  ({scaled.min():.2f} to {scaled.max():.2f})"
              f"                {np.median(ratios):.2f}" + ("" if len(runs) == len(k_runs[K]) else "  (some runs left out)"))
    all_runs = scan + [run for runs in k_runs.values() for run in runs]
    print()
    print(f"Grids used: N = {min(run['N'] for run in all_runs)} to {max(run['N'] for run in all_runs)};"
          f" largest edge probability in the runs kept: {max(run['edge'] for run in good(all_runs)):.1e}")
    print()
    print("Saved " + " and ".join(path.relative_to(REPO_ROOT).as_posix() for path in paths))
    print(f"Run time {elapsed:.0f} s")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Experiment 07 (bonus): scaling of xi with K and hbar.")
    parser.add_argument("--show", action="store_true", help="open the figure in a window after saving it")
    args = parser.parse_args(argv)

    plotting.use_backend(args.show)
    plotting.apply_style()
    start = time.perf_counter()

    all_K = sorted(set(K_VALUES) | {K_FIXED})
    D_values, D_errors = classical_D(all_K, N_PARTICLES, D_WINDOW)
    D, D_err = dict(zip(all_K, D_values)), dict(zip(all_K, D_errors))

    scan = [localization_run(K_FIXED, float(hbar), D[K_FIXED]) for hbar in HBAR_SCAN]
    print(f"  done: K = {K_FIXED:g}, {len(scan)} values of hbar  [{time.perf_counter() - start:.0f} s]", flush=True)
    from_scan = {round(run["hbar"], 2): run for run in scan}
    k_runs = {}
    for K in K_VALUES:
        k_runs[K] = [from_scan[round(h, 2)] if K == K_FIXED else localization_run(K, h, D[K]) for h in HBARS_FOR_K]
        print(f"  done: K = {K:g}  [{time.perf_counter() - start:.0f} s]", flush=True)
    print()

    fig = make_figure(scan, k_runs, D)
    paths = plotting.save_figure(fig, REPO_ROOT / "figures", FIGURE_STEM)
    print_summary(scan, k_runs, D, D_err, paths, time.perf_counter() - start)

    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
