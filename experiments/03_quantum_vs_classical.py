"""Experiment 03: quantum vs classical energy growth (the main result).

The classical ensemble (p = 0, theta uniform) spreads diffusively forever,
<p^2> ~ D n. The quantum rotor started in m = 0 follows the same growth at
first, then stops: dynamical localization. Two scales are predicted
(PLAN.md note 5): the localization length xi ~ D / hbar^2 momentum states,
which sets both the break time (in kicks) and the saturated energy
<p^2> ~ hbar^2 xi^2 ~ D^2 / hbar^2.

The break time t* is measured as in Izrailev (1990, p. 309): the time from
which the quantum energy stays more than 25% below the classical one.

Run from the repo root:
    python experiments/03_quantum_vs_classical.py           # save figures
    python experiments/03_quantum_vs_classical.py --show    # also open a window
"""

import argparse
import pathlib
import sys
import time

import matplotlib.pyplot as plt
import numpy as np

# "python experiments/<script>.py" puts experiments/ on the import path, not
# the repo root, so add the root to make the kicked_rotor package importable.
REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from kicked_rotor import classical, observables, plotting, quantum  # noqa: E402

FIGURE_STEM = "03_quantum_vs_classical"
K = 5.0
HBARS = (0.25, 0.5, 1.0)  # hbar / 4 pi irrational for all three: no quantum resonance
N = quantum.DEFAULT_N
N_KICKS = 2000
N_PARTICLES = 40_000
SLOPE_START = 200  # classical D from the slope of <p^2> over kicks 200..N_KICKS
BREAK_RATIO = 0.75  # Izrailev's 25% criterion
COLORS = {0.25: "C0", 0.5: "C1", 1.0: "C2"}


def quantum_run(K, hbar, N, n_kicks):
    """<p^2> after each kick and the largest edge probability, one FFT per kick."""
    m = quantum.integer_momenta(N).astype(float)
    outer = np.abs(m) > 3 * N / 8
    p2, edge = np.empty(n_kicks), 0.0
    for n, state in enumerate(quantum.iterate(quantum.momentum_eigenstate(N), K, hbar, n_kicks)):
        prob = np.abs(quantum.to_momentum(state)) ** 2
        p2[n] = hbar**2 * np.sum(m**2 * prob)
        edge = max(edge, prob[outer].sum())
    return p2, edge


def classical_run(K, n_particles, n_kicks):
    """<p^2> after each kick for p0 = 0, theta0 on a uniform grid, and the late-time slope D."""
    theta, p = classical.initial_ensemble(n_particles, sampling="grid")
    p2 = np.empty(n_kicks)
    for n in range(1, n_kicks + 1):
        theta, p = classical.standard_map(theta, p, K)
        p2[n - 1] = np.mean(p**2)
        if n == SLOPE_START:
            early = p**2
    diff = (p**2 - early) / (n_kicks - SLOPE_START)
    return p2, diff.mean(), diff.std() / np.sqrt(n_particles)


def break_time(p2_quantum, p2_classical, ratio):
    """First kick from which the quantum <p^2> stays below ratio * classical."""
    above = np.flatnonzero(p2_quantum >= ratio * p2_classical)
    return int(above[-1]) + 2 if above.size else 1  # kicks are numbered from 1


def make_figure(n, p2_cl, runs, D):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(6.5, 3.3), layout="constrained")
    ax.loglog(n, p2_cl, color="k", lw=1.4, label="classical")
    for hbar, run in runs.items():
        ax.loglog(n, run["p2"], color=COLORS[hbar], lw=1, label=rf"quantum, $\hbar = {hbar:g}$")
        t = run["t_break"]
        ax.plot(t, run["p2"][t - 1], "v", color=COLORS[hbar], ms=5, mec="k", mew=0.5)
    ax.set_xlim(1, n[-1])
    ax.set_xlabel("kicks $n$")
    ax.set_ylabel(r"$\langle p^2 \rangle$")
    ax.set_title(r"(a) $K = 5$: growth, then saturation")
    ax.legend(frameon=False, loc="upper left", fontsize=8)
    ax.text(0.97, 0.03, r"$\blacktriangledown$ break time $t^*$", transform=ax.transAxes, ha="right", fontsize=8)

    x_line = np.logspace(-3, 1.3, 50)
    ax2.loglog(x_line, x_line, ":", color="0.3", lw=1, label=r"diffusion, $\langle m^2\rangle = Dn/\hbar^2$")
    ax2.axhline(1.0, color="0.3", lw=1, ls="--", label=r"$\langle m^2\rangle = \xi^2$")
    for hbar, run in runs.items():
        xi = D / hbar**2
        ax2.loglog(n / xi, run["p2"] / hbar**2 / xi**2, color=COLORS[hbar], lw=1, label=rf"$\hbar = {hbar:g}$")
    ax2.set_xlim(1e-3, 20)
    ax2.set_ylim(1e-3, 3)
    ax2.set_xlabel(r"$n / \xi$")
    ax2.set_ylabel(r"$\langle m^2 \rangle / \xi^2$")
    ax2.set_title(r"(b) Scaled with $\xi = D/\hbar^2$")
    ax2.legend(frameon=False, loc="lower right", fontsize=7.5)
    return fig


def print_summary(p2_cl, D, D_err, runs, paths, elapsed):
    print("Experiment 03: quantum vs classical energy growth")
    print(f"K = {K:g}, quantum: m = 0, N = {N}; classical: p = 0, theta uniform ({N_PARTICLES} particles); {N_KICKS} kicks")
    print()
    print(f"Classical: <p^2> after {N_KICKS} kicks = {p2_cl[-1]:.0f}; D = {D:.2f} +- {D_err:.2f}"
          f" (D/(K^2/2) = {D / (K**2 / 2):.3f})")
    print()
    print("  hbar   xi = D/hbar^2   break time t*   t*/xi   <p^2> at n = 1000, 2000   saturation   D^2/hbar^2   ratio")
    for hbar, run in runs.items():
        xi = D / hbar**2
        saturation = run["p2"][N_KICKS // 2:].mean()
        print(f"  {hbar:4.2f}   {xi:12.1f}   {run['t_break']:13d}   {run['t_break'] / xi:5.2f}"
              f"   {run['p2'][999]:10.1f} {run['p2'][-1]:8.1f}   {saturation:10.1f}   {D**2 / hbar**2:10.1f}"
              f"   {saturation / (D**2 / hbar**2):5.2f}")
    print("  (saturation = mean of <p^2> over the second half of the run; ratio = saturation / (D^2/hbar^2))")
    print(f"  largest edge probability (|m| > 3N/8) in any run: {max(r['edge'] for r in runs.values()):.1e}")
    print()
    print("Izrailev (1990), Fig. 1, has the same case (K = 5, hbar = 0.25, 2000 kicks); read off the")
    print("figure, its energy levels off near <p^2> ~ 1400.")
    print("After the break the quantum energy still swings slowly. In a one-off 40000-kick run at hbar = 0.25,")
    print("<p^2> averaged over blocks of 2500 kicks ranged from 1100 to 1960 around a mean of about 1500")
    print("(17% rms, PLAN.md note 5), so the saturation values above are rough.")
    print()
    print("Saved " + " and ".join(path.relative_to(REPO_ROOT).as_posix() for path in paths))
    print(f"Run time {elapsed:.1f} s")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Experiment 03: quantum vs classical energy growth.")
    parser.add_argument("--show", action="store_true", help="open the figure in a window after saving it")
    args = parser.parse_args(argv)

    plotting.use_backend(args.show)
    plotting.apply_style()
    start = time.perf_counter()

    p2_cl, D, D_err = classical_run(K, N_PARTICLES, N_KICKS)
    runs = {}
    for hbar in HBARS:
        p2, edge = quantum_run(K, hbar, N, N_KICKS)
        runs[hbar] = {"p2": p2, "edge": edge, "t_break": break_time(p2, p2_cl, BREAK_RATIO)}

    n = np.arange(1, N_KICKS + 1)
    fig = make_figure(n, p2_cl, runs, D)
    paths = plotting.save_figure(fig, REPO_ROOT / "figures", FIGURE_STEM)
    print_summary(p2_cl, D, D_err, runs, paths, time.perf_counter() - start)

    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
