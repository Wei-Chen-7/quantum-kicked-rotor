"""Experiment 05: quantum resonance and antiresonance (validation check 2).

At hbar = 4 pi the free step exp(-i hbar m^2 / 2) = exp(-2 pi i m^2) is the
identity, so the kicks add up coherently: after n kicks the state is
exp(-i n K cos(theta) / hbar)|0>, its momentum amplitudes are Bessel functions,
c_m = (-i)^m J_m(n K / hbar), and <p^2>_n = K^2 n^2 / 2 grows quadratically
forever. At hbar = 2 pi the free step is (-1)^m, a shift of theta by pi, so
every second kick undoes the one before and <p^2> alternates K^2/2, 0, K^2/2, ...
(Casati et al. 1979 call this case anti-resonant).

Run from the repo root:
    python experiments/05_resonance_test.py           # save figures
    python experiments/05_resonance_test.py --show    # also open a window
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

from kicked_rotor import observables, plotting, quantum  # noqa: E402

FIGURE_STEM = "05_resonance_test"
EPS = np.finfo(float).eps
K = 5.0
N = quantum.DEFAULT_N
N_KICKS = 1000  # resonance run
N_DISTRIBUTION = 100  # kicks before the momentum distribution is compared with J_m^2
N_ANTI = 20  # antiresonance run


def resonance_run(K, N, n_kicks):
    """<p^2> after each kick at hbar = 4 pi."""
    hbar = 4 * np.pi
    psi0 = quantum.momentum_eigenstate(N)
    return np.array([observables.mean_p2(s, hbar) for s in quantum.iterate(psi0, K, hbar, n_kicks)])


def resonance_distribution(K, N, n_kicks):
    """Momentum distribution after n_kicks at hbar = 4 pi, and the exact J_m(nK/hbar)^2."""
    hbar = 4 * np.pi
    state = quantum.evolve(quantum.momentum_eigenstate(N), K, hbar, n_kicks)
    m, prob = observables.momentum_distribution(state)
    return m, prob, jv(m, n_kicks * K / hbar) ** 2


def antiresonance_run(K, N, n_kicks):
    """<p^2> after each kick at hbar = 2 pi, and the distance of each state from |0>."""
    hbar = 2 * np.pi
    psi0 = quantum.momentum_eigenstate(N)
    p2, distance = [], []
    for state in quantum.iterate(psi0, K, hbar, n_kicks):
        p2.append(observables.mean_p2(state, hbar))
        distance.append(np.linalg.norm(state - psi0))
    return np.array(p2), np.array(distance)


def make_figure(p2, dist_data, p2_anti):
    fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.6), layout="constrained")
    n = np.arange(1, p2.size + 1)
    exact = K**2 * n**2 / 2

    ax = axes[0, 0]
    show = np.unique(np.round(np.logspace(0, np.log10(n[-1]), 30)).astype(int))
    ax.loglog(n, exact, color="0.3", lw=1, label=r"exact $K^2 n^2/2$")
    ax.loglog(show, p2[show - 1], "o", ms=3.5, color="C0", label="simulation")
    ax.set_xlabel("kicks $n$")
    ax.set_ylabel(r"$\langle p^2 \rangle$")
    ax.set_title(r"(a) Resonance, $\hbar = 4\pi$")
    ax.legend(frameon=False)

    ax = axes[0, 1]
    rel = np.abs(p2 / exact - 1)
    positive = rel > 0
    ax.loglog(n[positive], rel[positive], ".", ms=2, color="C0", label="relative error")
    ax.loglog(n, n * np.log2(N) * EPS, "--", color="0.3", lw=1, label=r"allowance $n\,\log_2\! N\,\varepsilon$")
    ax.set_xlabel("kicks $n$")
    ax.set_ylabel(r"$|\langle p^2\rangle / (K^2n^2/2) - 1|$")
    ax.set_title("(b) Error stays at roundoff level")
    ax.legend(frameon=False, loc="upper left")

    ax = axes[1, 0]
    m, prob, exact_prob = dist_data
    window = np.abs(m) <= 70
    ax.semilogy(m[window], exact_prob[window], color="0.3", lw=1, label=r"$J_m^2(nK/\hbar)$")
    ax.semilogy(m[window], prob[window], "o", ms=2.5, color="C0", label="simulation")
    ax.set_ylim(1e-30, 1)
    ax.set_xlabel("$m$")
    ax.set_ylabel(r"$|\psi_m|^2$")
    ax.set_title(rf"(c) Distribution after $n = {N_DISTRIBUTION}$ kicks")
    ax.legend(frameon=False, loc="lower center")

    ax = axes[1, 1]
    n_anti = np.arange(0, p2_anti.size + 1)
    values = np.concatenate([[0.0], p2_anti])
    exact_anti = np.where(n_anti % 2 == 1, K**2 / 2, 0.0)
    ax.plot(n_anti, exact_anti, color="0.3", lw=1, drawstyle="steps-mid", label="exact")
    ax.plot(n_anti, values, "o", ms=3.5, color="C1", label="simulation")
    ax.set_xlabel("kicks $n$")
    ax.set_ylabel(r"$\langle p^2 \rangle$")
    ax.set_title(r"(d) Antiresonance, $\hbar = 2\pi$")
    ax.set_xticks(range(0, p2_anti.size + 1, 4))
    ax.set_ylim(-1, 17.5)
    ax.legend(frameon=False, loc="upper center", ncols=2)
    return fig


def print_summary(p2, dist_data, p2_anti, dist_anti, paths, elapsed):
    n = np.arange(1, p2.size + 1)
    rel = np.abs(p2 / (K**2 * n**2 / 2) - 1)
    m, prob, exact_prob = dist_data
    odd = np.arange(1, p2_anti.size + 1) % 2 == 1
    print("Experiment 05: quantum resonance and antiresonance (validation check 2)")
    print(f"K = {K:g}, N = {N} momentum states, start in m = 0")
    print()
    print(f"Resonance, hbar = 4 pi, {p2.size} kicks:")
    print(f"  <p^2> after {p2.size} kicks = {p2[-1]:.6f}; exact K^2 n^2 / 2 = {K**2 * p2.size**2 / 2:.6f}")
    print(f"  largest relative error over all kicks: {rel.max():.1e}"
          f" (roundoff allowance n log2(N) eps = {p2.size * np.log2(N) * EPS:.1e})")
    print(f"  distribution after {N_DISTRIBUTION} kicks vs J_m(nK/hbar)^2: largest difference {np.abs(prob - exact_prob).max():.1e}")
    print(f"Antiresonance, hbar = 2 pi, {p2_anti.size} kicks:")
    print(f"  odd kicks: <p^2> = K^2/2 = {K**2 / 2:g} to within {np.abs(p2_anti[odd] / (K**2 / 2) - 1).max():.1e} (relative)")
    print(f"  even kicks: back to m = 0, ||psi_n - psi_0|| <= {dist_anti[~odd].max():.1e}; largest <p^2> = {p2_anti[~odd].max():.1e}")
    print()
    print("Saved " + " and ".join(path.relative_to(REPO_ROOT).as_posix() for path in paths))
    print(f"Run time {elapsed:.1f} s")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Experiment 05: quantum resonance and antiresonance.")
    parser.add_argument("--show", action="store_true", help="open the figure in a window after saving it")
    args = parser.parse_args(argv)

    plotting.use_backend(args.show)
    plotting.apply_style()
    start = time.perf_counter()

    p2 = resonance_run(K, N, N_KICKS)
    dist_data = resonance_distribution(K, N, N_DISTRIBUTION)
    p2_anti, dist_anti = antiresonance_run(K, N, N_ANTI)
    fig = make_figure(p2, dist_data, p2_anti)
    paths = plotting.save_figure(fig, REPO_ROOT / "figures", FIGURE_STEM)
    print_summary(p2, dist_data, p2_anti, dist_anti, paths, time.perf_counter() - start)

    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
