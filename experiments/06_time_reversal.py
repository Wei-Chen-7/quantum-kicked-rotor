"""Experiment 06: time reversal (validation check 5).

Run n kicks forward, then n kicks backward with the inverse map, and measure
how far the system ends up from where it started.

Quantum: the evolution is linear and unitary, so the roundoff made in each
step is never amplified, and the return error grows only linearly with n.

Classical: each trajectory is chaotic, so the roundoff made in every step is
stretched by about e^lambda per kick. The return fails once eps * e^(lambda n)
reaches order one, after about ln(1/eps) / lambda = 37 kicks at K = 5.
(The classical phase-space density obeys the Liouville equation, which is also
linear and reversible; the irreversibility comes from following individual
chaotic trajectories.)

Run from the repo root:
    python experiments/06_time_reversal.py           # save figures
    python experiments/06_time_reversal.py --show    # also open a window
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

from kicked_rotor import classical, plotting, quantum  # noqa: E402

FIGURE_STEM = "06_time_reversal"
EPS = np.finfo(float).eps
K = 5.0
HBAR = 0.25
N = quantum.DEFAULT_N
N_ORBITS = 4000  # classical ensemble: p = 0, theta on a uniform grid
CLASSICAL_KICKS = np.arange(1, 71)
QUANTUM_KICKS = np.array([1, 2, 3, 5, 7, 10, 15, 20, 30, 40, 50, 60, 70, 100, 200, 500, 1000])


def classical_return_distances(K, theta0, p0, n_values):
    """Phase-space distance of each orbit from its start after n kicks forward and n back."""
    rows = []
    for n in n_values:
        theta, p = classical.evolve(theta0, p0, K, n)
        theta, p = classical.evolve(theta, p, K, n, backward=True)
        rows.append(np.hypot(classical.wrap_centered(theta - theta0), p - p0))
    return np.array(rows)  # shape (len(n_values), number of orbits)


def quantum_return_errors(K, hbar, N, n_values):
    """||psi_back - psi_0|| after n kicks forward and n back."""
    psi0 = quantum.momentum_eigenstate(N)
    errors = []
    for n in n_values:
        back = quantum.evolve(quantum.evolve(psi0, K, hbar, n), K, hbar, n, backward=True)
        errors.append(np.linalg.norm(back - psi0))
    return np.array(errors)


def growth_fit(n_values, median):
    """Fit median = A e^(rate n) where the median is between roundoff and saturation."""
    window = (median > 1e-13) & (median < 1e-3)
    rate, log_a = np.polyfit(n_values[window], np.log(median[window]), 1)
    return rate, np.exp(log_a), n_values[window]


def make_figure(distances, median, lam, q_errors):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(6.5, 3.1), layout="constrained", width_ratios=(1.35, 1))
    n = CLASSICAL_KICKS
    q10, q90 = np.percentile(distances, [10, 90], axis=1)
    floor = 1e-18  # zeros cannot be drawn on a log axis; this is below the plotted range
    ax.fill_between(n, np.maximum(q10, floor), np.maximum(q90, floor), color="C3", alpha=0.2, lw=0,
                    label="classical, 10% to 90% of orbits")
    show = median > 0
    ax.semilogy(n[show], median[show], "o-", ms=2.5, lw=0.8, color="C3", label="classical, median")
    ax.semilogy(n, EPS * np.exp(lam * n), "--", color="0.3", lw=1,
                label=rf"$\varepsilon\, e^{{\lambda n}}$, $\lambda = {lam:.2f}$")
    in_range = QUANTUM_KICKS <= n[-1]
    ax.semilogy(QUANTUM_KICKS[in_range], q_errors[in_range], "s", ms=3, color="C0", label="quantum")
    ax.axhline(1.0, color="0.5", lw=0.8, ls=":")
    ax.set_ylim(1e-17, 1e3)
    ax.set_xlim(0, n[-1])
    ax.set_xlabel("kicks forward, then the same number back")
    ax.set_ylabel("distance from the start")
    ax.set_title("(a) Classical return fails, quantum does not")
    ax.legend(frameon=False, loc="upper left", fontsize=8)

    ax2.loglog(QUANTUM_KICKS, q_errors, "s", ms=3, color="C0", label=r"$\|\psi_{\rm back} - \psi_0\|$")
    ax2.loglog(QUANTUM_KICKS, 2 * QUANTUM_KICKS * np.log2(N) * EPS, "--", color="0.3", lw=1,
               label=r"allowance $2n\,\log_2\! N\,\varepsilon$")
    ax2.set_ylim(1e-17, 1e-10)
    ax2.set_xlabel("kicks $n$")
    ax2.set_title("(b) Quantum error grows only linearly")
    ax2.legend(frameon=False, loc="upper left", fontsize=8)
    return fig


def print_summary(median, lam, rate, prefactor, window, q_errors, paths, elapsed):
    n = CLASSICAL_KICKS
    first_failure = int(n[np.argmax(median > 1.0)]) if np.any(median > 1.0) else None
    pick = lambda k: median[k - 1]
    print("Experiment 06: time reversal (validation check 5)")
    print("n kicks forward, then n kicks backward; distance from the starting point.")
    print()
    print(f"Classical: {N_ORBITS} orbits, K = {K:g}, starting at p = 0 with theta uniform")
    print(f"  Lyapunov exponent of these orbits: lambda = {lam:.3f} (median finite-time value, 2000 kicks)")
    print("  median return distance: " + ", ".join(f"n = {k}: {pick(k):.1e}" for k in (10, 20, 30, 40, 60)))
    print(f"  fit over n = {window[0]} to {window[-1]}: median = {prefactor:.1e} * exp({rate:.3f} n)"
          f"  (eps = {EPS:.1e}, lambda = {lam:.3f})")
    print(f"  the median first exceeds 1 at n = {first_failure}; ln(1/eps)/lambda = {np.log(1 / EPS) / lam:.1f}")
    print(f"Quantum: K = {K:g}, hbar = {HBAR:g}, N = {N}, start in m = 0")
    print("  ||psi_back - psi_0||: " + ", ".join(
        f"n = {k}: {q_errors[list(QUANTUM_KICKS).index(k)]:.1e}" for k in (10, 100, 1000)))
    print(f"  roundoff allowance at n = 1000: 2 n log2(N) eps = {2 * 1000 * np.log2(N) * EPS:.1e}")
    print()
    print("Why: the quantum evolution is unitary and linear, so roundoff is never amplified.")
    print("A chaotic classical orbit stretches the roundoff of every step by about e^lambda")
    print("per kick, so it is lost once eps * e^(lambda n) ~ 1.")
    print()
    print("Saved " + " and ".join(path.relative_to(REPO_ROOT).as_posix() for path in paths))
    print(f"Run time {elapsed:.1f} s")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Experiment 06: time reversal, quantum vs classical.")
    parser.add_argument("--show", action="store_true", help="open the figure in a window after saving it")
    args = parser.parse_args(argv)

    plotting.use_backend(args.show)
    plotting.apply_style()
    start = time.perf_counter()

    theta0, p0 = classical.initial_ensemble(N_ORBITS, sampling="grid")
    lam = float(np.median(classical.lyapunov_exponent(theta0, p0, K, 2000, n_transient=200)))
    distances = classical_return_distances(K, theta0, p0, CLASSICAL_KICKS)
    median = np.median(distances, axis=1)
    rate, prefactor, window = growth_fit(CLASSICAL_KICKS, median)
    q_errors = quantum_return_errors(K, HBAR, N, QUANTUM_KICKS)

    fig = make_figure(distances, median, lam, q_errors)
    paths = plotting.save_figure(fig, REPO_ROOT / "figures", FIGURE_STEM)
    print_summary(median, lam, rate, prefactor, window, q_errors, paths, time.perf_counter() - start)

    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
