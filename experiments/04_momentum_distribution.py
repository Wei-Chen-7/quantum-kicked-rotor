"""Experiment 04: the localized momentum distribution and its localization length.

After the break time the quantum momentum distribution stops spreading. Its
time average decays exponentially in |m|; in the project's convention
|psi_m|^2 ~ exp(-2|m|/xi). Izrailev (1990, Eqs. 2.3.13 and 2.3.14) gives the
full shape of the time-averaged distribution of a state started at m = 0,

    w_m = (1 + 2|m|/xi) exp(-2|m|/xi) / (2 xi),   with xi ~ D / hbar^2,

and Shepelyansky (1986) reports xi / (D/hbar^2) = 1.04 on average, with an
rms spread of 0.20 (PLAN.md note 5). The classical distribution at the same
time is a Gaussian that keeps widening.

The far tails fill in slowly, so the fit uses a late averaging window; an
early window is also reported to show the drift. (A one-off 40000-kick run at
hbar = 0.25 shows the fit leveling off at about 0.91 D/hbar^2 for windows from
10000 kicks on; see PLAN.md note 5.)

Also saves an animation (GIF) of the distribution evolving, quantum against
classical, at hbar = 0.25. Its frames are taken at kicks n = 1000 (k/80)^2,
so the classical width, which grows like sqrt(n), widens at a steady rate.

Run from the repo root:
    python experiments/04_momentum_distribution.py           # save figures and the GIF
    python experiments/04_momentum_distribution.py --show    # also open a window
"""

import argparse
import pathlib
import sys
import time

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation, PillowWriter

# "python experiments/<script>.py" puts experiments/ on the import path, not
# the repo root, so add the root to make the kicked_rotor package importable.
REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from kicked_rotor import classical, observables, plotting, quantum  # noqa: E402

FIGURE_STEM = "04_momentum_distribution"
K = 5.0
HBARS = (0.25, 0.5, 1.0)  # hbar / 4 pi irrational for all three: no quantum resonance
N = quantum.DEFAULT_N
N_KICKS = 5000
EARLY_WINDOW = (1000, 2000)  # averaging windows (first, last], in kicks
LATE_WINDOW = (2500, 5000)
N_PARTICLES = 40_000  # classical ensemble, the same as experiment 03
D_WINDOW = (200, 2000)  # classical D from the slope of <p^2> over these kicks, as in experiment 03
TAIL_LEVEL = 1e-14  # fit out to where the averaged distribution falls to this fraction of its peak
M_MIN = 5  # leave out the narrow peak at m = 0
MIN_COUNT = 10  # draw a classical histogram bin only if it holds at least this many particles
COLORS = {0.25: "C0", 0.5: "C1", 1.0: "C2"}

# Animation: hbar = 0.25, frames at kicks ANIM_KICKS * (k / ANIM_FRAMES)^2
ANIM_HBAR = 0.25
ANIM_KICKS = 1000
ANIM_FRAMES = 80
ANIM_PARTICLES = 200_000
ANIM_HALF_WIDTH = 1500  # plotted range of m
ANIM_BIN = 10  # classical histogram bin width, in m


def averaged_distribution(K, hbar, N, n_kicks, windows):
    """Time averages of |psi_m|^2 over the given (first, last] kick windows.

    Returns ``(m, averages, last, edge)``: the momenta sorted from -N/2 up,
    a dict of averages keyed by window, the distribution after the final kick,
    and the largest edge probability (|m| > 3N/8) over the run.
    """
    m = quantum.integer_momenta(N)
    outer = np.abs(m) > 3 * N / 8
    sums = {window: np.zeros(N) for window in windows}
    edge = 0.0
    for n, state in enumerate(quantum.iterate(quantum.momentum_eigenstate(N), K, hbar, n_kicks), start=1):
        prob = np.abs(quantum.to_momentum(state)) ** 2
        edge = max(edge, prob[outer].sum())
        for (first, last), total in sums.items():
            if first < n <= last:
                total += prob
    averages = {(first, last): np.fft.fftshift(total) / (last - first) for (first, last), total in sums.items()}
    return np.fft.fftshift(m), averages, np.fft.fftshift(prob), edge


def fit_xi(m, avg):
    """Izrailev-shape fits over [M_MIN, m_tail] and [M_MIN, m_tail/2], and sqrt<m^2>."""
    m_tail = np.abs(m[avg > TAIL_LEVEL * avg.max()]).max()
    xi_full, amp = observables.fit_averaged_localization_length(m, avg, M_MIN, m_tail)
    xi_half, _ = observables.fit_averaged_localization_length(m, avg, M_MIN, m_tail / 2)
    return {"xi": xi_full, "amp": amp, "xi_half": xi_half, "rms": np.sqrt(np.sum(m**2 * avg)), "m_tail": m_tail}


def classical_run(K, n_particles, n_kicks, d_window):
    """Momenta after n_kicks for p0 = 0, theta0 on a uniform grid, and D from the slope of <p^2> over d_window."""
    theta, p = classical.initial_ensemble(n_particles, sampling="grid")
    first, last = d_window
    for n in range(1, n_kicks + 1):
        theta, p = classical.standard_map(theta, p, K)
        if n == first:
            p2_first = p**2
        if n == last:
            growth = (p**2 - p2_first) / (last - first)
    return p, growth.mean(), growth.std() / np.sqrt(n_particles)


def classical_histogram(p, hbar, half_width, bin_width):
    """Classical probability per momentum state m = p/hbar, from a histogram of the ensemble.

    The bins are centered on m = 0, +-bin_width, ... Bins holding fewer than
    MIN_COUNT particles are returned as NaN so plots leave them out: their
    relative counting error, 1/sqrt(count), would be above about 30%, and a
    single stray particle would otherwise draw a false plateau.
    """
    n_side = int(np.ceil(half_width / bin_width))
    edges = (np.arange(-n_side, n_side + 2) - 0.5) * bin_width
    counts, _ = np.histogram(p / hbar, bins=edges)
    density = np.where(counts >= MIN_COUNT, counts / (p.size * bin_width), np.nan)
    return 0.5 * (edges[1:] + edges[:-1]), density


def izrailev_shape(m, xi, amp):
    x = np.abs(m)
    return amp * (1 + 2 * x / xi) * np.exp(-2 * x / xi)


def make_figure(results, p_classical, D):
    fig, axes = plt.subplots(1, 3, figsize=(6.5, 3.0), layout="constrained", sharey=True)
    for ax, (hbar, res) in zip(axes, results.items()):
        xi_pred = D / hbar**2
        m, avg, fit = res["m"], res["late"], res["fit_late"]
        half_width = 12 * xi_pred
        show = np.abs(m) <= half_width
        ax.semilogy(m[show], res["last"][show], color="0.8", lw=0.5, label=f"quantum, kick {N_KICKS}")
        ax.semilogy(m[show], avg[show], color=COLORS[hbar], lw=0.7, label="quantum, time average")
        ax.semilogy(m[show], izrailev_shape(m[show], fit["xi"], fit["amp"]), "k-", lw=0.8,
                    label=rf"fit, $\xi = {fit['xi']:.3g}$")
        ax.semilogy(m[show], izrailev_shape(m[show], xi_pred, 1 / (2 * xi_pred)), "k--", lw=0.8,
                    label=rf"$\xi = D/\hbar^2 = {xi_pred:.3g}$")
        centers, density = classical_histogram(p_classical, hbar, half_width, max(1, round(half_width / 60)))
        ax.semilogy(centers, density, color="C3", lw=0.8, drawstyle="steps-mid", label=f"classical, kick {N_KICKS}")
        ax.set_xlim(-half_width, half_width)
        ax.set_ylim(1e-20, 1)
        ax.set_xlabel("$m$")
        ax.set_title(rf"$\hbar = {hbar:g}$")
        ax.legend(frameon=False, fontsize=6.5, loc="lower center", handlelength=1.5)
    axes[0].set_ylabel(r"$|\psi_m|^2$")
    return fig


def make_animation(K, hbar, N, path):
    """GIF of the momentum distribution, quantum (m = 0 start) against classical (p = 0, theta uniform)."""
    frame_kicks = np.unique(np.rint(ANIM_KICKS * (np.arange(ANIM_FRAMES + 1) / ANIM_FRAMES) ** 2).astype(int))
    psi = quantum.momentum_eigenstate(N)
    kick, free = quantum.floquet_factors(N, K, hbar)
    theta, p = classical.initial_ensemble(ANIM_PARTICLES, sampling="grid")
    m = np.fft.fftshift(quantum.integer_momenta(N))
    show = np.abs(m) <= ANIM_HALF_WIDTH
    frames, n = [], 0
    for target in frame_kicks:
        while n < target:
            psi = quantum.floquet_step(psi, kick, free)
            theta, p = classical.standard_map(theta, p, K)
            n += 1
        prob = np.fft.fftshift(np.abs(quantum.to_momentum(psi)) ** 2)[show]
        frames.append((n, prob, classical_histogram(p, hbar, ANIM_HALF_WIDTH, ANIM_BIN)))

    fig, ax = plt.subplots(figsize=(5.2, 3.0), layout="constrained")
    (q_line,) = ax.semilogy(m[show], frames[0][1], color="C0", lw=0.7, label=rf"quantum, $\hbar = {hbar:g}$")
    (c_line,) = ax.semilogy([], [], color="C3", lw=1.0, drawstyle="steps-mid", label="classical")
    label = ax.text(0.02, 0.95, "", transform=ax.transAxes, va="top")
    ax.set_xlim(-ANIM_HALF_WIDTH, ANIM_HALF_WIDTH)
    ax.set_ylim(1e-16, 1.5)
    ax.set_xlabel("$m$")
    ax.set_ylabel(r"probability $|\psi_m|^2$")
    ax.set_title(rf"Momentum distribution, $K = {K:g}$")
    ax.legend(frameon=False, loc="upper right")

    def update(frame):
        n, prob, (centers, density) = frame
        q_line.set_ydata(np.maximum(prob, 1e-30))  # exact zeros (at n = 0) cannot go on a log scale
        c_line.set_data(centers, density)
        label.set_text(f"kick {n}")
        return q_line, c_line, label

    animation = FuncAnimation(fig, update, frames=frames, blit=False)
    animation.save(path, writer=PillowWriter(fps=12), dpi=100)
    plt.close(fig)
    return len(frames)


def print_summary(results, D, D_err, n_frames, gif_path, paths, elapsed):
    print("Experiment 04: the localized momentum distribution")
    print(f"K = {K:g}, N = {N}, start in m = 0, {N_KICKS} kicks; classical D = {D:.2f} +- {D_err:.2f}"
          f" (as in experiment 03)")
    print(f"Fit of (1 + 2|m|/xi) exp(-2|m|/xi) over {M_MIN} <= |m| <= m_tail, where the average falls to"
          f" {TAIL_LEVEL:g} of its peak")
    print()
    print("  hbar   D/hbar^2   window      m_tail   xi fit   xi fit (half range)   sqrt<m^2>   xi fit / (D/hbar^2)")
    for hbar, res in results.items():
        xi_pred = D / hbar**2
        for name, window in (("early", EARLY_WINDOW), ("late", LATE_WINDOW)):
            fit = res[f"fit_{name}"]
            print(f"  {hbar:4.2f}   {xi_pred:8.1f}   {window[0]:4d}-{window[1]:<5d}  {fit['m_tail']:6d}   {fit['xi']:6.1f}"
                  f"   {fit['xi_half']:19.1f}   {fit['rms']:9.1f}   {fit['xi'] / xi_pred:19.2f}")
    print(f"  largest edge probability (|m| > 3N/8) in any run: {max(r['edge'] for r in results.values()):.1e}")
    print()
    print("Shepelyansky (1986): xi / (D/hbar^2) = 1.04 on average, rms spread 0.20.")
    print("The far tail fills in slowly, so the early window gives a smaller xi. A one-off check at hbar = 0.25")
    print("out to 40000 kicks levels off at 0.91 (PLAN.md note 5), close to the late window here.")
    print("At hbar = 1, xi is only about 12 states and the averaged distribution is bumpy, so that fit is rough.")
    print()
    print("Saved " + " and ".join(path.relative_to(REPO_ROOT).as_posix() for path in paths))
    print(f"Saved {gif_path.relative_to(REPO_ROOT).as_posix()} ({n_frames} frames up to kick {ANIM_KICKS},"
          f" {gif_path.stat().st_size / 1e6:.1f} MB)")
    print(f"Run time {elapsed:.1f} s")


def main(argv=None):
    parser = argparse.ArgumentParser(description="Experiment 04: localized momentum distribution and xi.")
    parser.add_argument("--show", action="store_true", help="open the figure in a window after saving it")
    args = parser.parse_args(argv)

    plotting.use_backend(args.show)
    plotting.apply_style()
    start = time.perf_counter()

    p_classical, D, D_err = classical_run(K, N_PARTICLES, N_KICKS, D_WINDOW)
    results = {}
    for hbar in HBARS:
        m, averages, last, edge = averaged_distribution(K, hbar, N, N_KICKS, (EARLY_WINDOW, LATE_WINDOW))
        results[hbar] = {
            "m": m,
            "late": averages[LATE_WINDOW],
            "last": last,
            "edge": edge,
            "fit_early": fit_xi(m, averages[EARLY_WINDOW]),
            "fit_late": fit_xi(m, averages[LATE_WINDOW]),
        }

    fig = make_figure(results, p_classical, D)
    paths = plotting.save_figure(fig, REPO_ROOT / "figures", FIGURE_STEM)
    gif_path = REPO_ROOT / "figures" / f"{FIGURE_STEM}.gif"
    n_frames = make_animation(K, ANIM_HBAR, N, gif_path)
    print_summary(results, D, D_err, n_frames, gif_path, paths, time.perf_counter() - start)

    if args.show:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    main()
