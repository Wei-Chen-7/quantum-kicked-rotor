"""Observables of the quantum kicked rotor.

Every function takes a state psi stored on the angle grid as described in
quantum.py: sum_j |psi_j|^2 = 1, and the momentum amplitudes c_m come from
the unitary FFT.
"""

import numpy as np

from kicked_rotor import quantum


def norm(psi):
    """<psi|psi> = sum_j |psi_j|^2. It equals sum_m |c_m|^2, since the FFT is unitary."""
    return float(np.sum(np.abs(psi) ** 2))


def fidelity(psi, phi):
    """|<psi|phi>|^2: 1 for the same state (up to a phase), 0 for orthogonal states."""
    return float(np.abs(np.vdot(psi, phi)) ** 2)


def momentum_probabilities(psi):
    """Integer momenta m and probabilities P_m = |c_m|^2, both in FFT order."""
    return quantum.integer_momenta(psi.size), np.abs(quantum.to_momentum(psi)) ** 2


def momentum_distribution(psi):
    """``(m, P_m)`` sorted by m, from -N/2 up to N/2 - 1 (handy for plotting)."""
    m, prob = momentum_probabilities(psi)
    return np.fft.fftshift(m), np.fft.fftshift(prob)


def mean_p2(psi, hbar):
    """<p^2> = hbar^2 sum_m m^2 |c_m|^2."""
    m, prob = momentum_probabilities(psi)
    return float(hbar**2 * np.sum(m.astype(float) ** 2 * prob))


def edge_probability(psi, edge_fraction=3 / 8):
    """Total probability in the outer momentum states, |m| > edge_fraction * N.

    The FFT treats m as periodic, so probability that reaches |m| = N/2 would
    reappear on the other side of the grid. With the default edge_fraction =
    3/8 this is the probability in the outer quarter of the grid. It changes
    from kick to kick, so track its largest value over a run, not only the
    final one.
    """
    m, prob = momentum_probabilities(psi)
    return float(np.sum(prob[np.abs(m) > edge_fraction * psi.size]))


def fit_localization_length(m, prob, m_min, m_max):
    """Fit prob_m = A exp(-2|m| / xi) for m_min <= |m| <= m_max; return ``(xi, A)``.

    The fit is a straight line through ln(prob) against |m|, using both signs
    of m. This fixes the project's convention for the localization length:
    |psi_m|^2 ~ exp(-2|m|/xi), so |psi_m| ~ exp(-|m|/xi). Pass a time-averaged
    distribution, with m_min past the central peak and m_max before the tail
    reaches roundoff level or the grid edge.

    Caution: the time-averaged distribution from m = 0 is expected to be
    (1 + 2|m|/xi) exp(-2|m|/xi), not a pure exponential (PLAN.md note 5).
    A pure-exponential fit to that shape overestimates xi by about 30% over
    xi/2 <= |m| <= 3 xi and by 7% even over 4 xi <= |m| <= 10 xi.
    """
    m = np.asarray(m)
    prob = np.asarray(prob, dtype=float)
    size = np.abs(m)
    window = (size >= m_min) & (size <= m_max) & (prob > 0)
    if np.count_nonzero(window) < 2:
        raise ValueError("need at least two points with prob > 0 in the fit window")
    slope, intercept = np.polyfit(size[window].astype(float), np.log(prob[window]), 1)
    if slope >= 0:
        raise ValueError("the distribution does not decay in the fit window")
    return -2.0 / slope, float(np.exp(intercept))


def fit_averaged_localization_length(m, prob, m_min, m_max):
    """Fit prob_m = A (1 + 2|m|/xi) exp(-2|m|/xi) for m_min <= |m| <= m_max; return ``(xi, A)``.

    This is the shape Izrailev (1990, Eq. 2.3.13) derives for the time-averaged
    distribution of a state started at m = 0, built from exponentially
    localized eigenstates with random fluctuations. For a normalized
    distribution of this shape A = 1/(2 xi) and <m^2> = xi^2. The fit is a
    least-squares fit of ln(prob), using both signs of m.
    """
    from scipy.optimize import curve_fit

    m = np.asarray(m)
    prob = np.asarray(prob, dtype=float)
    size = np.abs(m).astype(float)
    window = (size >= m_min) & (size <= m_max) & (prob > 0)

    def log_model(x, xi, log_a):
        return log_a - 2 * x / xi + np.log1p(2 * x / xi)

    xi_start, a_start = fit_localization_length(m, prob, m_min, m_max)
    (xi, log_a), _ = curve_fit(log_model, size[window], np.log(prob[window]), p0=(xi_start, np.log(a_start)))
    return float(xi), float(np.exp(log_a))
