"""Check 4: the momentum grid is wide enough.

The FFT treats m as periodic, so probability that reaches |m| = N/2 would
wrap around. The test runs the hardest default case, K = 5 and hbar = 0.25
(the localization length grows like 1/hbar^2, so smaller hbar needs more
states), for 1000 kicks.

EDGE_LIMIT comes from PLAN.md. At hbar = 0.25 and N = 2^14, probability 1e-10
beyond |m| = 3N/8 could change <p^2> by at most hbar^2 (N/2)^2 * 1e-10 = 4e-4,
under 1e-6 of <p^2> (about 1150 after 1000 kicks).
"""

import numpy as np
from common import K, fft_roundoff

from kicked_rotor import observables, quantum

HBAR = 0.25
N_KICKS = 1000
EDGE_LIMIT = 1e-10


def run(N):
    """<p^2> after each kick, the largest edge probability, and the final distribution."""
    p2, edge = [], 0.0
    for state in quantum.iterate(quantum.momentum_eigenstate(N), K, HBAR, N_KICKS):
        p2.append(observables.mean_p2(state, HBAR))
        edge = max(edge, observables.edge_probability(state))
    m, prob = observables.momentum_distribution(state)
    return np.array(p2), edge, m, prob


def test_default_grid_is_wide_enough():
    N = quantum.DEFAULT_N
    p2, edge, m, prob = run(N)
    p2_big, edge_big, m_big, prob_big = run(2 * N)
    assert edge < EDGE_LIMIT, f"edge probability {edge:.1e}"

    # Doubling N changes nothing beyond machine precision.
    n = np.arange(1, N_KICKS + 1)
    assert np.all(np.abs(p2 / p2_big - 1) <= fft_roundoff(n, 2 * N))
    on_small_grid = np.isin(m_big, m)
    assert np.array_equal(m_big[on_small_grid], m)
    # An error delta in the state changes |c_m|^2 by at most about 2|delta|.
    assert np.max(np.abs(prob_big[on_small_grid] - prob)) <= 2 * fft_roundoff(N_KICKS, 2 * N)


def test_too_small_grid_is_caught():
    # Negative control: PLAN.md note 4 predicted that 2^12 states are too few
    # at hbar = 0.25. The edge criterion must flag it.
    _, edge, _, _ = run(2**12)
    assert edge > EDGE_LIMIT, f"edge probability {edge:.1e}"
