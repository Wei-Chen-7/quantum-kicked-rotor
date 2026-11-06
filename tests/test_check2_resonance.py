"""Check 2: exact results at the quantum resonance and antiresonance.

At hbar = 4 pi the free step exp(-i hbar m^2 / 2) = exp(-2 pi i m^2) is the
identity, so n kicks act like one kick of strength n K / hbar. Starting from
m = 0 this gives c_m = (-i)^m J_m(n K / hbar), hence <p^2>_n = K^2 n^2 / 2
(Casati et al. 1979, p. 349).

At hbar = 2 pi the free step is (-1)^m, a shift of theta by pi, so the second
kick undoes the first: the state returns exactly after every two kicks, for
any K and any initial state (the "anti-resonant" case of Casati et al. 1979).
"""

import numpy as np
from common import REFERENCE_TOL, K, fft_roundoff
from scipy.special import jv

from kicked_rotor import observables, quantum

N = quantum.DEFAULT_N


def test_resonance_energy_grows_as_n_squared():
    hbar, n_kicks = 4 * np.pi, 1000
    psi0 = quantum.momentum_eigenstate(N)
    p2 = np.array([observables.mean_p2(s, hbar) for s in quantum.iterate(psi0, K, hbar, n_kicks)])
    n = np.arange(1, n_kicks + 1)
    relative_error = np.abs(p2 / (K**2 * n**2 / 2) - 1)
    # Machine precision, growing with the number of kicks.
    assert np.all(relative_error <= fft_roundoff(n, N)), f"worst relative error {relative_error.max():.2e}"


def test_resonance_distribution_is_a_bessel_function():
    hbar, n_kicks = 4 * np.pi, 100
    state = quantum.evolve(quantum.momentum_eigenstate(N), K, hbar, n_kicks)
    m, prob = observables.momentum_probabilities(state)
    exact = jv(m, n_kicks * K / hbar) ** 2
    # An error delta in the state changes |c_m|^2 by at most about 2|delta|
    # (machine precision), plus the accuracy of scipy's Bessel functions.
    assert np.max(np.abs(prob - exact)) <= 2 * fft_roundoff(n_kicks, N) + REFERENCE_TOL


def test_antiresonance_alternates():
    hbar, n_kicks = 2 * np.pi, 20
    psi0 = quantum.momentum_eigenstate(N)
    for n, state in enumerate(quantum.iterate(psi0, K, hbar, n_kicks), start=1):
        if n % 2 == 1:
            # Odd n: one effective kick, so <p^2> = K^2/2 to machine precision.
            assert abs(observables.mean_p2(state, hbar) / (K**2 / 2) - 1) <= fft_roundoff(n, N)
        else:
            # Even n: back to the initial state, to machine precision.
            assert np.linalg.norm(state - psi0) <= fft_roundoff(n, N)
