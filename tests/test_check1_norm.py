"""Check 1: the norm stays at 1 to near machine precision."""

import numpy as np
import pytest
from common import K, fft_roundoff

from kicked_rotor import observables, quantum


@pytest.mark.parametrize("hbar", [0.25, 1.0])
def test_norm_stays_at_one(hbar):
    N, n_kicks = quantum.DEFAULT_N, 1000
    psi0 = quantum.momentum_eigenstate(N)
    errors = np.array([abs(observables.norm(s) - 1.0) for s in quantum.iterate(psi0, K, hbar, n_kicks)])
    # Machine precision: after n kicks the allowance is fft_roundoff(n, N).
    allowance = fft_roundoff(np.arange(1, n_kicks + 1), N)
    assert np.all(errors <= allowance), f"worst |norm - 1| = {errors.max():.2e}"


def test_inverse_step_keeps_the_norm():
    N, n_kicks = quantum.DEFAULT_N, 200
    psi = quantum.evolve(quantum.momentum_eigenstate(N), K, 0.5, n_kicks, backward=True)
    assert abs(observables.norm(psi) - 1.0) <= fft_roundoff(n_kicks, N)
