"""Check 5: time reversal.

The quantum evolution is linear and unitary, so running it backward returns
the initial state up to roundoff that grows only linearly with time. A
classical trajectory, run backward, amplifies roundoff by its Lyapunov
exponent lambda, so the return fails once eps * exp(lambda n) reaches order
one: after about ln(1/eps) / lambda = 36.0 / 0.98 = 37 kicks at K = 5
(PLAN.md note 6; lambda = 0.98 was measured in Step 1).
"""

import numpy as np
from common import EPS, K, fft_roundoff

from kicked_rotor import classical, quantum


def test_quantum_returns_after_1000_kicks():
    N, hbar, n_kicks = quantum.DEFAULT_N, 0.25, 1000
    psi0 = quantum.momentum_eigenstate(N)
    back = quantum.evolve(quantum.evolve(psi0, K, hbar, n_kicks), K, hbar, n_kicks, backward=True)
    # Machine precision: 2 * n_kicks periods in total.
    assert np.linalg.norm(back - psi0) <= fft_roundoff(2 * n_kicks, N)


def test_classical_inverse_undoes_one_kick():
    rng = np.random.default_rng(5)
    theta, p = rng.uniform(0, 2 * np.pi, 10**5), rng.uniform(-200, 200, 10**5)
    theta1, p1 = classical.standard_map(theta, p, K)
    theta_back, p_back = classical.inverse_standard_map(theta1, p1, K)
    # Machine precision: each operation rounds by at most half an ulp of its
    # result, so the angle comes back within a few eps * (|p1| + 2 pi), and
    # that angle error enters p through K sin(theta).
    tol_theta = 2 * EPS * (np.abs(p1) + 2 * np.pi)
    tol_p = K * tol_theta + 2 * EPS * (np.abs(p1) + K)
    assert np.all(np.abs(classical.wrap_centered(theta_back - theta)) <= tol_theta)
    assert np.all(np.abs(p_back - p) <= tol_p)


def median_return_distance(n_kicks, theta0, p0):
    theta, p = classical.evolve(theta0, p0, K, n_kicks)
    theta, p = classical.evolve(theta, p, K, n_kicks, backward=True)
    return np.median(np.hypot(classical.wrap_centered(theta - theta0), p - p0))


def test_classical_return_fails_after_the_lyapunov_time():
    theta0, p0 = classical.initial_ensemble(4000, sampling="grid")

    # Short runs come back. Even the most unstable orbit stretches errors by
    # at most the largest singular value of the one-kick Jacobian, about 8 at
    # K = 5, so after 10 kicks roundoff of 1e-14 grows to at most about 1e-5.
    assert median_return_distance(10, theta0, p0) < 1e-5

    # Long runs fail, and the first n where the median distance passes 1 must
    # match ln(1/eps_eff) / lambda: with eps_eff between eps/2 and 10 eps and
    # lambda between ln(K/2) = 0.92 and 1.05 this lies between 32 and 40.
    first_failure = next(n for n in range(1, 61) if median_return_distance(n, theta0, p0) > 1.0)
    assert 30 <= first_failure <= 45, f"return first failed at n = {first_failure}"
