"""Check 3: early times, against exact results (PLAN.md note 1).

Quantum (state m = 0) and classical (p0 = 0, theta0 uniform) agree exactly
after one kick. After two kicks each code is compared with its own exact
formula. After that, quantum and classical differ by an amount that grows
very fast with n (PLAN.md note 1), so at n = 3 the test checks the
correspondence limit instead: the gap shrinks like hbar^2 as hbar -> 0, and a
Richardson extrapolation recovers the classical value.
"""

import numpy as np
import pytest
from common import REFERENCE_TOL, K, fft_roundoff, summation_roundoff

from kicked_rotor import classical, observables, quantum

N = quantum.DEFAULT_N
K_VALUES = (2.0, 5.0, 9.5)
HBAR_VALUES = (0.25, 0.5, 1.0, 2.0, 3.0)


def quantum_p2(K, hbar, n_kicks):
    state = quantum.evolve(quantum.momentum_eigenstate(N), K, hbar, n_kicks)
    assert observables.edge_probability(state) < 1e-10  # the grid is wide enough
    return observables.mean_p2(state, hbar)


def classical_p2(K, n_kicks, M):
    theta, p = classical.initial_ensemble(M, sampling="grid")
    theta, p = classical.evolve(theta, p, K, n_kicks)
    return np.mean(p**2)


@pytest.mark.parametrize("hbar", HBAR_VALUES)
def test_one_kick_quantum(hbar):
    assert abs(quantum_p2(K, hbar, 1) / (K**2 / 2) - 1) <= fft_roundoff(1, N)


def test_one_kick_classical():
    # The midpoint rule averages sin^2 exactly for M >= 3, so only rounding is left.
    M = 4096
    assert abs(classical_p2(K, 1, M) / (K**2 / 2) - 1) <= summation_roundoff(M)


@pytest.mark.parametrize("K_value", K_VALUES)
@pytest.mark.parametrize("hbar", HBAR_VALUES)
def test_two_kicks_quantum_formula(K_value, hbar):
    relative_error = abs(quantum_p2(K_value, hbar, 2) / quantum.exact_p2_two_kicks(K_value, hbar) - 1)
    assert relative_error <= fft_roundoff(2, N) + REFERENCE_TOL


@pytest.mark.parametrize("K_value", K_VALUES)
def test_two_kicks_classical_formula_grid(K_value):
    # The average is an integral over theta0 of a smooth periodic function
    # whose Fourier coefficients decay faster than exponentially beyond about
    # e*K harmonics, so the midpoint rule with M = 4096 points is exact up to
    # rounding.
    M = 4096
    relative_error = abs(classical_p2(K_value, 2, M) / classical.exact_p2_two_kicks(K_value) - 1)
    assert relative_error <= summation_roundoff(M) + REFERENCE_TOL


@pytest.mark.parametrize("K_value", K_VALUES)
def test_two_kicks_classical_formula_random(K_value):
    # Statistical error: a random ensemble must agree within 4 standard errors
    # (a 6e-5 chance of failing for a fresh seed; the seed here is fixed).
    theta, p = classical.initial_ensemble(10**6, rng=2026)
    theta, p = classical.evolve(theta, p, K_value, 2)
    mean, standard_error = np.mean(p**2), np.std(p**2) / np.sqrt(p.size)
    assert abs(mean - classical.exact_p2_two_kicks(K_value)) <= 4 * standard_error


def test_three_kicks_converge_to_classical_like_hbar_squared():
    # Classical value at n = 3, converged in the number of grid points.
    c_coarse, c_fine = classical_p2(K, 3, 10**5), classical_p2(K, 3, 2 * 10**5)
    assert abs(c_fine / c_coarse - 1) <= summation_roundoff(2 * 10**5)
    c = c_fine

    q = {hbar: quantum_p2(K, hbar, 3) for hbar in (0.02, 0.01, 0.005)}

    # Known physical correction: <p^2> is even in hbar (hbar -> -hbar only
    # complex-conjugates the evolution), so q(hbar) = c + a hbar^2 + b hbar^4 + ...
    # The observed order from three halvings must be 2; leading errors of order
    # hbar or hbar^3 would give about 1 or 3, and the hbar^4 term shifts it by
    # only a few hundredths at these hbar.
    order = np.log2((q[0.02] - q[0.01]) / (q[0.01] - q[0.005]))
    assert 1.8 <= order <= 2.2, f"observed order {order:.3f}"

    # Richardson extrapolation removes the hbar^2 term and leaves O(hbar^4).
    # The difference between two successive extrapolations estimates that
    # leftover (it is about 15 times the error of the finer one).
    richardson_fine = q[0.005] + (q[0.005] - q[0.01]) / 3
    richardson_coarse = q[0.01] + (q[0.01] - q[0.02]) / 3
    assert abs(richardson_fine - c) <= abs(richardson_fine - richardson_coarse)
