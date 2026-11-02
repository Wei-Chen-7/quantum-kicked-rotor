"""Shared constants and the roundoff rule for the validation tests.

CLAUDE.md asks that every tolerance have a stated reason: machine precision,
the statistical error of an ensemble, or a known physical correction. The
machine-precision rule used throughout is defined here; the other tolerances
are explained where they are used.
"""

import numpy as np

EPS = np.finfo(float).eps  # 2.2e-16
K = 5.0  # the project's default kick strength

# The exact reference formulas use scipy's Bessel functions, which are accurate
# to about 1e-15; a few of them are combined with a little arithmetic.
REFERENCE_TOL = 1e-14


def fft_roundoff(n_kicks, N):
    """Machine-precision allowance after ``n_kicks`` periods on N states.

    Each period does one FFT and one inverse FFT. A length-N FFT has a relative
    rounding error of order eps * log2(N), and errors from successive kicks add
    at most linearly. So n_kicks * log2(N) * eps bounds the relative error of
    the state and of the observables computed from it. (In the measurements
    made for Step 3 the errors came out at about a tenth of this.)
    """
    return np.asarray(n_kicks) * np.log2(N) * EPS


def summation_roundoff(M):
    """Machine-precision allowance for an average over M numbers.

    NumPy sums with pairwise summation, whose relative error grows like
    eps * log2(M); a few more eps cover evaluating the summands.
    """
    return (np.log2(M) + 4) * EPS
