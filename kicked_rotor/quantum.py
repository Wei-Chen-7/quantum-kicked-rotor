"""Quantum kicked rotor: exact one-period (Floquet) evolution with the FFT.

Units match classical.py: the kick period is 1, the Hamiltonian is

    H = p^2/2 + K cos(theta) * sum_n delta(t - n),

and hbar is an effective Planck constant. The angular momentum is p = hbar m
with integer m, and <theta|m> = exp(i m theta) / sqrt(2 pi).

How a state is stored
---------------------
A state is a complex array psi of length N on the angle grid
theta_j = 2 pi j / N (j = 0, ..., N - 1), normalized so sum_j |psi_j|^2 = 1.
Its momentum amplitudes are the unitary FFT

    c_m = (1 / sqrt(N)) sum_j psi_j exp(-i m theta_j),

listed in the order of np.fft.fftfreq(N, d=1/N): m = 0, 1, ..., N/2 - 1,
-N/2, ..., -1. The FFT treats m as periodic mod N, so N must be large enough
that the state never gets near |m| = N/2 (observables.edge_probability checks
this).

One period
----------
From just before kick n to just before kick n + 1:

    1. kick:  psi(theta) -> exp(-i K cos(theta) / hbar) psi(theta)
    2. FFT to the momentum amplitudes c_m
    3. free rotation for one period:  c_m -> exp(-i hbar m^2 / 2) c_m
    4. inverse FFT back to psi(theta)

The kick is instantaneous, so these four steps are the exact evolution over one
period: there is no time step and no splitting error. The only approximation
is keeping N momentum states. The inverse period undoes the steps in reverse
order.

Run ``python -m kicked_rotor.quantum`` from the repo root for a short demo.
"""

import argparse
import time

import numpy as np

TWO_PI = 2.0 * np.pi

# Number of momentum states for the physics runs. Validation check 4 chose it:
# at K = 5 and hbar = 0.25 (the hardest case for hbar in [0.25, 1], since the
# localization length grows like 1/hbar^2), 1000 kicks keep the probability
# beyond |m| = 3N/8 below 1e-26, and doubling N changes <p^2> only at roundoff
# level. 2^13 also passes; 2^12 does not.
DEFAULT_N = 2**14


def angles(N):
    """Angle grid theta_j = 2 pi j / N, for j = 0, ..., N - 1."""
    return TWO_PI * np.arange(N) / N


def integer_momenta(N):
    """Integer momenta m in FFT order, from np.fft.fftfreq(N, d=1/N). p = hbar m."""
    return np.rint(np.fft.fftfreq(N, d=1.0 / N)).astype(np.int64)


def momentum_eigenstate(N, m0=0):
    """The state |m0> on the grid: psi_j = exp(i m0 theta_j) / sqrt(N)."""
    if not -(N // 2) <= m0 <= (N - 1) // 2:
        raise ValueError(f"m0 = {m0} is outside the grid of N = {N} momentum states")
    return np.exp(1j * m0 * angles(N)) / np.sqrt(N)


def to_momentum(psi):
    """Momentum amplitudes c_m (FFT order) of the angle-grid state psi."""
    return np.fft.fft(psi, norm="ortho")


def to_angle(c):
    """Angle-grid state psi from momentum amplitudes c_m (FFT order)."""
    return np.fft.ifft(c, norm="ortho")


def floquet_factors(N, K, hbar):
    """The two phase factors of one period, returned as ``(kick, free)``.

    kick[j] = exp(-i K cos(theta_j) / hbar) multiplies psi on the angle grid.
    free[k] = exp(-i hbar m_k^2 / 2) multiplies the momentum amplitudes.
    """
    kick = np.exp(-1j * (K / hbar) * np.cos(angles(N)))
    m = integer_momenta(N).astype(float)
    free = np.exp(-0.5j * hbar * m**2)
    return kick, free


def floquet_step(psi, kick, free):
    """One period: kick, FFT, free rotation, inverse FFT.

    NumPy's fft (no scaling) and ifft (scaling 1/N) cancel, so the norm of psi
    is unchanged.
    """
    return np.fft.ifft(free * np.fft.fft(kick * psi))


def inverse_floquet_step(psi, kick, free):
    """Undo one period: FFT, undo the free rotation, inverse FFT, undo the kick."""
    return np.conj(kick) * np.fft.ifft(np.conj(free) * np.fft.fft(psi))


def evolve(psi, K, hbar, n_kicks, backward=False):
    """Apply ``n_kicks`` periods to psi and return the final state.

    With ``backward=True`` the periods are undone instead, so
    ``evolve(evolve(psi, K, hbar, n), K, hbar, n, backward=True)`` gives back
    psi up to roundoff.
    """
    psi = np.asarray(psi, dtype=complex)
    kick, free = floquet_factors(psi.size, K, hbar)
    step = inverse_floquet_step if backward else floquet_step
    for _ in range(n_kicks):
        psi = step(psi, kick, free)
    return psi


def iterate(psi, K, hbar, n_kicks, backward=False):
    """Yield the state after each of ``n_kicks`` periods (n = 1, 2, ..., n_kicks).

    Handy for recording observables along the way, for example
        [observables.mean_p2(state, hbar) for state in iterate(psi0, K, hbar, 100)]
    """
    psi = np.asarray(psi, dtype=complex)
    kick, free = floquet_factors(psi.size, K, hbar)
    step = inverse_floquet_step if backward else floquet_step
    for _ in range(n_kicks):
        psi = step(psi, kick, free)
        yield psi


def exact_p2_two_kicks(K, hbar):
    """Exact <p^2> after two kicks, starting in m = 0 (derived in PLAN.md note 1).

    <p^2>_2 = K^2 [1 + cos(hbar/2) (J0(K_q) - J2(K_q)) - J2(K_2)/2], with
    K_q = (2K/hbar) sin(hbar/2) and K_2 = (2K/hbar) sin(hbar). After one kick
    the exact value is simply K^2/2.
    """
    from scipy.special import jv

    k_q = 2 * K / hbar * np.sin(hbar / 2)
    k_2 = 2 * K / hbar * np.sin(hbar)
    return K**2 * (1 + np.cos(hbar / 2) * (jv(0, k_q) - jv(2, k_q)) - 0.5 * jv(2, k_2))


def _exact_p2(n, K, hbar):
    """Exact <p^2> after n = 1 or 2 kicks from m = 0, or None for other n."""
    if n == 1:
        return K**2 / 2
    if n == 2:
        return exact_p2_two_kicks(K, hbar)
    return None


def main(argv=None):
    """Demo: norm and <p^2> for the first few kicks, starting in m = 0."""
    parser = argparse.ArgumentParser(description="Demo: the first kicks of the quantum kicked rotor, from m = 0.")
    parser.add_argument("--K", type=float, default=5.0, help="kick strength (default 5)")
    parser.add_argument("--hbar", type=float, default=0.5, help="effective Planck constant (default 0.5)")
    parser.add_argument("--N", type=int, default=DEFAULT_N, help="number of momentum states (default 2^14)")
    parser.add_argument("--kicks", type=int, default=10, help="number of kicks (default 10)")
    args = parser.parse_args(argv)

    # Imported here because observables itself imports this module.
    from kicked_rotor import observables

    start = time.perf_counter()
    psi = momentum_eigenstate(args.N, 0)
    print(f"Quantum kicked rotor: the first {args.kicks} kicks, starting in m = 0")
    print(f"K = {args.K:g}, hbar = {args.hbar:g}, N = {args.N} momentum states")
    print()
    print("   n    norm - 1        <p^2>    exact <p^2>   relative difference")
    worst_norm = 0.0
    for n, state in enumerate(iterate(psi, args.K, args.hbar, args.kicks), start=1):
        norm_error = observables.norm(state) - 1.0
        worst_norm = max(worst_norm, abs(norm_error))
        p2 = observables.mean_p2(state, args.hbar)
        exact = _exact_p2(n, args.K, args.hbar)
        line = f"{n:4d}   {norm_error:+.1e}   {p2:12.6f}"
        if exact is not None:
            line += f"   {exact:12.6f}   {p2 / exact - 1:+.1e}"
        print(line)
    print()
    print("Exact values: K^2/2 after one kick, and the two-kick formula in PLAN.md note 1.")
    print(f"Largest |norm - 1|: {worst_norm:.1e}. Run time {time.perf_counter() - start:.2f} s")


if __name__ == "__main__":
    main()
