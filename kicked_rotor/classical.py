"""Classical kicked rotor: the Chirikov standard map.

Units: the kick period is 1 and the Hamiltonian is

    H = p^2/2 + K cos(theta) * sum_n delta(t - n).

The kick changes p by -dH/dtheta = K sin(theta); between kicks the rotor turns
freely by p. So one period, from just before kick n to just before kick n + 1, is

    p_{n+1}     = p_n + K sin(theta_n)
    theta_{n+1} = theta_n + p_{n+1}    (mod 2 pi)

The angle is always kept in [0, 2 pi). The momentum is never wrapped in the
dynamics, because its unbounded growth is what we want to measure. It is folded
into [-pi, pi) only for drawing phase portraits (see ``wrap_centered``). That
folding is allowed because the map is unchanged by p -> p + 2 pi.

Every function takes NumPy arrays of any shape (one entry per orbit) and works
elementwise, so a whole ensemble moves at once. There is no hidden state:
everything a function needs is passed in as an argument.
"""

import numpy as np

TWO_PI = 2.0 * np.pi


def wrap_angle(theta):
    """Wrap angles into [0, 2 pi).

    ``np.mod`` alone can return exactly 2 pi for a tiny negative input, because
    2 pi - 1e-17 rounds to 2 pi. Those values are set to 0 so the result is
    always in [0, 2 pi).
    """
    wrapped = np.mod(theta, TWO_PI)
    return np.where(wrapped >= TWO_PI, 0.0, wrapped)


def wrap_centered(x):
    """Wrap values into [-pi, pi).

    Used to fold the momentum into one cell of the torus for phase portraits.
    """
    return wrap_angle(np.asarray(x) + np.pi) - np.pi


def standard_map(theta, p, K):
    """Advance one kick period. Returns ``(theta_next, p_next)``."""
    p_next = p + K * np.sin(theta)
    theta_next = wrap_angle(theta + p_next)
    return theta_next, p_next


def inverse_standard_map(theta, p, K):
    """Go back one kick period, undoing ``standard_map``.

    Solving the map for the earlier state gives
        theta_n = theta_{n+1} - p_{n+1}   (mod 2 pi)
        p_n     = p_{n+1} - K sin(theta_n)
    Returns ``(theta_prev, p_prev)``.
    """
    theta_prev = wrap_angle(theta - p)
    p_prev = p - K * np.sin(theta_prev)
    return theta_prev, p_prev


def _as_state(theta, p):
    """Float copies of (theta, p), broadcast to a common shape, theta wrapped."""
    theta = wrap_angle(np.asarray(theta, dtype=float))
    p = np.asarray(p, dtype=float)
    theta, p = np.broadcast_arrays(theta, p)
    return theta.copy(), p.copy()


def evolve(theta, p, K, n_kicks, backward=False):
    """Apply ``n_kicks`` periods and return the final ``(theta, p)``.

    With ``backward=True`` the inverse map is applied instead, so
    ``evolve(*evolve(theta, p, K, n), K, n, backward=True)`` returns to the
    start up to roundoff.
    """
    step = inverse_standard_map if backward else standard_map
    theta, p = _as_state(theta, p)
    for _ in range(n_kicks):
        theta, p = step(theta, p, K)
    return theta, p


def trajectory(theta0, p0, K, n_kicks):
    """Full history of each orbit.

    Returns ``(theta, p)``, each with shape ``(n_kicks + 1,) + shape``, where
    ``shape`` is the broadcast shape of the inputs. Row ``n`` holds the state
    just before kick ``n`` (row 0 is the initial state). The momentum is not
    wrapped.
    """
    theta, p = _as_state(theta0, p0)
    theta_hist = np.empty((n_kicks + 1,) + theta.shape)
    p_hist = np.empty_like(theta_hist)
    theta_hist[0], p_hist[0] = theta, p
    for n in range(n_kicks):
        theta, p = standard_map(theta, p, K)
        theta_hist[n + 1], p_hist[n + 1] = theta, p
    return theta_hist, p_hist


def symmetry_line_starts(n_per_line, offset=1e-3):
    """Starting points spread along the four symmetry lines of the map.

    The standard map is a product of two involutions (maps that are their own
    inverse): T = I2 I1, with
        I1(theta, p) = (-theta, p + K sin(theta)),   I2(theta, p) = (p - theta, p).
    I1 leaves the lines theta = 0 and theta = pi fixed, and I2 leaves the lines
    theta = p/2 and theta = p/2 + pi fixed. Symmetric periodic orbits, and the
    islands around them, cross these lines, so starts spread along them reach
    the main island chains as well as the chaotic regions in between. (For
    example, at K = 5 the largest islands surround a period-2 orbit on the
    line theta = p/2 + pi, which the lines theta = 0 and pi miss.)

    Along each line the momenta are the midpoints of ``n_per_line`` equal steps
    across [-pi, pi), shifted by ``offset`` so no start sits exactly on the
    fixed point at the origin.

    Returns ``(theta0, p0)``, each of length ``4 * n_per_line``, ordered by
    line: theta = pi, theta = 0, theta = p/2, theta = p/2 + pi.
    """
    p_line = -np.pi + (np.arange(n_per_line) + 0.5) * TWO_PI / n_per_line + offset
    theta0 = np.concatenate(
        [np.full(n_per_line, np.pi), np.zeros(n_per_line), p_line / 2, p_line / 2 + np.pi]
    )
    p0 = np.tile(p_line, 4)
    return wrap_angle(theta0), p0


def initial_ensemble(n_particles, p0=0.0, sampling="random", rng=None):
    """Classical ensemble matched to the quantum initial state m = 0.

    Every momentum equals ``p0`` and the angles are uniform on [0, 2 pi).

    sampling="random"
        Independent uniform angles. ``rng`` must be a seed (int) or a
        ``numpy.random.Generator``, so every run is reproducible.
    sampling="grid"
        The midpoint grid theta_j = 2 pi (j + 1/2) / n_particles. It is
        deterministic, and averaging over it is the trapezoid rule, which is
        extremely accurate for smooth periodic functions of theta (so for
        averages at early times).

    Returns ``(theta, p)``, each of length ``n_particles``.
    """
    if sampling == "grid":
        theta = TWO_PI * (np.arange(n_particles) + 0.5) / n_particles
    elif sampling == "random":
        if rng is None:
            raise ValueError("random sampling needs rng=<seed or Generator> so runs are reproducible")
        theta = np.random.default_rng(rng).uniform(0.0, TWO_PI, n_particles)
    else:
        raise ValueError(f"unknown sampling {sampling!r}; use 'random' or 'grid'")
    p = np.full(n_particles, float(p0))
    return theta, p


def torus_ensemble(n_particles, rng):
    """Points spread uniformly over the torus: theta and p both uniform on [0, 2 pi).

    The uniform measure on the torus is invariant under the map, so averages
    over this ensemble have no start-up transient. That makes it the natural
    ensemble for measuring the diffusion rate D = <(p_n - p_0)^2> / n. (Points
    that land in islands of regular motion stay there and do not diffuse.)
    ``rng`` must be a seed (int) or a ``numpy.random.Generator``.
    """
    if rng is None:
        raise ValueError("torus_ensemble needs rng=<seed or Generator> so runs are reproducible")
    rng = np.random.default_rng(rng)
    return rng.uniform(0.0, TWO_PI, n_particles), rng.uniform(0.0, TWO_PI, n_particles)


def exact_p2_two_kicks(K):
    """Exact <p^2> after two kicks for the ensemble p0 = 0, theta0 uniform.

    <p^2>_2 = K^2 [1 + J0(K) - J2(K) - J2(2K)/2] (derived in PLAN.md note 1).
    After one kick the exact value is simply K^2/2.
    """
    from scipy.special import jv

    return K**2 * (1 + jv(0, K) - jv(2, K) - 0.5 * jv(2, 2 * K))


def lyapunov_exponent(theta0, p0, K, n_kicks, n_transient=0):
    """Largest Lyapunov exponent of each orbit, measured over a finite run.

    A small displacement (dtheta, dp) is carried along the orbit by the
    linearized map
        dp'     = dp + K cos(theta) dtheta
        dtheta' = dtheta + dp'
    and rescaled to length 1 after every kick. The exponent is the average of
    log(growth factor) per kick.

    The first ``n_transient`` kicks are run but not counted. This lets the
    displacement line up with the most unstable direction. It also helps tell
    regular from chaotic orbits: on a regular orbit displacements grow only
    linearly, so the log growth from kick n to kick 2n is about ln 2, giving
    roughly ln(2) / n_kicks when ``n_transient == n_kicks``. On a chaotic orbit
    the result approaches the true exponent lambda > 0.

    Returns an array with the broadcast shape of ``theta0`` and ``p0``.
    """
    if n_kicks < 1:
        raise ValueError("n_kicks must be at least 1")
    theta, p = _as_state(theta0, p0)
    d_theta = np.full(theta.shape, np.sqrt(0.5))
    d_p = np.full(theta.shape, np.sqrt(0.5))
    log_growth = np.zeros(theta.shape)
    for n in range(n_transient + n_kicks):
        d_p = d_p + K * np.cos(theta) * d_theta
        d_theta = d_theta + d_p
        theta, p = standard_map(theta, p, K)
        norm = np.hypot(d_theta, d_p)
        d_theta = d_theta / norm
        d_p = d_p / norm
        if n >= n_transient:
            log_growth += np.log(norm)
    return log_growth / n_kicks
