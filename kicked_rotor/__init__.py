"""Classical and quantum kicked rotor.

Modules:
    classical    standard map, its inverse, ensembles, orbits and phase-portrait
                 starts, Lyapunov exponent
    quantum      exact one-period (Floquet) evolution with the FFT, its inverse,
                 grids and momentum eigenstates (run it with -m for a demo)
    observables  norm, fidelity, momentum distribution, <p^2>, edge
                 probability, fit of the localization length
    plotting     shared figure style and helpers

Import the module you need, for example ``from kicked_rotor import classical``.
Nothing is imported here, so loading the physics code never loads Matplotlib.
"""
