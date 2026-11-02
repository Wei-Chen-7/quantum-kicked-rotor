# Plan: quantum kicked rotor

Final project for PHY-315 (Quantum Mechanics). This file is the full spec. `CLAUDE.md` has the working rules.

## Decisions so far

- Write-up: a full LaTeX report, plus a Beamer deck built from it. (This fills the "[format]" blank in the original brief.)
- The code lives in this repo and is built in Claude Code cloud sessions. Wei runs it on Windows from a local clone.
- Python on Windows (decided 2026-11-01, revised 2026-11-02): Wei left the choice to Claude. He works in PowerShell, where `git` and `python` work but `conda` is not on the PATH, so the run instructions use a virtual environment in `.venv` inside the clone (`python -m venv .venv`). Commands call `.\.venv\Scripts\python` directly, so nothing needs activating and no PowerShell execution-policy change is needed. The README has the commands; VS Code can use the same `.venv`.
- The brief's opening step (ask for blanks, then give a one-paragraph summary) already happened in the planning chat. Start directly with Step 1.
- Proposed step map (Wei can change it):
  1. Classical code: `kicked_rotor/classical.py` (map, inverse map, ensembles, phase-portrait data, Lyapunov exponent) and the plotting helpers it needs. Demo: experiment 01.
  2. Quantum code: `kicked_rotor/quantum.py` and `kicked_rotor/observables.py`. Demo: a short run that prints the norm and $\langle p^2\rangle$ for the first few kicks.
  3. Validation: `tests/` for checks 1 to 5, using the exact results in the notes below. Figures: experiments 05 and 06.
  4. Physics experiments 02, 03 and 04, plus one or two animations. Experiment 07 only if time allows.
  5. Write-up: report first, then slides.

## Wei's brief

Condensed from his original message, with every requirement kept.

### About the project

- Wei is a physics and math double major. Explain things in plain, conversational language, use LaTeX for math, and avoid em dashes. Deliver write-ups as LaTeX files.
- The professor approved this topic. Wei's first idea (split-step Fourier for wave packets and tunneling) was rejected because it's covered in class, so don't drift back into that material.
- The code may later be extended for a numerical analysis project. Keep it modular so that analysis can be added later. Don't build numerical-analysis-only parts unless Wei asks.

### The approved pitch

Study the quantum kicked rotor, a particle on a ring that gets a sharp kick once every period. Classically, strong kicks make the motion chaotic, and the momentum keeps spreading out like a random walk. Quantum mechanically, the spreading follows the classical behavior at first and then stops. This effect is called dynamical localization, and it's closely related to Anderson localization. Simulate both versions, using the fast Fourier transform to switch between angle and momentum space for the quantum part, and compare how the energy grows over time. Check the code against known exact results, like the quantum resonance case where the energy grows quadratically.

### Physics setup

- $H = p^2/2 + K\cos\theta\,\sum_n \delta(t - n)$, period 1, scaled units with an effective $\hbar$.
- Classical standard map (one period): $p_{n+1} = p_n + K\sin\theta_n$, $\theta_{n+1} = \theta_n + p_{n+1}$. Wrap $\theta$ mod $2\pi$, but never wrap $p$.
- Quantum: $p = \hbar m$ with integer $m$. One period is:
  1. kick: multiply $\psi(\theta)$ by $e^{-iK\cos\theta/\hbar}$
  2. FFT to momentum space
  3. free: multiply $\psi_m$ by $e^{-i\hbar m^2/2}$
  4. inverse FFT

  This is exact, with no splitting error, because the kick is instantaneous.
- Matched initial conditions: the quantum state starts in $m = 0$; the classical ensemble starts with $p = 0$ and $\theta$ uniform on $[0, 2\pi)$.

### Key physics to understand and explain

- Classical chaos for $K$ above about 1, and diffusion $\langle p^2\rangle \approx Dn$ with $D \approx K^2/2$ for large $K$, plus Bessel-function corrections (show them).
- Dynamical localization: quantum $\langle p^2\rangle$ follows the classical growth up to a break time, then saturates. The momentum distribution becomes exponential, $|\psi_m|^2 \sim e^{-2|m|/\xi}$. Compare $\xi$ with the standard estimate $\xi \sim D/\hbar^2$ (check the exact prefactor in the literature).
- Floquet theory and quasienergies, and why a discrete quasienergy spectrum stops the growth.
- The Anderson localization analogy (Fishman, Grempel, Prange): explain it as an analogy, without a full derivation.
- Quantum-classical correspondence: when quantum follows classical and when it breaks away.
- Experiments: dynamical localization was observed with cold atoms (Raizen group, mid 1990s).
- Key references to verify and cite properly: Casati, Chirikov, Izrailev, Ford (1979); Fishman, Grempel, Prange (1982); Moore et al. (1995). Don't invent citations; say so if one can't be confirmed.

### Pitfalls

- Quantum resonance: if $\hbar$ is a rational multiple of $4\pi$, the energy grows forever. Avoid it for localization runs, but use $\hbar = 4\pi$ as a test.
- Classical accelerator modes: keep $K$ away from multiples of $2\pi$ (exact windows in note 7).
- The momentum grid must be wide enough that $\psi_m$ never reaches the edges, because the FFT treats momentum as periodic.
- Suggested starting point: $K = 5$, $\hbar$ between 0.25 and 1 (smaller $\hbar$ means a longer classical-like phase before saturation), $N = 2^{12}$ momentum states, a few hundred to a thousand kicks.

### Code

Python 3, NumPy, Matplotlib; SciPy only if needed.

```
kicked_rotor/
  classical.py    # standard map, ensembles, phase portraits, Lyapunov exponent
  quantum.py      # grids (np.fft.fftfreq(N, d=1/N) for integer m), Floquet step,
                  # inverse step for time reversal
  observables.py  # norm, <p^2>, momentum distribution, fit of xi
  plotting.py     # figures and animations
experiments/
  01_classical_phase_portraits.py   # K = 0.5, 0.97, 2, 5
  02_classical_diffusion.py         # <p^2> vs n against D n
  03_quantum_vs_classical.py        # main result: growth then saturation
  04_momentum_distribution.py       # log plot, straight line, fit xi
  05_resonance_test.py              # hbar = 4*pi
  06_time_reversal.py               # forward n kicks, then backward
  07_bonus_xi_scaling.py            # xi vs K and hbar (only if time)
tests/   # pytest for the validation checks
```

Every experiment saves figures to `figures/` and prints a short summary.

### Validation (must pass before physics runs)

1. The norm stays at 1 to near machine precision.
2. Resonance: $\hbar = 4\pi$ from $m = 0$ gives exactly $\langle p^2\rangle = K^2 n^2/2$.
3. At early times, quantum $\langle p^2\rangle$ matches the classical ensemble.
4. Doubling $N$ changes nothing, and the edge probability stays below about $10^{-10}$.
5. Time reversal: the quantum state returns to the initial state almost exactly. The classical run fails after some tens of kicks because chaos amplifies roundoff. Explain why.

### Deliverables

- Clean figures: phase portraits, quantum vs classical $\langle p^2\rangle$, log-scale momentum distribution, resonance test, time-reversal comparison.
- One or two animations of the momentum distribution evolving over kicks.
- A LaTeX report and a Beamer deck that explain the physics clearly for a QM class.

### Working style

- Build in order: classical code, quantum code, validation, physics experiments, write-up.
- After each step, stop, show Wei how to run it and what he should see, and wait for him before moving on.
- If a check fails, help debug instead of loosening tolerances.

## Notes from planning

These were checked numerically in a scratch script on 2026-11-01 (not part of this repo), unless marked otherwise.

### 1. Exact early-time results (use these for check 3)

The quantum state starts in $|m{=}0\rangle$; the classical ensemble starts at $p_0 = 0$ with $\theta_0$ uniform.

After one kick, both give exactly $\langle p^2\rangle_1 = K^2/2$.

After two kicks:

$$\langle p^2\rangle_2^{\mathrm{cl}} = K^2\left[1 + J_0(K) - J_2(K) - \tfrac{1}{2}J_2(2K)\right]$$

$$\langle p^2\rangle_2^{\mathrm{qm}} = K^2\left[1 + \cos\tfrac{\hbar}{2}\,\bigl(J_0(K_q) - J_2(K_q)\bigr) - \tfrac{1}{2}J_2(K_2)\right],\qquad K_q = \tfrac{2K}{\hbar}\sin\tfrac{\hbar}{2},\quad K_2 = \tfrac{2K}{\hbar}\sin\hbar$$

Derivation sketch (Heisenberg picture): write one period as $U = F\hat K$, with the kick $\hat K = e^{-iK\cos\theta/\hbar}$ first and then the free step $F$. Then $p_2 = p_0 + K\sin\theta_0 + K\sin\theta_1$ and $p_0|0\rangle = 0$. The free step gives $F^\dagger e^{i\theta}F = e^{i\theta}e^{ip}e^{i\hbar/2}$, and $e^{ip}$ shifts $\theta$ by $\hbar$, so $\hat K^\dagger e^{ip}\hat K|0\rangle = e^{iK_q\sin(\theta+\hbar/2)}|0\rangle$. The Jacobi-Anger expansion then gives the cross term $\langle\{\sin\theta_0,\sin\theta_1\}\rangle = \cos\tfrac{\hbar}{2}\,[J_0(K_q) - J_2(K_q)]$, and the same steps with $e^{2i\theta}$ give $\langle\sin^2\theta_1\rangle = \tfrac12 - \tfrac12 J_2(K_2)$. The classical result is the same calculation with $\theta_1 = \theta_0 + K\sin\theta_0$.

What was checked:

- An FFT simulation matches the quantum formula to about $10^{-15}$ (relative) for $K \in \{2, 5, 9.5\}$ and $\hbar \in \{0.25, 0.5, 1, 2, 3\}$.
- The classical formula is reproduced to all printed digits by a uniform grid $\theta_j = 2\pi(j+\tfrac12)/M$ with $M = 4096$ (the trapezoid rule on a smooth periodic function), and within $0.01\sigma$ by $10^6$ random samples.
- Limits: $\hbar \to 0$ gives the classical result, $\hbar = 4\pi$ gives $2K^2$ (matching $K^2n^2/2$), and $\hbar = 2\pi$ gives $0$.

So quantum and classical agree exactly only at $n = 1$. At $n = 2$ they differ at order $\hbar^2$. Relative gap (quantum minus classical) at $n = 2$, $K = 5$:

| $\hbar$ | 0.05 | 0.1 | 0.25 | 0.5 | 1.0 |
|---|---|---|---|---|---|
| gap | -0.05% | -0.18% | -1.03% | -2.77% | +5.5% |

How to use this in check 3: test each code against its own exact formula (quantum to machine precision, classical to its statistical error). For the next few kicks, compare quantum and classical with a tolerance built from the ensemble's standard error plus the expected size of the $O(\hbar^2)$ gap, or run that comparison at a small $\hbar$ such as 0.05. Don't pick a percentage by eye.

### 2. Resonance and antiresonance

- $\hbar = 4\pi$, $K = 5$, $N = 2^{12}$, 1000 kicks: the largest relative error against $K^2n^2/2$ was $2.2\times10^{-13}$, using the plain free phase $e^{-i\hbar m^2/2}$. A tolerance like $10^{-10}$ is safe and still catches any real bug. The resonant state spreads to $|m| \approx nK/\hbar$ (about 400 after 1000 kicks), so keep $nK/\hbar$ well below $N/2$.
- Antiresonance (a proposed extra test for check 2; ask Wei first): at $\hbar = 2\pi$ the free step is $(-1)^m$, which shifts $\theta$ by $\pi$. Consecutive kicks then cancel, and $\langle p^2\rangle$ alternates $K^2/2, 0, K^2/2, 0, \dots$ (confirmed numerically). It's also an exception to "rational multiples of $4\pi$ grow forever".

### 3. Classical diffusion and the Bessel corrections

Let $C_j = \langle\sin\theta_0\sin\theta_j\rangle$, averaged over the uniform (invariant) measure on the torus. Then $D = \lim_{n\to\infty}\langle(p_n-p_0)^2\rangle/n = K^2\bigl[C_0 + 2\sum_{j\ge1}C_j\bigr]$. Derived in planning: $C_0 = \tfrac12$, $C_1 = 0$, $C_2 = -\tfrac12 J_2(K)$, $C_3 = \tfrac12\bigl[J_3^2(K) - J_1^2(K)\bigr]$. Measured at $K = 5$: $C_2 = -0.0229$ (predicted $-0.0233$) and $C_3 = 0.0139$ (predicted $0.0129$).

These terms give the familiar corrected form

$$D \approx \frac{K^2}{2}\left[1 - 2J_2(K) - 2J_1^2(K) + 2J_2^2(K) + 2J_3^2(K)\right],$$

recalled from memory as Rechester and White (1980). Verify the exact form and the citation before using it.

Measured with a uniform torus ensemble ($M = 10^6$, $n = 1000$, statistical error about $\pm0.0014$):

| $K$ | measured $D/(K^2/2)$ | formula above | $1 - 2J_2(K)$ |
|---|---|---|---|
| 5 | 0.990 | 0.963 | 0.907 |
| 10 | 0.648 | 0.623 | 0.491 |
| 15 | 0.928 | 0.911 | 0.917 |
| 22 | 0.821 | 0.762 | 0.737 |
| 30 | 0.864 | 0.861 | 0.843 |

The truncated formula is an asymptotic expansion. At moderate $K$ the higher correlations aren't small (at $K = 5$, $C_4 \approx 0.021$ and $C_5 \approx -0.035$), so it misses by a few percent. The 7% miss at $K = 22$ isn't understood yet; look into it in experiment 02. With Wei's initial condition $p_0 = 0$, the first kicks are correlated differently (see the $n = 2$ formula), so fit $D$ from the slope at late times.

### 4. Grid size and edge probability

At $K = 5$, $D \approx 0.99\,K^2/2 \approx 12.4$, so $D/\hbar^2 \approx 198$ at $\hbar = 0.25$ and about 50 at $\hbar = 0.5$. If $\xi \approx D/\hbar^2$, then at $\hbar = 0.25$ the tail at the grid edge $|m| = N/2$ is about $e^{-2\cdot2048/198} \approx 10^{-9}$ of the peak for $N = 2^{12}$, which fails a $10^{-10}$ edge test. With $N = 2^{13}$ it's about $10^{-18}$. A thousand kicks at $N = 2^{14}$ take under a second, so be generous with $N$ and let check 4 confirm it. Define "edge probability" precisely in the code and the report (for example, the total probability with $|m| > 3N/8$). Check its maximum over time, not just its final value, because the instantaneous distribution fluctuates around its time average.

### 5. Localization length convention

Papers differ by factors of 2, depending on whether they describe the long-time-averaged distribution of a packet started at $m = 0$ or individual Floquet eigenstates, and whether they use the decay of $|\psi|$ or $|\psi|^2$. In Wei's convention ($|\psi_m|^2 \sim e^{-2|m|/\xi}$ for the time-averaged packet), the two candidates are $\xi \approx D/\hbar^2$ and $\xi \approx D/(2\hbar^2)$. From memory (unverified): eigenstates localize with a length of about $D/(2\hbar^2)$ in units of $m$, and the averaged packet about twice that. Settle this from the original papers before experiment 04, and state the convention in the report. Fit $\xi$ to a time average over many kicks after saturation, on a log scale, leaving out the central peak.

### 6. Time reversal (predictions, not yet tested)

- Classical: roundoff ($\epsilon \sim 10^{-16}$ to $10^{-15}$) grows like $e^{\lambda n}$, with $\lambda \approx \ln(K/2) \approx 0.92$ at $K = 5$. The return should fail once the forward run is longer than about $\ln(1/\epsilon)/\lambda \approx 35$ to 40 kicks. (Step 1 measured $\lambda = 0.98$ at $K = 5$, the median finite-time exponent of the chaotic orbits in experiment 01, which narrows this to about 35 to 38 kicks.)
- Quantum: the evolution is linear and unitary, so roundoff never grows exponentially, and the return overlap should stay close to machine precision.
- A good point for the report: the classical phase-space density (Liouville evolution) is also linear and reversible. The irreversibility comes from following individual chaotic trajectories, whose errors grow exponentially.

### 7. Accelerator modes

Stable accelerator modes exist for $2\pi l \le K \le \sqrt{(2\pi l)^2 + 16}$: about $[6.28, 7.45]$ for $l = 1$, $[12.57, 13.19]$ for $l = 2$, and $[18.85, 19.27]$ for $l = 3$. Stay out of these windows, with some margin, in any scan over $K$. $K = 5$ is well clear.

### 8. Keeping clear of the rejected topic

The kick, FFT, free phase, inverse FFT sequence looks like split-step Fourier, but here it's exact. In the write-up, cover it in one short paragraph as a tool, and put the weight on chaos, Floquet theory and localization.

## References to verify

All of these are from memory, and none has been checked yet. Verify each one before citing it.

- G. Casati, B. V. Chirikov, F. M. Izrailev, J. Ford, "Stochastic behavior of a quantum pendulum under a periodic perturbation," in *Stochastic Behavior in Classical and Quantum Hamiltonian Systems*, Lecture Notes in Physics 93 (Springer, 1979), p. 334.
- S. Fishman, D. R. Grempel, R. E. Prange, "Chaos, quantum recurrences, and Anderson localization," Phys. Rev. Lett. 49, 509 (1982).
- F. L. Moore, J. C. Robinson, C. F. Bharucha, B. Sundaram, M. G. Raizen, "Atom optics realization of the quantum delta-kicked rotor," Phys. Rev. Lett. 75, 4598 (1995).
- Classical background: B. V. Chirikov, Phys. Rep. 52, 263 (1979); J. M. Greene, J. Math. Phys. 20, 1183 (1979), for the critical $K_c \approx 0.9716$; A. B. Rechester and R. B. White, Phys. Rev. Lett. 44, 1586 (1980), for the Bessel corrections.
- Localization length and resonances: D. L. Shepelyansky, Phys. Rev. Lett. 56, 677 (1986); F. M. Izrailev, Phys. Rep. 196, 299 (1990); F. M. Izrailev and D. L. Shepelyansky, Theor. Math. Phys. 43, 553 (1980).

## Open questions for Wei

- Should the $\hbar = 2\pi$ antiresonance be added to check 2?
- Report length, slide count or talk length, and the due date (not given yet).
