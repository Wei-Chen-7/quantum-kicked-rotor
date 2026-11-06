# Plan: quantum kicked rotor

Final project for PHY-315 (Quantum Mechanics). This file is the full spec. `CLAUDE.md` has the working rules.

## Decisions so far

- Write-up: a full LaTeX report, plus a Beamer deck built from it. (This fills the "[format]" blank in the original brief.)
- The code lives in this repo and is built in Claude Code cloud sessions. Wei runs it on Windows from a local clone.
- Validation design (Wei agreed 2026-11-02): check 2 includes the $\hbar = 2\pi$ antiresonance; check 3 is exact at $n = 1, 2$ plus an $\hbar^2$-convergence test at $n = 3$ (note 1); check 4 chose $N = 2^{14}$ (note 4). The tests live in `tests/` (one file per check, shared roundoff rule in `tests/common.py`). Each of six deliberately planted bugs (wrong free phase, wrong momentum ordering, wrong classical update, wrong inverse map, broken edge probability, kick strength off by 0.1%) made at least one test fail. Flipping the sign of the kick or free phase changes no observable we compute (a symmetry), so no test can see it.
- Step 4 design (2026-11-02): $D$ for the quantum comparisons is measured, not taken from a formula. Experiments 03 and 04 both use the matched ensemble ($p = 0$, 40000 grid angles) and the slope of $\langle p^2\rangle$ over kicks 200 to 2000. The $D(K)$ scan in experiment 02 uses the invariant torus ensemble (`classical.torus_ensemble`), which has no start-up transient. The break time uses Izrailev's 25% rule. $\xi$ is fitted with Izrailev's shape (`observables.fit_averaged_localization_length`) to a late time average, kicks 2500 to 5000. The animation (from experiment 04, $\hbar = 0.25$) takes frames at kicks $n = 1000(k/80)^2$, so the classical width, which grows like $\sqrt n$, widens at a steady rate on screen. It leaves out classical histogram bins with fewer than 10 particles, so stray particles don't draw a false plateau.
- Experiment 07 design (2026-11-02, built at Wei's request before the write-up). Panels (a) and (b): $K = 5$, single runs at $\hbar = 0.20, 0.22, \dots, 1.00$. Panel (c): $K = 4.5, 5, 9.75, 10.5, 11$, the median over $\hbar = 0.4, 0.5, \dots, 0.8$. Each run lasts $\max(5000, 50\,D_f/\hbar^2)$ kicks and fits Izrailev's shape to the time average over its second half, where $D_f$ is the five-term formula: the measured $D$ differs slightly between computers (note 6), and sizing runs from it would make the quantum results machine-dependent too. The grid starts at $N \ge 40\,D_f/\hbar^2$ and doubles until the edge probability is below $10^{-10}$ and the fit range lies inside $|m| < 3N/8$, up to $2^{16}$; a run that still fails is reported as not localized and left out. $D$ is the matched-ensemble value at each $K$, and $K = 4$, 5.5 and 9.25 are excluded (note 3). A first design averaged five runs at $\hbar(1 + 0.01j)$, $j = -2, \dots, 2$; it was dropped because nearby $\hbar$ values give correlated results, and because $0.7 \times 1.02 = 0.714$ happens to lie within $10^{-7}$ of the resonance $\hbar/4\pi = 5/88$ (note 5).
- Write-up (Step 5, 2026-11-03): the report is a single-column `article` (11 pt, letter, 1 inch margins, so the 6.5-inch figures fit at full width) in Times (`txfonts`, close to the STIX fonts in the figures), split into `report/sections/*.tex` with BibTeX (`unsrtnat`) and DOIs; it includes the figures directly from `figures/`. The deck is Beamer 16:9 at 10 pt with the metropolis inner, outer and color themes (its font theme needs XeLaTeX, so it is left out). The animation is embedded with the `animate` package from palette-PNG frames cut losslessly from the GIF by `slides/make_animation_frames.py`, with the last frame as the poster for viewers that can't play it. Both PDFs are committed so Wei can read them without TeX. Both build with no LaTeX warnings.
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

How to use this in check 3: test each code against its own exact formula (quantum to machine precision, classical to its statistical error). Beyond $n = 2$, read the next paragraph before comparing quantum with classical. Don't pick a percentage by eye.

**Step 2 finding (2026-11-02): after $n = 2$ the gap grows very fast.** Gap (quantum minus classical, relative) at $K = 5$, with the classical values from the midpoint grid (converged to all digits: $\langle p^2\rangle_3 = 32.25286$, $\langle p^2\rangle_4 = 40.07047$). The quantum value at $n = 3$, $\hbar = 0.05$ was confirmed to 13 digits by an independent Bessel-matrix evolution, so this is physics, not a bug:

| $\hbar$ | 0.5 | 0.25 | 0.1 | 0.05 | 0.02 | 0.01 | 0.005 |
|---|---|---|---|---|---|---|---|
| gap at $n = 3$ | +0.5% | -6.4% | -10.8% | -6.3% | -1.4% | -0.35% | -0.09% |
| gap at $n = 4$ | +49% | +29% | +6.6% | +15% | +14% | +11% | -1.4% |

The gap does vanish like $\hbar^2$ once $\hbar$ is small enough, but the $\hbar^2$ coefficient grows very fast with $n$: about $-0.2$ at $n = 2$, $-35$ at $n = 3$, and roughly $-550$ at $n = 4$. So at $\hbar = 0.05$ the $n = 3$ gap is already 6%, and the earlier suggestion to compare at $\hbar = 0.05$ fails. Izrailev (1990), pp. 313 to 314, gives the reason: full quantum-classical correspondence lasts only up to the Ehrenfest time $t_E \sim \ln(1/\hbar)$ (his Eq. 2.2.20), a few kicks, and for $t_E \lesssim t \lesssim t_D$ "only rough statistical properties like diffusion remain." Decided with Wei (2026-11-02): check 3 tests $n = 1$ and $n = 2$ exactly, and at $n = 3$ tests the correspondence limit instead. Since $\langle p^2\rangle$ is even in $\hbar$, $q(\hbar) = c + a\hbar^2 + b\hbar^4 + \dots$; from $\hbar = 0.02, 0.01, 0.005$ the observed order is 1.97 (the test asks for 1.8 to 2.2), and the Richardson value $q(0.005) + [q(0.005) - q(0.01)]/3$ lands within $1.2\times10^{-4}$ of the classical $32.25286$, inside its own error estimate of $2.3\times10^{-3}$. Comparing the diffusion rates over many kicks belongs to experiment 03.

### 2. Resonance and antiresonance

- $\hbar = 4\pi$, $K = 5$, $N = 2^{12}$, 1000 kicks: the largest relative error against $K^2n^2/2$ was $2.2\times10^{-13}$, using the plain free phase $e^{-i\hbar m^2/2}$. A tolerance like $10^{-10}$ is safe and still catches any real bug. The resonant state spreads to $|m| \approx nK/\hbar$ (about 400 after 1000 kicks), so keep $nK/\hbar$ well below $N/2$.
- Antiresonance (part of check 2; Wei agreed on 2026-11-02): at $\hbar = 2\pi$ the free step is $(-1)^m$, which shifts $\theta$ by $\pi$. Consecutive kicks then cancel, and $\langle p^2\rangle$ alternates $K^2/2, 0, K^2/2, 0, \dots$ (confirmed numerically). It's also an exception to "rational multiples of $4\pi$ grow forever". Original source: Casati et al. (1979), p. 349, prove that at $\tau = 2\pi$ the solution returns exactly after two kicks, for any $k$ and any initial state, and call it "anti-resonant"; Moore et al. (1995) saw the sign-flipping early dynamics in experiment.

### 3. Classical diffusion and the Bessel corrections

Let $C_j = \langle\sin\theta_0\sin\theta_j\rangle$, averaged over the uniform (invariant) measure on the torus. Then $D = \lim_{n\to\infty}\langle(p_n-p_0)^2\rangle/n = K^2\bigl[C_0 + 2\sum_{j\ge1}C_j\bigr]$. Derived in planning: $C_0 = \tfrac12$, $C_1 = 0$, $C_2 = -\tfrac12 J_2(K)$, $C_3 = \tfrac12\bigl[J_3^2(K) - J_1^2(K)\bigr]$. Measured at $K = 5$: $C_2 = -0.0229$ (predicted $-0.0233$) and $C_3 = 0.0139$ (predicted $0.0129$).

These terms give a corrected $D$, but the sources keep different terms (checked 2026-11-02):

- **Rechester and White (1980), Eq. (20)** (read). With their small added diffusion $\sigma \to 0$ and converted to our $D = \langle(\Delta p)^2\rangle/n$: $D = \frac{K^2}{2}\left[1 - 2J_2 - 2J_1^2 + 2J_3^2\right]$. Their own $D$ is $\langle(\Delta v)^2\rangle/(2T)$ (their Eq. 15), half of ours, so the printed form is $D = \tfrac12\epsilon^2\left[\tfrac12 - J_2(\epsilon)e^{-\sigma} - J_1^2(\epsilon)e^{-\sigma} + J_3^2(\epsilon)e^{-3\sigma}\right]$. This equals $K^2[C_0 + 2C_2 + 2C_3]$ from the correlations above.
- **Izrailev (1990), Eq. (2.1.9)** (read): $D = \frac{K^2}{2}\{1 - 2J_2(K)[1 - J_2(K)]\}$ for $K \gtrsim 4.5$, citing Rechester and White (1980) and Rechester, Rosenbluth and White (1981). For $K_c < K < 4.5$ he gives $D \approx 0.30(K - K_c)^3$.
- **Rechester, Rosenbluth and White (1981), Eq. (26)** (read): $D/D_{QL} = 1 - 2J_2e^{-\sigma} - 2J_1^2e^{-\sigma} + 2J_3^2e^{-3\sigma} + 2J_2^2e^{-2\sigma}$ plus sums of products of three Bessel functions (some with arguments $|m+1|\epsilon$). The paper says this "extends Eq. (20) of Ref. 1, which included only the first four terms," and its Fig. 2 caption credits the $J_2^2$ term to a third-order path. Same conventions as the 1980 paper ($D_{QL} = \epsilon^2/4$, $D = \langle(\Delta v)^2\rangle/(2T)$), with the map written as $v_{t+1} = v_t - \epsilon\sin x_t$, a sign flip that doesn't change $D$. So **the planning chat's version**, $\frac{K^2}{2}\left[1 - 2J_2 - 2J_1^2 + 2J_2^2 + 2J_3^2\right]$, is exactly the first five terms of this Eq. (26) with $\sigma \to 0$; cite it that way. The three-Bessel terms are hard to read in the scan; transcribe them from the paper if experiment 02 needs them (for example, to chase the miss at $K = 22$).

Measured with a uniform torus ensemble ($M = 10^6$, $n = 1000$, statistical error about $\pm0.0014$), against each version:

| $K$ | measured $D/(K^2/2)$ | Rechester-White 1980 | Izrailev 1990 | planning version | $1 - 2J_2(K)$ |
|---|---|---|---|---|---|
| 5 | 0.990 | 0.958 | 0.911 | 0.963 | 0.907 |
| 10 | 0.648 | 0.494 | 0.620 | 0.623 | 0.491 |
| 15 | 0.928 | 0.908 | 0.920 | 0.911 | 0.917 |
| 22 | 0.821 | 0.727 | 0.772 | 0.762 | 0.737 |
| 30 | 0.864 | 0.848 | 0.855 | 0.861 | 0.843 |

The $J_2^2$ term matters most near $K = 10$, where the 1980 formula misses by 24%. Each truncated formula is an asymptotic expansion. At moderate $K$ the higher correlations aren't small (at $K = 5$, $C_4 \approx 0.021$ and $C_5 \approx -0.035$), so it misses by a few percent. With Wei's initial condition $p_0 = 0$, the first kicks are correlated differently (see the $n = 2$ formula), so fit $D$ from the slope at late times.

Measured in Step 4 (experiment 02, 2026-11-02):

- Matched ensemble ($p = 0$, 40000 grid angles, slope over kicks 200 to 1000): $D/(K^2/2) = 0.995 \pm 0.009$ at $K = 5$ (five-term formula 0.963) and $0.638 \pm 0.006$ at $K = 10$ (0.623).
- Scan over 103 values of $K$ from 4.5 to 30 (torus ensemble, 10000 points per $K$, slope over kicks 100 to 400, typical error 0.019): away from the accelerator windows (more than 1.2 from their edges) the five-term formula misses by 0.041 rms. The misses beyond 4 standard errors are at $K = 8.75$ ($-0.115$), 9.25 ($+0.153$), 9.5 ($+0.070$), 10.75 ($-0.070$) and 22 ($+0.062$). They are real, not heavy-tailed statistics (kurtosis about 3, and $D$ is stable in time). Just outside the windows $D$ is strongly enhanced, up to $D/(K^2/2) = 3.25$ at $K = 7.5$.
- The miss at $K = 22$ is explained by the correlations (400000 torus points, $C_j$ averaged over 4 time origins, error 0.0004). The Bessel results hold: $C_2 = -0.0658$ (formula $-0.0657$) and $C_3 = -0.0024$ ($-0.0025$), and the sum up to lag 3 gives 0.728, the same as the 1980 formula (0.727). The excess comes from $C_4 = +0.030$ (adding 0.12) and a slowly decaying tail whose sign repeats with period 4 in the lag. The running sum $1 + 4\sum_{i \le j}C_i$ is 0.849, 0.833 and 0.823 at $j = 4$, 8 and 16, matching the measured 0.821. The five-term formula has lags up to 3 exactly but only the $J_2^2$ part of lag 4, so it misses wherever lag 4 and beyond are large; the same pattern appears near $K = 9.25$. A possible cause at $K = 22$ is stickiness near the only island a search found there, around a period-4 orbit ($\theta \approx 4.56 \to 1.73 \to 1.78 \to 4.50$) covering about 0.002% of the torus. That part is a guess, not proven.
- Run-to-run check (2026-11-02, after Wei's Windows run flagged a different list: 8.75, 9, 9.25, 9.5, 10, 10.75, 22 and 28.25). Five seeds of the scan spread by 0.94 times the stated errors, so the error bars are honest. $K = 8.75$, 9.25 and 10.75 crossed the 4-error line in all five runs and 9.5 in four. $K = 9$, 22 and 28.25 have mean misses of $+0.056$, $+0.056$ and $+0.059$, right at the line, so they cross it in some runs only. Mean misses over the five runs: $-0.113$ (8.75), $+0.164$ (9.25), $+0.064$ (9.5), $+0.023$ (10), $-0.063$ (10.75).
- Matched ensemble ($p = 0$) against torus ensemble (2026-11-02, 40000 points each, slope over kicks 200 to 2000), checked for experiment 07. They agree within errors at $K = 4.5$, 5 and every tested $K$ from 8.75 to 11.25 except 9.25, where the matched ensemble gives $D/(K^2/2) = 0.718$ and the torus 0.800. Both are steady in time (4000 kicks), so this is not slow convergence: 9.2% of the line $p = 0$ lies in an island and never diffuses (its kurtosis 3.27 matches an 8% non-diffusing fraction), against 0.9% of the torus. That explains the gap between the two ensembles, but not the formula's miss at 9.25: the torus value, which the island barely affects, is itself 0.17 above the formula. At $K = 4$ the two differ too (12.5% and 6.8% of points stuck), and at $K = 5.5$, just below the $l = 1$ accelerator window, the diffusion is anomalous (kurtosis about 300, $D/(K^2/2) = 5.5$ matched and 1.5 torus).

### 4. Grid size and edge probability

At $K = 5$, $D \approx 0.99\,K^2/2 \approx 12.4$, so $D/\hbar^2 \approx 198$ at $\hbar = 0.25$ and about 50 at $\hbar = 0.5$. If $\xi \approx D/\hbar^2$, then at $\hbar = 0.25$ the tail at the grid edge $|m| = N/2$ is about $e^{-2\cdot2048/198} \approx 10^{-9}$ of the peak for $N = 2^{12}$, which fails a $10^{-10}$ edge test. With $N = 2^{13}$ it's about $10^{-18}$. A thousand kicks at $N = 2^{14}$ take under a second, so be generous with $N$ and let check 4 confirm it. Define "edge probability" precisely in the code and the report (for example, the total probability with $|m| > 3N/8$). Check its maximum over time, not just its final value, because the instantaneous distribution fluctuates around its time average.

Measured in Step 3 ($K = 5$, $\hbar = 0.25$, 1000 kicks, largest edge probability over the run): $N = 2^{11}$: $1.5\times10^{-3}$; $2^{12}$: $3.1\times10^{-8}$ (worse than the estimate above, and $\langle p^2\rangle$ shifts by $10^{-8}$ when $N$ doubles); $2^{13}$: $2.8\times10^{-19}$; $2^{14}$: $2.6\times10^{-27}$, with doubling changing $\langle p^2\rangle$ by $10^{-13}$, i.e. roundoff. At $\hbar = 0.5$ and 1 even $2^{11}$ passes. Chosen default: `quantum.DEFAULT_N = 2**14` (1000 kicks take about 1 s). $\langle p^2\rangle$ after 1000 kicks: 1149.9 ($\hbar = 0.25$), 274.7 ($\hbar = 0.5$), 154.07 ($\hbar = 1$).

### 5. Localization length convention (settled from the papers, 2026-11-02)

Wei's convention: $|\psi_m|^2 \sim e^{-2|m|/\xi}$ for the time-averaged distribution of a packet started at $m = 0$. Here $D$ is the classical momentum diffusion rate $\langle(\Delta p)^2\rangle/n$, including its Bessel oscillations, so $D/\hbar^2$ is the diffusion rate in $m$. (Both papers set $\hbar = 1$ and vary the kick period, $T$ or $\tau$. Their period plays the role of our $\hbar$, their $K = kT$ is our $K$, and their $D_n$ is our $D/\hbar^2$.)

- **Floquet eigenstates.** Shepelyansky (1986) defines $l$ by the decay of the amplitude, $u_n \propto e^{-|n|/l}$, so $|u_n|^2 \propto e^{-2|n|/l}$, and finds $l = D/(2\hbar^2)$. The factor $\tfrac12$ is exact for the Lloyd model; quantum standard map numerics give $\langle l/D_n\rangle = 0.57$ with rms spread 0.11 for $5 \le k \le 75$, $1.5 \le K \le 29$. Izrailev (1990), Eq. (2.3.9): $l_\infty \approx D_n/2$, accurate to better than 10% in that range.
- **Time-averaged packet from $m = 0$.** Shepelyansky (1986) writes $\bar f_n \propto e^{-2|n|/l_s}$ and reports $l_s \approx 2l$ (from Chirikov and Shepelyansky's data: $\langle l_s/D_n\rangle = 1.04$ with rms spread 0.20), the difference coming from strong fluctuations of the eigenfunctions. Izrailev (1990), Eqs. (2.3.13) and (2.3.14): $\bar w_n = \frac{1}{2l_s}e^{-2|n|/l_s}\left(1 + 2|n|/l_s\right)$, with $\langle n^2\rangle = l_s^2$ and $l_s \approx 2l_\infty \approx D_n$, checked for $k = 5$ to 120, $\tau \le 1$, $l_s = 9$ to 180.

So in Wei's convention $\xi = l_s \approx D/\hbar^2$: about 198 at $K = 5$, $\hbar = 0.25$ (with $D \approx 12.4$), and about 50 at $\hbar = 0.5$. Eigenstates are localized about twice as tightly. Consequences for experiment 04:

- Fit Izrailev's shape, not a pure exponential. On that shape a pure-exponential fit overestimates $\xi$ by 31% over $\xi/2 \le |m| \le 3\xi$, 21% over $[\xi, 4\xi]$, 13% over $[2\xi, 6\xi]$, and still 7% over $[4\xi, 10\xi]$ (computed 2026-11-02). Fit a time average over many kicks after saturation.
- $\langle m^2\rangle = \xi^2$ for that shape, so the saturated energy is $\langle p^2\rangle_\infty \approx \hbar^2\xi^2 \approx D^2/\hbar^2$ (Izrailev's Eq. 2.3.15 writes $E_\infty \approx \hbar^2 l_s^2/2$).
- At larger $\hbar$, Izrailev (p. 321) says the quantum diffusion rate follows Eq. (2.1.9) with $K$ replaced by $K_q = 2k\sin(\tau/2) = (2K/\hbar)\sin(\hbar/2)$, the same $K_q$ as in note 1.
- Reference curve: Izrailev's Fig. 1 (after Chirikov, Izrailev and Shepelyansky 1981) shows the energy for exactly our case ($K = 5$, $\hbar = 0.25$, start in $m = 0$) over 2000 kicks. Read off the figure, it levels off near $\langle p^2\rangle \approx 1400$, below the $\hbar^2\xi^2 \approx 2500$ from $\xi = D/\hbar^2$; compare in experiment 03.
- Watch the convention when quoting experiments: Moore et al. (1995) define $\xi$ as the $1/e$ point of the probability, $P \propto e^{-|n|/\xi}$, which is half of Wei's $\xi$ (see References).

Measured in Step 4 (experiments 03 and 04, 2026-11-02, $K = 5$, $N = 2^{14}$). Both use the same measured $D = 12.51 \pm 0.10$, that is $D/(K^2/2) = 1.001$ (matched ensemble, 40000 grid angles, slope over kicks 200 to 2000).

- Break time by Izrailev's 25% rule: $t^* = 20$, 12 and 4 kicks at $\hbar = 0.25$, 0.5 and 1, which is 0.10, 0.24 and 0.32 of $\xi = D/\hbar^2$. Rescaled as $\langle m^2\rangle/\xi^2$ against $n/\xi$, the three curves roughly collapse.
- $\langle p^2\rangle$ after 2000 kicks: 1315, 307 and 78; the mean over kicks 1000 to 2000 is 1305, 300 and 119, or 0.52, 0.48 and 0.76 of $D^2/\hbar^2$. After the break the energy keeps swinging slowly. In a one-off 40000-kick run at $\hbar = 0.25$, the means over blocks of 2500 kicks ranged from 1100 to 1960 around 1495 (17% rms), with no upward trend. So the Izrailev Fig. 1 read-off (about 1400) agrees within these swings.
- Fitted $\xi$ (Izrailev's shape over $5 \le |m| \le m_{\rm tail}$, where the time average falls to $10^{-14}$ of its peak), as a fraction of $D/\hbar^2$: 0.80 and 0.87 at $\hbar = 0.25$ (averaged over kicks 1000 to 2000 and 2500 to 5000), and 0.84 and 0.90 at $\hbar = 0.5$. The far tail fills in slowly. In the 40000-kick run the fit levels off at 0.90, 0.92 and 0.91 (windows 5000 to 10000, 10000 to 20000 and 20000 to 40000), while the average at $|m| \approx 3000$ (about $15\xi$) rises from $10^{-17.7}$ (kicks 1000 to 2000) to about $10^{-14.6}$ from kick 5000 on. A 10000-kick check at $\hbar = 0.5$ gave 0.92 for kicks 5000 to 10000. Converged values of about 0.91 and 0.92 are within the papers' $1.04 \pm 0.20$. At $\hbar = 1$ the fit gives 0.98, but it depends strongly on the range (12.2 over the full range, 17 to 18.5 over half of it), so it is only rough.
- The time average has a narrow central peak that Izrailev's shape doesn't have. $|m| < 5$ holds 9% of the probability at $\hbar = 0.25$ (3% under the fitted shape) and 19% at $\hbar = 0.5$ (6%); at $\hbar = 0.5$ the fitted shape carries only 64% of the probability. This, together with the data sitting a little below the fitted shape at intermediate $|m|$, is why $\sqrt{\langle m^2\rangle}$ (154 and 34) comes out below the fitted $\xi$ (175 and 45). Not explained yet. It may be the enhanced return probability expected with time-reversal symmetry, but that hasn't been checked against a source, so don't put it in the report as fact.
- The instantaneous distribution fluctuates strongly around the time average. At $\hbar = 0.25$, about 1.2% of the probability builds up near $|m| \approx 680$ between kicks 1000 and 1700 and then drains away; it shows as side peaks in the last frames of the animation. Averaged over kicks 2500 to 5000 that band is unremarkable (0.57%, against 0.85% in the band just inside it). It is a slow beat between Floquet eigenstates with nearly equal quasienergies, the same kind of fluctuation that makes $\langle p^2\rangle$ swing.

Measured in experiment 07 (bonus, 2026-11-02):

- Over 41 runs at $K = 5$ with $\hbar$ from 0.20 to 1.00, $\xi\hbar^2/D$ has median 1.19, rms spread 0.25 and range 0.83 to 2.17 (Shepelyansky: mean 1.04, rms spread 0.20). So $\xi \approx 1.2\,D/\hbar^2$ holds on average while $\xi$ itself changes by a factor of about 27.
- The ratio wanders in broad bands rather than scattering independently: medians 0.89 for $\hbar < 0.3$, 1.28 for 0.3 to 0.5, 1.24 for 0.5 to 0.75, 1.09 for 0.75 to 1. In an earlier scan with steps of 0.01 (time average over kicks 5000 to 10000), neighbors 0.01 apart had a lag-1 autocorrelation of 0.57. So nearby $\hbar$ values are not independent samples. The cause of the bands is not known. The quantum correction $K \to K_q$ changes $D$ by only about 2% at $\hbar = 0.35$, far too little. The three values in experiment 04 (0.87, 0.90 and 0.98 at $\hbar = 0.25$, 0.5 and 1) lie on the low side of this distribution: $\hbar = 0.25$ is in the low band, and $\hbar = 0.5$ sits in a local dip (0.94 and 0.90 at $\hbar = 0.48$ and 0.50).
- Over $K$ (medians over $\hbar = 0.4$ to 0.8), $\xi\hbar^2/(K^2/2) = 1.15$, 1.17, 0.60, 0.81, 0.86 at $K = 4.5$, 5, 9.75, 10.5, 11, against the classical $D/(K^2/2) = 0.86$, 1.00, 0.60, 0.71, 0.82. The localization length follows the Bessel dip of $D(K)$; the medians of $\xi\hbar^2/D$ are 1.34, 1.17, 1.01, 1.14, 1.06.
- Near-resonant $\hbar$ can fail to localize. At $\hbar = 0.714$, $\hbar/4\pi$ is within about $10^{-7}$ of $5/88$, so the free phase repeats with period 88 out to $|m| \sim 1000$. In 10000 kicks the state reached the edge of a 2048-state grid (edge probability $2\times10^{-5}$), and the fit gave a meaningless $4.8\,D/\hbar^2$. Experiment 07's grid-doubling rule reports such a run as not localized instead.

### 6. Time reversal (predictions; the quantum part checked in Step 2)

- Classical: roundoff ($\epsilon \sim 10^{-16}$ to $10^{-15}$) grows like $e^{\lambda n}$, with $\lambda \approx \ln(K/2) \approx 0.92$ at $K = 5$. The return should fail once the forward run is longer than about $\ln(1/\epsilon)/\lambda \approx 35$ to 40 kicks. (Step 1 measured $\lambda = 0.98$ at $K = 5$, the median finite-time exponent of the chaotic orbits in experiment 01, which narrows this to about 35 to 38 kicks.) Izrailev (1990), p. 315, citing Shepelyansky, Physica D 8, 208 (1983) and Chirikov, Izrailev and Shepelyansky, Physica D 33, 77 (1988), reports that at $K = 5$ with roundoff around $10^{-12}$ the classical map retraces its path for only about 10 kicks, shorter than $\ln(10^{12})/0.98 \approx 28$. Experiment 06 (Step 3) confirms the prediction for our roundoff: over 4000 orbits ($p_0 = 0$, $\theta_0$ uniform) the median return distance follows $1.5\times10^{-16}e^{0.972n}$ (with $\lambda = 0.981$ measured on the same orbits) and first exceeds 1 at $n = 38$, against $\ln(1/\varepsilon)/\lambda = 36.7$. The gap to Izrailev's 10 kicks is still unexplained; his failure criterion may differ.
- Quantum: the evolution is linear and unitary, so roundoff never grows exponentially, and the return overlap should stay close to machine precision. Step 2 check: $K = 5$, $\hbar = 0.25$, $N = 2^{14}$, 1000 kicks forward and 1000 back returned to the initial state within $3\times10^{-15}$ in every component. The norm drifts steadily upward by about $3\times10^{-16}$ per kick ($2.9\times10^{-13}$ after 1000 kicks), so the fidelity came out $1 + 6\times10^{-13}$; compare states directly or normalize first.
- Seen for free in Step 4 (2026-11-02): Wei's Windows run reproduced every quantum number to all printed digits, but every classical one differed, for example $D = 12.44 \pm 0.10$ there against $12.51 \pm 0.10$ here. The only operation in the classical map that IEEE arithmetic doesn't pin down is the sine, whose last bit can differ between NumPy builds (here `np.sin` matched `math.sin` on 200000 test values), and chaos amplifies that. Shifting every starting angle by one ulp here moves $D$ by about the stated error (spread 0.078 over 9 shifts), and 32 independent random-angle ensembles spread by 0.103 against a stated 0.099 (their mean, $D = 12.560 \pm 0.018$, is the most precise value we have for this ensemble and window). In experiment 02 the running sums of the lag correlations at $K = 22$ agree between the machines through lag 8 and differ at lag 16, as expected from $\lambda \approx \ln(K/2) \approx 2.4$ per kick there (roundoff reaches order 1 after about 15 kicks). A good point for the report: same code, two machines, identical quantum results, different classical orbits.
- A good point for the report: the classical phase-space density (Liouville evolution) is also linear and reversible. The irreversibility comes from following individual chaotic trajectories, whose errors grow exponentially.

### 7. Accelerator modes

Stable accelerator modes exist for $2\pi l \le K \le \sqrt{(2\pi l)^2 + 16}$: about $[6.28, 7.45]$ for $l = 1$, $[12.57, 13.19]$ for $l = 2$, and $[18.85, 19.27]$ for $l = 3$. Stay out of these windows, with some margin, in any scan over $K$. $K = 5$ is well clear. Source: Chirikov (1979), Eq. (5.54), $(2\pi n)^2 < K^2 < (2\pi n)^2 + 16$; inside a window $I = I_0 + 2\pi n t$ grows ballistically (Eq. 5.55), and the window width is about $8/K$. It also follows in two lines: a fixed point on the torus needs $K\sin\theta^* = 2\pi l$, and stability needs $|2 + K\cos\theta^*| < 2$. Chirikov (p. 332) adds period-2 accelerator orbits, $K(\sin\theta_1 + \sin\theta_2) = 2\pi l$, which occur at many more values of $K$ the larger $K$ is.

### 8. Keeping clear of the rejected topic

The kick, FFT, free phase, inverse FFT sequence looks like split-step Fourier, but here it's exact. In the write-up, cover it in one short paragraph as a tool, and put the weight on chaos, Floquet theory and localization.

## References

Status as of 2026-11-02. "Read" means the paper itself was read (Wei supplied the PDFs of the papers marked "read"; they are not in the repo). "Listed" means only reference lists or a bibliographic record were checked, which confirms where the paper is but not what it says. Verify the rest before citing.

- G. Casati, B. V. Chirikov, F. M. Izrailev, J. Ford, "Stochastic behavior of a quantum pendulum under a periodic perturbation," in *Stochastic Behavior in Classical and Quantum Hamiltonian Systems*, edited by G. Casati and J. Ford, Lecture Notes in Physics 93 (Springer, Berlin, 1979), pp. 334-352, DOI 10.1007/BFb0021757. Read: title page, abstract and p. 349. The chapter's own byline is "G. Casati; B. V. Chirikov and F. M. Izraelev; Joseph Ford" (so this order is right; the lists in Shepelyansky 1986 and Izrailev 1990 swap Ford and Izrailev). Abstract: the first numerical study of the quantum kicked rotor, finding diffusive-like linear energy growth and a Gaussian momentum distribution as in the classical motion, plus "several surprising distinctions." Page 349: at $\tau = 4\pi$ the free step is the identity, so $\psi(t) = \psi(0)e^{itk\cos\theta}$ and $\langle E\rangle \propto t^2$ (their Eqs. 37 and 38; also seen at $\tau = 4\pi/m$ for $m = 3, 4, 8, 32$); at $\tau = 2\pi$ the solution is strictly periodic with period 2 ("anti-resonant"); and, after the Como conference, the first report of a break time $t_B$ after which diffusion slows, with $t_B \propto k$ at fixed $K$ (run with $k = 40$, $\tau = 1/8$).
- S. Fishman, D. R. Grempel, R. E. Prange, "Chaos, quantum recurrences, and Anderson localization," Phys. Rev. Lett. 49, 509 (1982), pp. 509-512, DOI 10.1103/PhysRevLett.49.509. Read. The mapping: with $U(\theta) = -\tan[kV(\theta)/2]$ and $T_m = \tan(E_m/2)$, $E_m = \omega - \tau m^2/2$, quasienergy states obey $T_m u_m + \sum_{r \ne 0} U_r u_{m+r} = \epsilon u_m$ (their Eq. 3), a 1D tight-binding (Anderson) model with $\epsilon = -U_0$. For irrational $\tau/4\pi$ the $T_m$ act as pseudorandom diagonal disorder with a Cauchy distribution (the Lloyd model), so the states are localized; for rational $\tau/4\pi$ they are periodic, giving band states and quantum resonance. (Izrailev's list misprints the page as 508.)
- F. L. Moore, J. C. Robinson, C. F. Bharucha, B. Sundaram, M. G. Raizen, "Atom optics realization of the quantum $\delta$-kicked rotor," Phys. Rev. Lett. 75, 4598 (1995), pp. 4598-4601, issue 25, published 18 December 1995, DOI 10.1103/PhysRevLett.75.4598. Read. Ultracold sodium atoms in a standing wave of near-resonant light pulsed with period $T$; scaled kick strength $\kappa = 11.6$ (10% uncertainty) and scaled Planck constant $\bar k = 8\omega_r T = 2.0$. The momentum spread grows diffusively at the classical rate $\kappa^2/2$ until a break time of about 8 kicks, then stays exponentially localized (checked out to 50 kicks). Their localization length $\xi = \kappa^2/4\bar k^2 \approx 8.3$ (units of $2\hbar k_L$) is "the $1/e$ point of the momentum distribution", that is, $P \propto e^{-|n|/\xi}$: half of $\xi$ in Wei's convention, which converts to $\kappa^2/(2\bar k^2) = D/\hbar^2$, consistent with note 5. Quantum resonances appear when $\bar k/2$ is a rational multiple of $2\pi$; for odd multiples of $\pi$ they see the sign flip between kicks in the first 2 to 3 kicks (our antiresonance). At resonance their distributions are narrow, not quadratically growing, because the initial state is a Gaussian packet, not a single plane wave.
- B. V. Chirikov, "A universal instability of many-dimensional oscillator systems," Phys. Rep. 52, 263 (1979), pp. 263 to 379, DOI 10.1016/0370-1573(79)90023-1 (pages and DOI confirmed 2026-11-03 from the INSPIRE-HEP record). Listed by Rechester and White (1980), Rechester, Rosenbluth and White (1981) and Izrailev (1990); Shepelyansky's 1986 list misprints the page as 623. Read: sections 4.1 to 4.2 (pp. 286 to 292) and 5 (pp. 310 to 336), from excerpts Wei made; the title is the PDF's running head. Facts for the report: the kicked pendulum $H = J^2/2 + k\sum_n\cos(\theta - n\tau)$ with $k = K/(2\pi)^2$ has resonances $I_r = 2\pi n$ of half-width $2\sqrt K$ (Eqs. 4.18, 4.19), and the overlap criterion gives $K_T = \pi^2/4 \approx 2.47$, about 2.5 times the numerical border $K \approx 1$ (p. 291); a footnote on p. 321 cites Greene's 0.971635. The main fixed point is unstable for $K > 4$ and the period-2 orbit for $|K| > 2$ (p. 320). KS entropy $h \approx \ln(K/2)$ for $K \gg 1$ (Eq. 5.41, from averaging $\ln|K\cos\theta|$); his Table 5.1 (p. 322) gives $h/\ln(K/2) = 1.027$ at $K = 6.21$ and 1.20 at $K = 4$ (our $0.98/0.916 = 1.07$ at $K = 5$ fits between). Random phases give $D = \langle(\Delta I)^2\rangle/t = K^2/2$ (Eq. 5.46, our normalization); the measured rate varies with $K$ with period about $2\pi$, with large values near $K \approx 2\pi n$ (Fig. 5.6). Accelerator modes: see note 7.
- J. M. Greene, "A method for determining a stochastic transition," J. Math. Phys. 20, 1183 (1979), DOI 10.1063/1.524170, for the critical $K_c \approx 0.9716$. Listed by Rechester and White (1980) and Izrailev (1990, as "Green, G."); title and DOI from an index record. Wei's PDF is the earlier report version (Princeton Plasma Physics Laboratory, November 1978, same title), which was read for the map and the critical value: Greene's map is $r_{n+1} = r_n - \frac{k}{2\pi}\sin 2\pi\theta_n$, $\theta_{n+1} = \theta_n + r_{n+1}$ on the unit torus, which is ours with $p = 2\pi r$ and $K = k$ (the sign of the kick is a shift of $\theta$ by $\pi$), and it gives $k_c = 0.971635\ldots$. Cite the journal version; its end page is still unconfirmed (a web search on 2026-11-03 found the DOI again, in a Bohrium record, but no page range), so the report cites J. Math. Phys. 20, 1183 (1979). The value $K_{cr} = 0.971635$ also appears in Shepelyansky (1986), Fig. 5 caption.
- A. B. Rechester and R. B. White, "Calculation of turbulent diffusion for the Chirikov-Taylor model," Phys. Rev. Lett. 44, 1586 (1980), pp. 1586 to 1589, DOI 10.1103/PhysRevLett.44.1586 (confirmed 2026-11-03 from the Scilit and MIT DSpace records). Read; see note 3.
- A. B. Rechester, M. N. Rosenbluth, R. B. White, "Fourier-space paths applied to the calculation of diffusion for the Chirikov-Taylor model," Phys. Rev. A 23, 2664 (1981), pp. 2664 to 2672, DOI 10.1103/PhysRevA.23.2664 (confirmed 2026-11-03 from the APS and Scilit records). Read; source of the $J_2^2$ term in note 3.
- D. L. Shepelyansky, "Localization of quasienergy eigenfunctions in action space," Phys. Rev. Lett. 56, 677 (1986), pp. 677 to 680, DOI 10.1103/PhysRevLett.56.677 (confirmed 2026-11-03 from the APS record and the author's publication list). Read; see note 5.
- F. M. Izrailev, "Simple models of quantum chaos: spectrum and eigenfunctions," Phys. Rep. 196, 299 (1990), pp. 299 to 392, DOI 10.1016/0370-1573(90)90067-C (confirmed 2026-11-03 from the Scilit record). Read pp. 299 to 301 and 304 to 322, and the reference list (pp. 388 to 392); see notes 1, 3, 5 and 6.
- F. M. Izrailev and D. L. Shepelyansky, "Quantum resonance for a rotator in a nonlinear periodic field," Theor. Math. Phys. 43, 553 (1980), DOI 10.1007/BF01029131 [Teor. Mat. Fiz. 43, 417 (1980)]. Listed by Izrailev (1990); title and DOI from the Math-Net.Ru record. (Fishman, Grempel and Prange's list swaps the two page numbers.)
- F. Saif, "Classical and quantum chaos in atom optics," Phys. Rep. 419, 207 (2005), DOI 10.1016/j.physrep.2005.07.002. Wei has the PDF; only the title page and contents were read. A possible secondary source on the cold-atom experiments (its contents list sections on the kicked rotor model and on dynamical localization of atoms).

## Open questions for Wei

- Papers still wanted: none. The journal version of Greene (1979) would only add its end page.
- Course policy (Wei decides): if the course asks students to declare tools they used, `report/report.tex` has a commented placeholder for an acknowledgment just before the bibliography.
- Answered 2026-11-02: about one month left; report about 20+ pages (full theory, derivations such as the Bessel corrections and the exact two-kick results, appendices); talk about 20+ minutes (about 20 slides); TeX installed in the cloud session; the Beamer deck embeds the animation frames with the `animate` package; experiment 07 is built before the write-up, with a stop for Wei to check it.
