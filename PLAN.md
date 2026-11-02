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

How to use this in check 3: test each code against its own exact formula (quantum to machine precision, classical to its statistical error). Beyond $n = 2$, read the next paragraph before comparing quantum with classical. Don't pick a percentage by eye.

**Step 2 finding (2026-11-02): after $n = 2$ the gap grows very fast.** Gap (quantum minus classical, relative) at $K = 5$, with the classical values from the midpoint grid (converged to all digits: $\langle p^2\rangle_3 = 32.25286$, $\langle p^2\rangle_4 = 40.07047$). The quantum value at $n = 3$, $\hbar = 0.05$ was confirmed to 13 digits by an independent Bessel-matrix evolution, so this is physics, not a bug:

| $\hbar$ | 0.5 | 0.25 | 0.1 | 0.05 | 0.02 | 0.01 | 0.005 |
|---|---|---|---|---|---|---|---|
| gap at $n = 3$ | +0.5% | -6.4% | -10.8% | -6.3% | -1.4% | -0.35% | -0.09% |
| gap at $n = 4$ | +49% | +29% | +6.6% | +15% | +14% | +11% | -1.4% |

The gap does vanish like $\hbar^2$ once $\hbar$ is small enough, but the $\hbar^2$ coefficient grows very fast with $n$: about $-0.2$ at $n = 2$, $-35$ at $n = 3$, and roughly $-550$ at $n = 4$. So at $\hbar = 0.05$ the $n = 3$ gap is already 6%, and the earlier suggestion to compare at $\hbar = 0.05$ fails. Izrailev (1990), pp. 313 to 314, gives the reason: full quantum-classical correspondence lasts only up to the Ehrenfest time $t_E \sim \ln(1/\hbar)$ (his Eq. 2.2.20), a few kicks, and for $t_E \lesssim t \lesssim t_D$ "only rough statistical properties like diffusion remain." Options for check 3 beyond $n = 2$ (decide with Wei): compare at very small $\hbar$ (about 0.005 for $n = 3$) with a tolerance from the measured $\hbar^2$ scaling, or compare the average growth rate over many kicks instead of kick-by-kick values.

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

The $J_2^2$ term matters most near $K = 10$, where the 1980 formula misses by 24%. Each truncated formula is an asymptotic expansion. At moderate $K$ the higher correlations aren't small (at $K = 5$, $C_4 \approx 0.021$ and $C_5 \approx -0.035$), so it misses by a few percent. The 7% miss at $K = 22$ isn't understood yet; look into it in experiment 02. With Wei's initial condition $p_0 = 0$, the first kicks are correlated differently (see the $n = 2$ formula), so fit $D$ from the slope at late times.

### 4. Grid size and edge probability

At $K = 5$, $D \approx 0.99\,K^2/2 \approx 12.4$, so $D/\hbar^2 \approx 198$ at $\hbar = 0.25$ and about 50 at $\hbar = 0.5$. If $\xi \approx D/\hbar^2$, then at $\hbar = 0.25$ the tail at the grid edge $|m| = N/2$ is about $e^{-2\cdot2048/198} \approx 10^{-9}$ of the peak for $N = 2^{12}$, which fails a $10^{-10}$ edge test. With $N = 2^{13}$ it's about $10^{-18}$. A thousand kicks at $N = 2^{14}$ take under a second, so be generous with $N$ and let check 4 confirm it. Define "edge probability" precisely in the code and the report (for example, the total probability with $|m| > 3N/8$). Check its maximum over time, not just its final value, because the instantaneous distribution fluctuates around its time average.

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

### 6. Time reversal (predictions; the quantum part checked in Step 2)

- Classical: roundoff ($\epsilon \sim 10^{-16}$ to $10^{-15}$) grows like $e^{\lambda n}$, with $\lambda \approx \ln(K/2) \approx 0.92$ at $K = 5$. The return should fail once the forward run is longer than about $\ln(1/\epsilon)/\lambda \approx 35$ to 40 kicks. (Step 1 measured $\lambda = 0.98$ at $K = 5$, the median finite-time exponent of the chaotic orbits in experiment 01, which narrows this to about 35 to 38 kicks.) Izrailev (1990), p. 315, citing Shepelyansky, Physica D 8, 208 (1983) and Chirikov, Izrailev and Shepelyansky, Physica D 33, 77 (1988), reports that at $K = 5$ with roundoff around $10^{-12}$ the classical map retraces its path for only about 10 kicks, shorter than $\ln(10^{12})/0.98 \approx 28$; find out why in experiment 06.
- Quantum: the evolution is linear and unitary, so roundoff never grows exponentially, and the return overlap should stay close to machine precision. Step 2 check: $K = 5$, $\hbar = 0.25$, $N = 2^{14}$, 1000 kicks forward and 1000 back returned to the initial state within $3\times10^{-15}$ in every component. The norm drifts steadily upward by about $3\times10^{-16}$ per kick ($2.9\times10^{-13}$ after 1000 kicks), so the fidelity came out $1 + 6\times10^{-13}$; compare states directly or normalize first.
- A good point for the report: the classical phase-space density (Liouville evolution) is also linear and reversible. The irreversibility comes from following individual chaotic trajectories, whose errors grow exponentially.

### 7. Accelerator modes

Stable accelerator modes exist for $2\pi l \le K \le \sqrt{(2\pi l)^2 + 16}$: about $[6.28, 7.45]$ for $l = 1$, $[12.57, 13.19]$ for $l = 2$, and $[18.85, 19.27]$ for $l = 3$. Stay out of these windows, with some margin, in any scan over $K$. $K = 5$ is well clear.

### 8. Keeping clear of the rejected topic

The kick, FFT, free phase, inverse FFT sequence looks like split-step Fourier, but here it's exact. In the write-up, cover it in one short paragraph as a tool, and put the weight on chaos, Floquet theory and localization.

## References

Status as of 2026-11-02. "Read" means the paper itself was read (Wei supplied the PDFs of the papers marked "read"; they are not in the repo). "Listed" means only reference lists or a bibliographic record were checked, which confirms where the paper is but not what it says. Verify the rest before citing.

- G. Casati, B. V. Chirikov, F. M. Izrailev, J. Ford, "Stochastic behavior of a quantum pendulum under a periodic perturbation," in *Stochastic Behavior in Classical and Quantum Hamiltonian Systems*, edited by G. Casati and J. Ford, Lecture Notes in Physics 93 (Springer, Berlin, 1979), pp. 334-352, DOI 10.1007/BFb0021757. Read: title page, abstract and p. 349. The chapter's own byline is "G. Casati; B. V. Chirikov and F. M. Izraelev; Joseph Ford" (so this order is right; the lists in Shepelyansky 1986 and Izrailev 1990 swap Ford and Izrailev). Abstract: the first numerical study of the quantum kicked rotor, finding diffusive-like linear energy growth and a Gaussian momentum distribution as in the classical motion, plus "several surprising distinctions." Page 349: at $\tau = 4\pi$ the free step is the identity, so $\psi(t) = \psi(0)e^{itk\cos\theta}$ and $\langle E\rangle \propto t^2$ (their Eqs. 37 and 38; also seen at $\tau = 4\pi/m$ for $m = 3, 4, 8, 32$); at $\tau = 2\pi$ the solution is strictly periodic with period 2 ("anti-resonant"); and, after the Como conference, the first report of a break time $t_B$ after which diffusion slows, with $t_B \propto k$ at fixed $K$ (run with $k = 40$, $\tau = 1/8$).
- S. Fishman, D. R. Grempel, R. E. Prange, "Chaos, quantum recurrences, and Anderson localization," Phys. Rev. Lett. 49, 509 (1982), pp. 509-512, DOI 10.1103/PhysRevLett.49.509. Read. The mapping: with $U(\theta) = -\tan[kV(\theta)/2]$ and $T_m = \tan(E_m/2)$, $E_m = \omega - \tau m^2/2$, quasienergy states obey $T_m u_m + \sum_{r \ne 0} U_r u_{m+r} = \epsilon u_m$ (their Eq. 3), a 1D tight-binding (Anderson) model with $\epsilon = -U_0$. For irrational $\tau/4\pi$ the $T_m$ act as pseudorandom diagonal disorder with a Cauchy distribution (the Lloyd model), so the states are localized; for rational $\tau/4\pi$ they are periodic, giving band states and quantum resonance. (Izrailev's list misprints the page as 508.)
- F. L. Moore, J. C. Robinson, C. F. Bharucha, B. Sundaram, M. G. Raizen, "Atom optics realization of the quantum $\delta$-kicked rotor," Phys. Rev. Lett. 75, 4598 (1995), pp. 4598-4601, issue 25, published 18 December 1995, DOI 10.1103/PhysRevLett.75.4598. Read. Ultracold sodium atoms in a standing wave of near-resonant light pulsed with period $T$; scaled kick strength $\kappa = 11.6$ (10% uncertainty) and scaled Planck constant $\bar k = 8\omega_r T = 2.0$. The momentum spread grows diffusively at the classical rate $\kappa^2/2$ until a break time of about 8 kicks, then stays exponentially localized (checked out to 50 kicks). Their localization length $\xi = \kappa^2/4\bar k^2 \approx 8.3$ (units of $2\hbar k_L$) is "the $1/e$ point of the momentum distribution", that is, $P \propto e^{-|n|/\xi}$: half of $\xi$ in Wei's convention, which converts to $\kappa^2/(2\bar k^2) = D/\hbar^2$, consistent with note 5. Quantum resonances appear when $\bar k/2$ is a rational multiple of $2\pi$; for odd multiples of $\pi$ they see the sign flip between kicks in the first 2 to 3 kicks (our antiresonance). At resonance their distributions are narrow, not quadratically growing, because the initial state is a Gaussian packet, not a single plane wave.
- B. V. Chirikov, Phys. Rep. 52, 263 (1979). Listed by Rechester and White (1980), Rechester, Rosenbluth and White (1981) and Izrailev (1990). (Shepelyansky's 1986 list misprints the page as 623.)
- J. M. Greene, "A method for determining a stochastic transition," J. Math. Phys. 20, 1183 (1979), DOI 10.1063/1.524170, for the critical $K_c \approx 0.9716$. Listed by Rechester and White (1980) and Izrailev (1990, as "Green, G."); title and DOI from an index record. Wei's PDF is the earlier report version (Princeton Plasma Physics Laboratory, November 1978, same title), which was read for the map and the critical value: Greene's map is $r_{n+1} = r_n - \frac{k}{2\pi}\sin 2\pi\theta_n$, $\theta_{n+1} = \theta_n + r_{n+1}$ on the unit torus, which is ours with $p = 2\pi r$ and $K = k$ (the sign of the kick is a shift of $\theta$ by $\pi$), and it gives $k_c = 0.971635\ldots$. Cite the journal version; its page numbers weren't checked. The value $K_{cr} = 0.971635$ also appears in Shepelyansky (1986), Fig. 5 caption.
- A. B. Rechester and R. B. White, "Calculation of turbulent diffusion for the Chirikov-Taylor model," Phys. Rev. Lett. 44, 1586 (1980). Read; see note 3.
- A. B. Rechester, M. N. Rosenbluth, R. B. White, "Fourier-space paths applied to the calculation of diffusion for the Chirikov-Taylor model," Phys. Rev. A 23, 2664 (1981). Read; source of the $J_2^2$ term in note 3.
- D. L. Shepelyansky, "Localization of quasienergy eigenfunctions in action space," Phys. Rev. Lett. 56, 677 (1986). Read; see note 5.
- F. M. Izrailev, "Simple models of quantum chaos: spectrum and eigenfunctions," Phys. Rep. 196, 299 (1990). Read pp. 299 to 301 and 304 to 322, and the reference list (pp. 388 to 392); see notes 1, 3, 5 and 6.
- F. M. Izrailev and D. L. Shepelyansky, "Quantum resonance for a rotator in a nonlinear periodic field," Theor. Math. Phys. 43, 553 (1980), DOI 10.1007/BF01029131 [Teor. Mat. Fiz. 43, 417 (1980)]. Listed by Izrailev (1990); title and DOI from the Math-Net.Ru record. (Fishman, Grempel and Prange's list swaps the two page numbers.)
- F. Saif, "Classical and quantum chaos in atom optics," Phys. Rep. 419, 207 (2005), DOI 10.1016/j.physrep.2005.07.002. Wei has the PDF; only the title page and contents were read. A possible secondary source on the cold-atom experiments (its contents list sections on the kicked rotor model and on dynamical localization of atoms).

## Open questions for Wei

- Should the $\hbar = 2\pi$ antiresonance be added to check 2?
- Check 3 beyond $n = 2$: compare at very small $\hbar$, or compare average growth rates (note 1)?
- Papers still wanted for the write-up: Chirikov (1979), and the journal version of Greene (1979) if a page range is needed.
- Report length, slide count or talk length, and the due date (not given yet).
