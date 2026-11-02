# CLAUDE.md

Quantum kicked rotor: Wei Chen's final project for PHY-315 (Quantum Mechanics). Wei is a physics and math double major. Read `PLAN.md` before starting any step. It has his full brief, the decisions made so far, and notes from a planning chat that were checked numerically.

## How Wei wants to work

- Build in this order, one step at a time (the step map is in `PLAN.md`):
  1. classical code
  2. quantum code
  3. validation tests
  4. physics experiments
  5. write-up: a LaTeX report, plus a Beamer deck built from it
- After each step, run everything yourself and look at every figure you made. Then commit, push, and stop. Tell Wei the branch name, the exact commands to run the step on his Windows computer, and what he should see (key numbers, and what each figure looks like). Wait for his go-ahead before starting the next step.
- If a check fails, debug it. Never loosen a tolerance to make a test pass. Every tolerance needs a stated reason: machine precision, the ensemble's statistical error, or a known physical correction.
- Explain things in plain, conversational language. Use LaTeX notation for math. Never use em dashes (in chat, code comments, docs, or LaTeX). Be concise, direct and technically precise.
- Don't invent citations. Check every reference against a real source, and say plainly when one can't be confirmed.
- Don't add experiments or features beyond `PLAN.md` without asking.

## Scope rules

- Wei's first idea (split-step Fourier for wave packets and tunneling) was rejected because it's covered in class. Don't drift into wave-packet propagation or tunneling. The FFT kick/free step here is exact, not a splitting approximation. Present it briefly as a tool and keep the focus on chaos, Floquet theory and dynamical localization.
- The code may later be extended for a numerical analysis project. Keep it modular and parameter-driven (plain functions, no hidden globals) so that analysis can be added later, but don't build numerical-analysis-only parts unless Wei asks.

## Physics conventions

Scaled units, kick period 1, effective Planck constant written $\hbar$.

- $H = p^2/2 + K\cos\theta\,\sum_n \delta(t-n)$.
- Classical standard map: $p_{n+1} = p_n + K\sin\theta_n$, $\theta_{n+1} = \theta_n + p_{n+1}$. Wrap $\theta$ mod $2\pi$. Never wrap $p$ in the dynamics (wrap it only to draw phase portraits).
- Quantum: $p = \hbar m$ with integer $m$ from `np.fft.fftfreq(N, d=1/N)`. One period: multiply $\psi(\theta)$ by $e^{-iK\cos\theta/\hbar}$, FFT, multiply $\psi_m$ by $e^{-i\hbar m^2/2}$, inverse FFT. The inverse step undoes these in reverse order.
- Initial states: quantum $m = 0$; classical ensemble with $p = 0$ and $\theta$ uniform on $[0, 2\pi)$.
- Defaults: $K = 5$, $\hbar$ between 0.25 and 1, a few hundred to a thousand kicks. Let validation check 4 choose $N$; $2^{12}$ is probably too small at $\hbar = 0.25$ (see `PLAN.md`).
- Avoid quantum resonances ($\hbar/4\pi$ rational) in localization runs; $\hbar = 4\pi$ is the resonance test. Keep $K$ out of the accelerator-mode windows $2\pi l \le K \le \sqrt{(2\pi l)^2 + 16}$.

## Code layout and commands

```
kicked_rotor/   classical.py, quantum.py, observables.py, plotting.py
experiments/    01_classical_phase_portraits.py ... 07_bonus_xi_scaling.py
tests/          pytest validation checks
figures/        every experiment saves its figures here
report/         LaTeX report
slides/         Beamer deck
```

- Dependencies are in `requirements.txt` (NumPy, SciPy, Matplotlib, pytest). In cloud sessions, the SessionStart hook in `.claude/settings.json` runs `scripts/install_pkgs.sh` to install them.
- From the repo root, `python -m pytest` runs the tests and `python experiments/01_classical_phase_portraits.py` runs an experiment. Both must also work on Wei's Windows computer from the repo root.
- Every experiment saves its figures to `figures/` and prints a short summary. Scripts never block on a plot window; they open one only when `--show` is passed. Save figures as PNG and PDF. Make animations as GIFs with Matplotlib's `PillowWriter`, so no ffmpeg is needed. Fix random seeds.
- Commit the final figures so the LaTeX can include them. Keep animations small.
- TeX isn't preinstalled in the cloud image. At the write-up step, install it with apt (for example `texlive-latex-recommended texlive-latex-extra texlive-science latexmk`), or ask Wei whether he'd rather compile locally.
- Wei runs the code on Windows from a clone outside Box (Box syncing and `.git` don't mix). If he hasn't said how he runs Python (plain terminal, Anaconda Prompt, or VS Code), ask before writing run instructions.

## Status

- [x] Planning, 2026-11-01 (claude.ai chat): brief, decisions and verified notes in `PLAN.md`
- [x] Step 1: classical code (demo: experiment 01), 2026-11-01
- [x] Step 2: quantum code (demo: `python -m kicked_rotor.quantum`), 2026-11-02
- [x] Step 3: validation tests (figures: experiments 05 and 06), 2026-11-02
- [ ] Step 4: physics experiments 02, 03, 04 and animations (07 only if time allows)
- [ ] Step 5: LaTeX report and Beamer deck

Update this list at the end of each step.
