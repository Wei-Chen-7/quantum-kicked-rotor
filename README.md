# Quantum kicked rotor

Classical and quantum simulations of the kicked rotor, a particle on a ring that gets a sharp kick once per period. With strong kicks, the classical momentum spreads forever like a random walk. The quantum spread follows the classical growth for a while and then stops. This effect, dynamical localization, is the subject of this final project for an undergraduate quantum mechanics course.

Work in progress. The plan is in [PLAN.md](PLAN.md).

## How to run

Keep the clone outside Box (Box syncing and `.git` don't mix).

### One-time setup (Windows PowerShell)

Clone the repo, make a virtual environment in a `.venv` folder inside it, and install the packages from `requirements.txt` into that environment:

```powershell
git clone https://github.com/Wei-Chen-7/quantum-kicked-rotor.git
cd quantum-kicked-rotor
python -m venv .venv
.\.venv\Scripts\python -m pip install -r requirements.txt
```

Until the code is merged into `main`, also check out the work branch after cloning, for example `git checkout claude/wonderful-allen-wl1uau`.

### Running an experiment

From the folder of your clone:

```powershell
git pull
.\.venv\Scripts\python experiments\01_classical_phase_portraits.py
```

`.\.venv\Scripts\python` is the project's own Python, so there is nothing to activate. Each experiment prints a short summary and saves its figures to `figures/` as PNG and PDF. Add `--show` to also open the figures in a window.

A quick demo of the quantum code prints the norm and $\langle p^2\rangle$ for the first few kicks, next to the exact values after one and two kicks:

```powershell
.\.venv\Scripts\python -m kicked_rotor.quantum
```

Options such as `--hbar 0.25 --kicks 20` change the parameters (`--help` lists them).

### Small differences between computers

The quantum numbers come out the same on every computer, to all printed digits. The classical ones differ slightly. Two computers can round the sine of a number differently in the last bit, and chaos multiplies such a difference by about $e^{\lambda}$ per kick ($\lambda \approx 1$ at $K = 5$). After about 37 kicks every orbit is different, which is the effect experiment 06 shows. Ensemble averages then agree only within their printed error bars.

So running an experiment rewrites the figure files with slightly different versions. The report uses the committed ones; to get them back (do this before `git pull`):

```powershell
git restore figures
```

### Validation tests

```powershell
.\.venv\Scripts\python -m pytest
```

This runs the five validation checks from `PLAN.md` (norm, resonance and antiresonance, early times, grid size, time reversal): 39 tests in about 10 seconds. Every tolerance in `tests/` states its reason.

In VS Code, open the clone's folder and choose the `.venv` interpreter (Command Palette, **Python: Select Interpreter**).

On macOS or Linux, use `python3 -m venv .venv` and `.venv/bin/python` in place of `.\.venv\Scripts\python`.

## Experiments

| Script | What it shows |
|---|---|
| `experiments/01_classical_phase_portraits.py` | Phase portraits of the classical standard map for $K = 0.5, 0.97, 2, 5$, with chaotic and regular orbits told apart by their Lyapunov exponent |
| `experiments/02_classical_diffusion.py` | Classical diffusion: $\langle p^2\rangle$ grows like $Dn$; $D$ against $K$ from 4.5 to 30 compared with the Bessel-function formulas, with the accelerator-mode windows marked; the lag correlations behind the remaining misses (about 30 s) |
| `experiments/03_quantum_vs_classical.py` | The main result: the quantum $\langle p^2\rangle$ follows the classical growth, then stops, for $\hbar = 0.25, 0.5, 1$; break times, and the curves rescaled with $\xi = D/\hbar^2$ |
| `experiments/04_momentum_distribution.py` | The localized momentum distribution on a log scale against the classical Gaussian, with fitted localization lengths compared with $D/\hbar^2$. Also saves the animation `figures/04_momentum_distribution.gif` (about 50 s) |
| `experiments/05_resonance_test.py` | Validation check 2 as a figure: quadratic energy growth at the quantum resonance $\hbar = 4\pi$, with the Bessel-function momentum distribution, and the exact return every two kicks at the antiresonance $\hbar = 2\pi$ |
| `experiments/06_time_reversal.py` | Validation check 5 as a figure: run forward, then backward; the classical return is lost after about 37 kicks, while the quantum state comes back to roundoff level |
| `experiments/07_bonus_xi_scaling.py` | Bonus: the localization length against $\hbar$ (41 values from 0.2 to 1 at $K = 5$) and against $K$ (4.5 to 11, medians over five values of $\hbar$), compared with $D/\hbar^2$ (about 3 minutes) |
