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

In VS Code, open the clone's folder and choose the `.venv` interpreter (Command Palette, **Python: Select Interpreter**).

On macOS or Linux, use `python3 -m venv .venv` and `.venv/bin/python` in place of `.\.venv\Scripts\python`.

## Experiments

| Script | What it shows |
|---|---|
| `experiments/01_classical_phase_portraits.py` | Phase portraits of the classical standard map for $K = 0.5, 0.97, 2, 5$, with chaotic and regular orbits told apart by their Lyapunov exponent |
