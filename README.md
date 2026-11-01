# Quantum kicked rotor

Classical and quantum simulations of the kicked rotor, a particle on a ring that gets a sharp kick once per period. With strong kicks, the classical momentum spreads forever like a random walk. The quantum spread follows the classical growth for a while and then stops. This effect, dynamical localization, is the subject of this final project for an undergraduate quantum mechanics course.

Work in progress. The plan is in [PLAN.md](PLAN.md).

## How to run

Keep the clone outside Box (Box syncing and `.git` don't mix). Run every command from the repo root.

### One-time setup (Windows, Anaconda Prompt)

Make a separate conda environment for the project. The first command takes Python from the free conda-forge channel only; pip then installs the packages listed in `requirements.txt`. Replace `C:\path\to\quantum-kicked-rotor` with the folder of your clone.

```bat
conda create -n kicked-rotor --override-channels -c conda-forge python=3.13 pip
conda activate kicked-rotor
cd /d C:\path\to\quantum-kicked-rotor
python -m pip install -r requirements.txt
```

### Each session

```bat
conda activate kicked-rotor
cd /d C:\path\to\quantum-kicked-rotor
python experiments\01_classical_phase_portraits.py
```

Each experiment prints a short summary and saves its figures to `figures/` as PNG and PDF. Add `--show` to also open the figures in a window.

The same environment works in VS Code: open the folder, run **Python: Select Interpreter** from the Command Palette, pick `kicked-rotor`, and use the same `python` commands in a new terminal.

On other systems, `python -m pip install -r requirements.txt` followed by `python experiments/01_classical_phase_portraits.py` does the same.

## Experiments

| Script | What it shows |
|---|---|
| `experiments/01_classical_phase_portraits.py` | Phase portraits of the classical standard map for $K = 0.5, 0.97, 2, 5$, with chaotic and regular orbits told apart by their Lyapunov exponent |
