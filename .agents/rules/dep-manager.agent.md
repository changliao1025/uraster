# Rule: Scientific Dependencies & Environment Management

Guidelines for managing dependencies in the pymath mathematical and physical testbed.

## Environment & Tool Execution
- **Resolved Interpreter:** Always execute package tools via the resolved Python interpreter from `.vscode/settings.json`:
  `<python_interpreter> -m pip install <package>`
- **No Blind Editable Installs:** Do NOT run `pip install -e .` unless a formal `pyproject.toml` or `setup.py` exists in the repository root.

## Scientific Stack Guidelines
- **Core Stack Awareness:** `pymath` relies on standard scientific and data analysis packages (e.g., `numpy`, `scipy`, `matplotlib`, `pandas`, `scikit-learn`, `pymc`, `torch`, `chaospy`).
- **Targeted Installations:** When a script requires an additional mathematical or scientific library, verify if it is already present before installing.
- **Fail-Fast Guard:** If a dependency installation fails or produces dependency conflicts with the `pymath` Conda environment, STOP immediately. Output the exact package manager log rather than retrying in an automated loop.
