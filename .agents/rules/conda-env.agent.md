# Rule: Python Environment & Conda Resolution

Guidelines for resolving and verifying Python environments in the pymath workspace.

## Environment & Interpreter Resolution
- **Primary Resolution:** Read `.vscode/settings.json` for `"python.defaultInterpreterPath"`.
- **Target Conda Environment:** This workspace uses the Conda environment `pymath` (typically located at `/Users/changliao/miniforge3/envs/pymath/bin/python`).
- **Invoke via Interpreter:** Always invoke Python scripts and CLI tools directly through the resolved interpreter:
  `<python_interpreter> <path_to_script.py>`
  `<python_interpreter> -m <tool> ...`
- **Never Run Bare `python` or `pip`:** Subshells may default to the system Python (`/usr/bin/python3`) lacking scientific libraries (NumPy, SciPy, Matplotlib, PyTorch).

## Constraints
- **Token Guard:** Never execute verbose discovery commands such as `conda config --show` or `conda info` (they waste tens of thousands of tokens).
- **Environment Stability:** Never create, destroy, or run unbounded updates on the `pymath` Conda environment without explicit user instructions.
