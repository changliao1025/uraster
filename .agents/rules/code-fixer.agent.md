# Rule: Numerical & Code Debugging

Guidelines for diagnosing and resolving numerical, mathematical, and script runtime issues in pymath.

## Python Execution & Verification
- **Dynamic Interpreter Resolution:** Always resolve the active interpreter from `.vscode/settings.json` (`python.defaultInterpreterPath`).
- **Direct Script Verification:** Validate fixes by executing the specific script directly:
  `<python_interpreter> <path_to_script.py>`
  Check console outputs, matrix dimensions, convergence metrics, and generated figures.
- **Headless Plotting:** For scripts that render figures, ensure non-blocking / headless execution where appropriate (e.g., `plt.savefig(...)` instead of blocking `plt.show()` in unattended subshells).

## Numerical Bug Fixing Constraints
- **Numerical Failure Awareness:** Watch for common numerical edge cases:
  - Matrix shape and broadcasting mismatches.
  - Singular or ill-conditioned matrices in inversions or decompositions (use pseudo-inverse or add small regularization $\epsilon$ if appropriate).
  - Division by zero, numerical underflow/overflow, and `NaN` / `Inf` propagation.
  - Convergence failures in iterative solvers, MCMC chains, or optimization routines.
- **Minimal Blast Radius:** Limit modifications strictly to the specific script or function being repaired (maximum 2 files). Do not reformat unrelated demos.
- **Fail-Fast Guard (Max 2 Attempts):** If a numerical bug or script failure cannot be resolved after 2 attempts, STOP immediately. Report the exact traceback, array shapes, and values, and ask the user for mathematical clarification.
