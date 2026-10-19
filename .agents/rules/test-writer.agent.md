# Rule: Mathematical & Numerical Verification

Guidelines for authoring and repairing verification routines and numerical tests in pymath.

## Test Execution & Environment
- **Dynamic Interpreter Resolution:** Resolve the active interpreter from `.vscode/settings.json` (`python.defaultInterpreterPath`).
- **Execution:**
  - For standalone test/verification scripts: `<python_interpreter> <path_to_script>`
  - For pytest files (if configured): `<python_interpreter> -m pytest <path_to_test_file> -v`
- **Targeted Execution Only:** Run ONLY the specific test or script relevant to the task.

## Numerical & Scientific Testing Guidelines
- **Floating-Point Tolerances:** NEVER use exact equality (`==`) on floating-point results. Always use tolerance-based assertions:
  - `np.testing.assert_allclose(actual, expected, rtol=1e-5, atol=1e-8)`
  - `np.allclose(...)`
- **Seed Determinism:** Any test involving stochastic sampling (Monte Carlo, MCMC, random initializations, synthetic datasets) MUST use a fixed random seed.
- **Physical & Mathematical Invariants:** Prioritize testing mathematical properties and conservation laws over pure software mocks:
  - Matrix properties: Orthogonality ($U^T U = I$), symmetry ($A = A^T$), positive definiteness.
  - Probability properties: Probability distributions sum or integrate to $1$, variances are non-negative.
  - Physical quantities: Energy conservation, monotonicity, boundary conditions.
- **Fail-Fast Guard (Max 2 Attempts):** If a mathematical assertion fails after 2 repair attempts, STOP immediately. Report expected vs. actual values, array shapes, and tolerance settings to the user.