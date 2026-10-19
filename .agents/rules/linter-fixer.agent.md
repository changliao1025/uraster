# Rule: Code Quality & Scientific Formatting

Guidelines for code style, documentation, and linting across mathematical and physical scripts.

## Scope & Target Bounds
- **Targeted Scope Only:** Inspect and format ONLY the specific files requested by the user or touched during the current task. Never run repository-wide lint or format sweeps.
- **Preserve Academic & Tutorial Context:** Preserve all comments explaining mathematical derivations, formulas, and source references (e.g., links to papers, Colab notebooks, or tutorial articles).

## Mathematical Naming & Style Conventions
- **Standard Math Notation Allowed:** Do NOT rename conventional mathematical symbols to conform to generic style rules:
  - Uppercase letters for matrices and datasets (e.g., `X`, `A`, `B`, `W`).
  - Standard short vector/tensor names (e.g., `y`, `b`, `u, s, vh`).
  - Greek letter transliterations (e.g., `alpha`, `beta`, `mu`, `sigma`, `theta`).
- **Array Dimension Documentation:** In function docstrings or comments, document the expected tensor/array shapes (e.g., `X: (N, D) array`, `y: (N,) array`).

## Tool Execution & Fail-Fast Guard
- **Execute via Resolved Python:**
  `<python_interpreter> -m ruff check <target_file>`
- **Missing Tool Guard:** If a linter tool is not installed in the `pymath` environment, report that it is missing and complete your turn immediately. Do NOT attempt automated installation or retry in a loop.
- **Single Pass Only:** Never iterate continuously trying to achieve "zero warnings". Apply clean fixes in one pass, report any remaining warnings, and stop.