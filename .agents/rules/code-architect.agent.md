# Rule: Math & Physics Testbed Architecture & Script Organization

Guidelines for structuring mathematical routines, physical simulations, and exploratory demo scripts.

## Scope & Token Efficiency
- **Targeted Scope:** Focus analysis strictly on the specific mathematical topic, script, or question requested. Do not perform repo-wide sweeps.
- **Context Guard:** Never read more than 5 full files into context. Use targeted symbol navigation or grep.

## Topic-Based Organization
Organize code cleanly by mathematical and physical domains:
1. **Domain Folders:** Keep scripts categorized by topic (e.g., `linear_algebra/`, `statistics/`, `hilbert_space/`, `deep_learning/`, `numeric/`, `geometry/`).
2. **Self-Contained Demos:** Demo scripts should be runnable as standalone scripts (`python path/to/script.py`) with runnable entry points or standard `if __name__ == "__main__":` blocks.
3. **Data & Artifact Separation:**
   - Raw datasets belong in `pymath/data/`.
   - Output plots, figures, and charts must be saved to `figures/<topic>/`.
   - Never hardcode absolute or brittle working-directory paths (e.g., `os.getcwd() + '/pymath/data/...'`). Always resolve paths relative to the file using `pathlib.Path(__file__).resolve().parent`.

## Mathematical Reproducibility
- **Explicit Random Seeds:** Any script using pseudo-random generation (NumPy, SciPy, PyTorch, random) must set an explicit seed (e.g., `rng = np.random.default_rng(42)`) to ensure deterministic figures and calculations.
- **Invariant Preservation:** When refactoring formulas, preserve mathematical equivalence and document mathematical equations in comments.
