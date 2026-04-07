## 2026-04-07 - [Optimization of Stan ODE Models]
**Learning:** Algebraic simplification of analytic ODE solutions significantly reduces computational overhead. Specifically, reducing `exp()` calls from 6 to 3 in the main logic path and hoisting clone-level parameter exponentiation out of the observation loops provided a ~1.6x speedup in the `m1.stan` model benchmark. Additionally, adding a small floor (e.g., `1e-9`) to return values ensures numerical stability when these values are passed to `log()` functions.

**Action:** Always look for algebraic simplifications in analytic solutions and lift expensive, redundant operations like `exp()` out of loops when the parameters are constant across multiple iterations. Use `fmin` and `fmax` for vectorized stability checks.
