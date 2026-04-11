## 2026-04-11 - Optimizing ODE Analytic Solutions in Stan
**Learning:** Significant performance gains (up to 2x speedup) can be achieved in Stan models by:
1. Simplifying complex analytic ODE solutions to minimize transcendental function calls (`exp()`).
2. Hoisting exponentiation of clone-level parameters out of observation loops.
3. Vectorizing error calculations using `fmin()` and `log()` instead of scalar loops with conditional logic.
**Action:** Always check for repeated `exp()` calls in observation loops and identify if parameters can be pre-calculated at a higher level (e.g., clone-level instead of observation-level).
