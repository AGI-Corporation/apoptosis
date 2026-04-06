## 2026-04-06 - Optimized ODE analytic solution and loop hoisting in Stan models
**Learning:** Algebraic simplification of analytic solutions in Stan (reducing redundant `exp()` calls) and hoisting parameter transformations out of observation loops provide significant performance gains by reducing the computational load on the automatic differentiation stack.
**Action:** Always check for redundant transcendental function calls in loops and use piece-wise simplification for complex analytic solutions.
