## 2026-03-28 - Stan Model Optimization: Algebraic Simplification and Loop Hoisting
**Learning:** In Bayesian models with analytical solutions, reducing the number of expensive transcendental function calls (like `exp()`) and hoisting redundant transformations out of observation loops (`N \approx 900`) to parameter-level loops (`C \approx 26`) yields significant performance gains.
**Action:** Always check Stan models for `exp()`, `log()`, or `pow()` calls inside the `model` or `transformed parameters` blocks that depend only on index-mapped parameters; pre-calculate these at the higher level.
