## 2026-04-16 - Stan Model Optimization
**Learning:** Hoisting redundant exponentiations and simplifying analytic ODE solutions in Stan results in a significant performance boost (~2x speedup in this case, from 9.52s to 4.65s). Vectorizing parameters where possible and minimizing `exp()` calls in tight observation loops are critical patterns for Stan performance.
**Action:** Always check for redundant indexing or function calls (like `exp()`) inside `for` loops in Stan `transformed parameters` and `model` blocks.
