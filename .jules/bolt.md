## 2025-05-14 - Stan Model Optimization: Analytic Simplification and Loop Hoisting

**Learning:** Combining algebraic simplification of analytic ODE solutions with loop hoisting and vectorization in Stan models leads to significant performance gains.
1. Reducing `exp()` calls from 6 to 3 in the analytic solution for cell density (`yt`) directly impacts the inner loop.
2. Hoisting exponentiation of clone-level parameters (`kq`, `td`, `kd`) out of the observation loop ((N)$) into a clone-level vector ((C)$) avoids redundant expensive calculations, as  \gg C$.
3. Vectorizing error calculations (`x_small`) using `fmin()` reduces the size of the expression graph for Stan's auto-differentiation engine.

**Action:** Always check for analytical simplifications in ODE solutions and identify clone-level or design-level parameters that can be pre-calculated outside of observation loops in Stan's `transformed parameters` or `generated quantities` blocks.
