## 2026-04-20 - Loop Hoisting and Analytic Simplification in Stan
**Learning:** Hoisting exponentiations of clone-level parameters out of observation loops ((N)$ to (C)$) provides a significant performance boost in Stan models. Additionally, algebraic simplification of piecewise ODE solutions and the use of `expm1()` improves both speed and numerical stability.
**Action:** Always check for redundant exponentiations or indexing inside observation loops in Stan models and pre-calculate them in the `transformed parameters` block.
