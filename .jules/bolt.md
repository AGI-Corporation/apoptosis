## 2026-04-02 - Stan Analytic Solution and Vectorization
**Learning:** Significant performance gains in Stan can be achieved by algebraically simplifying analytic solutions to minimize 'exp()' calls and lifting expensive exponentiated parameter transformations out of observation loops.
**Action:** Always check if complex analytic solutions can be simplified and pre-calculate any clone-level or replicate-level parameter transformations in 'transformed parameters' instead of doing them repeatedly in loops.
