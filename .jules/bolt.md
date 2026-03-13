## 2025-05-14 - Redundant Transcendental Operations in Stan Loops
**Learning:** In Stan models, transcendental functions like `exp()` are computationally expensive, especially when used inside observation loops ($N$). In this codebase, parameters like `kq`, `td`, and `kd` are constant per clone ($C=27$), but were being re-exponentiated for every observation ($N=895$) in every leapfrog step. Additionally, the analytical ODE solution in `functions.stan` had a highly redundant form (7 `exp` calls) that could be mathematically simplified.

**Action:** Move per-clone parameter transformations outside of observation loops in Stan files. Simplify analytical formulas to minimize the number of `exp()` calls, reducing the size of the autodiff graph and improving sampling speed.
