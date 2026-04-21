## 2026-04-21 - Stan Model Optimization
**Learning:** Hoisting `exp()` calls from observation-level loops to clone-level vectors significantly reduces the computational load on the Stan auto-diff stack. Combining this with algebraic simplification of analytic ODE solutions (e.g., using `expm1()` and reducing total exponentiations) provides a measurable speedup.
**Action:** Always check for redundant exponentiations or parameter transformations inside loops in Stan models.

**Learning:** Modern `cmdstanpy` (>= 1.2.0) and Stan (2.33+) require explicit `array[]` syntax and have removed `jsondump`.
**Action:** Use `write_stan_json` and update array syntax when working with newer Stan versions.
