## 2026-05-10 - Algebraic Simplification and Loop Hoisting in Stan
**Learning:** Combining algebraic simplification of analytic ODE solutions with loop hoisting (pre-calculating clone-level parameters and ODE coefficients) provides a massive performance boost, measured at ~1.8x speedup on minimal benchmarks. Reducing the number of transcendental function calls (`exp()`) and using `expm1()` for numerical stability are key.
**Action:** Always check if complex analytic solutions can be simplified and ensure all invariant parameters are pre-calculated outside the innermost observation loops in Stan models.

## 2026-05-10 - CmdStan 2.33+ Syntax and Configuration
**Learning:** CmdStan 2.33.0+ removes support for legacy array declarations (`type var[size];`). Also, `cmdstanpy>=1.2.0` changed some utility functions (e.g., `jsondump` is gone, use `write_stan_json`) and requires `LD_LIBRARY_PATH` to include the TBB library for executing compiled binaries.
**Action:** Use `array[size] type var;` syntax and ensure `cmdstanpy` and environment variables are correctly configured for newer CmdStan versions.
