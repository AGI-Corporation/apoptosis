## 2026-04-24 - Stan ODE Analytic Solution & Loop Hoisting
**Learning:** Significant performance gains (~1.8x to 4.5x) can be achieved in Stan by (1) simplifying complex analytic ODE solutions to reduce transcendental function calls and (2) hoisting invariant exponentiations (`exp()`) out of large observation loops into `transformed parameters` at the clone/replicate level.
**Action:** Always check if expensive operations in Stan loops can be pre-calculated at a higher hierarchical level and simplify analytic expressions using `expm1()` where applicable.

## 2026-04-24 - CmdStanPy 1.3+ Compatibility
**Learning:** `jsondump` was removed from `cmdstanpy.utils` in recent versions; use `write_stan_json` instead. Also, `fixed_param=True` with `iter_warmup=0` requires `adapt_engaged=False` to pass validation.
**Action:** Use `write_stan_json` for data serialization and ensure `adapt_engaged=False` for fixed parameter runs.
