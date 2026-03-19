## 2025-05-15 - Algebraic Simplification of ODE Analytic Solutions
**Learning:** Algebraically simplifying the analytic solution of the cell density ODE reduced the number of `exp()` calls from 6 to 3 per observation. Combined with pre-calculating clone-level parameters and vectorizing the error calculation, this significantly reduces the computational load on the Stan auto-diff stack.
**Action:** Always check if complex ODE solutions can be simplified to minimize transcendental function calls and lift parameter-specific transformations out of observation loops.

## 2025-05-15 - Stan Variable Scope and Reusability
**Learning:** Variables declared inside a nested local block `{}` in the `transformed parameters` section are not accessible in `generated quantities`. Re-calculating them in `generated quantities` is redundant and slows down the model.
**Action:** Declare performance-critical intermediate variables at the top level of the `transformed parameters` block to ensure they are available for reuse in `generated quantities`.
