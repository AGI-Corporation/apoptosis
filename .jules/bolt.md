## 2026-03-31 - Stan Loop Optimization & Analytic Simplification
**Learning:** Combining loop hoisting (pre-calculating clone-level parameters) with algebraic simplification of analytic solutions (reducing `exp()` calls) provides a significant performance boost in Stan models. Vectorizing the error calculation further improves efficiency by avoiding scalar loops and conditional branching.
**Action:** Always look for opportunities to lift expensive operations out of observation loops and simplify mathematical expressions in Stan functions. Use `fmin` and `log` for vectorizing piecewise linear or logarithmic functions.

## 2026-03-31 - Stan Scoping in Transformed Parameters
**Learning:** Variables declared at the top level of the `transformed parameters` block are automatically in scope for `generated quantities`. Re-declaring them in `generated quantities` causes a parser error.
**Action:** Check if variables needed in `generated quantities` are already available from `transformed parameters`.

## 2026-03-31 - Stan Modern Array Syntax
**Learning:** Modern CmdStan (2.33+) requires the `array[]` keyword for array declarations. Old syntax (e.g., `int design[C]`) causes compilation errors.
**Action:** Use `stanc --print-canonical` to upgrade Stan files to modern syntax. Be careful as it expands `#include` directives; manually restore them if modularity is desired.
