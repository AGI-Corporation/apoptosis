## 2026-05-08 - Stan Model Optimization - Loop Hoisting & Analytic Simplification
**Learning:** In Stan models, transcendental functions like `exp()` are computationally expensive. Moving these calls out of the main observation loop (hoisting) into `transformed parameters` significantly reduces total execution time, especially when parameters are constant for a group (e.g., clone-level effects). Additionally, simplifying complex analytic solutions to ODEs can yield substantial performance gains and improve numerical stability (e.g., using `expm1`).
**Action:** Always check for redundant calculations in Stan's observation loops. Hoist parameter transformations to the highest possible scope. Use `expm1` for stability in growth models.

## 2026-05-08 - CmdStan 2.33+ Syntax Migration
**Learning:** Modern CmdStan (2.33+) requires the new `array` syntax (e.g., `array[N] int y;`). Using `stanc --print-canonical` to migrate can be problematic as it expands `#include` directives, breaking modularity.
**Action:** Use regex-based migration for legacy Stan arrays to preserve `#include` directives and maintain model structure.
