## 2026-04-01 - Optimizing Analytic Cell Density Solution in Stan

**Learning:** Reducing expensive transcendental function calls like `exp()` and hoisting parameter transformations out of observation loops provide significant performance gains in Stan models. Algebraic simplification of piecewise analytic solutions can reduce `exp()` calls from 6 per observation down to 1 (for $t < t_d$) or 3 (for $t \ge t_d$).

**Action:** Always check for redundant `exp()` or `log()` calls inside observation loops. Pre-calculate parameter vectors at the clone or design level ($O(C)$ or $O(D)$) before entering the observation loop ($O(N)$). Vectorize error calculations and intermediate variables where possible using Stan's built-in functions like `fmin()` and `log()`.

## 2026-04-01 - Modern Stan Array Syntax

**Learning:** CmdStan 2.33.0+ removed the legacy `type name[]` array syntax. Using `array[] type name` is now required and prevents compilation errors.

**Action:** Update all Stan files to the modern `array` syntax when working with recent CmdStan versions.
