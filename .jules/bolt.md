## 2026-04-27 - Optimized Apoptosis models (m1, m2, null)

**Learning:** Significant performance gains in Stan can be achieved by combining algebraic simplification of analytic solutions (reducing `exp()` calls) with loop hoisting (pre-calculating clone-level parameters outside the observation loop). Modern array syntax (`array[N] int var;`) is required for CmdStan 2.33.1+.

**Action:** Always check if expensive transcendental functions (`exp`, `log`, etc.) or redundant lookups can be lifted out of observation-level loops. Use modern array syntax and `expm1()` for numerical stability.
