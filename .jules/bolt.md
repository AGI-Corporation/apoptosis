## 2026-05-09 - Stan Performance Optimization

**Learning:**
1. Algebraic simplification of analytic ODE solutions can yield significant speedups (approx. 1.8x) by reducing the number of transcendental function calls (`exp`, `expm1`).
2. Loop hoisting of clone-level parameters from (N)$ to (C)$ (where  \gg C$) further reduces overhead.
3. Upgrading to modern CmdStan (2.33+) requires migrating legacy array syntax (`int var[N]` to `array[N] int var`). Automated migration using regex is safer than `stanc --print-canonical` as it preserves modular `#include` directives.
4. Using `expm1(x)` instead of `exp(x) - 1` is more numerically stable for small `x`, which is critical in growth models.

**Action:**
Always check for redundant `exp()` or `log()` calls inside observation loops in Stan models. Prefer analytic solutions over numerical ODE solvers when available, and simplify them algebraically before implementation.
