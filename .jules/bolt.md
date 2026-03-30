## 2026-03-30 - Optimized Analytic ODE Solution in Stan
**Learning:** Algebraically simplifying the analytic solution for cell density reduced `exp()` calls from ~7 per evaluation to a maximum of 3. Combined with loop hoisting for clone-level parameters, this yielded a 1.38x speedup in sampling.
**Action:** Always check if complex analytic solutions in Stan `functions` blocks can be further simplified or if redundant transcendental functions can be lifted out of loops.

## 2026-03-30 - Stan 2.33+ Syntax Compatibility
**Learning:** Modern CmdStan (2.33.0+) strictly enforces the `array[]` keyword syntax. Old syntax (`int design[C]`) causes compilation failures.
**Action:** Use `stanc --print-canonical` to automatically upgrade legacy Stan files when moving to newer environments.
