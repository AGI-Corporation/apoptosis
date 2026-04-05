## 2026-04-05 - Stan Model Optimization
**Learning:** Algebraic simplification of analytic ODE solutions combined with hoisting exponentiated parameters from observation loops to clone-level vectors significantly reduces computational overhead (measured ~1.6x-1.9x speedup). Using a small floor (1e-9) on densities improves numerical stability without affecting physical realism.
**Action:** Always check for opportunities to pre-calculate shared parameters and simplify math to minimize expensive `exp()` or `log()` calls in Stan loops.
