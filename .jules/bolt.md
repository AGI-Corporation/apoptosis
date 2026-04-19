## 2026-04-19 - Stan model optimization
**Learning:** Stan models with analytic solutions can be significantly optimized by algebraic simplification (reducing expensive `exp()` calls) and loop hoisting of exponentiated parameters. Modern array syntax is also preferred for compatibility with newer CmdStan.
**Action:** Apply algebraic simplification to `yt` function and hoist `exp()` calls out of observation loops in all Stan models.
