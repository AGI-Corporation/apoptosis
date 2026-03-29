## 2024-05-24 - Stan Analytic Optimization
**Learning:** Algebraic simplification of analytic solutions (reducing `exp()` calls) combined with loop hoisting (pre-calculating clone-level parameters) provides a significant performance boost, measured at ~1.9x speedup (reducing `m1.stan` runtime from 5.32s to 2.79s for 100 iterations).
**Action:** Always look for opportunities to simplify complex mathematical expressions in Stan and lift expensive operations out of observation loops.

## 2024-05-24 - Stan Array Syntax
**Learning:** Modern Stan (2.33+) requires the `array[N] type name` syntax instead of the legacy `type name[N]` syntax. Using `stanc --print-canonical` is an efficient way to upgrade old models.
**Action:** Use the new array syntax for all Stan models to ensure compatibility with recent CmdStan versions.
