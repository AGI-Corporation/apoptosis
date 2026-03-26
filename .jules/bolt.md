## 2025-05-15 - Algebraic Simplification of ODE Analytic Solutions
**Learning:** Consolidating multiple analytic sub-functions (Rt, Qat, Qct) into a single optimized function (yt) reduced the number of expensive `exp()` calls from 6 to at most 3 per observation. In a Python benchmark, this resulted in a 2.4x speedup (2.0s -> 0.8s for 1M calls).
**Action:** Always look for algebraic simplifications in Stan functions that reduce the number of calls to transcendental functions (`exp`, `log`, `pow`).

## 2025-05-15 - Vectorization of Parameters in Transformed Parameters Block
**Learning:** Lifting `exp()` operations on clone-level parameters out of the observation loop (N observations) into a vector operation on clone parameters (C clones) reduces the number of operations by -C$. In Stan, this also reduces the size of the auto-differentiation graph.
**Action:** Pre-calculate exponentiated or transformed parameters outside of loops whenever they only depend on a lower-dimension index (like clone or design).

## 2025-05-15 - Stan 2.33.0+ Array Syntax Migration
**Learning:** Modern Stan (2.33.0+) has deprecated the old array syntax (`int design[C]`) in favor of `array[C] int design`. Using the old syntax results in compilation errors.
**Action:** Use `stanc --print-canonical` to automatically migrate old Stan files to modern syntax when working with updated CmdStan versions.
