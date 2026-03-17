## 2026-03-17 - Algebraic Simplification of ODE Solutions
**Learning:** Algebraic simplification of analytic ODE solutions can significantly reduce the number of expensive transcendental functions (like `exp()`) called per observation. In this case, reducing `exp()` calls from 6 to 3 for $t \ge td$.
**Action:** Always check if multiple `exp()` calls can be combined using exponent rules (e.g., $e^a \cdot e^b = e^{a+b}$) before implementing in Stan.

## 2026-03-17 - Lifting Operations out of Observation Loops
**Learning:** Calling `exp()` on clone-level parameters inside an observation loop is redundant and expensive. Stan's autodiff engine has to track these operations for every observation.
**Action:** Pre-calculate `exp()` of parameters at the hierarchy level where they are defined (e.g., clone level) before entering observation-level loops.
