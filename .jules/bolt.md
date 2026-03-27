## 2025-05-14 - Stan Optimization: Algebraic Simplification and Loop Lifting
**Learning:** Significant performance gains in Stan models can be achieved by algebraically simplifying analytic solutions to reduce the number of transcendental function calls (like `exp()`). Furthermore, lifting redundant exponentiations of group-level parameters out of observation-level loops (where $N \gg C$) saves a significant number of computations per HMC leapfrog step.
**Action:** Always check if group-level parameters can be pre-calculated before entering observation loops, and look for mathematical simplifications in analytic ODE solutions.

## 2025-05-14 - Stan 2.33+ Array Syntax and Include Behavior
**Learning:** Stan 2.33.0+ removed the old array syntax (e.g., `int arr[N]`) in favor of `array[N] int arr`. Using `stanc --print-canonical` automatically upgrades the syntax but also inlines any `#include` directives, which can lead to code duplication if not careful.
**Action:** When upgrading Stan files with `stanc --print-canonical`, be aware that includes will be expanded. If maintainability via includes is desired, manually restore the `#include` directive after the syntax upgrade.
