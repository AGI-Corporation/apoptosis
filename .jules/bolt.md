## 2025-05-15 - Algebraic Simplification and Loop Hoisting in Stan
**Learning:** Re-calculating `exp(log_parameter)` inside observation loops ($O(N)$) is a significant bottleneck when $N$ is large. Pre-calculating these into clone-level vectors ($O(C)$) in the `transformed parameters` block provides a measurable speedup. Additionally, simplifying complex analytic solutions to reduce the number of transcendental function calls (`exp()`, `expm1()`) directly impacts sampling efficiency.
**Action:** Always check for redundant exponentiations or complex arithmetic that can be hoisted out of loops or simplified algebraically in Stan models.

## 2025-05-15 - Stan Modern Array Syntax Migration
**Learning:** Modern Stan (2.33+) has removed legacy array syntax (`int var[N];`). Using a regex-based migration tool is more modular than `stanc --print-canonical` because it preserves `#include` directives.
**Action:** Use the validated regex-based migration script for updating legacy Stan codebases while maintaining modularity.
