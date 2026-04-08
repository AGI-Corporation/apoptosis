## 2026-04-08 - Hoisting clone-level operations in Stan

**Learning:** In hierarchical Stan models where parameters are group-level (e.g., per clone) but used in observation-level loops, performing expensive operations like `exp()` or redundant indexing inside the $O(N)$ loop is a significant bottleneck. Hoisting these to $O(C)$ in the `transformed parameters` block, where $C$ is the number of groups, drastically reduces the computational load and the size of the autodiff stack.

**Action:** Always look for group-level parameter transformations and hoist them out of observation-level loops. Use `fmin`/`fmax` to avoid branching when possible.
