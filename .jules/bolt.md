## 2025-05-22 - [Optimizing ODE solutions in Stan]
**Learning:** Algebraic simplification of analytic ODE solutions can significantly reduce computational cost. Reducing `exp()` calls from 6 to 1-3 per observation, combined with hoisting exponentiation out of $O(N)$ loops, yielded a ~78% performance improvement in the core cell density calculation.
**Action:** Always look for mathematical simplifications and hoist parameter transformations (like `exp()`) out of observation loops in Stan models.

## 2025-05-22 - [Binary Artifacts in PRs]
**Learning:** Tools like CmdStan generate large compiled binaries and header files. Including these in a PR violates repository hygiene and causes review failures.
**Action:** Explicitly clean up all build artifacts, benchmark scripts, and header files (`.hpp`) before submission.
