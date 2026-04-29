## 2025-05-14 - Performance Bottlenecks in Stan Models
**Learning:** In Stan, redundant calculations of `exp()` inside large observation loops (e.g.,  \approx 1000$) can be significantly optimized by hoisting them to clone-level vectors (e.g.,  \approx 10$) and by algebraically simplifying the analytic ODE solutions.
**Action:** Always look for parameter transformations that can be moved outside observation loops and simplify mathematical expressions to minimize transcendental function calls.
