## 2025-05-14 - Initializing Bolt Journal
**Learning:** Found several performance opportunities in Stan models and functions. Specifically, the analytic solution in `functions.stan` is not simplified, leading to redundant `exp()` calls (6 vs 3 for  \ge t_d$).
**Action:** I will focus on optimizing the `yt` function in `functions.stan` by simplifying it algebraically to reduce the number of `exp()` calls.
