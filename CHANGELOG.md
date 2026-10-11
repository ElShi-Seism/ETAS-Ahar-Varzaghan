# Changelog

## [1.0.0] - 2026-10-10

Initial public release. It supersedes an internal analysis script for the same
case study. Verified by a complete run with seed 20250923 (ETAS 0.7.2,
changepoint 2.3, R 4.x); two runs gave identical results.

### Fixed

- **Uncontrolled randomness.** `ETAS::catalog()` with its default
  `roundoff = TRUE` adds uniform noise to epicentre coordinates before the
  spatial window is applied, so the number of target events varied between
  unseeded runs (249-252). `set.seed(20250923)` is now called before every
  `catalog()` call, and `ROUNDOFF` can be switched off for sensitivity checks.
- **Wrong residual field.** The earlier script used `resid.etas(fit)$tres`, a
  binned diagnostic residual on a fixed grid. The pipeline uses `$tau` (per-event
  transformed time, non-decreasing, one value per target event) and `$U` for the
  Kolmogorov-Smirnov test.
- **Wrong breakpoint-to-date conversion.** The earlier script used
  `(index / total_events) * total_days`, which assumes evenly spaced events. The
  pipeline reads the recorded elapsed times of the target events.
- **Segment counts ignored the spatial window.** Events outside the study window
  were counted, which inflated the post-break count (44 instead of 27). Counts
  now apply the window and the magnitude threshold.
- **Sparse segment.** A full 8-parameter ETAS fit on a segment with few events is
  not attempted (threshold `MIN_EVENTS_FOR_ETAS = 50`); the observed count after
  the break is compared with the count expected under the stationary model.
- **Hard-coded study length** (4490 days) is now computed from the study dates.
- **Degrees of freedom** of the likelihood-ratio test exclude `q`, which is not
  estimated under `mver = 2`.

### Investigated

- **Weak identifiability of D and q under `mver = 1`.** Different runs gave very
  different (D, q) pairs, and this run gave q = 8.8 and D = 0.0012. A large q with
  small D is the power-law kernel approaching its Gaussian limit, where the
  likelihood is nearly flat. `mver = 2` (Gaussian kernel, no q) gave stable
  estimates from three random starting points: A, c, alpha, p, D and gamma agreed
  to about four significant digits, and mu to within about 2%.
- **Package version.** The original analysis cited ETAS 0.5.1. A source comparison
  of the C/C++ kernel and optimiser code between 0.5.1 and 0.7.2 found no
  functional difference, so the version does not explain the issues above.
  `mver = 2` was verified to be present in 0.7.2 and absent in 0.5.1; other
  versions were not tested.

### Not implemented

- Profile likelihood or bootstrap for `q`.
- A formal test of non-stationarity that accounts for the breakpoint having been
  chosen from the same data.
