# Changelog

## [1.1.0] - 2026-10-19

Bootstrap rerun with the corrected simulator (179 replicates), independent likelihood, step-change model,
completeness check and jitter check (see the sections below). The main pipeline script is unchanged.

## [1.0.0] - 2026-10-10

Initial public release. It supersedes an internal analysis script for the same
case study. Verified by a complete run with seed 20250923 (ETAS 0.7.2,
changepoint 2.3, R 4.x); two runs gave identical results.

### Added

- **Step 4b, Monte Carlo calibration of the data-driven break.** An exact scan
  likelihood-ratio statistic for a change in the mean of exponential gaps, with the
  full search over split points repeated on each of `N_SIM` (default 9999) replicates
  drawn under the time-rescaling null; a Monte Carlo check of the maximum departure of
  N(tau) from the diagonal; and the null frequency of a PELT break. It gives a
  valid p-value for the selected break (post-selection is accounted for) but does
  not propagate parameter-estimation uncertainty (see README, Limitations).

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
- **Degrees of freedom** of the likelihood-ratio test count the 7 parameters of the Gaussian-kernel model.

- **Observed count in the post-break comparison.** The expected count
  (`tau[n] - tau[idx]`, 44.8) is a compensator increment over (t_idx, t_n]; the
  observed count over the same interval is `n - idx` = 27 (the earlier code
  subtracted one more, giving 26 and p = 0.0017). The corrected lower-tail Poisson
  probability is 0.0029 (descriptive).

- **Step 4c** (serial independence of the transformed gaps) and **Step 4d** (branching ratio
  with delta-method standard error), and the helper file `R/etas_helpers.R`.
- **Parametric bootstrap with refitting**: `R/etas_sim_gauss.R`, `R/etas_parametric_bootstrap.R`, `results/bootstrap_*.csv` (179 replicates).
- **`R/etas_sensitivity_Mc.R`**: refits and diagnostics for Mc = 2.5, 2.7, 2.9, 3.2.

### Corrected (bootstrap simulator)

- The simulator drew the number of background events with mean `integ0` instead of `mu * integ0`; fixed, and the bootstrap
  was rerun (179 replicates; the earlier 102-replicate results are superseded).
- Added `R/etas_ownlik.R` (independent likelihood, Hessian SE, profile likelihood), `R/etas_stepchange_bootstrap.R`
  (competing step-change model with bootstrap-calibrated LRT) and `R/etas_jitter_check.R`.

### Added (referee 4)

- `R/completeness_check.py` (magnitude distribution by period and around the break) and the source-level explanation of the package standard errors (README).

### Investigated

- **Stability of the Gaussian-kernel fit.** Three random starting points gave A, c, alpha, p, D and gamma agreeing to about four significant digits, and mu within about 2%.
- **Package version.** The original analysis cited ETAS 0.5.1. A source comparison
  of the C/C++ kernel and optimiser code between 0.5.1 and 0.7.2 found no
  functional difference, so the version does not explain the issues above.
  `mver = 2` was verified to be present in 0.7.2 and absent in 0.5.1; other
  versions were not tested.

### Corrected
- **Meaning of D.** In the Gaussian kernel of ETAS 0.7.2 (`ffun2`), `D*exp(gamma*m)` is a standard deviation, so D is a
  length in projected degrees (0.0088 deg, about 1 km), not a variance in deg^2. Verified by reproducing the package's own
  `lambd` values.

### Not implemented

- Bootstrap with more than 179 replicates, and any test of the spatial kernel itself.
