# ETAS Analysis Pipeline — Ahar-Varzaghan 2012 Doublet Sequence

Authors: Davood Kazemi Lafmejani (Islamic Azad University North Tehran Branch) and Elham Shabani (Institute of Geophysics, University of Tehran)

An analysis pipeline built on top of the [`ETAS`](https://cran.r-project.org/package=ETAS) R package
(Jalilian, 2019) to model the space-time seismicity of the 11 August 2012
Ahar-Varzaghan (NW Iran) earthquake doublet sequence, detect temporal
non-stationarity via changepoint analysis, and diagnose spatial-kernel
parameter identifiability.

**This repository is not a modification of the `ETAS` package.** It is an
independent analysis workflow that calls the package's public API. All
credit for the underlying statistical model and its C++/R implementation
belongs to the original authors (see [Citation](#citation)).

## What this pipeline does

1. Loads a raw nine-column earthquake catalog (`long, lat, year, month, day,
   mag, sec, hour, minute`) and validates its time fields.
2. Builds an `ETAS::catalog` object over the Ahar-Varzaghan study region.
3. Fits the space-time ETAS model using the Gaussian spatial kernel
   (`mver = 2`), which has no spatial decay exponent.
4. Runs a multi-start stability check to confirm convergence.
5. Computes per-event transformed-time residuals (`tau`, not the binned
   `tres` field) and a Kolmogorov-Smirnov goodness-of-fit test.
6. Detects structural breakpoints in the residual sequence with the PELT
   algorithm (`changepoint` package) and maps them back to calendar dates
   using the catalog's own recorded elapsed times.
7. Calibrates the data-driven break by Monte Carlo (Step 4b): an exact scan
   likelihood-ratio statistic for a change in mean of the transformed gaps, with the
   whole search for the best split repeated on each of 9999 replicates simulated under
   the time-rescaling null (i.i.d. Exp(1) gaps), plus a Monte Carlo check of the
   maximum departure of the cumulative count from the diagonal.
8. Tests serial independence of the transformed gaps (Step 4c: lag-1 autocorrelation,
   Ljung-Box at lags 5/10/20 and a runs test, each calibrated by Monte Carlo) and
   computes the branching ratio with a delta-method standard error (Step 4d).
9. Compares the periods before and after the breakpoint: a likelihood-ratio test
   when both have enough events, otherwise a descriptive comparison of observed
   and expected counts.

## Verified results (seed 20250923, ETAS 0.7.2, changepoint 2.3)

| Quantity | Value |
|---|---|
| Target events (M >= 3.2, study window) | 251 |
| Model | space-time ETAS, Gaussian spatial kernel (`mver = 2`) |
| Log-likelihood | 945.56 |
| Main parameters | mu 1.24, A 0.333, c 0.0427, alpha 1.560, p 1.180, D 0.00875, gamma 0.600 |
| Stability (3 random starts) | A, c, alpha, p, D, gamma agree to about 4 significant digits; mu varies within ~2% |
| Residual K-S test (`res$U`) | D = 0.084, p = 0.058 (borderline) |
| PELT breakpoint | event index 224 = day 217.7 = 2013-03-16 (same for minseglen 10-25) |
| Monte Carlo calibration of the break (9999 null replicates, minseglen 20) | scan LR = 10.82 at event 224; p = 0.022 two-sided, 0.010 one-sided (deficit); same split and p between 0.009 and 0.027 for minseglen 10-25 |
| PELT false-positive rate under the null | 3.7% (2000 replicates) |
| Early excess (max abs departure of N(tau) from tau) | 22.1 events, null 95% quantile 35.3, Monte Carlo p = 0.30 (not significant) |
| Serial independence of the gaps (Monte Carlo) | lag-1 r = 0.068 (p = 0.27); Ljung-Box lags 5/10/20 p = 0.23/0.30/0.55; runs test p = 0.048 |
| Branching ratio (delta method) | 1.57 +/- 0.40 without magnitude cutoff; 1.18 +/- 0.17 (95% interval 0.85-1.51) with cutoff at 6.5 |
| Events after the breakpoint | 27 observed vs 44.8 expected under the stationary model (lower-tail Poisson p = 0.0029, descriptive) |

Two complete runs with the same seed gave identical event counts and log-likelihood.

## Parametric bootstrap (refitting)

`R/etas_sim_gauss.R` simulates the Gaussian-kernel ETAS model (conditioned on the doublet and on the observed
events outside the window; non-doublet magnitudes truncated at 5.4; its intensity formula reproduces the
`lambd` values computed by the package; the background count is Poisson with mean `mu * integ0`).
`R/etas_parametric_bootstrap.R` simulates, refits (`mver = 2`) and recomputes the diagnostics:

```bash
Rscript R/etas_parametric_bootstrap.R run <worker> <n_workers> <B> <outdir>   # needs <outdir>/reference_fit.rds
Rscript R/etas_parametric_bootstrap.R summarise <outdir>
```

Run each worker only once per output directory (a worker resumes from its CSV, but two copies of the same worker
write duplicate rows; remove duplicates by `b` before summarising). About 1 minute per replicate and core; 179 replicates are in `results/`.
Results (95% range of the refits): mu 1.06-1.46, A 0.22-0.39, c 0.027-0.088, alpha 1.45-1.74, p 1.15-1.26,
`D` (a standard deviation, not a variance) 0.0073-0.0106 deg, gamma 0.49-0.67. The asymptotic 95% intervals
cover the generating value in only 25% (alpha), 32% (p) and 61% (gamma) of replicates.
Bootstrap p-values (smallest possible 0.0056): scan statistic 0.0056, K-S 0.022, early excess 0.044,
runs test 0.050, Ljung-Box (lag 10) 0.31. The sensitivity to the magnitude cutoff (5.4) of the simulator was not examined.

## Competing model with a step change, independent likelihood, jitter check

- `R/etas_ownlik.R` (and `etas_ownlik_cpp.R`) is an independent implementation of the package's conditional
  log-likelihood; it reproduces the package value (945.49) exactly. `R/own_se.R` uses it for Hessian standard
  errors and profile likelihoods. The Hessian standard errors (alpha 0.099, p 0.038, gamma 0.068, D 0.0012,
  c 0.018) agree with the bootstrap, not with the `asd` values returned by the package, whose transformation
  appears to be inverted (too small for alpha, p, gamma; too large for D, c). Profile intervals: p 1.13-1.25,
  D 0.0044-0.0175.
- `R/etas_stepchange_bootstrap.R`: M1 = the ETAS intensity multiplied by a factor r after a change time chosen
  by scanning (8 parameters) against the stationary model M0. Observed LR = 27.3 (r = 0.24, break at day 219.4);
  bootstrap null (100 replicates, whole procedure repeated) has 95% quantile 8.6 and maximum 11.9, p = 0.0099.
  The AIC difference (25) is not informative: 81% of null replicates also have a lower AIC for M1.
  With all parameters free, M1 shifts A (0.33 to 0.82) and p (1.18 to 1.06), so the ETAS parameters and the step are
  not separately identified.
- `R/etas_jitter_check.R`: roundoff = FALSE and four seeds with roundoff = TRUE (`results/jitter_check.csv`):
  248-251 target events; alpha 1.52-1.61, p 1.19-1.21, D 0.0080-0.0106, gamma 0.53-0.65, c 0.047-0.061.

## Completeness over time (`R/completeness_check.py`, `results/completeness_check.txt`)

Maximum-curvature Mc stays at 2.5-2.7 in every period, and the share of events with 2.5 <= M < 3.2 does not drop after
the break (0.80 against 0.71 before; against days 30-217: 0.80 against 0.78, Fisher p = 0.61). The b-value of M >= 3.2 is
0.85 before and 0.94 after (28 events; Kolmogorov-Smirnov p = 0.22). A change of completeness at the break is therefore
not suggested by the magnitude distribution; this does not exclude changes in locations or in the network.

## Why the package standard errors are wrong

In `etasfit.R` the optimiser works on tht = sqrt(theta), `H` is the inverse Hessian with respect to tht, and the
package sets `avcov = 1/4 diag(1/tht) H diag(1/tht)`. The delta method needs `4 diag(tht) H diag(tht)`. The reported
standard error is therefore the correct one multiplied by 1/(4 theta): for D (0.00875) a factor 28.6 (reported 0.035,
Hessian 0.0012), for alpha (1.56) a factor 0.16 (0.0159 against 0.0995), for p (1.18) a factor 0.21 (0.0081 against 0.0380).

## Sensitivity to the completeness magnitude

`R/etas_sensitivity_Mc.R` refits the Gaussian-kernel model for Mc = 2.5, 2.7, 2.9 and 3.2 (same window,
period and seed) and repeats the diagnostics; results are in `results/sensitivity_Mc.csv`
(about 9 minutes in total for N_SIM = 9999). Summary (scan-test p-values are two-sided; 1e-4 is the
smallest value possible with 9999 replicates):

| Mc | events | alpha | p | gamma | break (day) | scan p | early-excess p | Ljung-Box lag 10 (MC) p | branching ratio, cutoff |
|---|---|---|---|---|---|---|---|---|---|
| 2.5 | 904 | 1.551 | 1.204 | 0.535 | 205.4 | 1e-4 | 0.030 | 1e-4 | 1.26 +/- 0.12 |
| 2.7 | 659 | 1.488 | 1.224 | 0.524 | 219.5 | 2e-4 | 0.095 | 3e-4 | 1.13 +/- 0.10 |
| 2.9 | 446 | 1.560 | 1.197 | 0.554 | 219.5 | 0.011 | 0.21 | 0.89 | 1.15 +/- 0.14 |
| 3.2 | 251 | 1.560 | 1.180 | 0.600 | 217.7 | 0.022 | 0.30 | 0.30 | 1.18 +/- 0.17 |

## Methodological notes and pitfalls addressed

The pipeline was refined from an earlier analysis script for the same case
study. The points below are pitfalls that affect any similar ETAS workflow
(details in [`CHANGELOG.md`](CHANGELOG.md)):

- **Randomness was uncontrolled.** `ETAS::catalog()` (default `roundoff = TRUE`)
  adds uniform noise to epicentre coordinates before the spatial window is
  applied, so target-event counts near the window edge change between
  unseeded runs (249-252 were observed). Optimiser starting points add further
  variability. The seed is set before every `catalog()` call.
- **Wrong residual field.** `resid.etas()$tres` is a binned diagnostic on an
  artificial grid; the per-event transformed time is `$tau`.
- **Wrong breakpoint-to-date conversion.** A proportional index-to-date formula
  assumes evenly spaced events; the pipeline uses the recorded elapsed times.
- **Segment counts ignored the spatial window**, and a full ETAS fit was
  attempted on a sparse segment. Segments are now counted inside the window,
  and a segment with fewer than 50 events is compared with the stationary
  model's expected count instead of being fitted.

## Limitations

- The asymptotic standard errors returned by `ETAS` are unreliable here (see the bootstrap); D is a standard deviation in projected degrees (about 1 km at Mc), of the order of catalog location error, so it is not interpreted as a fault dimension.
- The Poisson comparison for the sparse segment is descriptive: the breakpoint
  was chosen from the same data, so the p-value is not a formal test.
- A formal piecewise-vs-stationary likelihood-ratio test is skipped when the
  second segment has fewer than 50 events.
- The Monte Carlo calibration in Step 4b holds the 7 estimated parameters fixed. The parametric
  bootstrap (below) removes this restriction, but it uses 179 replicates (smallest p-value 0.0056),
  tests recoverability within the Gaussian model only, and does not test the kernel itself.
- Results depend on the catalog and on the random coordinate jitter; set
  `ROUNDOFF <- FALSE` in the script to check sensitivity to it.

## Requirements

- R >= 4.0
- R packages: `ETAS` (verified with 0.7.2; `mver = 2` is not available in 0.5.1), `changepoint`
- See [`install.R`](install.R) for exact installation steps.

## Usage

```bash
Rscript install.R
Rscript R/etas_ahar_varzaghan_pipeline.R
```

A complete run takes about 6 minutes; Step 4b takes a few seconds (set `N_SIM <- 999` for a quick test, `N_SIM <- 0` to skip it).

The script expects a tab-separated catalog file named `ascii_correct.txt`
in the working directory (rename your catalog file accordingly). See [`data/README.md`](data/README.md) for the
expected format and data provenance.

## Repository structure

```
.
├── README.md, CHANGELOG.md, CITATION.cff, LICENSE (GPL-3), install.R, .gitignore
├── R/
│   ├── etas_ahar_varzaghan_pipeline.R   main pipeline (Steps 0-5)
│   ├── etas_helpers.R                   Monte Carlo calibration, independence tests, branching ratio
│   ├── etas_sensitivity_Mc.R            refits for Mc = 2.5, 2.7, 2.9, 3.2
│   ├── etas_sim_gauss.R                 simulator of the Gaussian-kernel ETAS model
│   ├── etas_parametric_bootstrap.R      parametric bootstrap with refitting
│   ├── etas_ownlik.R, etas_ownlik_cpp.R independent likelihood (R and Rcpp)
│   ├── own_se.R                         Hessian standard errors and profile likelihoods
│   ├── etas_stepchange_bootstrap.R      step-change model and bootstrap-calibrated LRT
│   ├── etas_jitter_check.R              roundoff / seed sensitivity
│   └── completeness_check.py            magnitude distribution by period
├── results/                             CSV and text outputs quoted in the paper
└── data/README.md
```

## Citation

If you use this pipeline, please cite both this repository and the
underlying `ETAS` package:

> Jalilian, A. (2019). ETAS: An R Package for Fitting the Space-Time ETAS
> Model to Earthquake Data. *Journal of Statistical Software*, 88.
> https://doi.org/10.18637/jss.v088.c01

See [`CITATION.cff`](CITATION.cff) for the citation of this repository.

## License

GPL-3.0-or-later (see [`LICENSE`](LICENSE)). The `ETAS` package that this
code calls is licensed GPL (>= 2), which is compatible. The IRSC catalog is
not included and is subject to its own terms.

## Disclaimer

This pipeline is provided for research transparency and reproducibility.
It has not been peer-reviewed as a standalone software tool. Users should
independently verify results before relying on them for publication.

