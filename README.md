# ETAS Analysis Pipeline — Ahar-Varzaghan 2012 Doublet Sequence

Authors: Davood Kazemi Lafmejani (Islamic Azad University, Tehran North Branch) and Elham Shabani (Institute of Geophysics, University of Tehran)

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
   (`mver = 2`), because the default power-law kernel (`mver = 1`) is weakly
   identified for this catalog (see [Methodological notes](#methodological-notes-and-pitfalls-addressed)).
4. Runs a multi-start stability check to confirm convergence.
5. Computes per-event transformed-time residuals (`tau`, not the binned
   `tres` field) and a Kolmogorov-Smirnov goodness-of-fit test.
6. Detects structural breakpoints in the residual sequence with the PELT
   algorithm (`changepoint` package) and maps them back to calendar dates
   using the catalog's own recorded elapsed times.
7. Compares the periods before and after the breakpoint: a likelihood-ratio test
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
| Events after the breakpoint | 27 observed vs 44.8 expected under the stationary model (lower-tail Poisson p = 0.0017, descriptive) |

Two complete runs with the same seed gave identical event counts and log-likelihood.
For comparison, the default power-law kernel (`mver = 1`) gave D = 0.0012 and q = 8.8
in this run: a very large q, close to the kernel's Gaussian limit, which is weakly identified.

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
- **Spatial-kernel weak identifiability** under `mver = 1` (see above) is
  documented, and `mver = 2` is the reported model.

## Limitations

- A profile-likelihood analysis of `q` is **not implemented**.
- The Poisson comparison for the sparse segment is descriptive: the breakpoint
  was chosen from the same data, so the p-value is not a formal test.
- A formal piecewise-vs-stationary likelihood-ratio test is skipped when the
  second segment has fewer than 50 events.
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

A complete run takes about 6 minutes.

The script expects a tab-separated catalog file named `ascii_correct.txt`
in the working directory (rename your catalog file accordingly). See [`data/README.md`](data/README.md) for the
expected format and data provenance.

## Repository structure

```
.
├── README.md
├── .gitignore
├── LICENSE            (GPL-3)
├── CITATION.cff
├── CHANGELOG.md
├── install.R
├── R/
│   └── etas_ahar_varzaghan_pipeline.R
└── data/
    └── README.md
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
