# ============================================================================
# Environment setup for the Ahar-Varzaghan ETAS pipeline
# ============================================================================
# The ETAS package version used here (0.7.2 verified; mver = 2 is absent in 0.5.1) is not
# always the latest version on CRAN. If install.packages("ETAS") does not
# provide mver = 2, install directly from the source repository instead:
#
#   https://github.com/cran/ETAS
#
# ----------------------------------------------------------------------------

required_pkgs <- c("ETAS", "changepoint", "Rcpp")  # Rcpp: only for R/etas_ownlik_cpp.R

installed <- rownames(installed.packages())
missing <- setdiff(required_pkgs, installed)

if (length(missing) > 0) {
  install.packages(missing, repos = "https://cloud.r-project.org")
}

etas_version <- as.character(packageVersion("ETAS"))
cat("ETAS package version installed:", etas_version, "\n")

if (utils::compareVersion(etas_version, "0.7.2") < 0) {
  warning(paste(
    "Installed ETAS version is older than 0.7.2 (the version this pipeline was verified with) and may not support",
    "mver = 2 (Gaussian spatial kernel), which this pipeline requires.",
    "Install from source instead:",
    "  # in a terminal:",
    "  git clone https://github.com/cran/ETAS.git",
    "  R CMD INSTALL ETAS"
  ))
}

cat("Setup complete.\n")
