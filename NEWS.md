# varGuidTS 0.1.9

* Revised the introductory vignette to combine the MathML-rendered model equations from version 0.1.8 with the simpler basic workflow, coefficient-summary workflow, and prediction examples from version 0.1.6.
* Kept `inst/doc/index.html` out of the source tarball while retaining `vignettes/intro.Rmd`, `inst/doc/intro.Rmd`, `inst/doc/intro.R`, and `inst/doc/intro.html`.

# varGuidTS 0.1.8

* Fixed the prebuilt vignette HTML so equations are rendered as HTML MathML instead of appearing as raw \[...\] text.
* Removed `inst/doc/index.html` from the source tarball as requested.
* Kept `vignettes/intro.Rmd`, `inst/doc/intro.Rmd`, `inst/doc/intro.R`, and `inst/doc/intro.html`.

# varGuidTS 0.1.7

* Expanded the vignette using the project README/model-specification material.
* Clarified model orders, data layout, shared versus subject-specific coefficients, and subject-specific thresholds.
* Clarified that `df_t` in `predict()` is the Student-t degrees of freedom, not a data frame.
* Computed Student-t exceedance scores using a standardized Student-t innovation distribution.
* Kept `inst/doc/index.html` as the prebuilt vignette index for CRAN-style source packages.

# varGuidTS 0.1.6

* Restored the introductory vignette for CRAN submission.
* Added prebuilt HTML vignette files under `inst/doc/` so the source package includes both the vignette source and rendered HTML.
* Kept the explicit CRAN-readable Maintainer field in `DESCRIPTION`.

# varGuidTS 0.1.5

* Updated the DESCRIPTION file for CRAN submission with explicit Author and Maintainer fields.
* Removed the package vignette directory and VignetteBuilder metadata to avoid prebuilt-vignette warnings for the CRAN source submission.
* Kept the public API focused on `lmvt()`, `predict.lmvt()`, `summary.lmvt()`, `coef.lmvt()`, and `simulate_scenario()`.

# varGuidTS 0.1.4

* Updated DESCRIPTION metadata for CRAN maintainer parsing using an explicit Maintainer field.
* Renamed the package from `vgrisk` to `varGuidTS`.
* Updated package metadata, documentation, README, examples, and tests
  to use the GitHub URL `https://github.com/zionwzz/variance-guided-risk-demo`.
* Updated authorship metadata: Zihao Wang is an author; Min Lu is author and
  maintainer.

# vgrisk 0.1.2

* Added `summary.lmvt()` and `print.summary.lmvt()` for model diagnostics,
  selected/nonzero coefficient counts, and a full coefficient summary table.
* Expanded `predict.lmvt()` documentation with examples for pooled thresholds,
  subject-specific thresholds, and Gaussian/Student-t exceedance probabilities.
* Fixed subject-aware lag lookup in the forward variance recursion used by
  `lmvt()` and `predict.lmvt()` so that ARCH/GARCH lags are resolved correctly
  for all subjects, not only the first subject.
