# varGuidTS

> Variance-guided time-series modeling for personalized temporal risk detection

`varGuidTS` fits balanced-panel autoregressive models with conditional
heteroscedasticity. The package implements a penalized ARX--GARCHX workflow
that simultaneously estimates:

- conditional mean dynamics with autoregressive and exogenous covariate effects;
- conditional variance dynamics with ARCH/GARCH and optional covariate effects;
- subject-specific baseline terms with shared population-level coefficients;
- exceedance-based risk scores defined as conditional threshold-exceedance probabilities.

The package is motivated by wearable sensor studies, but the same data structure
applies to other high-frequency or repeated-measures panel time series.

## Installation

```r
# From CRAN, once released
install.packages("varGuidTS")

# From GitHub
install.packages("remotes")
remotes::install_github("luminwin/varGuidTS")
```

Project URL: <https://github.com/zionwzz/variance-guided-risk-demo>

## Model orders

The main function `lmvt()` uses four order parameters:

| Argument | Meaning |
|---|---|
| `p` | autoregressive order in the mean, using lags of `y` |
| `q` | distributed-lag order for exogenous predictors `x` |
| `r` | ARCH order in the variance, using lags of squared residuals |
| `s_ord` | GARCH order in the variance, using lags of conditional variance |

When `use_x_in_variance = TRUE`, the same lag block of exogenous predictors is
used in both the mean and variance equations.

## Quick start

```r
library(varGuidTS)

set.seed(1)
sim <- simulate_scenario(scen = 2, S = 4, T = 60,
                         d_noise = 4, noise_kind = "ar1", seed = 1)

df <- sim[, c("s", "t", "y",
              grep("^X|^Noise|^Xbin", names(sim), value = TRUE))]

fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
            lambda_beta = 0.05, lambda_gamma = 0.05,
            use_x_in_variance = TRUE, maxit = 8)

fit
summary(fit, include_zero = FALSE)
```

## Exceedance-based risk scores

After fitting the model, `predict()` returns the estimated conditional mean,
conditional standard deviation, and optional exceedance-based risk scores.
A scalar threshold gives a pooled cutoff, while a named vector supplies
subject-specific cutoffs.

```r
# Pooled threshold
pred_global <- predict(fit, df, threshold = quantile(df$y, 0.90))
head(pred_global)

# Subject-specific thresholds
subject_cutoffs <- tapply(df$y, df$s, quantile, probs = 0.75, na.rm = TRUE)

pred_subject <- predict(fit, df, threshold = subject_cutoffs,
                        innov_g = TRUE, innov_t = TRUE, df_t = 8)
head(pred_subject)
```

Here `df_t = 8` means 8 degrees of freedom for the standardized Student-t
innovation distribution used to compute `risk_t`; it is unrelated to the data
frame named `df` in the example.

## Main functions

| Function | Purpose |
|---|---|
| `lmvt()` | Fit penalized panel ARX--GARCHX models |
| `predict()` | Estimate conditional mean, variance, and exceedance-based risk scores |
| `summary()` | Return model diagnostics and a coefficient summary table |
| `coef()` | Extract fitted coefficients |
| `simulate_scenario()` | Generate reproducible simulation scenarios |

## Interpretation

The exceedance-based risk score is a model-based temporal feature. It summarizes
how likely a future or held-out observation is to exceed a chosen threshold,
conditional on the estimated mean, variance, covariates, and lagged history. It
should not be interpreted as a clinically validated diagnostic probability unless
validated in a prospective clinical study.

## Reference

Wang, Z. and Lu, M. (2026). *Variance-Aware Penalized Panel Models for Temporal Risk Detection from Wearable Sensor Data* (submitted).

## Maintainer

Min Lu <luminwin@gmail.com>.

## License

MIT.
