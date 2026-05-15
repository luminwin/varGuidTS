\name{summary.lmvt}
\alias{summary.lmvt}
\title{Summarize a fitted \code{lmvt} object}
\description{
Produces model-level diagnostics and a coefficient summary table for a
fitted penalized panel ARX-GARCHX model. The table includes
subject-specific mean intercepts, autoregressive coefficients, mean
covariate coefficients, subject-specific variance baselines, ARCH/GARCH
coefficients, and variance covariate coefficients. Because \code{lmvt()}
fits a penalized quasi-likelihood model, the table reports estimates and
variable-selection indicators rather than classical standard errors or
p-values.
}
\usage{
\method{summary}{lmvt}(object, tol = 1e-8, include_zero = TRUE,
        sort_by = c("block", "abs_estimate", "none"), top = NULL, ...)
}
\arguments{
\item{object}{A fitted \code{lmvt} object.}
\item{tol}{Absolute value used to flag a coefficient as selected/nonzero.}
\item{include_zero}{Logical; if \code{FALSE}, return only coefficients with
\code{abs(estimate) > tol}.}
\item{sort_by}{Sorting rule for the coefficient table. Use \code{"block"}
for model order, \code{"abs_estimate"} for decreasing absolute estimate,
or \code{"none"} to keep construction order.}
\item{top}{Optional positive integer limiting the returned coefficient table
to the first \code{top} rows after filtering and sorting.}
\item{...}{Unused, kept for S3 compatibility.}
}
\value{
\code{summary.lmvt()} returns an object of class \code{"summary.lmvt"}, a
list containing:
\itemize{
  \item \code{model}: model-level diagnostics;
  \item \code{nonzero}: counts of selected/nonzero coefficients by parameter block;
  \item \code{table}: the coefficient summary table after filtering/sorting;
  \item \code{full_table}: the complete coefficient table before filtering.
}
}
\examples{
\donttest{
set.seed(12)
sim <- simulate_scenario(scen = 1, S = 4, T = 60, d_noise = 4,
                         noise_kind = "iid", seed = 12)
df  <- sim[, c("s", "t", "y", grep("^X|^Noise|^Xbin",
                                   names(sim), value = TRUE))]
fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
            lambda_beta = 0.05, lambda_gamma = 0.05, maxit = 8)

# Printed summary
summary(fit)

# Programmatic coefficient table
sm <- summary(fit, include_zero = FALSE, sort_by = "abs_estimate")
sm$table
}
}
\seealso{\code{\link{lmvt}}, \code{\link{coef.lmvt}}, \code{\link{predict.lmvt}}}
