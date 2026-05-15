\name{predict.lmvt}
\alias{predict.lmvt}
\title{Predict from a fitted \code{lmvt} object}
\description{
Computes the conditional mean \eqn{\hat\mu_{s,t}}, conditional
standard deviation \eqn{\hat\sigma_{s,t}}, and the exceedance-based
risk score
\eqn{\hat\pi_{s,t}(c) = P(Y_{s,t} > c \,|\, F_{t-1})} on a new (or
held-out) panel data frame.
}
\usage{
\method{predict}{lmvt}(object, newdata, threshold = NULL,
        innov_g = TRUE, innov_t = FALSE, df_t = 6, ...)
}
\arguments{
\item{object}{A fitted \code{lmvt} object.}
\item{newdata}{A data frame containing columns \code{s}, \code{t},
\code{y} and the same covariate columns used during fitting.}
\item{threshold}{Threshold \eqn{c}. Numeric scalar, vector of length
\code{nrow(newdata)}, or per-subject vector. Defaults to the pooled
90th percentile of \code{newdata$y}.}
\item{innov_g}{Logical; return Gaussian-innovation exceedance
probability \code{risk_g}.}
\item{innov_t}{Logical; return Student-t-innovation exceedance
probability \code{risk_t}.}
\item{df_t}{Degrees of freedom for the standardized Student-t innovation. Must be greater than 2 when \code{innov_t = TRUE}.}
\item{...}{Unused, kept for S3 compatibility.}
}
\value{
A data frame with columns \code{s}, \code{t}, \code{yhat = muhat},
\code{muhat}, \code{sigmahat}, \code{threshold_used}, and one or both
of \code{risk_g}, \code{risk_t}.
}
\examples{
\donttest{
set.seed(1)
sim <- simulate_scenario(scen = 2, S = 4, T = 70, d_noise = 4,
                         noise_kind = "iid", seed = 1)
df  <- sim[, c("s", "t", "y", grep("^X|^Noise|^Xbin",
                                   names(sim), value = TRUE))]
fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
            lambda_beta = 0.05, lambda_gamma = 0.05, maxit = 8)

# 1) Pooled/global threshold: one cutoff is applied to every row.
pr_global <- predict(fit, df, threshold = stats::quantile(df$y, 0.90))
head(pr_global)

# 2) Subject-specific thresholds: use a named vector indexed by subject ID.
subject_cutoffs <- tapply(df$y, df$s, stats::quantile, probs = 0.75,
                          na.rm = TRUE)
pr_subject <- predict(fit, df, threshold = subject_cutoffs)
head(pr_subject)

# 3) Return both Gaussian and standardized Student-t risk scores.
# Here df_t = 8 means 8 degrees of freedom for the Student-t
# innovation distribution; it is not the data frame named df.
pr_t <- predict(fit, df, threshold = subject_cutoffs,
                innov_g = TRUE, innov_t = TRUE, df_t = 8)
head(pr_t[, c("s", "t", "muhat", "sigmahat", "risk_g", "risk_t")])
}
}
\seealso{\code{\link{lmvt}}, \code{\link{summary.lmvt}}}
