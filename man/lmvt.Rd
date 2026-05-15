\name{lmvt}
\alias{lmvt}
\title{Fit a penalized panel ARX-GARCHX model}
\description{
Joint maximum-quasi-likelihood fit of the conditional mean and variance
across a balanced panel of subjects, with subject-specific intercepts,
shared population coefficients, and L1 penalties on the mean and
variance covariate vectors.
}
\usage{
lmvt(data,
     p = 1, q = 0,
     r = 1, s_ord = 1,
     lambda_beta = 0,
     lambda_gamma = 0,
     maxit = 200, tol = 1e-4,
     standardize_X = TRUE,
     use_x_in_variance = TRUE,
     phi_cap = 0.995,
     omega_min = 1e-6,
     omega_cap_mult = 10,
     gamma_steps = 8,
     gamma_step = 1e-3,
     verbose = FALSE)
}
\arguments{
\item{data}{A data frame containing columns \code{s} (subject ID),
\code{t} (integer time), \code{y} (outcome), and any number of
covariate columns.}
\item{p}{Autoregressive order in the mean.}
\item{q}{Distributed-lag order for covariates \eqn{x}.}
\item{r}{ARCH order in the variance.}
\item{s_ord}{GARCH order in the variance.}
\item{lambda_beta}{L1 penalty for the mean covariate coefficients.}
\item{lambda_gamma}{L1 penalty for the variance covariate coefficients.}
\item{maxit}{Maximum number of outer block-coordinate iterations.}
\item{tol}{Absolute tolerance on the QML objective for convergence.}
\item{standardize_X}{Center and scale covariates before fitting.}
\item{use_x_in_variance}{Include the lagged-X block in the variance equation.}
\item{phi_cap}{Upper bound on \eqn{\sum a + \sum b}.}
\item{omega_min}{Lower bound on subject-specific baseline variance.}
\item{omega_cap_mult}{Multiplier on \code{median(eps^2)} used as the
upper bound for \eqn{\omega_s} in the L-BFGS-B update.}
\item{gamma_steps}{Inner proximal-gradient steps for \eqn{\gamma} per outer iter.}
\item{gamma_step}{Step size for the proximal-gradient update of \eqn{\gamma}.}
\item{verbose}{Logical; print per-iteration diagnostics.}
}
\value{
An object of class \code{"lmvt"}: a list with elements
\code{alpha, theta, beta, omega, a, b, gamma, ids, xcols, x_center,
x_scale, orders, penalties, call, converged, iters, mean_times,
var_times, iter_times}.
}
\references{
Wang Z. and Lu M. (2026). Variance-Aware Penalized Panel Models for
Temporal Risk Detection from Wearable Sensor Data.
}
\examples{
set.seed(1)
sim <- simulate_scenario(scen = 1, S = 5, T = 80, d_noise = 5,
                         noise_kind = "iid", seed = 1)
df  <- sim[, c("s","t","y", grep("^X|^Noise|^Xbin",
                                 names(sim), value = TRUE))]
fit <- lmvt(df, p = 1, q = 0, r = 1, s_ord = 1,
            lambda_beta = 0.05, lambda_gamma = 0.05)
print(fit)
}
\seealso{\code{\link{predict.lmvt}}, \code{\link{summary.lmvt}}, \code{\link{simulate_scenario}}}
