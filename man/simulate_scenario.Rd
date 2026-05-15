\name{simulate_scenario}
\alias{simulate_scenario}
\title{Simulate a panel ARX-GARCHX dataset}
\description{
Simulates a balanced panel from one of the four canonical scenarios
in Wang and Lu (2026), Table I. Returns a long-format data frame with
the response \code{y}, oracle mean and variance, the covariate block,
optional regime indicator, and the oracle exceedance probability.
}
\usage{
simulate_scenario(scen = 1L, S = 8L, T = 240L, L = 20L,
                  d_noise = 20L, noise_kind = c("iid", "ar1"),
                  threshold_c = 1.0, seed = NULL,
                  floor_sigma2 = 1e-10)
}
\arguments{
\item{scen}{Integer 1:4 selecting one of the four scenarios.}
\item{S}{Number of subjects.}
\item{T}{Number of post-burn-in time points per subject.}
\item{L}{Mean run length for the piecewise regime in Scenario 4.}
\item{d_noise}{Number of nuisance covariates with zero true effect.}
\item{noise_kind}{Either \code{"iid"} or \code{"ar1"}.}
\item{threshold_c}{Threshold \eqn{c} at which the oracle exceedance
probability is computed.}
\item{seed}{Random seed.}
\item{floor_sigma2}{Lower clamp on simulated variance.}
}
\value{
A data frame with columns \code{s}, \code{t}, \code{y}, oracle
\code{mu} and \code{sigma}, covariate columns, optional \code{Xbin},
and the oracle exceedance probability \code{pi_true}. Generating
parameters are stored in \code{attr(., "params")}.
}
\examples{
set.seed(7)
sim <- simulate_scenario(scen = 2, S = 4, T = 50,
                         d_noise = 3, noise_kind = "iid", seed = 7)
head(sim)
}
\references{
Wang Z. and Lu M. (2026). Variance-Aware Penalized Panel Models for
Temporal Risk Detection from Wearable Sensor Data.
}
