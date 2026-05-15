\name{varGuidTS-package}
\alias{varGuidTS-package}
\alias{varGuidTS}
\docType{package}
\title{varGuidTS: Variance-Guided Time-Series Modeling for Temporal Risk Detection}
\description{
Fits a balanced panel ARX-GARCHX model with subject-specific
intercepts, shared population coefficients, and joint L1
penalization on the conditional mean and variance. Provides an
exceedance-based risk score
\eqn{\pi_{s,t}(c) = P(Y_{s,t} > c \,|\, F_{t-1})} for personalized
temporal risk assessment from wearable-sensor and other high-frequency
panel time-series data.
}
\author{
Maintainer: Min Lu \email{luminwin@gmail.com}\cr
Authors: Zihao Wang, Min Lu
}
\references{
Wang Z. and Lu M. (2026). Variance-Aware Penalized Panel Models for
Temporal Risk Detection from Wearable Sensor Data.
}
\seealso{
\code{\link{lmvt}}, \code{\link{predict.lmvt}},
\code{\link{summary.lmvt}}, \code{\link{simulate_scenario}}
}
\keyword{package}
