\name{coef.lmvt}
\alias{coef.lmvt}
\alias{print.lmvt}
\title{Extract or print coefficients from a fitted lmvt model}
\description{
\code{coef()} returns the list of fitted parameter blocks
(subject-specific intercepts, AR coefficients, mean covariate
coefficients, baseline variances, ARCH/GARCH parameters, and
variance covariate coefficients).

\code{print()} produces a one-screen summary of the fit.
}
\usage{
\method{coef}{lmvt}(object, ...)
\method{print}{lmvt}(x, ...)
}
\arguments{
\item{object, x}{A fitted \code{lmvt} object.}
\item{...}{Unused.}
}
\seealso{\code{\link{lmvt}}, \code{\link{summary.lmvt}}}
