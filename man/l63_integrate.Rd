% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lorenz63.R
\name{l63_integrate}
\alias{l63_integrate}
\title{Integrate state forward with lorenz63}
\usage{
l63_integrate(state, nsteps = 1, deltat = 1/75, params)
}
\arguments{
\item{deltat}{time interval for model integration}

\item{params}{named vector (sigma, rho, beta) for the dynamics}
}
\value{
state at the end of integration (either a vector or a matrix depending on initial state)
}
\description{
Time integration with lorenz63: implemented with a euler scheme in fortran90.
}
\examples{
state1 <- l63_integrate(rnorm(3), params=c(sigma=10, rho=28, beta=8/3), nsteps=10, deltat=1/75)
}
