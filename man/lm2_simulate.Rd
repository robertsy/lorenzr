% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lorenzM2.R
\name{lm2_simulate}
\alias{lm2_simulate}
\title{Simulate data from the lorenz MII model}
\usage{
lm2_simulate(duration, freq, ndim = 960, kint = 32, Forcing = 15,
  sig = sqrt(0.5), obs_type = "partial", R_sig = sig, nobs = ndim/20,
  deltat = 1/200)
}
\arguments{
\item{duration}{total time to integrate model (duration=freq*(nsteps-1)  (nsteps = duration/freq + 1)}

\item{freq}{frequency of time integration (freq=duration/(nsteps-1))}

\item{ndim}{dimension of system}

\item{kint}{parameter K in model specification, related to number of waves}

\item{Forcing}{external forcing on dynamics}

\item{sig}{standard deviation of error measurments}

\item{obs_type}{(all: all sites, odd: every other site, 
partial: regular observations at nobs sites, one: unique observation in the middle)}

\item{R_sig}{standard deviation of error measurements for assimilation (default=sig)}

\item{deltat}{time interval for model integration}
}
\value{
list with: state.ts matrix of nsteps x ndim
y.ts list of nsteps lists, each of them containing y, R, H, d and y.loc (sites which are observed)
f.propagate function to propagate the ensemble according to the model specification 
f.propagate <- function(state=initial state, nsteps=how many steps to propagate, ...)
ndim, duration and freq to keep track of some parameters
}
\description{
Simulate data from the lorenz MII model
}
\examples{
lm2_run <- lm2_simulate(10 * .2, .2, sig=3)
lorenz_plot(lm2_run$state.ts[11,], lm2_run$y.ts[[11]])
lorenz_heatmap(lm2_run) 
}
