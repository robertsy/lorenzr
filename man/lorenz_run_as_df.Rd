% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lorenz_all.R
\name{lorenz_run_as_df}
\alias{lorenz_run_as_df}
\title{Apply lorenz_as_df to a model run}
\usage{
lorenz_run_as_df(model_run, duration, freq, lorenz_run = NULL)
}
\arguments{
\item{model_run}{as returned from a cycling experiment}

\item{duration}{total time to integrate model (duration=freq*(nsteps-1)  (nsteps = duration/freq + 1)}

\item{freq}{frequency of time integration (freq=duration/(nsteps-1))}

\item{lorenz_run}{reference run as returned from l**_simulate}
}
\value{
tbl object with time evolution of ensemble and truth
}
\description{
used to take a model run lists, transform into a data frame that can be used to
evaluate forecasts in lorenz_da_eval
}
