## Lorenz 96 model



#' Simulate data from the lorenz96 model
#'
#' @param duration total time to integrate model (duration=freq*(nsteps-1)  (nsteps = duration/freq + 1)
#' @param freq frequency of time integration (freq=duration/(nsteps-1))
#' @param ndim dimension of system
#' @param Forcing external forcing on dynamics
#' @param sig standard deviation of error measurments
#' @param obs_type (all: all sites, odd: every other site, 
#' partial: regular observations at nobs sites, one: unique observation in the middle)
#' @param R_sig standard deviation of error measurements for assimilation (default=sig) 
#' @param deltat time interval for model integration
#' @param ... additional arguments passed to lorenz_observe (e.g. nobs in case obs_type=partial)
#' @return list with: state.ts matrix of nsteps x ndim
#' y.ts list of nsteps lists, each of them containing y, R, H, d and y.loc (sites which are observed)
#' f.propagate function to propagate the ensemble according to the model specification 
#' f.propagate <- function(state=initial state, nsteps=how many steps to propagate, ...)
#' ndim, duration and freq to keep track of some parameters
#' @examples 
#' l96_run <- l96_simulate(0.1*100, 0.1, ndim=40, Forcing=8, deltat=1/1000)
#' lorenz_plot(l96_run$state.ts[101,], l96_run$y.ts[[101]])
l96_simulate <- function(duration, freq, ndim=40, Forcing=8,
                         sig=sqrt(0.5), obs_type='odd', R_sig=sig, deltat=1/1000,...){
  ## initialize random state0 with a burning time:
  state <- l96_integrate( rnorm(ndim), 2/deltat, ndim=ndim, Forcing=Forcing, deltat=deltat)

  ## result in state.ts matrix, with row=time, col=state:
  state.ts <- matrix(NA, nrow=ceiling(duration/freq)+1, ncol=ndim)
  state.ts[1,] <- state

  ## integrate for time=duration and record every nsteps=freq:
  for (i in 2:((duration/freq)+1)){
    state.ts[i,] <- l96_integrate(state.ts[i-1,], nsteps=freq/deltat, ndim=ndim, Forcing=Forcing, deltat=deltat)
  }

  ## observations:
  ## y.ts is a list with y, H, R, etc.
  y.ts <- lapply(seq_len(nrow(state.ts)), function(i)
    lorenz_observe(state.ts[i,], ndim=ndim, sig=sig, obs_type = obs_type, R_sig=R_sig,...))

  ## propagate function closure:
  ## (include all the parameters, ku, kr, alpha, etc!)
  f.propagate <- function(state, nsteps, ...){
    l96_integrate(state, nsteps/deltat, ndim=ndim, Forcing = Forcing, deltat=deltat, ...)
  }

  return(list(state.ts=state.ts, y.ts=y.ts, f.propagate=f.propagate,
              ndim=ndim, duration=duration, freq=freq))
}


#' Simulate an initial ensemble with lorenz96
#' 
#' @description 
#' For each ensemble member, start with iid normal state and integrate for nsteps time steps (250)
#'
#' @param K ensemble size
#' @param l96_run output from l96_simulate
#' @return ndim x K matrix of initial ensemble 
#' @examples 
#' l96_run <- l96_simulate(0.1*100, 0.1, ndim=40, Forcing=8, deltat=1/1000)
#' ens0 <- l96_ens0(10, l96_run, 100)
#' lorenz_plot(l96_run$state.ts[101,], l96_run$y.ts[[101]], ens0)
l96_ens0 <-  function(K, l96_run, nsteps=250,...){
  ens0 <- replicate(K, rnorm(l96_run$ndim))
  ens0 <- l96_run$f.propagate(ens0, nsteps, ...)
  return(ens0)
}





#' Integrate state forward with lorenz96
#' 
#' @description 
#' Time integration with lorenz96: implemented with a euler scheme in fortran90.
#'
#' @param state initial state (can be an ensemble matrix of ndim x K)
#' @param nsteps how many steps ahead
#' @inheritParams l96_simulate
#' @return state at the end of integration (either a vector or a matrix depending on initial state)
#' @examples 
#' current_state <- rnorm(40)
#' par(mfrow=c(4,2), mai=c(0.2,0.4,.5,.2))
#' for (i in 1:8){
#'  plot.ts(current_state, main=paste('Time=',i), xlab='', ylab='')
#'  current_state <- l96_integrate( current_state, 500, deltat=1/1000)
#' }
l96_integrate <- function(state, nsteps=1, deltat=1/1000, ndim=40, Forcing=8, ...){

  ## set number of ensemble (1 if just the state):
  kens <- 1
  if (!is.null(dim(state))) {
    ## unstack ensemble matrix as long vector (ens_1, ens_2, ...)
    kens <- ncol(state)
    state <- c(state)
  }

  ## integrate nsteps in time:
  output <- .Fortran('l96integrate',
                     state = as.double(state),
                     kens=as.integer(kens),
                     ndim=as.integer(ndim),
                     forcing=as.double(Forcing),
                     nsteps=as.integer(nsteps),
                     deltat=as.double(deltat))


  ## if ensemble, unstack the long vector to matrix:
  ## (ens_1 | ens_2 | ...)
  if (kens > 1){
    state <- t(matrix(output$state, nrow=kens, byrow=TRUE))
  } else{
    state <- output$state
  }

  return(state)
}


