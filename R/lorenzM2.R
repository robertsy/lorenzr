## Lorenz M2 model



#' Simulate data from the lorenz MII model
#'
#' @param kint parameter K in model specification, related to number of waves
#' @inheritParams l96_simulate
#' @return list with: state.ts matrix of nsteps x ndim
#' y.ts list of nsteps lists, each of them containing y, R, H, d and y.loc (sites which are observed)
#' f.propagate function to propagate the ensemble according to the model specification 
#' f.propagate <- function(state=initial state, nsteps=how many steps to propagate, ...)
#' ndim, duration and freq to keep track of some parameters
#' @examples 
#' lm2_run <- lm2_simulate(10 * .2, .2, sig=3)
#' lorenz_plot(lm2_run$state.ts[11,], lm2_run$y.ts[[11]])
#' lorenz_heatmap(lm2_run) 
lm2_simulate <- function(duration, freq, ndim=960, kint=32, Forcing=15,
                         sig=sqrt(0.5), obs_type='partial', R_sig=sig, nobs=ndim/20, deltat=1/200){

  ## initialize random state0 with a burning time:
  state <- lm2_integrate( runif(ndim, -1, 1), 10/deltat, ndim=ndim, kint=kint, Forcing=Forcing, deltat=deltat)

  ## result in state.ts matrix, with row=time, col=state:
  state.ts <- matrix(NA, nrow=ceiling(duration/freq)+1, ncol=ndim)
  state.ts[1,] <- state

  ## integrate for time=duration and record every nsteps=freq:
  for (i in 2:((duration/freq)+1)){
    state.ts[i,] <- lm2_integrate(state.ts[i-1,], nsteps=freq/deltat, ndim=ndim, kint=kint, Forcing=Forcing, deltat=deltat)
  }

  ## observations:
  ## y.ts is a list with y, H, R, etc.
  y.ts <- lapply(seq_len(nrow(state.ts)), function(i)
    lorenz_observe(state.ts[i,], ndim=ndim, sig=sig, obs_type=obs_type, R_sig=R_sig, nobs=nobs))

  ## propagate function closure:
  ## (include all the parameters, ku, kr, alpha, etc!)
  f.propagate <- function(state, nsteps, ...){
    lm2_integrate(state, nsteps/deltat, ndim=ndim, kint=kint, Forcing=Forcing, deltat=deltat)
  }

  return(list(state.ts=state.ts, y.ts=y.ts, f.propagate=f.propagate,
              ndim=ndim, duration=duration, freq=freq))

}



#' Simulate an initial ensemble with lorenz MII
#' 
#' @description 
#' For each ensemble member, start with iid uniform(-1,1) state and integrate for duration
#'
#' @param K ensemble size
#' @param lm2_run output from lm2_simulate
#' @return ndim x K matrix of initial ensemble 
#' @examples 
#' lm2_run <- lm2_simulate(10 * .2, .2, sig=3)
#' ens0 <- lm2_ens0(10, lm2_run)
#' lorenz_plot(lm2_run$state.ts[11,], lm2_run$y.ts[[11]], ens0)
lm2_ens0 <-  function(K, lm2_run, duration=5,...){
  ens0 <- replicate(K, runif(lm2_run$ndim, -1,1))
  ens0 <- lm2_run$f.propagate(ens0, duration)
  return(ens0)
}



#' Integrate state forward with lorenz63
#' 
#' @description 
#' Time integration with lorenz MII: implemented with a euler scheme in fortran90.
#'
#' @inheritParams lm2_simulate
#' @inheritParams l96_simulate
#' @return state at the end of integration (either a vector or a matrix depending on initial state)
#' @examples 
#' current_state <- runif(960, -1, 1)
#' par(mfrow=c(5,1), mai=c(0.2,0.4,.5,.2))
#' for (i in 1:5){
#'  plot.ts(current_state, main=paste('Time=',i), xlab='', ylab='')
#'  current_state <- lm2_integrate( current_state, 500, deltat=1/200)
#' }
lm2_integrate <- function(state, nsteps=1, deltat=1/200, ndim=960, kint=32, Forcing=15){

  if ( kint %%2 != 0 ) warning('kint must be an even number!')

  ## set number of ensemble (1 if just the state):
  kens <- 1
  if (!is.null(dim(state))) {
    ## unstack ensemble matrix as long vector (ens_1, ens_2, ...)
    kens <- ncol(state)
    state <- c(state)
  }

  ## integrate nsteps in time:
  output <- .Fortran('lm2integrate',
                     state = as.double(state),
                     kens=as.integer(kens),
                     ndim=as.integer(ndim),
                     kint=as.integer(kint),
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
