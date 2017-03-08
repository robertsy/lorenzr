## Lorenz M2 model
#
# PUBLIC:  lm2_simulate, lm2_ens0
# PRIVATE: lm2_integrate

# ###################################################
# # ## PUBLIC:
lm2_simulate <- function(duration, freq, ndim=960, kint=15, Forcing=15,
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




lm2_ens0 <-  function(K, lm2_run, duration=5,...){
  ens0 <- replicate(K, runif(lm2_run$ndim, -1,1))
  ens0 <- lm2_run$f.propagate(ens0, duration)
  return(ens0)
}




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
