## Lorenz 96 model
#
# PUBLIC:  l96_simulate, l96_ens0
# PRIVATE: l96_integrate


# ###################################################
# # ## PUBLIC:
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
    l96_integrate(state, nsteps/deltat, ndim=ndim, Forcing = Forcing, deltat=deltat)
  }

  return(list(state.ts=state.ts, y.ts=y.ts, f.propagate=f.propagate,
              ndim=ndim, duration=duration, freq=freq))
}



l96_ens0 <-  function(K, l96_run, ...){
  ens0 <- replicate(K, rnorm(l96_run$ndim))
  ens0 <- l96_run$f.propagate(ens0, 250)
  return(ens0)
}





###################################################
## PRIVATE:
l96_integrate <- function(state, nsteps=1, deltat=1/1000, ndim=40, Forcing=8, forcingvec=rep(Forcing,ndim)){

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
                     deltat=as.double(deltat),
                     forcingvec=as.double(forcingvec))


  ## if ensemble, unstack the long vector to matrix:
  ## (ens_1 | ens_2 | ...)
  if (kens > 1){
    state <- t(matrix(output$state, nrow=kens, byrow=TRUE))
  } else{
    state <- output$state
  }

  return(state)
}


