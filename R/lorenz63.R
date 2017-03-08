

l63_simulate <- function(duration, freq, params=c(sigma=10, rho=28, beta=8/3),
                         sig=sqrt(0.5), obs_type='one', R_sig=sig, deltat=1/75,...){
  ## initialize random state0 with a burning time:
  state <- l63_integrate( rnorm(3), nsteps=2/deltat, deltat=deltat, params=params)

  ## result in state.ts matrix, with row=time, col=state:
  state.ts <- matrix(NA, nrow=ceiling(duration/freq)+1, ncol=3)
  state.ts[1,] <- state

  ## integrate for time=duration and record every nsteps=freq:
  for (i in 2:((duration/freq)+1)){
    state.ts[i,] <- l63_integrate(state.ts[i-1,], nsteps=freq/deltat, deltat=deltat, params=params)
  }

  ## observations:
  ## y.ts is a list with y, H, R, etc.
  y.ts <- lapply(seq_len(nrow(state.ts)), function(i)
    lorenz_observe(state.ts[i,], ndim=3, sig=sig, obs_type = obs_type, R_sig=R_sig,...))

  ## propagate function closure:
  ## (include all the parameters
  f.propagate <- function(state, nsteps, ...){
    l63_integrate(state, nsteps=nsteps/deltat, deltat=deltat,params=params)
  }

  return(list(state.ts=state.ts, y.ts=y.ts, f.propagate=f.propagate,
              ndim=3, duration=duration, freq=freq))
}



l63_ens0 <-  function(K, l63_run, ...){
  ens0 <- replicate(K, rnorm(l63_run$ndim))
  ens0 <- l63_run$f.propagate(ens0, 250)
  return(ens0)
}


l63_integrate <- function(state, nsteps=1, deltat=1/75,  params){
  ## set number of ensemble (1 if just the state):
  kens <- 1
  if (!is.null(dim(state))) {
    ## unstack ensemble matrix as long vector (ens_1, ens_2, ...)
    kens <- ncol(state)
    state <- c(state)
  }

  ## integrate nsteps in time:
  output <- .Fortran('l63integrate',
                     state = as.double(state),
                     kens=as.integer(kens),
                     params=as.double(params),
                     nsteps=as.integer(nsteps),
                     deltat=as.double(deltat)
                     )


  ## if ensemble, unstack the long vector to matrix:
  ## (ens_1 | ens_2 | ...)
  if (kens > 1){
    state <- t(matrix(output$state, nrow=kens, byrow=TRUE))
  } else{
    state <- output$state
  }

  return(state)
}

