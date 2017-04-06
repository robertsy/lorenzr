



#' Simulate data from the lorenz63 model
#'
#' @inheritParams l96_simulate
#' @param params named vector (sigma, rho, beta) for the dynamics
#' @examples 
#' l63_run <- l63_simulate(100*1/5, 1/5, deltat=1/75, sig=3)
#' lorenz_plot(l63_run$state.ts[101,], l63_run$y.ts[[101]])
#' l63_df <- lorenz_as_df(l63_run$state.ts) %>%
#'           rename(time=x, variable=ensemble)
#' l63_df$variable <- factor(l63_df$variable, levels = c('ens_1', 'ens_2', 'ens_3'), labels=c('x(1)', 'x(2)', 'x(3)'))
#' l63_y_df <-
#'   foreach(i=1:length(l63_run$y.ts), .combine='rbind')%do%{
#'     yy <- l63_run$y.ts[[i]]
#'     data_frame(variable=paste('x(',yy$y.loc,')',sep=''), value=yy$y, time=i)
#'   }
#' l63_df %>%
#'   ggplot(aes(x=time, y=value)) + geom_line() +
#'   geom_point(data=l63_y_df, aes(x=time, y=value), color='red')+
#'   facet_wrap(~variable, ncol=1) +
#'   theme_bw()
#' @return list with: state.ts matrix of nsteps x ndim
#' y.ts list of nsteps lists, each of them containing y, R, H, d and y.loc (sites which are observed)
#' f.propagate function to propagate the ensemble according to the model specification 
#' f.propagate <- function(state=initial state, nsteps=how many steps to propagate, ...)
#' ndim, duration and freq to keep track of some parameters
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




#' Simulate an initial ensemble
#' 
#' @description 
#' For each ensemble member, start with iid normal state and integrate for nsteps time steps (250)
#'
#' @param K ensemble size
#' @param l63_run output from l63_simulate
#' @return ndim x K matrix of initial ensemble 
#' @examples 
#' l63_run <- l63_simulate(100*1/5, 1/5, deltat=1/75, sig=3)
#' ens0 <- l63_ens0(20, l63_run)
#' lorenz_plot(l63_run$state.ts[101,], l63_run$y.ts[[101]], ens=ens0)
l63_ens0 <-  function(K, l63_run, nsteps=250, ...){
  ens0 <- replicate(K, rnorm(l63_run$ndim))
  ens0 <- l63_run$f.propagate(ens0, nsteps, ...)
  return(ens0)
}


#' Integrate state forward with lorenz63
#' 
#' @description 
#' Time integration with lorenz63: implemented with a euler scheme in fortran90.
#'
#' @inheritParams l63_simulate  
#' @inheritParams l96_simulate
#' @return state at the end of integration (either a vector or a matrix depending on initial state)
#' @examples 
#' state1 <- l63_integrate(rnorm(3), params=c(sigma=10, rho=28, beta=8/3), nsteps=10, deltat=1/75)
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

