## utils function, useful for both models.



#' plot state as a line
#'
#' @param lorenz_state vector of state
#' @param lorenz_obs list object with y and y.loc
#' @param ens ensemble matrix (can be passed as a list for multiple ensembles)
#' @param tit title of plot
#' @param ylim passed to plot()
#' @param xlab passed to plot()
#' @param ylab passed to plot()
#' @param ... extra parameters passed to plot()
lorenz_plot <- function(lorenz_state, lorenz_obs, ens=NULL, tit='', ylim, xlab='site', ylab='state',...){
  yrange <- range(c(lorenz_state, ens, lorenz_obs$y))
  if (!missing(ylim)) yrange <- ylim
  plot(lorenz_state, type='l', ylim=yrange, main=tit, lwd=2, xlab=xlab, ylab=ylab,...)
  points(lorenz_obs$y.loc, lorenz_obs$y, col='red', pch=3)
  if (!is.null(ens)){
    ## plot ensembles successively (with different color)
    if (!is.list(ens)) ens <- list(ens) ## if only one ensemble, make it a list of length 1
    colvec <- RColorBrewer::brewer.pal(max(length(ens),3), "Accent")
    for (ensi in 1:length(ens)){
      for (k in 1:ncol(ens[[ensi]])){
        lines(ens[[ensi]][,k], col=colvec[ensi])
      }
    }
  }
  points(lorenz_state, type='l', lwd=1, col='black') ## reprint state over
}




#' plot time evolution of the state as a heatmap 
#' 
#' @description 
#' x-axis=sites, y-axis=time, color=state
#'
#' @param l_run output of l**_simulate
#' @param vmean value of state around which to spread the color (e.g. 2.4=climatological mean of lorenz96)
#' @param ... extra parameters in case (not used)
lorenz_heatmap <- function(l_run, vmean=2.4, ...){
  lorenz_df <- reshape2::melt(l_run$state.ts, varnames=c('time', 'location'))
  lorenz_df %>%
    ggplot(aes(x=time, y=location, fill=value)) +
    geom_raster() +
    scale_fill_gradient2(midpoint = vmean) + ## 2.4: climatological mean (for l96)
    theme_bw()
}




#' Observation operator
#' 
#' @description 
#' used in l**_simulate to produce observations
#'
#' @param state vector of state
#' @inheritParams l96_simulate
lorenz_observe <- function(state, sig=sqrt(0.5), ndim, obs_type, R_sig=sig, nobs=ndim/10){
  ## observations parameters:
  H.full <- diag(ndim) #every location is observed

  stopifnot( obs_type %in% c('all', 'odd', 'partial', 'one'))
  y.loc <- switch (obs_type,
                   'all'= 1:ndim,
                   'odd'= seq(1,ndim, by=2),
                   'partial'=  {
                     yspace <- floor( ndim/(nobs+1) )
                     seq(yspace, ndim-yspace, by=yspace)},
                   'one'=      round(ndim/2)
  )

  d <- length(y.loc)
  H <- H.full[y.loc,,drop=FALSE]

  y <- rnorm(d, H%*%state, sig)

  R <- diag(d)*R_sig^2 #R diagonal
  R2 <- msqrt(R)

  return( list (y=y, R=R, H=H, d=d, y.loc=y.loc))
}


#' Put state in long format
#' 
#' @param state state vector
#' @return tbl object with value and x(=site) as variables
lorenz_as_df <- function(state){
  ## ensemble or state?
  is_ens <- is.matrix(state)

  if(is_ens) {
    state0 <- state[,1]
  } else {
    state0 <- state
  }

  ## domain information:
  ndim <- length(state0)

  if (is_ens){
    output <- data.frame(state)
    colnames(output) <- paste('ens_',1:ncol(state), sep='')
    output$x <- 1:ndim
    output <- gather(output, ensemble, value, -c(x))
  } else{
    output <- data.frame( value=state,  x=1:ndim)
  }

  return( tbl_df(output ) )
}



#' Apply lorenz_as_df to a model run
#' 
#' @description 
#' used to take a model run lists, transform into a data frame that can be used to
#' evaluate forecasts in lorenz_da_eval
#' 
#' @param model_run as returned from a cycling experiment
#' @param lorenz_run reference run as returned from l**_simulate
#' @inheritParams  l96_simulate
#' @return tbl object with time evolution of ensemble and truth
lorenz_run_as_df <- function(model_run, duration, freq, lorenz_run=NULL){
  time_vec <- seq(0,duration, by=freq)

  ens_as_df <- function(ens, ens_name){
    ens_all <- lapply( ens, lorenz_as_df )
    n_one_ens <- nrow(ens_all[[1]])
    ens_all <- bind_rows(ens_all)
    ens_all %>%
      mutate( time=rep(time_vec, each=n_one_ens ) ) %>%
      mutate( type=ens_name)
  }


  ## ENSEMBLE
  ## make each ensemble a df (with time as variable):
  ensB_df <- ens_as_df(model_run$ensB, 'forecast')
  ensA_df <- ens_as_df(model_run$ensA, 'analysis')

  ## bind together
  ens_all <- bind_rows(ensB_df, ensA_df)

  ## TRUE STATE
  ## make state a list:
  if(!missing(lorenz_run)){
    state_ls <- split(t(lorenz_run$state.ts), rep(1:nrow(lorenz_run$state.ts), each = ncol(lorenz_run$state.ts)))
    state_ts <- ens_as_df( state_ls , 'state') %>%
      rename(true_value=value)

    ## JOIN TOGETHER:
    ens_all <- ens_all %>% left_join( select(state_ts, -type) , by=c('x', 'time'))
  }
  return(ens_all)
}




#' Compute diagnostics of cycling experiment
#' 
#' @description 
#' compute error measures (bias, mse, rmse, crps) per time and 
#' per ensemble type (analysis, forecast)
#' 
#' @param model_run as returned from a cycling experiment
#' @param lorenz_run reference run as returned from l**_simulate
#' @return tbl object with bias, mse, rmse, nobs, crps for background/analysis ensembles
lorenz_da_eval <- function(lorenz_run, model_run){

  ## transform to nice df:
  ens_df <- lorenz_run_as_df(model_run, duration, freq, lorenz_run=lorenz_run)

  ## compute the mean ensemble at each location/time/etc:
  mean_df <- ens_df %>%
    group_by(time, type, x, true_value) %>%
    summarise( mean_value=mean(value))

  ## Compute RMSE and bias for the mean ensemble:
  diagnostic_df <- mean_df %>%
    group_by( time, type) %>%
    summarise( bias= mean(true_value-mean_value),                 ## BIAS
               mse =  mean((true_value-mean_value)^2),            ## MSE
               rmse= sqrt( mse ),                                 ## RMSE
               nobs=n()/K )

  ## First compute CRPS by location/time/etc
  ## then average the CRPS over location.
  crps_df <- ens_df %>%
    group_by( time, type, x) %>%
    summarise(crps=mean( f_crps(value, true_value))) %>%
    group_by( time, type) %>%
    summarise( crps=mean(crps))


  ## join both errors together:
  diagnostic_df <- diagnostic_df %>%
    left_join(crps_df, by=c("time", "type")) %>%
    gather( error_type, error, c(bias, mse, rmse, nobs, crps))

  return(diagnostic_df)
}




#' Plot time evolution of errors
#'
#' @param ens_eval data frame returned from lorenz_da_eval
lorenz_plot_eval <- function(ens_eval){
  g <-
    filter(ens_eval, error_type!='nobs') %>%
    ggplot(aes(x=time, y=error, linetype=type, color=method)) + geom_line() + facet_wrap(~error_type, scales = 'free')
  print(g)
}


#' Symmetric matrix square root of a positive definite matrix
#' 
#' @description 
#' use eigenvalue decompostion
#'
#' @param A a symmetric matrix
#' @return S such that S'S = A
msqrt <- function(A){
  if(ncol(A)==1) return(sqrt(A))
  e <- eigen(A)
  V <- e$vectors
  return(V %*% diag(sqrt(e$values)) %*% t(V))
}



