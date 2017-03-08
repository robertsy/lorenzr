## utils function, useful for both models.
#
# PUBLIC:  lorenz_plot, lorenz_heatmap, lorenz_observe, lorenz_as_df, lorenz_run_as_df, lorenz_da_eval, lorenz_plot_eval
# PRIVATE: msqrt

###################################################
## PUBLIC:




lorenz_plot <- function(lorenz_state, lorenz_obs, ens=NULL, tit='',ylim,...){
  yrange <- range(c(lorenz_state, ens, lorenz_obs$y))
  if (!missing(ylim)) yrange <- ylim

  plot(lorenz_state, type='l', ylim=yrange, main=tit, lwd=2,...)
  points(lorenz_obs$y.loc, lorenz_obs$y, col='red', pch=3)
  for (k in 1:ncol(ens)){
    lines(ens[,k], col='purple')
  }
  points(lorenz_state, type='l', lwd=2, col='black') ## reprint state over
}



lorenz_heatmap <- function(l96_run,...){
  lorenz_df <- reshape2::melt(l96_run$state.ts, varnames=c('time', 'location'))
  lorenz_df %>%
    ggplot(aes(x=time, y=location, fill=value)) +
    geom_raster() +
    scale_fill_gradient2(midpoint = 2.4)## 2.4: climatological mean
}


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




## apply lorenz_as_df to a model run (cycled DA)
## return a nice df object
## example:
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


lorenz_plot_eval <- function(ens_eval){
  g <-
    filter(ens_eval, error_type!='nobs') %>%
    ggplot(aes(x=time, y=error, linetype=type, color=method)) + geom_line() + facet_wrap(~error_type, scales = 'free')
}


###################################################
## PRIVATE:
msqrt <- function(A){
  if(ncol(A)==1) return(sqrt(A))
  e <- eigen(A)
  V <- e$vectors
  return(V %*% diag(sqrt(e$values)) %*% t(V))
}



