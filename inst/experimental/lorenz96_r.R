## L96 Model by Lorentz
## Ref: "Designing Chaotic Models", Lorenz 2004 (JAS)

#############################################################
### MODEL PARAMETERS:
n <- 40
Forcing <- 8
deltat <- 1/40

Forcing_vec <- Forcing+stats::filter(rnorm(n), 0.9)




#############################################################
###

L96.diff <- function(x){
  N <- length(x)
  x.diff <- numeric(N)
  for (i in 1:N){
    xx <- L96.neighbours(x, i)
    x.diff[i] <- -(xx[1] * xx[2]) + xx[2] * xx[4] - xx[3] + Forcing_vec[i]
  }
  return(x.diff)
}


runge.kutta <- function(x, f.diff, deltat, ...){
  k1 <- f.diff(x, ...)
  k2 <- f.diff(x + deltat/2 * k1, ...)
  k3 <- f.diff(x + deltat/2 * k2, ...)
  k4 <- f.diff(x + deltat * k3, ...)
  return(1/6 * (k1 + 2 * k2 + 2 * k3 + k4))
}



#############################################################
### RING MANIPULATIONS:

mod <- function(a,b){
  # a mod b
  # with special treatment for a==0
  c <- a %% b
  ifelse(c==0, b, c)
  #ifelse(a == 0, b, a %% b)
}


L96.block <- function(x, n, l=1, ind=(n-l):(n+l)){
  # n = center of block
  # 2*l + 1 = size of block
  N <- length(x)
  block.ind <- mod(ind, N)#mod((n-l):(n+l), N)
  return(x[block.ind])
}


## compute minimal distances (discrete) on a ring...
ring.dist <- function(i, j, N){
  abs(N/2 - abs(abs(i-j) - N/2))
}

ring.alldist <- function(N){
  ind <- 1:N
  A <- matrix(NA, nrow=N, ncol=N)
  for (i in ind){
    A[i,] <- ring.dist(i, ind, N)
  }
  return(A)
}

L96.neighbours <- function(x, n){
  L96.block(x, NA, ind=(n-2):(n+1))
}








#############################################################
### INTEGRATION:

L96.integrate <- function(x, duration, deltat){
  ## duration is in unit of time (5 days)
  ## deltat is for the increment
  for (t in 1:(duration/deltat)){
    x <- x + deltat * runge.kutta(x, L96.diff, deltat)#L96.diff(x)#
  }
  return(x)
}




## example of x.integrate: L96 model
## should integrate the whole matrix xa 1 time step ahead
## it uses L96.diff() and runge.kutta() from the environment
L96.ens.integrate <- function(xa, duration=1, deltat=deltat){
  require(foreach)
  #xb <- apply(xa, 1, L96.integrate, duration=duration, deltat=deltat)
  #tt1 <- L96.integrate(xa[1,], duration, deltat)
  #return(t(xb))
  foreach(i=1:nrow(xa), .combine='rbind') %dopar% L96.integrate(xa[i,], duration=duration, deltat=deltat)
}







L96.integrate.ts <- function(state, duration, freq, deltat, ind=NULL){
  state.ts <- matrix(NA, nrow=ceiling(duration/freq)+1, ncol=length(state))
  state.ts[1,] <- state
  for (i in 2:((duration/freq)+1)){
    state.ts[i,] <- L96.integrate(state.ts[i-1,], duration=freq, deltat)
  }
  if (!is.null(ind)) state.ts <- state.ts[,ind]
  return(state.ts)
}


L96.ens.integrate.ts <- function(xa, duration, freq, deltat, ind){
  require(foreach)
  foreach(i=1:nrow(xa), .combine='rbind') %do%
    L96.integrate.ts(xa[i,], duration, freq, deltat, ind)
}











#############################################################
### VIZZ

## function to plot nicely the state of L96 with ensemble spread
## can still be greatly improved, of course
L96.plot <- function(state, ens=NULL, obs=NULL){
  require(ggplot2)
  n <- length(state)
  x <- 1:n

  if (is.null(ens)){
    toplot  <- data.frame(x=x, value=state)
    if (!is.null(obs)){
      toplot <- cbind(toplot, obs)
    }
    g <- ggplot(data=toplot, aes(x=x, y=value)) + geom_line()
  }
  else{
    ens.mean  <- apply(ens, 2, mean)
    ens.quant <- apply(ens, 2, quantile, probs=c(0.05, 0.25, 0.75, 0.95))
    rownames(ens.quant) <- c('q05', 'q25', 'q75', 'q95')
    toplot <- data.frame(x=x, value=ens.mean, truth=state,
                         t(ens.quant))
    if (!is.null(obs)){
      toplot <- cbind(toplot, obs)
    }
    g <- ggplot(data=toplot, aes(x=x, y=value, ymin=q05, ymax=q95)) + geom_line()
    g <- g + geom_ribbon(aes(ymin = q05, ymax = q95, fill = "05%-95%"), alpha = .25)
    g <- g + geom_ribbon(aes(ymin = q25, ymax = q75, fill = "25%-75%"), alpha = .25)
    g <- g + geom_line(aes(x=x, y=truth, colour="truth"))
  }
  if (!is.null(obs)){
    g <- g + geom_point(aes(x=x, y=obs))
  }
  return(g)
}





#############################################################
### EVAL

L96.eval <- function(state, ens, crps.ind=c(1,20)){
  require(foreach)
  RMSE <- function(truth, x){
    model <- colMeans(x)
    sqrt(mean((truth-model)^2))
  }


  err <- RMSE(state, ens)

  crps <- foreach (j = crps.ind, .combine='c') %do% {
    f.crps(ens[,j], state[j])
  }
  names(crps) <- as.character(crps.ind) #paste('crps', crps.ind, sep='')
  return(c(err=err, crps=crps))
}


# gaussian mixture version:
L96.gmeval <- function(state, ens, crps.ind=c(1,20)){
  require(foreach)
  RMSE <- function(truth, x){
    model <- colMeans(x)
    sqrt(mean((truth-model)^2))
  }


  err <- RMSE(state, ens$xa)

  crps <- foreach (j = crps.ind, .combine='c') %do% {
    f.crps(ens$xa[,j], state[j])
  }

  crps.gm <- foreach (j = crps.ind, .combine='c') %do% {
    #f.crps.gm(ens$ma[,j], ens$Pa[j,j], state[j])
    #browser()
    if(is.null(ens$alpha)) ens$alpha <- rep(1/K, K)
    f.crps.gm.w(ens$ma[,j], ens$Pa[j,j], state[j], ens$alpha)
  }

  names(crps) <- names(crps.gm) <- as.character(crps.ind)
  return(c(err=err, crps=crps, crps.gm=crps.gm))
}






#################################################################################
#### MODEL SETUP
obs.full <- function(x){
  return(x)
}

obs.odd <- function(x){
  return(x[seq(1,length(x),by=2)])
}

obs.hole <- function(x, hole.size=5){
  if (hole.size%%2 == 0) hole.size <- hole.size + 1
  n.hole <- (hole.size-1)/2
  hole.ind <- (n/2 - n.hole):(n/2 + n.hole)
  return(x[!1:n %in% hole.ind])
}

create.H <- function(h.operator, all=FALSE, n=n, ...){
  require(foreach)
  e0 <- rep(0,n)
  H <- foreach(i = 1:n, .combine='cbind') %do% {
    ei <- e0
    ei[i] <- 1
    h.operator(ei, ...)
  }

  # create R and d
  if(all){
    print("d and R set in global env")
    d <<- nrow(H)
    R <<- 0.5 * diag(d)
  }
  return(H)
}



## needed that because global variables don't seem
## to work in parallel environment, for some obvious
## reasons...

create.context <- function(type, n=n, ...){
  if (type=='full'){
    obs.ind <- 1:n
    H <- create.H(obs.full, n=n, ...)
  } else if(type=='odd') {
    obs.ind <- seq(1, n, by=2)
    H <- create.H(obs.odd, n=n, ...)
  } else if(type=='hole'){
    obs.ind <- 1:n
    n.hole <- 3
    obs.ind <- obs.ind[!1:n %in% (n/2 - n.hole):(n/2 + n.hole)]
    H <- create.H(obs.hole, n=n, hole.size=2*n.hole+1, ...)
  } else {
    H <- NA
    print("not a valid obs type, choose full, odd or hole")
  }
  R <- 0.5 * diag(nrow(H))
  return(list(R=R, H=H, obs.ind=obs.ind))
}



## deprecated:
## instead use: create.H(obs.operator)
obs.setup <- function(type, sig=0.5){
  ## create d, H, and R in the global environment
  if (type == "full"){
    d <<- n
    obs.ind <<- 1:n
    H <<- diag(n)
    #####################################
  } else if (type == "odd"){
    d <- 20
    obs.ind <- seq(1, n, by=2)
    H <- matrix(0, nrow=d, ncol=n)
    for (j in 1:d){
      H[j,obs.ind[j]]  <- 1
    }

    H <<- H
    d <<- d
    obs.ind <<- obs.ind
    #####################################
  } else if (type == "hole"){
    obs.ind <- 1:n
    n.hole <- 3
    obs.ind <- obs.ind[!1:n %in% (n/2 - n.hole):(n/2 + n.hole)]
    d <- length(obs.ind)

    H <- matrix(0, nrow=d, ncol=n)
    for (j in 1:d){
      H[j,obs.ind[j]]  <- 1
    }

    H <<- H
    obs.ind <<- obs.ind
    d <<- d
    #####################################
  } else {
    print("not a valid obs type, choose full, odd or hole")
  }

  R <<- sig * diag(d)
  print("d, H, R and obs.ind set in global environment")
  return(type)
}



generate.xy <- function(tot.time, freq, myseed=1){
  source("dynamical_models/fL96.R")
  set.seed(myseed)
  x0 <- L96.integrate(rnorm(n), 2, deltat)
  ## Generate the process and the observations:
  x.true <<- fL96.integrate.ts(x0, tot.time, freq, deltat)
  y.true <<- allobs.H(x.true, H, R)
  print("x.true and y.true set in global environment")
  return(1)
}


generate.xy.o <- function(tot.time, freq, myseed=1, H=H, R=R, deltat=deltat){
  source("dynamical_models/fL96.R")
  set.seed(myseed)
  x0 <- L96.integrate(rnorm(n), 2, deltat)
  ## Generate the process and the observations:
  x.true <- fL96.integrate.ts(x0, tot.time, freq, deltat)
  y.true <- allobs.H(x.true, H, R)
  print("x.true and y.true set in global environment")
  return(list(x.true=x.true, y.true=y.true))
}


init.ens <- function(K, myseed=1){
  set.seed(myseed)
  t(replicate(K, fL96.integrate(rnorm(n), 1, deltat)))
}

