

## integrate nsteps in time:
set.seed(1)

ndim <- 960
deltat <- 1/200
state <- runif(ndim, -1, 1)
nsteps <-1000

kint <- 32

state1 <- lm2_integrate(state, nsteps, deltat=deltat, kint=kint, Forcing=15)
plot.ts(state1)

## ensemble ( to test)
ens0 <- cbind(state, state) #cbind( runif(ndim, -1, 1), runif(ndim, -1, 1))
ens1 <-    lm2_integrate(ens0, nsteps, deltat=deltat, kint=kint, Forcing=15)



### Cycles:
freq <- .2
duration <- 10*freq
lm2_run <- lm2_simulate(duration, freq, ndim=ndim, kint=kint, Forcing=15, deltat=deltat)

## between two time steps:
i <- 1
plot.ts(lm2_run$state.ts[i,])
plot.ts(lm2_run$state.ts[i+1,])



K <- 20
ens0 <- lm2_ens0(K, lm2_run)

lorenz_plot(lm2_run$state.ts[1,], lm2_run$y.ts[[1]], ens0)

lorenz_heatmap(lm2_run)

