## run basic Lorenz 96 model:

set.seed(1)
ndim <- 40
kens <- 1
Forcing <- 10#8
x0 <- rnorm(ndim)
deltat <- 1/40
deltat <- 1/500
nsteps <- 10/deltat

state <- x0
state1 <- l96_integrate(x0, nsteps, deltat=deltat)
# state1 <- L96.integrate(x0, nsteps, deltat=deltat)
plot.ts(state1)


## Forcing vector:
state <- x0
Forcing_vec <- Forcing+stats::filter(rnorm(ndim), 0.9)
Forcing_vec <- c(rep(5,20), rep(8,20))

state1 <- l96_integrate(x0, nsteps, deltat=deltat, forcingvec=Forcing_vec)
# state1 <- L96.integrate(x0, nsteps, deltat=deltat)
plot.ts(state1)



## ensemble:
ens0 <- cbind(state, state)
ens1 <-    l96_integrate(ens0, nsteps, deltat=deltat)
plot.ts(ens1[,1])


## df format
state_df <- lorenz_as_df(state1)
ens_df <- lorenz_as_df(ens0)


set.seed(1)
deltat <- 1/1000
freq <- 0.4
duration <- 200*freq
l96_run <- l96_simulate(duration, freq, ndim=ndim, Forcing=Forcing, deltat=deltat)
#plot.ts(l96_run$state.ts[duration/freq+1,])
lorenz_plot(l96_run$state.ts[duration/freq+1,], l96_run$y.ts[[duration/freq+1]])

## heatmap:
lorenz_heatmap(l96_run)


## between two time steps:
i <-4
plot.ts(l96_run$state.ts[i,])
plot.ts(l96_run$state.ts[i+1,])



K <- 20
ens0 <- l96_ens0(K, l96_run)


lorenz_plot(l96_run$state.ts[1,], l96_run$y.ts[[1]], ens0)





#
#   source('~/workspace/COMPLETED/DA/dynamical_models/L96.R')

# x1 <- L96.integrate(x0, nsteps*deltat, deltat)
#
#   all.equal(output$state, x1)




