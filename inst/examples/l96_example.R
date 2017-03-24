## run basic Lorenz 96 model:

set.seed(1)
deltat <- 1/1000
freq <- 0.4
duration <- 200*freq
l96_run <- l96_simulate(duration, freq, ndim=ndim, Forcing=Forcing, deltat=deltat)
lorenz_plot(l96_run$state.ts[duration/freq+1,], l96_run$y.ts[[duration/freq+1]])

## heatmap:
lorenz_heatmap(l96_run)


## between two time steps:
i <- 4
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




