## run basic Lorenz 96 model:

set.seed(1)
deltat <- 1/1000
freq <- 0.4
duration <- 200*freq
ndim <- 40
Forcing <- 8
l96_run <- l96_simulate(0.1*100, 0.1, ndim=40, Forcing=8, deltat=1/1000)
lorenz_plot(l96_run$state.ts[101,], l96_run$y.ts[[101]])

## heatmap:
lorenz_heatmap(l96_run)


## between two time steps:
i <- 4
plot.ts(l96_run$state.ts[i,])
plot.ts(l96_run$state.ts[i+1,])



K <- 20
ens0 <- l96_ens0(K, l96_run)


lorenz_plot(l96_run$state.ts[1,], l96_run$y.ts[[1]], ens0)





