## Assimilation example with Lorenz96

library(assimilr)



## Generate the data:
set.seed(1)
ndim <- 40
Forcing <- 10
deltat <- 1/1000
freq <- 0.1
duration <- 200*freq
l96_run <- l96_simulate(duration, freq, ndim=ndim, Forcing=Forcing, deltat=deltat)


## Generate initial ensemble of size K:
K <- 100
ens0 <- l96_ens0(K, l96_run)

## Plot state and ensemble at time tt:
tt <- 1
lorenz_plot(l96_run$state.ts[tt,], l96_run$y.ts[[tt]], ens0)



#### EnKF #######
l <- 10
taper <- ring_GC(ndim, l)
enkf_run <- da_cycle(ens0, l96_run, EnKPF, gam.fix=1, taper=taper)
tt <- duration/freq +1
lorenz_plot(l96_run$state.ts[tt,], l96_run$y.ts[[tt]], enkf_run$ensA[[tt]], tit='EnKF')



#### PF #######
pf_run <- da_cycle(ens0, l96_run, EnKPF,gam.fix=0)
lorenz_plot(l96_run$state.ts[tt,], l96_run$y.ts[[tt]], pf_run$ensA[[tt]], tit='PF')



#### EnKPF #######
enkpf_run <- da_cycle(ens0, l96_run, EnKPF,e.0=0.5, e.1=0.5, taper=taper)
lorenz_plot(l96_run$state.ts[tt,], l96_run$y.ts[[tt]], enkpf_run$ensA[[tt]], tit='EnKPF')



#### Local methods ##############################

#### LEnKF #######
lenkf_run <- da_cycle(ens0, l96_run, naive_LEnKPF, l=l,gam.fix=1, taper=taper,
                      ndim=ndim, get_neighbours=ring_neighbours)
lorenz_plot(l96_run$state.ts[tt,], l96_run$y.ts[[tt]], lenkf_run$ensA[[tt]], tit='LEnKF')


#### naive-LEnKPF #######
naive_run <- da_cycle(ens0, l96_run, naive_LEnKPF, l=l,e.0=0.5, e.1=0.5, taper=taper,
                      ndim=ndim, get_neighbours=ring_neighbours)
lorenz_plot(l96_run$state.ts[tt,], l96_run$y.ts[[tt]], naive_run$ensA[[tt]], tit='naive-LEnKPF')


#### block-LEnKPF #######
block_run <- da_cycle(ens0, l96_run, block_LEnKPF, l=l, e.0=0.5, e.1=0.5, taper=taper, block_size=l/2,
                      ndim=ndim, get_partition=ring_partition)
lorenz_plot(l96_run$state.ts[tt,], l96_run$y.ts[[tt]], block_run$ensA[[tt]], tit='block-LEnKPF')




