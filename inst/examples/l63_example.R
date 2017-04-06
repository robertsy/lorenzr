


set.seed(1)
kens <- 10
params <- c(sigma=10, rho=28, beta=8/3)
x0 <- rnorm(3)
deltat <- 1/75
nsteps <- 10/deltat

state <- x0
state1 <- l63_integrate(x0, params=params,nsteps=nsteps, deltat=deltat)
plot.ts(state1)


## df format
state_df <- lorenz_as_df(state1)
ens_df <- lorenz_as_df(ens0)


## simulate
set.seed(1)
deltat <- 1/75
freq <- 1/5
duration <- 100*freq
l63_run <- l63_simulate(duration, freq, deltat=deltat, sig = 3)

# lorenz_heatmap(l63_run)


## ensemble:
ens0 <- l63_ens0(kens, l63_run)

l63_df <- lorenz_as_df(l63_run$state.ts) %>%
  rename(time=x, variable=ensemble)
l63_df$variable <- factor(l63_df$variable, levels = c('ens_1', 'ens_2', 'ens_3'), labels=c('x(1)', 'x(2)', 'x(3)'))


l63_y_df <-
  foreach(i=1:length(l63_run$y.ts), .combine='rbind')%do%{
    yy <- l63_run$y.ts[[i]]
    data_frame(variable=paste('x(',yy$y.loc,')',sep=''), value=yy$y, time=i)
  }

l63_df %>%
  ggplot(aes(x=time, y=value)) + geom_line() +
  geom_point(data=l63_y_df, aes(x=time, y=value), color='red')+
  facet_wrap(~variable, ncol=1) +
  theme_bw()

