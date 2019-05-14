###        ODE in STAN           ###
#     12.11.2018-14.05.2019        #    
#           M.J.Counotte           #
####################################

rm(list=ls())

library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')
##############
setwd("data")
nicaI<-read.csv(file="nicaI.csv", header=TRUE)

inits <-c(290497,709503) # total number Ageclass1,
params<-c("beta1", "beta2", "gamma")#,"gamma")
N<-c(3740,1074)#,230,250) #sample size per age group
positives<-c(1347,601)#,148,150)
stan_d = list(n_obs = length(nicaI$new),
              n_params = length(params),
              n_difeq = 7,
              y = nicaI$new,
              initialsize = inits,
              ts = 1:100,
              tmax =100,
              agegroups = length(N),
              t0 = 0,
              N=N, # at time 100 sampled
              positives=positives,
              inference=1,
              prior_beta1=5,
              prior_beta2=3.3,
              prior_gamma=c(5,1),
              prior_rho=c(1,400),
              prior_init_inf=0.01
)

initf2 <- function() { # initial values for model with 4 betas, gamma, rho, sigma [add init_inf prior]
  list(params=c(rexp(1,stan_d$prior_beta1),
                rexp(1,stan_d$prior_beta2),
                rgamma(1,stan_d$prior_gamma[1],stan_d$prior_gamma[2])),
      rho=rbeta(1,stan_d$prior_rho[1],stan_d$prior_rho[2]), 
      init_inf=rexp(1,stan_d$prior_init_inf))
}

# Test / debug the model:
test1 = stan("ODE_2classes_inits.stan",
             data = stan_d, init=initf2,
             chains = 1, iter = 10)

# Prior predictive check
stan_d$inference=0

prior_pc=stan("ODE_2classes_inits.stan", init=initf2,
              data = stan_d, chains = 1, iter = 1000)
print(prior_pc,pars=c("params","rho","init_inf"))

summary(prior_pc,pars=c("pred_Y"))[1] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  add_rownames() %>%
  mutate(date=1:38) %>%
  ggplot(.) +
  geom_ribbon(aes(x=date,ymax=summary.97.5.,ymin=summary.2.5.),alpha=0.5) +
  geom_line(aes(x=date,y=summary.mean))+
  geom_point(data=nicaI,aes(x=wk,y=stan_d$y, color="Observed counts", fill="Observed counts"))


summary(prior_pc,pars=c("pred_p"))[1] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  add_rownames() %>%
  mutate(age=1:2) %>%
  ggplot(.) +
  geom_point(aes(x=factor(age),y=summary.mean/N)) +
  geom_errorbar(aes(x=c(1:2), ymin=summary.2.5./N, ymax=summary.97.5./N), width=0.5)+
  geom_point(aes(x=c(1:2),y=positives/N, color=factor(c(1:2))))+
  geom_errorbar(aes(x=c(1:2), ymin=c(0.34,0.50), ymax=c(0.38,0.58), color=factor(c(1:2)), width=0.5))


stan_d$inference=1
mod_2classes = stan("ODE_2classes_inits.stan", init=initf2,
                   data = stan_d, chains = 4, iter = 2000)

print(mod_2classes,pars=c("params","rho","init_inf"))

p1<-summary(mod_2classes,pars=c("pred_Y"))[1] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  add_rownames() %>%
  mutate(date=1:38) %>%
  ggplot(.) +
  geom_ribbon(aes(x=date,ymax=summary.97.5.,ymin=summary.2.5.),alpha=0.5) +
  geom_line(aes(x=date,y=summary.mean)) +
  geom_point(data=nicaI,aes(x=wk,y=new, color="Observed counts"))

p2<-summary(mod_2classes,pars=c("pred_p"))[1] %>%
  as.data.frame(.) %>%
  tbl_df() %>%
  add_rownames() %>%
  mutate(age=1:2) %>%
  ggplot(.) +
  geom_point(aes(x=age,y=summary.mean/N)) +
  geom_errorbar(aes(x=age, ymin=summary.2.5./N, ymax=summary.97.5./N), width=0.5)+
  geom_point(aes(x=c(1:2),y=positives/N, color=factor(c(1:2))))
  #geom_errorbar(aes(x=c(1:2), ymin=c(0.34,0.51), ymax=c(0.38,0.58), color=factor(c(1:2)), width=0.5))

p3<-stan_trace(mod_2classes,pars=c("params","rho","init_inf"),inc_warmup = TRUE)


library(gridExtra)
pdf("../../writing/figures/figure_fit2ageclasses.pdf", height=10, width=10)
grid.arrange(p1, p2, p3, nrow = 2)
dev.off()
