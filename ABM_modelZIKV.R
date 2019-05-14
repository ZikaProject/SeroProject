###        ABM in STAN           ###
#     20.12.2018-14.05.2019        #    
#           M.J.Counotte           #
####################################



rm(list=ls())
library(deSolve)
library(dplyr)
library(ggplot2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

##############
setwd("data")



propstart <-read.csv(file="popprop5.csv", header=TRUE) #proportion per age group.
ps_mat=cbind(age=rep(0:99),prop1=rep(propstart$prop5/5,rep(5,20)))

dr<-read.csv(file="death_rate.csv",header=TRUE)
dr0<-dr$rate_yr[1] # rate for age <1
dr<-dr[2:19,]

model <- lm(log(dr$rate_yr+0.0001) ~ poly(dr$age,3, raw=TRUE))

a<-model$coefficients[1]
b<-model$coefficients[2]
c<-model$coefficients[3]
d<-model$coefficients[4]

mort_rates<-unname(c(dr0,a,b,c,d)) # pass to data in STAN: [1]: <1 yr mortality, [2:5] polynomial coefficients


N0<-10000 #initial population size
init_pop<-function(n){ # return initial population of size N0 with distribution according to popprop5.csv
  tmp<-round(ps_mat[,2]/100*n,0)
  # correct for rounding:
  if(sum(tmp)!=n){
    tmp[1]<-tmp[1]+n-sum(tmp)
  }
  return(rep(1:100,tmp))
}
pop<-init_pop(N0) #inital ages for N0 individuals distributed according to ps_mat


#pop<-rep(1:100,(N0/100)) #inital ages uniform distributed.


initf <- function() { # initial values for model with 2 betas, gamma: 
  list(params=c(runif(1,0.13,0.21),runif(1,0.28,0.36),runif(1,4,6)))
}




load(file = "parameters.Rda")
# get arguments
runN<-as.numeric(commandArgs(trailingOnly=TRUE))
#runN<-1
#scenarios<-expand.grid(n=1:50, introduction_n=c(1,5,10))
scenarios<-expand.grid(n=1:80, introduction_n=c(1,5,10))

introduction_wk<-52*scenarios$n[runN]
introduction_n<-scenarios$introduction_n[runN]

wks<-introduction_wk+1820
stan_d = list(N0=N0, # 'active' (alive) individuals in population
              nend=N0*3, # maximum allowed individuals
              units=wks, #allow halve a year for the outbreak.
              age_init=c(pop, rep(0,N0*2)), # initial age for nend 
              birthrate=2.2, # life time birth rate
              mort_rates=mort_rates, # mortality rates based on fitted curve >1, fixed value<1
              introduction_wk=introduction_wk, # week at which introduction is introduced
              introduction_n=introduction_n,
              beta1=beta1,
              beta2=beta2,
              gamma=gamma,
              prop1=prop_end1,
              prop2=prop_end2)# number of infected individuals introduced

cat("Scenario:", runN, " Introduction wk:", introduction_wk)

 # model<-stan_model(file = "ABM_ZIKV.stan", model_name = "ABM Stan1", allow_undefined = TRUE,
 #           includes = paste0('\n#include "',
 #                              file.path(getwd(), 'count.hpp'), '"\n')) #count hpp keeps track of iteration number

load(file = "model.Rda")#load precompiled model
 
test1 = sampling(model,
             data = stan_d, 
             chains = 1, 
             warmup=0,
             iter = 1000,
             init=initf,
             pars = c("c_inf", "prop_c_inf", "c_inf1", "c_inf2", "c_inf3", "c_inf4", "alive_start", "S_start",
                      "S_start1","S_start2","S_start3","S_start4","R_start","R_start1","R_start2","R_start3","R_start4"))


#require(reshape2)
extr<-extract(test1, c("c_inf", "prop_c_inf", "c_inf1", "c_inf2", "c_inf3", "c_inf4", "alive_start", "S_start",
                 "S_start1","S_start2","S_start3","S_start4","R_start","R_start1","R_start2","R_start3","R_start4"))
n<-rep(introduction_n,length(extr[[1]]))
year<-rep(scenarios$n[runN],length(extr[[1]]))
df <- cbind(as.data.frame(extr),n,year)
save(df, file=paste0("results",runN,".Rda"))
#q()
