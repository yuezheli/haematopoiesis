rm(list = ls())

library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(nloptr)
library(numDeriv)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


# observed data from Shahaf et al., 2016
alldata = read.csv('data/BrdUfraction_control_obs.csv', header = TRUE) %>% mutate(dvtype = 0) %>% 
          mutate(ID = 1) %>% mutate(DV = BrdUfrac) %>% mutate(time = time_days * 4) %>% arrange(time)

for(i in 1:length(alldata[,1]))
{
    if (alldata$Bcelltype[i] == 'immatureB') {alldata$dvtype[i] = 1}
    if (alldata$Bcelltype[i] == 'maturecircul') {alldata$dvtype[i] = 2}
    if (alldata$Bcelltype[i] == 'transitionalB') {alldata$dvtype[i] = 3}
    if (alldata$Bcelltype[i] == 'splenicmatureB') {alldata$dvtype[i] = 4}
}

alldata = alldata %>% select(c(ID, time, dvtype, DV))

sim_pre <- mread("Bcell_mus") %>% mrgsim(end = 600) %>% as_tibble() %>% tail(n=1)

yobs <- alldata[["DV"]]

## set up objective function
OF <- function(pars, obs = alldata, dvcol = "DV", pred = F)
{
    pars <- lapply(pars,exp)  #Get out of log domain for MLE
    pars <- as.list(pars)
    names(pars) <- names(theta)
    
    out <- mread("Bcell_mus3") %>% init(UBoe = sim_pre$Boe) %>% init(UBi = sim_pre$Bi) %>% 
        init(UBt = sim_pre$Bt) %>% init(UBMspl = sim_pre$BMspl) %>% init(UBMrec = sim_pre$BMrec) %>%
        #param(pars) %>% mrgsim_d(end = -1, data = obs, output="df")
        param(pars) %>% mrgsim_d(end = -1, data = distinct(obs, time, .keep_all = TRUE), output="df") %>%
        #select(time, BrdUimmature, BrdUmaturerec, BrdUtransit, BrdUsplenicmature) 
        select(time, DV)

    #out3 = cbind(out2[1], stack(out2[2:5]))
         
    
    if(pred){ return(out)}
    else{
        #y_hat <- out3$values
        #llike <- dnorm(yobs, y_hat, pars$sigma, log = TRUE)
        llike <- dnorm(yobs, out[[dvcol]], pars$sigma, log = TRUE)
        return(-1*sum(llike, na.rm=TRUE))
    }
    
}

theta <- log(c(gamma= 0.3, delta_oe =0.5, mu_t = 0.03, delta_t = 0.03, phi_BM = 0.94, phi_s = 0.03, mu_re = 0.008, epsilon_spl = 0.008, sigma=1))


fit <- newuoa(theta, OF) 

p <- as.list(exp(fit$par))  #get the parameters on the linear scale
names(p) <- names(theta)
p

