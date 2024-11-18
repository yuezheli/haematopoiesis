# this is a file testing whether changing in conditioning intensity alters steady state

# this use all the parameters from Ribeil et al., 2017

rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

modSCD <- mread("fullmodel2") 

preconditioning_scan <- function(conditioningintensity = 0.9, exvivoCD34 = 2.8 * 1e8, mod = modSCD)
{
  # set the initial condition of the system preconditioning
  # the initial values are taken from SCD subject that runs to a steady state
  
  sim0 <- mod %>% param(totalCD34infused = 0) %>% 
    init(LT0 = 1) %>% param(ksynalpha = 6e-7) %>%
    mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble() 
  
  data = tail(sim0, n = 1)
  
  Data <- data[,-(1:2),drop=FALSE]  
  
  # steady state of sickle cell patients
  mod2 <- mod %>% init(Data)
  
  # myeloablative preconditioning
  ratioleft = 1-conditioningintensity
  
  mod3 <- mod2 %>% init(LT0 = Data$LT0 * ratioleft) %>% 
    init(ST01 = Data$ST01 * ratioleft) %>% init(ST02 = Data$ST02 * ratioleft) %>%
    init(MPP01 = Data$MPP01 * ratioleft) %>% init(MPP02 = Data$MPP02 * ratioleft) %>%
    init(CMP01 = Data$CMP01 * ratioleft) %>% init(CMP02 = Data$CMP02 * ratioleft) 
  
  sim <- mod3 %>% param(totalCD34infused = exvivoCD34) %>%
    param(ksynalpha = 6e-7) %>% mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble() %>% 
    select(RBC0, RBC1, RBC2, RBC3, RBC4, RET0, RET1, RET2, RET3, RET4, LT0, LT1, totalHb, totalHbA, totalHbS, totalHbF) 
  
  steadystate <- tail(sim, n = 1) %>% as.data.frame()
  
  steadystate['conditioning_intensity'] = conditioningintensity
  
  return(steadystate)
  
}

## keep the infusion the same, but the condition level different

init_cond <- c(1:19)/20

conditioningscan = preconditioning_scan(0) # this is assuming no conditioning exist

for(i in 1:length(init_cond))
{
  tmp = preconditioning_scan(init_cond[i])
  
  conditioningscan = rbind(conditioningscan, tmp)
  
  rm(tmp)
}

# post-processing 
conditioningscan['transducedRBC'] = conditioningscan['RBC1'] + conditioningscan['RBC2'] + conditioningscan['RBC3']  + conditioningscan['RBC4']  
conditioningscan['transducedRET'] = conditioningscan['RET1'] + conditioningscan['RET2'] + conditioningscan['RET3']  + conditioningscan['RET4']


RBCcount <- ggplot(data = conditioningscan) + 
  geom_point(aes(x = conditioning_intensity, y = RBC0/ (5e6), color = 'endogenous RBC')) + 
  geom_point(aes(x = conditioning_intensity, y = transducedRBC/ (5e6), color = 'transduced RBC')) + 
  geom_point(aes(x = conditioning_intensity, y = (RBC0 + transducedRBC) / (5e6), color = 'total RBC')) + 
  labs(y = 'RBC count (#/uL)', x = 'conditioning intensity', color = '') + 
  scale_y_continuous(limits = c(1, 4.5e6), breaks = c(1e6, 2e6, 3e6, 4e6), 
                     labels = c('1M', '2M', '3M', '4M'))  + theme_bw() + theme(legend.position = "bottom")

RETcount <- ggplot(data = conditioningscan) + 
  geom_point(aes(x = conditioning_intensity, y = RET0/ (5e6), color = 'endogenous RET')) + 
  geom_point(aes(x = conditioning_intensity, y = transducedRET/ (5e6), color = 'transduced RET')) + 
  geom_point(aes(x = conditioning_intensity, y = (RET0 + transducedRET) / (5e6), color = 'total RET')) + 
  labs(y = 'RET count (#/uL)', x = 'conditioning intensity', color = '') + theme_bw() + theme(legend.position = "bottom")


Hbconc <- ggplot(data = conditioningscan) + 
  geom_point(aes(x = conditioning_intensity, y = totalHb, color = 'Hb')) + 
  geom_point(aes(x = conditioning_intensity, y = totalHbS, color = 'HbS')) + 
  geom_point(aes(x = conditioning_intensity, y = totalHbA, color = 'HbA')) +
  geom_point(aes(x = conditioning_intensity, y = totalHbF, color = 'HbF')) + 
  labs(y = 'Hb conc (g/dL)', x = 'conditioning intensity', color = '') + theme_bw() + theme(legend.position = "bottom")


grid.arrange(RBCcount, RETcount, Hbconc, ncol = 3)

## dynamics of RBC and Hb recovery,  different conditioning

preconditioning_dynamics <- function(conditioningintensity = 0.9, exvivoCD34 = 2.8 * 1e8, mod = modSCD)
{
  # set the initial condition of the system preconditioning
  # the initial values are taken from SCD subject that runs to a steady state
  
  sim0 <- mod %>% param(totalCD34infused = 0) %>% 
    init(LT0 = 1) %>% param(ksynalpha = 6e-7) %>%
    mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble() 
  
  data = tail(sim0, n = 1)
  
  Data <- data[,-(1:2),drop=FALSE]  
  
  # steady state of sickle cell patients
  mod2 <- mod %>% init(Data)
  
  # myeloablative preconditioning
  ratioleft = 1-conditioningintensity
  
  mod3 <- mod2 %>% init(LT0 = Data$LT0 * ratioleft) %>% 
    init(ST01 = Data$ST01 * ratioleft) %>% init(ST02 = Data$ST02 * ratioleft) %>%
    init(MPP01 = Data$MPP01 * ratioleft) %>% init(MPP02 = Data$MPP02 * ratioleft) %>%
    init(CMP01 = Data$CMP01 * ratioleft) %>% init(CMP02 = Data$CMP02 * ratioleft) 
  
  sim <- mod3 %>% param(totalCD34infused = exvivoCD34) %>%
    param(ksynalpha = 6e-7) %>% mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble() %>% 
    select(time, RBCconc, totalHb, totalHbA, totalHbS) %>% as.data.frame()
  
  return(sim)
}

conditioningintensity_point9 <- preconditioning_dynamics(0.9)
conditioningintensity_point5 <- preconditioning_dynamics(0.5)
conditioningintensity_point1 <- preconditioning_dynamics(0.1)

# plot RBC & Hb recovery

RBCrecoverydynamics = ggplot() + 
  geom_line(data = conditioningintensity_point9, aes(time/30, RBCconc, color = '90%')) + 
  geom_line(data = conditioningintensity_point5, aes(time/30, RBCconc, color = '50%')) + 
  geom_line(data = conditioningintensity_point1, aes(time/30, RBCconc, color = '10%')) + 
  labs(y = 'RBC count (#/uL)', x = 'time (month)', color = 'preconditioning intensity') +
  scale_y_continuous(limits = c(1e6, 4e6),
                     breaks = c(1e6, 2e6, 3e6, 4e6), 
                     labels = c('1M', '2M', '3M', '4M')) + 
  theme_bw() + theme(legend.position = "bottom")


Hbrecoverydynamics = ggplot() + 
  geom_line(data = conditioningintensity_point9, aes(time/30, totalHb, color = '90%')) + 
  geom_line(data = conditioningintensity_point5, aes(time/30, totalHb, color = '50%')) + 
  geom_line(data = conditioningintensity_point1, aes(time/30, totalHb, color = '10%')) + 
  labs(y = 'Hb concentration (g/dL)', x = 'time (month)', color = 'preconditioning intensity') + theme_bw() + theme(legend.position = "bottom")

grid.arrange(RBCrecoverydynamics, Hbrecoverydynamics, ncol = 2)
