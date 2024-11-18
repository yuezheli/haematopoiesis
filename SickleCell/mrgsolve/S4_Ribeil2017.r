# the goal of this file is to simulate the result of Figure S4 in the appendix. 

rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

mod <- mread("fullmodel2") 

# set the initial condition of the system preconditioning
# assuming 90% endogenous progenitor pool is lost during myeloablative preconditioning
# the initial values are taken from SCD subject that runs to a steady state

sim0 <- mod %>% param(totalCD34infused = 0) %>% 
  init(LT0 = 1) %>% param(ksynalpha = 6e-7) %>%
  mrgsim(end = 30 * 18) %>% as_tibble() 

data = tail(sim0, n = 1)

Data <- data[,-(1:2),drop=FALSE]  

# steady state of sickle cell patients
mod2 <- mod %>% init(Data)

# myeloablative preconditioning
ratioleft = 1-0.9

mod3 <- mod2 %>% init(LT0 = Data$LT0 * ratioleft) %>% 
  init(ST01 = Data$ST01 * ratioleft) %>% init(ST02 = Data$ST02 * ratioleft) %>%
  init(MPP01 = Data$MPP01 * ratioleft) %>% init(MPP02 = Data$MPP02 * ratioleft) %>%
  init(CMP01 = Data$CMP01 * ratioleft) %>% init(CMP02 = Data$CMP02 * ratioleft) 


# simulation after the patient is given 2.8e8 (5.6e6 per kg)
# change the LT-HSC prliferation by assuming only 10% of LT-HSC are actively proliferating
sim <- mod3 %>% # param(rLT = 1/(7*9)) %>% 
  param(ksynalpha = 6e-7) %>% mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble() 

# read in the observed data from Ribeil et al., 2017
RBC_ribeil <- read.csv('../data/RBC_Ribeil.csv', header = T)
Hb_ribeil <- read.csv('../data/Hb_Ribeil.csv', header = T)

hb_hbs <- Hb_ribeil[Hb_ribeil['Hb'] == 'HbS',]
hb_hbaT87Q <- Hb_ribeil[Hb_ribeil['Hb'] == 'HbA_T87Q',]


# plot what is available in Figure S4
RBCnum = ggplot() + 
  geom_line(data = sim, aes(x = time/30, y = RBCconc, color = 'total RBC')) +
  geom_line(data = sim, aes(x = time/30, y = RBC0/5e6, color = 'endogenous RBC')) +
  geom_line(data = sim, aes(x = time/30, y = (RBC1 + RBC2 + RBC3 + RBC4)/5e6, color = 'transduced RBC')) +
  geom_point(data = RBC_ribeil, aes(x = month, y = RBC_in_millions * 1e6, color = 'Ribeil et al., 2017')) + 
  labs(y = 'RBC count (#/uL)', x = 'time (months)', color = '') + 
  scale_x_continuous(breaks = c(0,3,6,9,12,15)) + 
  scale_y_continuous(limits = c(1, 4.5e6), breaks = c(1e6, 2e6, 3e6, 4e6), 
                     labels = c('1M', '2M', '3M', '4M')) + theme(legend.position = "bottom") + theme_bw()

HbAconc = ggplot() + 
  geom_line(data = sim, aes(x = time/30, y = totalHbS, color = 'total HbS')) + 
  geom_line(data = sim, aes(x = time/30, y = totalHbA, color = 'total HbA_T87Q')) + 
  geom_point(data = hb_hbs, aes(x = month, y = Hbconc, color = 'HbS, Ribeil et al., 2017')) + 
  geom_point(data = hb_hbaT87Q, aes(x = month, y = Hbconc, color = 'HbA_T87Q, Ribeil et al., 2017')) + 
  scale_x_continuous(breaks = c(0,3,6,9,12,15)) + 
  labs(y = 'hemoglobin conc (g/dL)', x = 'time (months)', color = '') + theme(legend.position = "bottom") + theme_bw()

# grid.arrange(RBCnum, HbAconc, ncol = 1)


# RET data reported in Ribeil et al., 2017
ret_ribeil <- data.frame(c(3,6,9,12,15), c(259000, 132000, 131000, 143000, 141000)) # unit in uL-1
names(ret_ribeil) <- c('month', 'RET')

RETnum = ggplot() + 
  geom_line(data = sim, aes(x = time/30, y = RET/(5e6), color = 'total RET')) + 
  geom_line(data = sim, aes(x = time/30, y = RET0/(5e6), color = 'endogenous RET')) + 
  geom_line(data = sim, aes(x = time/30, y = (RET1 + RET2 + RET3 + RET4)/(5e6), color = 'transduced RET')) + 
  geom_point(data = ret_ribeil, aes(x = month, y = RET, color = 'RET, Ribeil et al., 2017')) + 
  scale_x_continuous(breaks = c(0,3,6,9,12,15)) + 
  labs(y = 'RET count (#/uL)', x = 'time (months)', color = '') + theme(legend.position = "bottom") + theme_bw()


# check ratio between endogenous and transduced HSC
HSCnum = ggplot(data = sim) + 
  geom_line(aes(x = time/30, y = LT0, color = 'endogenous LT-HSC')) + 
  geom_line(aes(x = time/30, y = LT1, color = 'transduced LT-HSC')) + 
  scale_x_continuous(breaks = c(0,3,6,9,12,15)) + 
  labs(y = 'LT-HSC count (#)', x = 'time (months)', color = '') + theme(legend.position = "bottom") + theme_bw()

# Hb conc in RBC data reported in Ribeil et al., 2017
hbinrbc_ribeil <- data.frame(c(3,6,9,12,15), c(34, 35, 36, 35, 35)) # unit in g.dL-1
names(hbinrbc_ribeil) <- c('month', 'HbinRBC')

HbinRBCconc = ggplot() + 
  geom_line(data = sim, aes(x = time/30, y = HbinRBC/10, color = 'Hb conc in RBC')) + 
  geom_point(data = hbinrbc_ribeil, aes(x = month, y = HbinRBC, color = 'Hb conc in RBC, Ribeil et al., 2017')) + 
  scale_x_continuous(breaks = c(0,3,6,9,12,15)) + 
  labs(y = 'Hb in RBC (g/dL)', x = 'time (months)', color = '') + theme(legend.position = "bottom") + theme_bw() 


grid.arrange(RBCnum, HbAconc, RETnum, HbinRBCconc, ncol = 1)
