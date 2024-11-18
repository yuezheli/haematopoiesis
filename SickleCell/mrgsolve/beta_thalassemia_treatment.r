rm(list = ls())
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

mod <- mread("thalassemia") 

modB <- mod %>% param(totalCD34infused = 0) %>% param(ksynalpha = 6e-7) %>% # these are for SCD parameter change
  param(HbS_saturation = 0.74)  %>% # change the blood carrying capacity of mutated HbA (in SCD, it is HbS)
  init(LT0 = 1) %>% # set the system with only 1 stem cell
  param(tauRBCsickle  = 36) %>% # change the avg RBC life span 
  param(E_synbeta_effect = 0.5) 

sim0 <- modB %>% mrgsim(end = 30 * 18) %>% as_tibble() 

data = tail(sim0, n = 1)

Data <- data[,-(1:2),drop=FALSE]  

# steady state of sickle cell patients
mod2 <- modB %>% init(Data)

# myeloablative preconditioning
ratioleft = 1-0.9

mod3 <- mod2 %>% init(LT0 = Data$LT0 * ratioleft) %>% 
  init(ST01 = Data$ST01 * ratioleft) %>% init(ST02 = Data$ST02 * ratioleft) %>%
  init(MPP01 = Data$MPP01 * ratioleft) %>% init(MPP02 = Data$MPP02 * ratioleft) %>%
  init(CMP01 = Data$CMP01 * ratioleft) %>% init(CMP02 = Data$CMP02 * ratioleft)

## -------------------- use the default model for simulation, VCN = 0.5 -------------------- ##
# use HGB204, non-beta0-beta0 type, assume body weight = 71kg
# sim_hgb204 <- mod3 %>% param(totalCD34infused = 7.1e6 * 71) %>% mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble()
# use HGB205, non-beta0-beta0 type, assume body weight = 71kg
# sim_hgb205 <- mod3 %>% param(totalCD34infused = 12e6 * 71) %>% mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble() 

## -------------------- HGB204 VCN = 0.7, HGB205 VCN = 1.3 -------------------- ##
sim_hgb204 <- mod3 %>% param(totalCD34infused = 7.1e6 * 71) %>% 
  param( ratiosynbetanew = (0.7*2)/(1+0.7*2)  ) %>% mrgsim(delta = 0.5, end = 30 * 36) %>% as_tibble()

sim_hgb205 <- mod3 %>% param(totalCD34infused = 12e6 * 71) %>% 
  param( ratiosynbetanew = (1.3*2)/(1+1.3*2)  ) %>% mrgsim(delta = 0.5, end = 30 * 36) %>% as_tibble() 


# read in observed data from clinical trial
hgb204 <- read.csv('../data/HGB204_HbAT87Q_nonbeta0_Thompson.csv', header = TRUE )
hgb205 <- read.csv('../data/HGB205_HbAT87Q_nonbeta0_Thompson.csv', header = TRUE )


HbA_204 = ggplot() + 
  geom_line(data = sim_hgb204, aes(x = time/30, y = totalHbA, color = 'simulated HbA_T87Q')) + 
  geom_point(data = hgb204, aes(x = month, y = HbAT87Q, color = 'observed HbA_T87Q in Thompton et al., 2018')) + 
  labs(y = 'hemoglobin conc (g/dL)', x = 'time (months)', color = ' ') + theme_bw() + 
  scale_x_continuous(breaks = c(0,3,6,9,12,15, 18, 21, 24, 30, 36)) + 
  theme(legend.position = "right") + ggtitle('HGB 204')

HbA_205 = ggplot() + 
  geom_line(data = sim_hgb205, aes(x = time/30, y = totalHbA, color ='simulated HbA_T87Q')) + 
  geom_point(data = hgb205, aes(x = month, y = HbAT87Q, color = 'observed HbA_T87Q in Thompton et al., 2018')) + 
  labs(y = 'hemoglobin conc (g/dL)', x = 'time (months)', color = ' ') + theme_bw() + 
  scale_x_continuous(breaks = c(0,3,6,9,12,15, 18, 21, 24, 36)) + 
  theme(legend.position = "right")  + ggtitle('HGB 205')


grid.arrange(HbA_204, HbA_205, ncol = 1)
