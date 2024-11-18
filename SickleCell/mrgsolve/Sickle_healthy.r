# the goal of this file is to create simulation in the condition similar to Figure 4 in Zheng et al., 2021
# https://ascpt.onlinelibrary.wiley.com/doi/full/10.1002/psp4.12638

rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

# read in the model
# mod <- mread("fullmodel") 
mod <- mread("fullmodel2") 

##------------------- sickle cell anemia simulation -------------------##

sim0 <- mod %>% param(totalCD34infused = 0) %>% 
  init(LT0 = 1) %>% param(ksynalpha = 6e-7) %>%
  mrgsim(end = 30 * 18) %>% as_tibble() 

notreatmentbloodcell = ggplot() + 
  geom_line(data = sim0, aes(x = time/30, y = RBC, color = 'RBC')) + 
  geom_line(data = sim0, aes(x = time/30, y = RET, color = 'RET')) + 
  geom_line(data = sim0, aes(x = time/30, y = LTHSC, color = 'LTHSC')) + 
  geom_line(data = sim0, aes(x = time/30, y = STHSC, color = 'STHSC')) + 
  geom_line(data = sim0, aes(x = time/30, y = MPP, color = 'MPP')) + 
  geom_line(data = sim0, aes(x = time/30, y = CMP, color = 'CMP')) + 
  geom_line(data = sim0, aes(x = time/30, y = BFUE, color = 'BFU-E')) + 
  geom_line(data = sim0, aes(x = time/30, y = CFUE, color = 'CFU-E')) + 
  labs(y = 'cell count', x = 'time (months)', color = '') + 
  scale_y_continuous(trans='log10', limits = c(1, 1e14), 
                     breaks = c(1, 1e3, 1e6, 1e9, 1e12), 
                     labels = c('1', '1K', '1M', '1B', '1T')) + theme(legend.position = "bottom") + 
  ggtitle('sickle cell patient')


notreatmenthemoglobin = ggplot() + 
  geom_line(data = sim0, aes(x = time/30, y = totalHb, color = 'total Hb')) + 
  geom_line(data = sim0, aes(x = time/30, y = totalHbA, color = 'HbA')) + 
  geom_line(data = sim0, aes(x = time/30, y = totalHbS, color = 'HbS')) + 
  geom_line(data = sim0, aes(x = time/30, y = totalHbF, color = 'HbF')) + 
  geom_line(data = sim0, aes(x = time/30, y = totalHbA2, color = 'HbA2')) + 
  labs(y = 'hemoglobin concentration (g/dL)', x = 'time (months)', color = 'SCD') 

notreatmentO2 = ggplot() + 
  geom_line(data = sim0, aes(x = time/30, y = bloodO2, color = 'vO2')) + 
  labs(y = 'blood oxygen (ml/dL)', x = 'time (months)', color = 'SCD') 

##------------------- healthy people simulation -------------------##

sim1 <- mod %>% param(totalCD34infused = 0) %>% 
  param(Kdalphabeta_sickle = 1e-3) %>% # change the dissociation from sickle Hb to normal HbA
  init(LT0 = 1) %>%
  param(tauRBCsickle = 120) %>%  # this parameter change simulate healthy people situation
  param(ksynalpha = 6e-7) %>%  #decrease hemoglobin synthesis rate
  param(HbS_saturation = 0.74) %>% # change HbS to HbA parameters
  mrgsim(end = 30 * 18) %>% as_tibble() 


healthybloodcell = ggplot() + 
  geom_line(data = sim1, aes(x = time/30, y = RBC, color = 'RBC')) + 
  geom_line(data = sim1, aes(x = time/30, y = RET, color = 'RET')) + 
  geom_line(data = sim1, aes(x = time/30, y = LTHSC, color = 'LTHSC')) + 
  geom_line(data = sim1, aes(x = time/30, y = STHSC, color = 'STHSC')) + 
  geom_line(data = sim1, aes(x = time/30, y = MPP, color = 'MPP')) + 
  geom_line(data = sim1, aes(x = time/30, y = CMP, color = 'CMP')) + 
  geom_line(data = sim1, aes(x = time/30, y = BFUE, color = 'BFU-E')) + 
  geom_line(data = sim1, aes(x = time/30, y = CFUE, color = 'CFU-E')) + 
  labs(y = 'cell count', x = 'time (months)', color = '') + 
  scale_y_continuous(trans='log10', limits = c(1, 1e14), 
                     breaks = c(1, 1e3, 1e6, 1e9, 1e12), 
                     labels = c('1', '1K', '1M', '1B', '1T')) + theme(legend.position = "bottom") + 
  ggtitle('healthy subject')


healthyhemoglobin = ggplot() + 
  geom_line(data = sim1, aes(x = time/30, y = totalHb, color = 'total Hb')) + 
#  geom_line(data = sim1, aes(x = time/30, y = totalHbA, color = 'HbA')) + 
  geom_line(data = sim1, aes(x = time/30, y = totalHbS, color = 'HbA')) + 
  geom_line(data = sim1, aes(x = time/30, y = totalHbF, color = 'HbF')) + 
  geom_line(data = sim1, aes(x = time/30, y = totalHbA2, color = 'HbA2')) + 
  labs(y = 'hemoglobin concentration (g/dL)', x = 'time (months)', color = 'healthy') 

healthyO2 = ggplot() + 
  geom_line(data = sim1, aes(x = time/30, y = bloodO2, color = 'vO2')) + 
  labs(y = 'blood oxygen (ml/dL)', x = 'time (months)', color = 'healthy') 


##-------------------- display data --------------------##

grid.arrange(notreatmentbloodcell, healthybloodcell,
             notreatmenthemoglobin, healthyhemoglobin, 
             notreatmentO2, healthyO2,
             ncol = 2)

print(sim0[5401, 183:199])
print(sim1[5401, 183:199])

grid.arrange(notreatmentbloodcell, notreatmenthemoglobin, notreatmentO2, 
             healthybloodcell, healthyhemoglobin, healthyO2,
             ncol = 3)


##-------------------- save data --------------------##
sickle_RBC <- sim0 %>% select(time, LT0, ST01, ST02) %>% filter(time < 24)
write_csv(sickle_RBC, file = '../data/sickle_simul.csv')
