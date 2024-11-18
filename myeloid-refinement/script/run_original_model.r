# date: 10/1/24
# author: Yuezhe Li 
# purpose of this code: to run the original model for ADA-SCID gene therapy simulation, as presented in the poster 

rm(list=ls())  #Clear out existing objects
gc()

library(here)
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(ggtext)

patient1 <- read.csv('data/Aiuti2009.csv', header = TRUE) 

adahomo <- mread("model/human_erythroid_lymphoid_myeloid")

# to simulate a normal person:
# mod <- mread("model/human_erythroid_lymphoid_myeloid") %>% param(epsilon_spl = 0.032, delta_dp = 1e-6)

# to simulate a normal person without transplant 
# mod <- mread("model/human_erythroid_lymphoid_myeloid_simplified") 

# pre-equlibrium
pre_sim <- adahomo  %>% init(LT0 = 1000) %>%  mrgsim(end = 365 * 7) %>% as_tibble() %>% tail(n = 1) %>% select(-c(ID, time))

# set up preconditioning
ratioleft = 1-0.9
mod2 <- adahomo %>% init(pre_sim) %>% 
  init(pre_sim[,1:22] * ratioleft) %>% # progenitor in BM conditioning
  init(pre_sim[,67:72] * ratioleft) %>% # CLP and propreB conditioning
  init(pre_sim[,73:78] * ratioleft) %>% # splenic B cell conditioning
  init(pre_sim[,79:116] * ratioleft) %>% # thymic endures similar conditioning
  init(pre_sim[,117:130] * ratioleft) # conditioning on T cell outside thymus and GMs

## set up initial condition for testing
totalcd34 = 65e6
lthscinit = 1e-6 * totalcd34; 
sthscinit = 1e-5 * totalcd34; 
mppinit = 1e-3 * totalcd34;
cmpinit = 1e-2 * totalcd34;
clpinit = 1e-4 * totalcd34;
gmpinit = 1e-2 * totalcd34;
boeinit = 5e-2 * totalcd34;
biinit = 3e-2 * totalcd34; 
sim00 <- mod2 %>% init(LT1 = lthscinit, ST12 = sthscinit, MPP12 = mppinit, CMP12 = cmpinit, GMP12 = gmpinit, 
                       CLP1 = clpinit, Boe1 = boeinit, Bi1 = biinit) %>% 
  mrgsim(end = 365 * 7) %>% as_tibble() 

cd15plot = ggplot(data = sim00) + 
  geom_line(aes(x = time/365, y = GM1/(GM0 + GM1) * 100, color = 'simul')) + 
  geom_point(data = patient1 %>% filter(celltype == "CD15") ,aes(x = time, y = cellpercentage, color = 'obs'), size = 3) + 
  scale_y_continuous(trans='log10',  limits = c(1, 100), breaks = c(1, 10, 100), labels = c('1', '10','100')) + 
  scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7), labels = c('0','1', '2', '3', '4', '5', '6', '7')) + 
  labs(y = 'vector-positive cells (%)', x = 'time (year)', color = '') + theme_bw() + theme(legend.position = c(0.4, 0.2)) + 
  ggtitle('CD15+ GM')

cd19plot = ggplot(data = sim00) + 
  geom_line(aes(x = time/365, y = transB * 100, color = 'simul')) + 
  geom_point(data = patient1 %>% filter(celltype == "CD19") ,aes(x = time, y = cellpercentage, color = 'obs'), size = 3) + 
  scale_y_continuous(trans='log10',  limits = c(1, 100), breaks = c(1, 10, 100), labels = c('1', '10','100')) + 
  scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7), labels = c('0','1', '2', '3', '4', '5', '6', '7')) + 
  labs(y = 'vector-positive cells (%)', x = 'time (year)', color = '') + theme_bw() + theme(legend.position = c(0.5, 0.2)) + 
  ggtitle('CD19+ B Cells')

cd3plot = ggplot(data = sim00) + 
  geom_line(aes(x = time/365, y = transT * 100, color = 'simul')) + 
  geom_point(data = patient1 %>% filter(celltype == "CD3") ,aes(x = time, y = cellpercentage, color = 'obs'), size = 3) + 
  scale_y_continuous(trans='log10',  limits = c(1, 100), breaks = c(1, 10, 100), labels = c('1', '10','100')) + 
  scale_x_continuous(limits = c(0, 7), breaks = c(0,1,2,3,4,5,6,7), labels = c('0','1', '2', '3', '4', '5', '6', '7')) + 
  labs(y = 'vector-positive cells (%)', x = 'time (year)', color = '') + theme_bw() +  theme(legend.position = c(0.5, 0.2)) +
  ggtitle('CD3+ T Cells')

grid.arrange(cd15plot, cd19plot, cd3plot, ncol = 3)
