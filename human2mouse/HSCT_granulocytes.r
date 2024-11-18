rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script

library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

source('transplant_functions.r')

modo <- mread("erythrocytes_Hb_lymphoid_myeloid") 

hsct_gm <- read.csv('MouseData/Boyer2019_1A_HSC.csv', header = TRUE) %>% filter(cell_type == "GM")
mppt_gm <- read.csv('MouseData/Boyer2019_1B_MPP.csv', header = TRUE) %>% filter(cell_type == "GM")
gmpt_gm <- read.csv('MouseData/Boyer2019_1B_MPP.csv', header = TRUE) %>% filter(cell_type == "GM")


hsct_1 = conditioning_proliferation_hsct(progenitor_conditioning_strength = 0.1, simendtime = 100)
hsct_5 = conditioning_proliferation_hsct(progenitor_conditioning_strength = 0.5,simendtime = 100)

mppt1 = conditioning_mppt(conditioning_strength = 0.1, latempp = 0, simendtime = 50)
mppt5 = conditioning_mppt(conditioning_strength = 0.5, latempp = 0, simendtime = 50)

gmpt1 = conditioning_gmpt(conditioning_strength = 0.1, lategmp = 0, simendtime = 50)
gmpt5 = conditioning_gmpt(conditioning_strength = 0.5, lategmp = 0, simendtime = 50)

gm_hsct <- ggplot() + 
  geom_line(data = hsct_1, aes(x = time, y = GM1/(GM0 + GM1) * 100, color = 'conditioning = 10%')) + 
  geom_line(data = hsct_5, aes(x = time, y = GM1/(GM0 + GM1) * 100, color = 'conditioning = 50%' )) +  
  geom_point(data = hsct_gm, aes(x = time, y = cell_fraction, color = 'obs')) + 
  labs(y = 'granulocyte percentage (%)', x = 'time (days)', color = 'GM') + theme_bw() + theme(legend.position='bottom') +
  ggtitle('HSC transplant')  

gm_mppt <- ggplot() + 
  geom_line(data = mppt1, aes(x = time, y = GM1/(GM0 + GM1) * 100, color = 'conditioning = 10%')) + 
  geom_line(data = mppt5, aes(x = time, y = GM1/(GM0 + GM1) * 100, color = 'conditioning = 50%' )) +  
  geom_point(data = mppt_gm, aes(x = time, y =  cell_fraction, color = 'obs')) + xlim(0, 50) + 
  labs(y = 'granulocyte percentage (%)', x = 'time (days)', color = 'GM') + theme_bw() + theme(legend.position='bottom') +
  ggtitle('MPP transplant')  


gm_gmpt <- ggplot() + 
  geom_line(data = gmpt1, aes(x = time, y = GM1/(GM0 + GM1) * 100, color = 'conditioning = 10%')) + 
  geom_line(data = gmpt5, aes(x = time, y = GM1/(GM0 + GM1) * 100, color = 'conditioning = 50%' )) +  
  geom_point(data = gmpt_gm, aes(x = time, y =  cell_fraction, color = 'obs')) +  xlim(0, 50) + 
  labs(y = 'granulocyte percentage (%)', x = 'time (days)', color = 'GM') + theme_bw() + theme(legend.position='bottom') +
  ggtitle('GMP transplant')  


## ---------------------------------------------- ##
#png('img/HSCT_withHb_GM.png', width = 15, height = 4, units = 'in', res = 300)
grid.arrange(gm_hsct, gm_mppt, gm_gmpt, ncol = 3)
#dev.off()