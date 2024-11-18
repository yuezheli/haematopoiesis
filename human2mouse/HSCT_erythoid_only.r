rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

source('transplant_functions.r')

# read data from Boyer et al., 2019
hsct_rbc <- read.csv('MouseData/Boyer2019_1M_HSC.csv', header = TRUE) %>% filter(cell_type == "RBC")
hsct_rbc[,4] <- read.csv('MouseData/Boyer2019_1M_HSC.csv', header = TRUE) %>% filter(cell_type == "RBClower") %>% select(cell_count) 
colnames(hsct_rbc)[4] <- "RBC_lower"
hsct_rbc$RBC_upper <- hsct_rbc$cell_count + (hsct_rbc$cell_count - hsct_rbc$RBC_lower)

mppt_rbc <- read.csv('MouseData/Boyer2019_1N_MPP.csv', header = TRUE) %>% filter(cell_type == "RBC")
mppt_rbc[,4] <- read.csv('MouseData/Boyer2019_1N_MPP.csv', header = TRUE) %>% filter(cell_type == "RBClower") %>% select(cell_count) 
colnames(mppt_rbc)[4] <- "RBC_lower"
mppt_rbc$RBC_upper <- mppt_rbc$cell_count + (mppt_rbc$cell_count - mppt_rbc$RBC_lower)


modo <- mread("erythrocytes_Hb") %>% param(aST = 32, aBFUE = 16) 

# for HSCT, conditioning strength & proliferation rate
hsct1_2point5 <- conditioning_proliferation_hsct(0.1, hscprof = 1/(2.5 * 7))
hsct5_2point5 <- conditioning_proliferation_hsct(0.5, hscprof = 1/(2.5 * 7))
hsct9_2point5 <- conditioning_proliferation_hsct(0.9, hscprof = 1/(2.5 * 7))
hsct9_1point5 <- conditioning_proliferation_hsct(0.9, hscprof = 1/(0.5 * 7))

# for MPPT, conditioning and ratio between MPP1 and MPP2
mppt1_0 = conditioning_mppt(0.1, 0)
mppt5_0 = conditioning_mppt(0.5, 0)
mppt5_1 = conditioning_mppt(0.5, 0.1)
mppt5_6 = conditioning_mppt(0.5, 0.6)
mppt9_0 = conditioning_mppt(0.9, 0)


rbc_hsct <- ggplot() + 
  geom_line(data = hsct1_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%' )) + 
  geom_line(data = hsct5_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
  geom_line(data = hsct9_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 90%' )) + 
  #  geom_line(data = hsct9_1point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 90%, HSC doubling time = 0.5 weeks' )) + 
  geom_errorbar(data=hsct_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
  geom_point(data = hsct_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
  labs(y = 'RBC count (#/uL)', x = 'time (days)', color = '') +
  theme_bw() + ggtitle('HSC transplant') + 
  scale_y_continuous(trans='log10',  limits = c(1, 5e6), breaks = c(1, 1e3, 1e6), 
                     labels = c('1', '1K', '1M')) 

rbc_mppt <- ggplot() + 
  geom_line(data = mppt1_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%')) + 
  geom_line(data = mppt5_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
  #  geom_line(data = mppt5_1, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%, late stage MPP = 10%' )) + 
  #  geom_line(data = mppt5_6, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%, late stage MPP = 60%' )) + 
  geom_line(data = mppt9_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 90%')) + 
  geom_errorbar(data=mppt_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
  geom_point(data = mppt_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
  labs(y = 'RBC count (#/uL)', x = 'time (days)', color = '') + theme_bw() + 
  ggtitle('MPP transplant')  


grid.arrange(rbc_hsct, rbc_mppt, ncol = 1)

## save this graph as a validation
png('img/HSCT_withHb_erythocytes.png', width = 8, height = 4, units = 'in', res = 300)
grid.arrange(rbc_hsct, rbc_mppt, ncol = 1)
dev.off()


hsct2 <- read.csv('MouseData/Boyer2019_1A_HSC.csv', header = TRUE) %>% filter(cell_type == "RBC")
mppt2 <- read.csv('MouseData/Boyer2019_1B_MPP.csv', header = TRUE) %>% filter(cell_type == "RBC")


rbc_hsct2 <- ggplot() + 
  geom_line(data = hsct1_2point5, aes(x = time, y = RBC1/(RBC0 + RBC1), color = 'conditioning = 10%' )) + 
  geom_line(data = hsct5_2point5, aes(x = time, y = RBC1/(RBC0 + RBC1), color = 'conditioning = 50%' )) + 
  geom_line(data = hsct9_2point5, aes(x = time, y = RBC1/(RBC0 + RBC1), color = 'conditioning = 90%' )) + 
  geom_point(data = hsct2, aes(x = time, y = cell_fraction/100, color = 'obs')) + 
  labs(y = 'RBC fraction', x = 'time (days)', color = '') +
  theme_bw() + ggtitle('HSC transplant') 

rbc_mppt2 <- ggplot() + 
  geom_line(data = mppt1_0, aes(x = time, y = RBC1/(RBC0 + RBC1), color = 'conditioning = 10%')) + 
  geom_line(data = mppt5_0, aes(x = time, y = RBC1/(RBC0 + RBC1), color = 'conditioning = 50%' )) + 
  #  geom_line(data = mppt5_1, aes(x = time, y = RBC1/(RBC0 + RBC1), color = 'conditioning = 50%, late stage MPP = 10%' )) + 
  #  geom_line(data = mppt5_6, aes(x = time, y = RBC1/(RBC0 + RBC1), color = 'conditioning = 50%, late stage MPP = 60%' )) + 
  geom_line(data = mppt9_0, aes(x = time, y = RBC1/(RBC0 + RBC1), color = 'conditioning = 90%')) + 
  geom_point(data = mppt2, aes(x = time, y = cell_fraction/100, color = 'obs')) + 
  labs(y = 'RBC fraction', x = 'time (days)', color = '') + ggtitle('MPP transplant')  


grid.arrange(rbc_hsct2, rbc_mppt2, ncol = 1)

