rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(mrggsave)

source('transplant_functions.r')
source('LoadObs.r')
source('Depleted_B.r')


# generate figure with multiple conditioning strength
## note this figure is only used in the readme in this folder, not the readme at repo level
if(FALSE)
{
modo <- mread("erythrocytes_Hb_lymphoid_myeloid") %>% param(depletedparam)
# modo <- mread("erythrocytes_Hb_lymphoid_myeloid") %>% param(depletedparam) %>% param(alpha_muN = 0.5, pN = 0.5, delta = 0.001) 
}

#-------------------------- Steady State --------------------------#
sim0 <- modo %>% init(LT0 = 1000) %>% mrgsim(end = 1500) %>% filter(time == 1500) %>% as.tibble()

#-------------------------- HSCT --------------------------#
hsct1_2point5 <- conditioning_proliferation_hsct(0.1, hscprof = 1/(2.5 * 7))
hsct2_2point5 <- conditioning_proliferation_hsct(0.2, hscprof = 1/(2.5 * 7))
hsct3_2point5 <- conditioning_proliferation_hsct(0.3, hscprof = 1/(2.5 * 7))
hsct4_2point5 <- conditioning_proliferation_hsct(0.4, hscprof = 1/(2.5 * 7))
hsct5_2point5 <- conditioning_proliferation_hsct(0.5, hscprof = 1/(2.5 * 7))
hsct6_2point5 <- conditioning_proliferation_hsct(0.6, hscprof = 1/(2.5 * 7))
hsct7_2point5 <- conditioning_proliferation_hsct(0.7, hscprof = 1/(2.5 * 7))
hsct8_2point5 <- conditioning_proliferation_hsct(0.8, hscprof = 1/(2.5 * 7))
hsct9_2point5 <- conditioning_proliferation_hsct(0.9, hscprof = 1/(2.5 * 7))

rbc_hsct <- ggplot() + 
  geom_line(data = hsct1_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%' )) + 
  geom_line(data = hsct2_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 20%' )) +
  geom_line(data = hsct3_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 30%' )) +
  geom_line(data = hsct4_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 40%' )) +
  geom_line(data = hsct5_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
  geom_line(data = hsct6_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 60%' )) + 
  geom_line(data = hsct7_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 70%' )) + 
  geom_line(data = hsct8_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 80%' )) + 
  geom_line(data = hsct9_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 90%' )) + 
  geom_errorbar(data=hsct_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
  geom_point(data = hsct_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
  labs(y = 'donor RBC (#/uL)', x = 'time (days)', color = '') +
  theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") +
  scale_y_continuous(trans='log10',  limits = c(1, 1e7), breaks = c(1, 1e3, 1e6), labels = c('1', '1K', '1M')) 

b_hsct <- ggplot() + 
  geom_line(data = hsct1_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 10%' )) + 
  geom_line(data = hsct2_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 20%' )) +
  geom_line(data = hsct3_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 30%' )) +
  geom_line(data = hsct4_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 40%' )) +
  geom_line(data = hsct5_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) + 
  geom_line(data = hsct6_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 60%' )) + 
  geom_line(data = hsct7_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 70%' )) + 
  geom_line(data = hsct8_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 80%' )) + 
  geom_line(data = hsct9_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 90%' )) + 
  geom_errorbar(data=hsct_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
  geom_point(data = hsct_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
  labs(y = 'donor B cell (%)', x = 'time (days)', color = '') +
  theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") 

t_hsct <- ggplot() + 
  geom_line(data = hsct1_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 10%' )) + 
  geom_line(data = hsct2_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 20%' )) +
  geom_line(data = hsct3_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 30%' )) +
  geom_line(data = hsct4_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 40%' )) +
  geom_line(data = hsct5_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 50%' )) + 
  geom_line(data = hsct6_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 60%' )) + 
  geom_line(data = hsct7_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 70%' )) + 
  geom_line(data = hsct8_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 80%' )) + 
  geom_line(data = hsct9_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 90%' )) + 
  geom_errorbar(data=hsct_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
  geom_point(data = hsct_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
  labs(y = 'donor T cell (%)', x = 'time (days)', color = '') +
  theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") 

#-------------------------- MPPT --------------------------#
mppt1_0 = conditioning_mppt(0.1, 0)
mppt2_0 = conditioning_mppt(0.2, 0)
mppt3_0 = conditioning_mppt(0.3, 0)
mppt4_0 = conditioning_mppt(0.4, 0)
mppt5_0 = conditioning_mppt(0.5, 0)
mppt6_0 = conditioning_mppt(0.6, 0)
mppt7_0 = conditioning_mppt(0.7, 0)
mppt8_0 = conditioning_mppt(0.8, 0)
mppt9_0 = conditioning_mppt(0.9, 0)

rbc_mppt <- ggplot() + 
  geom_line(data = mppt1_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%')) +
  geom_line(data = mppt2_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 20%' )) + 
  geom_line(data = mppt3_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 30%' )) + 
  geom_line(data = mppt4_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 40%' )) + 
  geom_line(data = mppt5_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
  geom_line(data = mppt6_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 60%' )) + 
  geom_line(data = mppt7_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 70%' )) + 
  geom_line(data = mppt8_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 80%' )) + 
  geom_line(data = mppt9_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 90%')) + 
  geom_errorbar(data=mppt_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
  geom_point(data = mppt_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
  labs(y = 'donor RBC (#/uL)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
  ggtitle('MPP transplant')  

b_mppt <- ggplot() + 
  geom_line(data = mppt1_0, aes(x = time, y = transB * 100, color = 'conditioning = 10%')) +
  geom_line(data = mppt2_0, aes(x = time, y = transB * 100, color = 'conditioning = 20%' )) + 
  geom_line(data = mppt3_0, aes(x = time, y = transB * 100, color = 'conditioning = 30%' )) + 
  geom_line(data = mppt4_0, aes(x = time, y = transB * 100, color = 'conditioning = 40%' )) + 
  geom_line(data = mppt5_0, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) + 
  geom_line(data = mppt6_0, aes(x = time, y = transB * 100, color = 'conditioning = 60%' )) + 
  geom_line(data = mppt7_0, aes(x = time, y = transB * 100, color = 'conditioning = 70%' )) + 
  geom_line(data = mppt8_0, aes(x = time, y = transB * 100, color = 'conditioning = 80%' )) + 
  geom_line(data = mppt9_0, aes(x = time, y = transB * 100, color = 'conditioning = 90%')) + 
  geom_errorbar(data=mppt_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
  geom_point(data = mppt_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
  labs(y = 'donor B cell (%)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
  ggtitle('MPP transplant')  

t_mppt <- ggplot() + 
  geom_line(data = mppt1_0, aes(x = time, y = transT * 100, color = 'conditioning = 10%')) +
  geom_line(data = mppt2_0, aes(x = time, y = transT * 100, color = 'conditioning = 20%' )) + 
  geom_line(data = mppt3_0, aes(x = time, y = transT * 100, color = 'conditioning = 30%' )) + 
  geom_line(data = mppt4_0, aes(x = time, y = transT * 100, color = 'conditioning = 40%' )) + 
  geom_line(data = mppt5_0, aes(x = time, y = transT * 100, color = 'conditioning = 50%' )) + 
  geom_line(data = mppt6_0, aes(x = time, y = transT * 100, color = 'conditioning = 60%' )) + 
  geom_line(data = mppt7_0, aes(x = time, y = transT * 100, color = 'conditioning = 70%' )) + 
  geom_line(data = mppt8_0, aes(x = time, y = transT * 100, color = 'conditioning = 80%' )) + 
  geom_line(data = mppt9_0, aes(x = time, y = transT * 100, color = 'conditioning = 90%')) + 
  geom_errorbar(data=mppt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
  geom_point(data = mppt_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
  labs(y = 'donor T cell (%)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
  ggtitle('MPP transplant')  

#-------------------------- CLPT --------------------------#
clpt_1 = conditioning_clpt(conditioning_strength = 0.1, t_conditioning_strength = 0.5)
clpt_2 = conditioning_clpt(conditioning_strength = 0.2, t_conditioning_strength = 0.5)
clpt_3 = conditioning_clpt(conditioning_strength = 0.3, t_conditioning_strength = 0.5)
clpt_4 = conditioning_clpt(conditioning_strength = 0.4, t_conditioning_strength = 0.5)
clpt_5 = conditioning_clpt(conditioning_strength = 0.5, t_conditioning_strength = 0.5)
clpt_6 = conditioning_clpt(conditioning_strength = 0.6, t_conditioning_strength = 0.5)
clpt_7 = conditioning_clpt(conditioning_strength = 0.7, t_conditioning_strength = 0.5)
clpt_8 = conditioning_clpt(conditioning_strength = 0.8, t_conditioning_strength = 0.5)
clpt_9 = conditioning_clpt(conditioning_strength = 0.9, t_conditioning_strength = 0.5)

b_clpt <- ggplot() + 
  geom_line(data = clpt_1, aes(x = time, y = transB * 100, color = 'conditioning = 10%')) + 
  geom_line(data = clpt_2, aes(x = time, y = transB * 100, color = 'conditioning = 20%' )) +  
  geom_line(data = clpt_3, aes(x = time, y = transB * 100, color = 'conditioning = 30%' )) +  
  geom_line(data = clpt_4, aes(x = time, y = transB * 100, color = 'conditioning = 40%' )) +  
  geom_line(data = clpt_5, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) +  
  geom_line(data = clpt_6, aes(x = time, y = transB * 100, color = 'conditioning = 60%' )) +  
  geom_line(data = clpt_7, aes(x = time, y = transB * 100, color = 'conditioning = 70%' )) +  
  geom_line(data = clpt_8, aes(x = time, y = transB * 100, color = 'conditioning = 80%' )) +  
  geom_line(data = clpt_9, aes(x = time, y = transB * 100, color = 'conditioning = 90%' )) +  
  geom_point(data = clpt_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
  geom_errorbar(data=clpt_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
  labs(y = 'donor B (%)', x = 'time (days)', color = 'B cells') + theme_bw() + theme(legend.position="bottom") +
  ggtitle('CLP transplant')  

t_clpt <- ggplot() + 
  geom_line(data = clpt_1, aes(x = time, y = transT * 100, color = 'conditioning = 10%')) + 
  geom_line(data = clpt_2, aes(x = time, y = transT * 100, color = 'conditioning = 20%')) + 
  geom_line(data = clpt_3, aes(x = time, y = transT * 100, color = 'conditioning = 30%')) + 
  geom_line(data = clpt_4, aes(x = time, y = transT * 100, color = 'conditioning = 40%')) + 
  geom_line(data = clpt_5, aes(x = time, y = transT * 100, color = 'conditioning = 50%')) + 
  geom_line(data = clpt_6, aes(x = time, y = transT * 100, color = 'conditioning = 60%')) + 
  geom_line(data = clpt_7, aes(x = time, y = transT * 100, color = 'conditioning = 70%')) + 
  geom_line(data = clpt_8, aes(x = time, y = transT * 100, color = 'conditioning = 80%')) + 
  geom_line(data = clpt_9, aes(x = time, y = transT * 100, color = 'conditioning = 90%')) + 
  geom_errorbar(data=clpt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
  geom_point(data = clpt_T, aes(x = time, y = cell_fraction, color = 'obs')) +
  labs(y = 'donor T (%)', x = 'time (days)', color = 'naive T') + theme_bw() + theme(legend.position="bottom") +
  ggtitle('CLP transplant') 

#-------------------------- save the plots --------------------------#
if(FALSE)
  png('img/HSCT_MPPT_CLPT_fullmodel.png', width = 16, height = 12, units = 'in', res = 300)
  grid.arrange(rbc_hsct, rbc_mppt, 
               b_hsct, b_mppt,
               t_hsct, t_mppt, 
               b_clpt, t_clpt, ncol = 3)
  dev.off()
}

grid.arrange(rbc_hsct, rbc_mppt, 
             b_hsct, b_mppt,
             t_hsct, t_mppt, 
             b_clpt, t_clpt, ncol = 3)
## save the file to pdf
#mrggsave(list(rbc_hsct, rbc_mppt, b_hsct, b_mppt,t_hsct, t_mppt, b_clpt, t_clpt), script = 'HSCT_fullmodel.r', dir = 'img/', stem = paste("HSCT_MPPT_CLPT"), width = 8, height = 4)
}


# add this section to create a more simplified validational image for the summary of the repo
if(FALSE)
{
  modo <- mread("erythrocytes_Hb_lymphoid_myeloid") %>% param(depletedparam)
  sim0 <- modo %>% init(LT0 = 1000) %>% mrgsim(end = 1500) %>% filter(time == 1500) %>% as.tibble()
  
  # reduce the number of simulations to run
  hsct1_2point5 <- conditioning_proliferation_hsct(0.1, hscprof = 1/(2.5 * 7))
  hsct5_2point5 <- conditioning_proliferation_hsct(0.5, hscprof = 1/(2.5 * 7))
  
  mppt1_0 = conditioning_mppt(0.1, 0)
  mppt5_0 = conditioning_mppt(0.5, 0)
  
  clpt_1 = conditioning_clpt(conditioning_strength = 0.1, t_conditioning_strength = 0.5)
  clpt_5 = conditioning_clpt(conditioning_strength = 0.5, t_conditioning_strength = 0.5)
  
  # plot simpliefied figures
  rbc_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
    geom_point(data = hsct_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
    labs(y = 'donor RBC (#/uL)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") +
    scale_y_continuous(trans='log10',  limits = c(1, 1e7), breaks = c(1, 1e3, 1e6), labels = c('1', '1K', '1M')) 
  
  b_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    geom_point(data = hsct_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor B cell (%)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") 
  
  t_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = hsct_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor T cell (%)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") 
  
  
  rbc_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
    geom_point(data = mppt_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
    labs(y = 'donor RBC (#/uL)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  b_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = transB * 100, color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    geom_point(data = mppt_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor B cell (%)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  t_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = transT * 100, color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = transT * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = mppt_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor T cell (%)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  
  b_clpt <- ggplot() + 
    geom_line(data = clpt_1, aes(x = time, y = transB * 100, color = 'conditioning = 10%')) + 
    geom_line(data = clpt_5, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) +    
    geom_point(data = clpt_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    geom_errorbar(data=clpt_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    labs(y = 'donor B (%)', x = 'time (days)', color = 'B cells') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('CLP transplant')  
  
  t_clpt <- ggplot() + 
    geom_line(data = clpt_1, aes(x = time, y = transT * 100, color = 'conditioning = 10%')) + 
    geom_line(data = clpt_5, aes(x = time, y = transT * 100, color = 'conditioning = 50%')) + 
    geom_errorbar(data=clpt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = clpt_T, aes(x = time, y = cell_fraction, color = 'obs')) +
    labs(y = 'donor T (%)', x = 'time (days)', color = 'naive T') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('CLP transplant') 
  
  
  # save the picture
  png('img/boyer_full_unadjusted.png', width = 1500, height = 800, res = 100)
  grid.arrange(b_hsct, t_hsct, # rbc_hsct, 
               b_mppt, t_mppt, # rbc_mppt, 
               b_clpt, t_clpt, 
               ncol = 2)
  dev.off()
}

# after DN/ DP dynamics is adjusted model is adjusted
if(FALSE)
{
  modo <- mread("erythrocytes_Hb_lymphoid_myeloid") %>% param(depletedparam) %>% param(alpha_muN = 0.5, pN = 0.5, delta = 0.001) 
  sim0 <- modo %>% init(LT0 = 1000) %>% mrgsim(end = 1500) %>% filter(time == 1500) %>% as.tibble()
  
  # reduce the number of simulations to run
  hsct1_2point5 <- conditioning_proliferation_hsct(0.1, hscprof = 1/(2.5 * 7))
  hsct5_2point5 <- conditioning_proliferation_hsct(0.5, hscprof = 1/(2.5 * 7))
  
  mppt1_0 = conditioning_mppt(0.1, 0)
  mppt5_0 = conditioning_mppt(0.5, 0)
  
  clpt_1 = conditioning_clpt(conditioning_strength = 0.1, t_conditioning_strength = 0.5)
  clpt_5 = conditioning_clpt(conditioning_strength = 0.5, t_conditioning_strength = 0.5)
  
  # plot simpliefied figures
  rbc_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
    geom_point(data = hsct_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
    labs(y = 'donor RBC (#/uL)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") +
    scale_y_continuous(trans='log10',  limits = c(1, 1e7), breaks = c(1, 1e3, 1e6), labels = c('1', '1K', '1M')) 
  
  b_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    geom_point(data = hsct_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor B cell (%)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") 
  
  t_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = hsct_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor T cell (%)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") 
  
  
  rbc_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
    geom_point(data = mppt_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
    labs(y = 'donor RBC (#/uL)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  b_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = transB * 100, color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    geom_point(data = mppt_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor B cell (%)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  t_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = transT * 100, color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = transT * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = mppt_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor T cell (%)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  
  b_clpt <- ggplot() + 
    geom_line(data = clpt_1, aes(x = time, y = transB * 100, color = 'conditioning = 10%')) + 
    geom_line(data = clpt_5, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) +    
    geom_point(data = clpt_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    geom_errorbar(data=clpt_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    labs(y = 'donor B (%)', x = 'time (days)', color = 'B cells') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('CLP transplant')  
  
  t_clpt <- ggplot() + 
    geom_line(data = clpt_1, aes(x = time, y = transT * 100, color = 'conditioning = 10%')) + 
    geom_line(data = clpt_5, aes(x = time, y = transT * 100, color = 'conditioning = 50%')) + 
    geom_errorbar(data=clpt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = clpt_T, aes(x = time, y = cell_fraction, color = 'obs')) +
    labs(y = 'donor T (%)', x = 'time (days)', color = 'naive T') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('CLP transplant') 
  
  
  # save the picture
  png('img/boyer_full_adjusted.png', width = 1500, height = 800, res = 100)
  grid.arrange(b_hsct, t_hsct, # rbc_hsct, 
               b_mppt, t_mppt, # rbc_mppt, 
               b_clpt, t_clpt, 
               ncol = 2)
  dev.off()
}

# final model 
if(TRUE)
{
  modo <- mread("erythrocytes_Hb_lymphoid_myeloid") %>% param(depletedparam) %>% param(alpha_muN = 0.5, pN = 0.5, delta = 0.001, alphaCLP2DN = 7.5e-4) 
  sim0 <- modo %>% init(LT0 = 1000) %>% mrgsim(end = 1500) %>% filter(time == 1500) %>% as.tibble()
  
  # reduce the number of simulations to run
  hsct1_2point5 <- conditioning_proliferation_hsct(0.1, hscprof = 1/(2.5 * 7))
  hsct5_2point5 <- conditioning_proliferation_hsct(0.5, hscprof = 1/(2.5 * 7))
  
  mppt1_0 = conditioning_mppt(0.1, 0)
  mppt5_0 = conditioning_mppt(0.5, 0)
  
  clpt_1 = conditioning_clpt(conditioning_strength = 0.1, t_conditioning_strength = 0.5)
  clpt_5 = conditioning_clpt(conditioning_strength = 0.5, t_conditioning_strength = 0.5)
  
  # plot simpliefied figures
  rbc_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
    geom_point(data = hsct_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
    labs(y = 'donor RBC (#/uL)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") +
    scale_y_continuous(trans='log10',  limits = c(1, 1e7), breaks = c(1, 1e3, 1e6), labels = c('1', '1K', '1M')) 
  
  b_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    geom_point(data = hsct_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor B cell (%)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") 
  
  t_hsct <- ggplot() + 
    geom_line(data = hsct1_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 10%' )) + 
    geom_line(data = hsct5_2point5, aes(x = time, y = transT * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=hsct_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = hsct_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor T cell (%)', x = 'time (days)', color = '') +
    theme_bw() + ggtitle('HSC transplant') + theme(legend.position="bottom") 
  
  
  rbc_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = RBC1/(2e3), color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_rbc, aes(x=time, ymin=RBC_lower,ymax=RBC_upper), alpha = 0.5) +
    geom_point(data = mppt_rbc, aes(x = time, y = cell_count, color = 'obs')) + 
    labs(y = 'donor RBC (#/uL)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  b_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = transB * 100, color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    geom_point(data = mppt_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor B cell (%)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  t_mppt <- ggplot() + 
    geom_line(data = mppt1_0, aes(x = time, y = transT * 100, color = 'conditioning = 10%')) +
    geom_line(data = mppt5_0, aes(x = time, y = transT * 100, color = 'conditioning = 50%' )) + 
    geom_errorbar(data=mppt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = mppt_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
    labs(y = 'donor T cell (%)', x = 'time (days)', color = '') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('MPP transplant')  
  
  
  b_clpt <- ggplot() + 
    geom_line(data = clpt_1, aes(x = time, y = transB * 100, color = 'conditioning = 10%')) + 
    geom_line(data = clpt_5, aes(x = time, y = transB * 100, color = 'conditioning = 50%' )) +    
    geom_point(data = clpt_B, aes(x = time, y = cell_fraction, color = 'obs')) + 
    geom_errorbar(data=clpt_B, aes(x=time, ymin=B_lower,ymax=B_upper), alpha = 0.5) +
    labs(y = 'donor B (%)', x = 'time (days)', color = 'B cells') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('CLP transplant')  
  
  t_clpt <- ggplot() + 
    geom_line(data = clpt_1, aes(x = time, y = transT * 100, color = 'conditioning = 10%')) + 
    geom_line(data = clpt_5, aes(x = time, y = transT * 100, color = 'conditioning = 50%')) + 
    geom_errorbar(data=clpt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
    geom_point(data = clpt_T, aes(x = time, y = cell_fraction, color = 'obs')) +
    labs(y = 'donor T (%)', x = 'time (days)', color = 'naive T') + theme_bw() + theme(legend.position="bottom") +
    ggtitle('CLP transplant') 
  
  
  # save the picture
  png('img/boyer_full_final.png', width = 1500, height = 800, res = 100)
  grid.arrange(b_hsct, t_hsct, # rbc_hsct, 
               b_mppt, t_mppt, # rbc_mppt, 
               b_clpt, t_clpt, 
               ncol = 2)
  dev.off()
}