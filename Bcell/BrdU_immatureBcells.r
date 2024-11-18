rm(list = ls())

library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


source('DepletionParam.r')

# simulate to a B cell steady state
sim_pre <- mread("Bcell_mus") %>% mrgsim(end = 600) %>% as_tibble() %>% tail(n=1)

# BrdU injection, control
mod <- mread("Bcell_mus3") %>% 
        init(UBoe = sim_pre$Boe) %>% init(UBi = sim_pre$Bi) %>% 
        init(UBt = sim_pre$Bt) %>% init(UBMspl = sim_pre$BMspl) %>% init(UBMrec = sim_pre$BMrec) %>% 
        mrgsim(end = 32) %>% as_tibble() 

# BrdU injection, depleted; conditioning strength scan
mod_point7 <- mread("Bcell_mus3") %>%  param(depletedparam) %>% 
        init(UBoe = sim_pre$Boe ) %>% init(UBi = sim_pre$Bi * 0.3) %>% 
        init(UBt = sim_pre$Bt * 0.6) %>% init(UBMspl = sim_pre$BMspl * 0.2) %>% init(UBMrec = 100) %>% mrgsim(end = 32) %>% as_tibble()

mod_point9 <- mread("Bcell_mus3") %>%  param(depletedparam) %>% init(UBoe = sim_pre$Boe ) %>% init(UBi = sim_pre$Bi * 0.1) %>% 
        init(UBt = sim_pre$Bt * 0.6) %>% init(UBMspl = sim_pre$BMspl * 0.2) %>% init(UBMrec = 100) %>% mrgsim(end = 32) %>% as_tibble()

mod_point1 <- mread("Bcell_mus3") %>%  param(depletedparam) %>% init(UBoe = sim_pre$Boe ) %>% init(UBi = sim_pre$Bi * 0.9) %>% 
        init(UBt = sim_pre$Bt * 0.6) %>% init(UBMspl = sim_pre$BMspl * 0.2) %>% init(UBMrec = 100) %>% mrgsim(end = 32) %>% as_tibble() 

Bi_conditiong <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUimmature, color = 'ctrl simul')) + 
        geom_line(data = mod_point7, aes(time/4, BrdUimmature, color = 'depl simul, condition 30%')) +
        geom_line(data = mod_point9, aes(time/4, BrdUimmature, color = 'depl simul, condition 10%')) +
        geom_line(data = mod_point1, aes(time/4, BrdUimmature, color = 'depl simul, condition 90%')) +
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,1) +
        theme_bw() + theme(legend.position = "bottom") + ggtitle("immature B cells")

Brec_conditiong <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUmaturerec, color = 'ctrl simul')) + 
        geom_line(data = mod_point7, aes(time/4, BrdUmaturerec, color = 'depl simul, condition 30%')) +
        geom_line(data = mod_point9, aes(time/4, BrdUmaturerec, color = 'depl simul, condition 10%')) +
        geom_line(data = mod_point1, aes(time/4, BrdUmaturerec, color = 'depl simul, condition 90%')) +
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,1) +
        theme_bw() + theme(legend.position = "bottom") + ggtitle("circulating B cells")

Bt_conditiong <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUtransit, color = 'ctrl simul')) + 
        geom_line(data = mod_point7, aes(time/4, BrdUtransit, color = 'depl simul, condition 30%')) +
        geom_line(data = mod_point9, aes(time/4, BrdUtransit, color = 'depl simul, condition 10%')) +
        geom_line(data = mod_point1, aes(time/4, BrdUtransit, color = 'depl simul, condition 90%')) +
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,1) +
        theme_bw() + theme(legend.position = "bottom") + ggtitle("transitioinal B cells")

BMspl_conditiong <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUsplenicmature, color = 'ctrl simul')) + 
        geom_line(data = mod_point7, aes(time/4, BrdUsplenicmature, color = 'depl simul, condition 30%')) +
        geom_line(data = mod_point9, aes(time/4, BrdUsplenicmature, color = 'depl simul, condition 10%')) +
        geom_line(data = mod_point1, aes(time/4, BrdUsplenicmature, color = 'depl simul, condition 90%')) +
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,1) +
        theme_bw() + theme(legend.position = "bottom") + ggtitle("splenic B cells")

grid.arrange(Bi_conditiong, Brec_conditiong, Bt_conditiong, BMspl_conditiong)


# BrdU injection, depleted; gamma, delta_oe, delta_i_re, mu_i scan

mod_default <- mread("Bcell_mus3") %>%  param(depletedparam) %>% 
        init(UBoe = sim_pre$Boe ) %>% init(UBi = sim_pre$Bi * 0.3) %>% 
        init(UBt = sim_pre$Bt * 0.6) %>% init(UBMspl = sim_pre$BMspl * 0.2) %>% 
        init(UBMrec = 100)


exidata <- expand.idata(gamma = c(0.2, 0.38, 0.9))
mod_default %>%idata_set(exidata) %>% mrgsim(end = 32) %>% plot(BrdUimmature + BrdUmaturerec + BrdUtransit~time) 


exidata2 <- expand.idata(delta_oe = c(0.125 ,0.25, 0.375), delta_i_re = c(0.075, 0.15, 0.22), mu_i = c(0.26, 0.53, 0.77), delta_i_t = c(0.1, 0.2, 0.3))
#exidata2 <- expand.idata(delta_i_re = c(0.075, 0.15, 0.22), mu_i = c(0.26, 0.53, 0.77), delta_i_t = c(0.1, 0.2, 0.3))
mod_default %>%idata_set(exidata2) %>% mrgsim(end = 32) %>% plot(BrdUimmature + BrdUmaturerec + BrdUtransit~time) 

# BrdU injection; gamma scan for the control case
exidata3 <- expand.idata(gamma = c(0.2, 0.38, 0.9))
mread("Bcell_mus3") %>% init(UBoe = sim_pre$Boe) %>% init(UBi = sim_pre$Bi) %>% 
        init(UBt = sim_pre$Bt) %>% init(UBMspl = sim_pre$BMspl) %>% init(UBMrec = sim_pre$BMrec) %>% 
        idata_set(exidata3) %>% mrgsim(end = 32) %>% plot(BrdUimmature + BrdUmaturerec + BrdUtransit~time) 

