rm(list = ls())

library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library("rstudioapi")
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script
theme_set(theme_bw())

source('DepletionParam.r')

# simulated data from Shahaf et al., 2016 (control)
immature <- read.csv('data/BrdUfraction_control_simul.csv', header = TRUE) %>% filter(Bcelltype == "immatureB")
brec <- read.csv('data/BrdUfraction_control_simul.csv', header = TRUE) %>% filter(Bcelltype == "maturecircul")
transb <- read.csv('data/BrdUfraction_control_simul.csv', header = TRUE) %>% filter(Bcelltype == "transitionalB")
spleenb <- read.csv('data/BrdUfraction_control_simul.csv', header = TRUE) %>% filter(Bcelltype == "splenicmatureB")

# observed data from Shahaf et al., 2016
immature2 <- read.csv('data/BrdUfraction_control_obs.csv', header = TRUE) %>% filter(Bcelltype == "immatureB")
brec2 <- read.csv('data/BrdUfraction_control_obs.csv', header = TRUE) %>% filter(Bcelltype == "maturecircul")
transb2 <- read.csv('data/BrdUfraction_control_obs.csv', header = TRUE) %>% filter(Bcelltype == "transitionalB")
spleenb2 <- read.csv('data/BrdUfraction_control_obs.csv', header = TRUE) %>% filter(Bcelltype == "splenicmatureB")

# simulated data from Shahaf et al., 2016 (depleted)
immature3 <- read.csv('data/BrdUfraction_depleted_simul.csv', header = TRUE) %>% filter(Bcelltype == "immatureB")
brec3 <- read.csv('data/BrdUfraction_depleted_simul.csv', header = TRUE) %>% filter(Bcelltype == "maturecircul")
transb3 <- read.csv('data/BrdUfraction_depleted_simul.csv', header = TRUE) %>% filter(Bcelltype == "transitionalB")
spleenb3 <- read.csv('data/BrdUfraction_depleted_simul.csv', header = TRUE) %>% filter(Bcelltype == "splenicmatureB")

# simulate to a B cell steady state
sim_pre <- mread("Bcell_mus") %>% mrgsim(end = 600) %>% as_tibble() %>% tail(n=1)

sim2_pre <- mread("Bcell_mus") %>% param(depletedparam) %>% mrgsim(end = 600) %>% as_tibble() %>% tail(n=1)


###----------------- test 1 -----------------###

# BrdU injection, control
mod <- mread("Bcell_mus3") %>% 
        init(UBoe = sim_pre$Boe) %>% init(UBi = sim_pre$Bi) %>% 
        init(UBt = sim_pre$Bt) %>% init(UBMspl = sim_pre$BMspl) %>% 
        init(UBMrec = sim_pre$BMrec) %>% 
        mrgsim(end = 32) %>% as_tibble() %>% select(-ID)

# BrdU injection, depleted
mod3 <- mread("Bcell_mus3") %>% param(depletedparam) %>%
        init(UBoe = sim_pre$Boe ) %>% init(UBi = sim_pre$Bi * 0.3) %>% 
        init(UBt = sim_pre$Bt * 0.6) %>% init(UBMspl = sim_pre$BMspl * 0.2) %>% 
        init(UBMrec = 100) %>% 
        mrgsim(end = 32) %>% as_tibble() %>% select(-ID)



Bi <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUimmature, color = 'ctrl simul')) + 
        geom_point(data = immature, aes(time_days, BrdUfrac, color = 'ref ctrl simul'), size = 3) + 
#        geom_point(data = immature2, aes(time_days, BrdUfrac, color = 'ref ctrl obs'), alpha = 0.5) +
        geom_line(data = mod3, aes(time/4, BrdUimmature, color = 'depl simul')) + 
        geom_point(data = immature3, aes(time_days, BrdUfrac, color = 'ref depl simul'), size = 3, shape=17) + 
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,1) +
        theme_bw() + theme(legend.position = "bottom") + ggtitle("immature B cells")

Brec <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUmaturerec, color = 'ctrl simul')) +
        geom_point(data = brec, aes(time_days, BrdUfrac, color = 'ref ctrl simul'), size = 3) + 
#        geom_point(data = brec2, aes(time_days, BrdUfrac, color = 'ref ctrl obs'), alpha = 0.5) + 
        geom_line(data = mod3, aes(time/4, BrdUmaturerec, color = 'depl simul')) +
        geom_point(data = brec3, aes(time_days, BrdUfrac, color = 'ref depl simul'), size = 3, shape=17) + 
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,0.8) +
        theme_bw() + theme(legend.position = "bottom") + ggtitle('circulating mature B cells')

Bt <-  ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUtransit, color = 'ctrl simul')) + 
        geom_point(data = transb, aes(time_days, BrdUfrac, color = 'ref ctrl simul'), size = 3) + 
#        geom_point(data = transb2, aes(time_days, BrdUfrac, color = 'ref ctrl obs'), alpha = 0.5) + 
        geom_line(data = mod3, aes(time/4, BrdUtransit, color = 'depl simul')) + 
        geom_point(data = transb3, aes(time_days, BrdUfrac, color = 'ref depl simul'), size = 3, shape=17) + 
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,0.8) +
        theme_bw() + theme(legend.position = "bottom") + ggtitle('transitional B cells')

BMspl <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUsplenicmature, color = 'ctrl simul')) +
        geom_point(data = spleenb, aes(time_days, BrdUfrac, color = 'ref ctrl simul'), size = 3) + 
#        geom_point(data = spleenb2, aes(time_days, BrdUfrac, color = 'ref ctrl obs'), alpha = 0.5) + 
        geom_line(data = mod3, aes(time/4, BrdUsplenicmature, color = 'depl simul')) +
        geom_point(data = spleenb3, aes(time_days, BrdUfrac, color = 'ref depl simul'), size = 3, shape=17) + 
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,0.4) +
        theme_bw() + theme(legend.position = "bottom") + ggtitle('splenic mature B cells')

grid.arrange(Bi, Brec, Bt, BMspl, ncol = 2)

###----------------- test 2 -----------------###

# BrdU injection, depleted
mod3 <- mread("Bcell_mus3") %>% param(depletedparam) %>%
        init(UBoe = sim2_pre$Boe ) %>% init(UBi = sim2_pre$Bi * 0.3) %>% 
        init(UBt = sim2_pre$Bt * 0.9) %>% init(UBMspl = sim2_pre$BMspl * 0.9) %>% 
        init(UBMrec = 100) %>% 
        mrgsim(end = 32) %>% as_tibble() %>% select(-ID)


Bi <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUimmature, color = 'ctrl simul')) + 
        geom_point(data = immature, aes(time_days, BrdUfrac, color = 'ref ctrl simul'), size = 3) + 
        geom_line(data = mod3, aes(time/4, BrdUimmature, color = 'depl simul')) + 
        geom_point(data = immature3, aes(time_days, BrdUfrac, color = 'ref depl simul'), size = 3, shape=17) + 
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,1) + theme(legend.position = "bottom") + ggtitle("immature B cells")

Brec <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUmaturerec, color = 'ctrl simul')) +
        geom_point(data = brec, aes(time_days, BrdUfrac, color = 'ref ctrl simul'), size = 3) + 
        geom_line(data = mod3, aes(time/4, BrdUmaturerec, color = 'depl simul')) +
        geom_point(data = brec3, aes(time_days, BrdUfrac, color = 'ref depl simul'), size = 3, shape=17) + 
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,0.8) + theme(legend.position = "bottom") + ggtitle('circulating mature B cells')

Bt <-  ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUtransit, color = 'ctrl simul')) + 
        geom_point(data = transb, aes(time_days, BrdUfrac, color = 'ref ctrl simul'), size = 3) + 
        geom_line(data = mod3, aes(time/4, BrdUtransit, color = 'depl simul')) + 
        geom_point(data = transb3, aes(time_days, BrdUfrac, color = 'ref depl simul'), size = 3, shape=17) + 
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,0.8) + theme(legend.position = "bottom") + ggtitle('transitional B cells')

BMspl <- ggplot() + 
        geom_line(data = mod, aes(time/4, BrdUsplenicmature, color = 'ctrl simul')) +
        geom_point(data = spleenb, aes(time_days, BrdUfrac, color = 'ref ctrl simul'), size = 3) + 
        geom_line(data = mod3, aes(time/4, BrdUsplenicmature, color = 'depl simul')) +
        geom_point(data = spleenb3, aes(time_days, BrdUfrac, color = 'ref depl simul'), size = 3, shape=17) + 
        labs(y = 'BrdU fraction', x = 'time (days)', color = '') + ylim(0,0.4) + theme(legend.position = "bottom") + ggtitle('splenic mature B cells')

grid.arrange(Bi, Brec, Bt, BMspl, ncol = 2)