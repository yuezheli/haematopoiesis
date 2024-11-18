# Goal: reproduce Fig4 in Thomas-Vaslin et al., 2008

rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


mod <- mread("Tcell")
# run the system to steady state
pre = mod %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1) %>% select(c(-ID, -time))

# GCV perturbation
post <- mod %>% init(pre * 0.3)%>% mrgsim(end = 70, delta = 1) %>% as_tibble() 
post$time = post$time + 7 # add to match time in the paper

# read in the observed data
obs = read.csv('ThomasVaslin2008.csv', header = TRUE)

# compare data
dnplot <- ggplot() + 
  geom_line(data = post, aes(x = time, y = DN, color = "simul"), size = 2) + 
  geom_point(data = filter(obs, cell_type == "DN"), aes(x = time, y = cell_count, color = "obs"), size = 4) + 
  labs(y = 'cell count (M)', x = 'time (days)', color = 'DN cells') + theme(legend.position = "bottom") + theme_bw()

dpplot <- ggplot() + 
  geom_line(data = post, aes(x = time, y = DP, color = "simul"), size = 2) + 
  geom_point(data = filter(obs, cell_type == "DP"), aes(x = time, y = cell_count, color = "obs"), size = 4) + 
  labs(y = 'cell count (M)', x = 'time (days)', color = 'DP cells') + theme(legend.position = "bottom") + theme_bw()

sp4plot <- ggplot() + 
  geom_line(data = post, aes(x = time, y = SP4, color = "simul"), size = 2) + 
  geom_point(data = filter(obs, cell_type == "SP4"), aes(x = time, y = cell_count, color = "obs"), size = 4) + 
  labs(y = 'cell count (M)', x = 'time (days)', color = 'SP4 cells') + theme(legend.position = "bottom") + theme_bw()

sp8plot <- ggplot() + 
  geom_line(data = post, aes(x = time, y = SP8, color = "simul"), size = 2) + 
  geom_point(data = filter(obs, cell_type == "SP8"), aes(x = time, y = cell_count, color = "obs"), size = 4) + 
  labs(y = 'cell count (M)', x = 'time (days)', color = 'SP8 cells') + theme(legend.position = "bottom") + theme_bw()


png('dynamics_verification.png', width = 10, height = 4, units = 'in', res = 300)
grid.arrange(dnplot, dpplot, sp4plot, sp8plot, ncol = 2)
dev.off()

