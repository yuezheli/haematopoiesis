rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(mrggsave)

source('Depleted_B.r')

basemodel <- mread("erythrocytes_Hb_lymphoid_myeloid") # %>% param(depletedparam)

ssbase <- basemodel %>% init(LT0 = 1000) %>% mrgsim(end = 1500) %>% filter(time == 1500) %>% as.tibble()

ssmodel <- basemodel %>% init(ssbase)

# simulation on depleting only B cell outside bone marrow

## 100% effective depletion
sim11 <- ssmodel %>% init(Bt0 = 0, BMspl0 = 0, BMrec0 = 0) %>% mrgsim(end = 350) %>% as.tibble() %>% select(time, Bt0, BMspl0, BMrec0)

## 90% effective depletion
sim12 <- ssmodel %>% init(Bt0 = ssbase$Bt0 * 0.1, BMspl0 = ssbase$BMspl0 * 0.1, BMrec0 = ssbase$BMrec0 * 0.1) %>% mrgsim(end = 350) %>% as.tibble() %>% select(time, Bt0, BMspl0, BMrec0)

## 50% effective depletion
sim13 <- ssmodel %>% init(Bt0 = ssbase$Bt0 * 0.5, BMspl0 = ssbase$BMspl0 * 0.5, BMrec0 = ssbase$BMrec0 * 0.5) %>% mrgsim(end = 350) %>% as.tibble() %>% select(time, Bt0, BMspl0, BMrec0)


ggplot() + 
  geom_line(data = sim11, aes(x = time/7, y = BMrec0, color = '100%')) + 
  geom_line(data = sim12, aes(x = time/7, y = BMrec0, color = '90%')) + 
  geom_line(data = sim13, aes(x = time/7, y = BMrec0, color = '50%')) + 
  labs(y = 'B cells count in blood', x = 'time (weeks)', color = 'depletion percentage') + theme_bw()

