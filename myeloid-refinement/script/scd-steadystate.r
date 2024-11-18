# date: 10/8/24
# author: Yuezhe Li 
# purpose of this code: steady state of SCD, published from Zheng et al., 2021


rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)

mod <- mread("model/zheng2021-cleanup") 

sim0 <- mod %>% param(totalCD34infused = 0) %>% 
  init(LT0 = 1) %>% param(ksynalpha = 6e-7) %>%
  mrgsim(end = 30 * 18) %>% as_tibble() 

data = tail(sim0, n = 1)

data$HbSinRET
data$HbSinRBC
