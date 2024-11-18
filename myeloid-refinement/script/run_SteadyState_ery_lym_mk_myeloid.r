# date: 10/1/24
# author: Yuezhe Li 
# purpose of this code: to test integration of megakaryocytes and steady state of the model 

rm(list=ls())  #Clear out existing objects
gc()

library(here)
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(ggtext)

# previous model for health people 
mod <- mread("model/human_erythroid_lymphoid_myeloid") %>% param(epsilon_spl = 0.032, delta_dp = 1e-6)

# steady state of the orginal model 
ss <- mod %>% init(LT0 = 1000) %>%  
  param(betaGM = 2.5) %>% # adjusted based on new information on half life; https://pubmed.ncbi.nlm.nih.gov/28303306/
  mrgsim(end = 365 * 7) %>% as_tibble() %>% tail(n = 1) %>% select(-c(ID, time)) %>% 
  mutate(MPP = MPP01 + MPP02) %>% 
  mutate(CMP = CMP01 + CMP02) %>%
  mutate(GMP = GMP01 + GMP02) %>%
  select(MPP, CMP, GMP, GMconc, RETconc, RBCconc, Bconc, Tconc, totalHb, HbinRBC)

# run simplified model 
mod_simp <- mread("model/human_erythroid_lymphoid_myeloid_simplified") 

ss_simp <- mod_simp %>% init(LT0 = 1000) %>%  mrgsim(end = 365 * 7) %>% as_tibble() %>% 
  mutate(MPP = MPP01 + MPP02) %>% 
  mutate(CMP = CMP01 + CMP02) %>%
  mutate(GMP = GMP01 + GMP02) %>%
  select(time, MPP, CMP, GMP, GMconc, RETconc, RBCconc, Bconc, Tconc, totalHb, HbinRBC)

tail(ss_simp, 1)

# new model after adding megakaryocytes
mod_new <- mread("model/human_ery_lymp_mk_myeloid") 

# simulate to steady state 
ss_new <- mod_new  %>% init(LT0 = 1000) %>%  mrgsim(end = 365 * 7) %>% as_tibble() %>% 
  mutate(MPP = MPP01 + MPP02) %>% 
  mutate(CMP = CMP01 + CMP02) %>%
  mutate(GMP = GMP01 + GMP02) %>% 
  mutate(BFUMK = BFUMK01 + BFUMK02) %>%
  mutate(CFUMK = CFUMK01 + CFUMK02) 

tail(ss_new %>% select(time, MPP, CMP, CFUMK, GMP, MK0, RETconc, RBCconc, PLTconc, TPO, Bconc, Tconc, totalHb, HbinRBC), 1)

tail(ss_new %>% select(Tconc, DN, DP, thymic_output, nT_lymphnodes, nT_lymphnodes, cd4rec0, cd8rec0, naiveTprof), 1)

# new model with neutrophil 
mod_neut <- mread("model/human_ery_lymp_mk_neutrophil") 

ss_neut <- mod_neut  %>% init(LT0 = 1275) %>% mrgsim(end = 365 * 7) %>% as_tibble() %>% 
  mutate(myelocyte = myelocyte01 + myelocyte02) %>% 
  mutate(MPP = MPP01 + MPP02) %>% 
  mutate(CMP = CMP01 + CMP02) %>%
  mutate(BFUMK = BFUMK01 + BFUMK02) %>%
  mutate(CFUMK = CFUMK01 + CFUMK02)

tail(ss_neut %>% select(time, MPP, CMP, CFUMK, GMP, MK0, RETconc, RBCconc, PLTconc, TPO, Bconc, Tconc, totalHb, HbinRBC), 1)

tail(ss_neut %>% select(Tconc, DN, DP, thymic_output, nT_lymphnodes, nT_lymphnodes, cd4rec0, cd8rec0, naiveTprof), 1)

     