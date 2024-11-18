# date: 10/8/24
# author: Yuezhe Li 
# purpose of this code: to test G-CSF in the model 

rm(list=ls())  #Clear out existing objects
gc()

library(here)
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(ggtext)

# pulse labeling data from Fig 2 of Price et al., 1996; https://pubmed.ncbi.nlm.nih.gov/8704192/
price1996 <- read.csv('data/Price1996.csv', header = TRUE) 

mod <- mread("model/human_mep_lymp_neut_gcsf") 

# simulation of pulse labeling 
pulse_label_exp <- function(gcsf_sc, mod = mod){
  pre_sim <- mod  %>% init(LT0 = 1000) %>% mrgsim(end = 100) %>% as_tibble() %>% select(-c(ID, time)) %>% tail(n = 1)
  mod_d1 <- mod %>% init(pre_sim) %>% init(gcsf_depot = gcsf_sc) %>% mrgsim(end = 1) %>% as_tibble() %>% select(-c(ID, time)) %>% tail(n = 1)
  mod_d2 <- mod %>% init(mod_d1) %>% init(gcsf_depot = gcsf_sc) %>% mrgsim(end = 1) %>% as_tibble() %>% select(-c(ID, time)) %>% tail(n = 1)
  mod_d3 <- mod %>% init(mod_d2) %>% init(gcsf_depot = gcsf_sc) %>% mrgsim(end = 1) %>% as_tibble() %>% select(-c(ID, time)) %>% tail(n = 1)
  mod_d4 <- mod %>% init(mod_d3) %>% init(gcsf_depot = gcsf_sc) %>% mrgsim(end = 1) %>% as_tibble() %>% select(-c(ID, time)) %>% tail(n = 1)
  mod_d5 <- mod %>% init(mod_d4) %>% init(gcsf_depot = gcsf_sc) %>% mrgsim(end = 1) %>% as_tibble() %>% select(-c(ID, time)) %>% tail(n = 1)
  
  mod_d6 <- mod %>% init(mod_d5) %>% init(gcsf_depot = gcsf_sc) %>%
    param(kCMP2GMP = 0) %>% 
    # assumed infused back cells is 10% or the original, estimated on plasma volume (Price et al., 1966)
    init(metamyelocyte = 0, bandcell = 0, neutrophil_pool = 0, neutrophil_blood = 0, 
         GMP = 1.1 * mod_d5$GMP, 
         myeloblast = 1.1 * mod_d5$myeloblast, 
         promyelocyte = 0, myelocyte01 = 0, myelocyte02 = 0)
  
  pulse_labeling <- mod_d6 %>% mrgsim(end = 8, delta = 0.1) %>% as_tibble() %>% 
    mutate(PMN = bandcell + neutrophil_blood) %>%
    mutate(relPMN = PMN/max(PMN) * 100) %>% 
    select(time, relPMN, bandcell, neutrophil_blood, Neuconc, gcsf)
  
  return(pulse_labeling)
}

pulse_label_0 <-  pulse_label_exp(0, mod) 
pulse_label_1point6e3 <-  pulse_label_exp(1.6E3, mod) 
pulse_label_16e3 <- pulse_label_exp(16E3, mod) 

pulse_labeling <- rbind(
  pulse_label_0 %>% mutate(G_CSF_ug = 0), 
  pulse_label_1point6e3 %>% mutate(G_CSF_ug = 30),
  pulse_label_16e3 %>% mutate(G_CSF_ug = 300)
)

p_pulse_labeling <- ggplot(data = pulse_labeling, aes(x = time, y = relPMN, col = as.factor(G_CSF_ug), group = as.factor(G_CSF_ug) )) + 
  geom_line(alpha = 0.7) + 
  geom_point(data = price1996, aes(x = Time_day, y = rel_Blood_PMN), alpha = 0.4) + 
  theme_bw() + theme(legend.position = "bottom") + 
  labs(x = "Time (day)", y = "Normalized cell count", col = "G-CSF dose (ug)", caption = "Generated by script/run_SteadyState_neutrophil.r") 

p_pulse_labeling

ggsave("deliv/figure/dynamics_neutrophil_pulselabeling.png", plot = p_pulse_labeling, width = 10, height = 10, units = "cm")

# simulate PK of G-CSF administration 
pre_sim <- mod  %>% init(LT0 = 1275) %>% mrgsim(end = 100) %>% as_tibble() %>% select(-c(ID, time)) %>% tail(n = 1)
print(pre_sim["gcsf"])
print(pre_sim["Neuconc"])

mod2 <- mod %>% init(pre_sim) %>% init(gcsf_depot = 16E3)  # 300ug G-CSF in pM

gcsf_dose <- mod2 %>% mrgsim(delta = 0.01, end = 6) %>% as_tibble() 

## observed data from de Haas et al., 1994; https://pubmed.ncbi.nlm.nih.gov/7524751/
deHaas1994_gcsf = data.frame(
  time_hr = c(0, 1, 4, 8, 12, 24, 1, 2, 4, 12, 24), 
  gcsf_ngmL = c(0.1, 9.8, 44.8, 12.5, 2.7, 0.2, 12.1, 19.8, 25.1, 2.4, 0.7)
) %>% 
  mutate(gcsf_pM = gcsf_ngmL*1E-6/18.8E3*1E12)

deHaas1994_neut = data.frame(
  time_hr = c(1, 2, 4, 8, 12, 24, 48, 72, 96, 144), 
  neut_1E9L = c(2.96, 5.16, 11.61, 18.98, 21.56, 17.26, 6.67, 4.09, 3.28, 2.85)
)

p_gcsf_neut <- ggplot() + 
  geom_line(data = gcsf_dose, aes(x = time * 24, y = Neuconc, col = "Simulated")) + 
  geom_point(data = deHaas1994_neut, aes(x = time_hr, y = neut_1E9L * 1E3, col = "de Haas et al., 1994"), alpha = 0.4) + 
  theme_bw() + theme(legend.position = "bottom") + 
  labs(x = "Time (hour)", y = "blood neutrophil count (#/uL)", col = "", caption = "Generated by script/run_neut_gcsf.r") 

p_gcsf_neut

ggsave("deliv/figure/dynamics_neutrophil.png", plot = p_gcsf_neut, width = 10, height = 10, units = "cm")

p_gcsf_dose <- ggplot() + 
  geom_line(data = gcsf_dose, aes(x = time * 24, y = gcsf, col = "Simulated")) + 
  geom_point(data = deHaas1994_gcsf, aes(x = time_hr, y = gcsf_pM, col = "de Haas et al., 1994"), alpha = 0.4) + 
  theme_bw() + theme(legend.position = "bottom") + 
  labs(x = "Time (hour)", y = "serum G-CSF (pM)", col = "", caption = "Generated by script/run_neut_gcsf.r") 

p_gcsf_dose

ggsave("deliv/figure/dynamics_gcsf.png", plot = p_gcsf_dose, width = 10, height = 10, units = "cm")
