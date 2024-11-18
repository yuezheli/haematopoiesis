rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

mod <- mread("fullmodel2") 

transplant_simul <- function(hscprof = 1/(2.5 * 7), modori = mod, ksyn = 6e-7, 
                             CD34infused = 2.8 * 1e8, ratioleft = 0.1)
{
  # default values in the model
  # hscprof: proliferation rate of LT-HSC, data from mouse
  # ksyn: synthesis rate of alpha globin; tuned; see other script for tuning data
  # CD34infused: total CD34+ infused, treatment data from Ribeil et al., 2017
  # ratioleft: myeloablative preconditioning; assuming only 10% of the endogenous progenitor cells left
  sim0 <- modori %>% param(totalCD34infused = 0) %>% param(rLT = hscprof) %>%
    init(LT0 = 1) %>% param(ksynalpha = ksyn) %>%
    mrgsim(delta = 1 ,end = 30 * 400) %>% as_tibble() 
  
  sim00 <- modori %>% param(totalCD34infused = 0) %>%
    init(LT0 = 1) %>% param(ksynalpha = ksyn) %>%
    mrgsim(delta = 1 ,end = 30 * 15) %>% as_tibble() 
  
  data = tail(sim00, n = 1)
  
  Data <- data[,-(1:2),drop=FALSE]  
  
  # steady state of sickle cell patients
  mod2 <- modori %>% init(Data) %>% param(rLT = hscprof) %>% param(ksynalpha = ksyn) 
  
  # myeloablative preconditioning
  # ratioleft = 1-0.9
  
  mod3 <- mod2 %>% init(LT0 = Data$LT0 * ratioleft) %>% 
    init(ST01 = Data$ST01 * ratioleft) %>% init(ST02 = Data$ST02 * ratioleft) %>%
    init(MPP01 = Data$MPP01 * ratioleft) %>% init(MPP02 = Data$MPP02 * ratioleft) %>%
    init(CMP01 = Data$CMP01 * ratioleft) %>% init(CMP02 = Data$CMP02 * ratioleft)
  
  
  # simulation after the patient is given 2.8e8 (5.6e6 per kg)
  sim <- mod3 %>% param(totalCD34infused = CD34infused) %>% 
                  mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble() 
  
  preinfusion <- sim0 %>% select(time, RBCconc, totalHbA, totalHbS, totalHb, HbinRBC)
  postinfusion <- sim %>% select(time, RBCconc, totalHbA, totalHbS, totalHb, HbinRBC)
  
  return(list(preinfusion = preinfusion, postinfusion = postinfusion))
}

# scanning for different proliferation of HSC
# values from Catlin et al., 2011
# https://ashpublications.org/blood/article/117/17/4460/20856/The-replication-rate-of-human-hematopoietic-stem
# 4 values: human, NHP, cat, mouse
hscprof <- c( 1/(40*7), 1/(25*7), 1/(7*8.3), 1/(7*2.5), 1/(7*9) )

human <- transplant_simul(hscprof = hscprof[1])
nhp <- transplant_simul(hscprof = hscprof[2])
cat <- transplant_simul(hscprof = hscprof[3])
mouse <- transplant_simul(hscprof = hscprof[4])
adjusted_human <- transplant_simul(hscprof = hscprof[5])

# compare simulation result before infusion
pre_RBC <- ggplot() + ggtitle('before infusion RBC') +  
  geom_line(data = human$preinfusion, aes(x = time/30, y = RBCconc, color = "human")) + 
  geom_line(data = nhp$preinfusion, aes(x = time/30, y = RBCconc, color = "NHP")) +
  geom_line(data = cat$preinfusion, aes(x = time/30, y = RBCconc, color = "cat")) +
  geom_line(data = mouse$preinfusion, aes(x = time/30, y = RBCconc, color = "mouse")) +
  labs(y = 'RBC cell count (#/uL)', x = 'time (months)', color = '') + theme(legend.position = "bottom") + theme_bw()

pre_Hb <- ggplot() + ggtitle('before infusion Hb') + 
  geom_line(data = human$preinfusion, aes(x = time/30, y = totalHb, color = "human")) + 
  geom_line(data = nhp$preinfusion, aes(x = time/30, y = totalHb, color = "NHP")) +
  geom_line(data = cat$preinfusion, aes(x = time/30, y = totalHb, color = "cat")) +
  geom_line(data = mouse$preinfusion, aes(x = time/30, y = totalHb, color = "mouse")) +
  labs(y = 'blood Hb (g/dL)', x = 'time (months)', color = '') + theme(legend.position = "bottom") + theme_bw()

# compare simulation result after infusion
post_RBC <- ggplot() + ggtitle('after infusion RBC') + 
  geom_line(data = human$postinfusion, aes(x = time/30, y = RBCconc, color = "human")) + 
  geom_line(data = nhp$postinfusion, aes(x = time/30, y = RBCconc, color = "NHP")) +
  geom_line(data = cat$postinfusion, aes(x = time/30, y = RBCconc, color = "cat")) +
  geom_line(data = mouse$postinfusion, aes(x = time/30, y = RBCconc, color = "mouse")) +
  geom_line(data = adjusted_human$postinfusion, aes(x = time/30, y = RBCconc, color = "human, adjusted")) +
  labs(y = 'RBC cell count (#/uL)', x = 'time (months)', color = '') + theme(legend.position = "bottom") + theme_bw()

post_Hb <- ggplot() + ggtitle('after infusion Hb') + 
  geom_line(data = human$postinfusion, aes(x = time/30, y = totalHb, color = "human")) + 
  geom_line(data = nhp$postinfusion, aes(x = time/30, y = totalHb, color = "NHP")) +
  geom_line(data = cat$postinfusion, aes(x = time/30, y = totalHb, color = "cat")) +
  geom_line(data = mouse$postinfusion, aes(x = time/30, y = totalHb, color = "mouse")) +
  geom_line(data = adjusted_human$postinfusion, aes(x = time/30, y = totalHb, color = "human, adjusted")) +
  labs(y = 'blood Hb (g/dL)', x = 'time (months)', color = '') + theme(legend.position = "bottom") + theme_bw()

# plot all the figures together
grid.arrange(pre_RBC, pre_Hb, post_RBC, post_Hb, ncol = 2)
