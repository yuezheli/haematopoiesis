# the goal of this script is to find parameters that make the simulation result more realistic
# this script use a healthy person for testing
# note that in this case, HbS is actually HbA

rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)

mod <- mread("fullmodel2") 

# sensitivity analysis of alpha globin synthesis rate
# for parameter tuning

synalphaimpact <- function(ksyn, Kdalphabeta = 1e-3, Hb_saturation = 0.74, tauRBC = 120)
{
  # all default parameters are for healthy people
  # if simulation for SCD is wanted, set Kdalphabeta_sickle = 1e-2, Hb_saturation = 0.68, tauRBC = 12
  
  sim1 <- mod %>% param(totalCD34infused = 0) %>%   
    param(Kdalphabeta_sickle = Kdalphabeta) %>% param(HbS_saturation = Hb_saturation) %>% param(tauRBCsickle = tauRBC) %>% # change HbS to HbA parameters
    init(LT0 = 1) %>%
    param(ksynalpha = ksyn) %>%  #decrease hemoglobin synthesis rate
    mrgsim(end = 30 * 18) %>% as_tibble() 
  
  # return the data at steady state
  data = tail(sim1, n = 1)
  
  output = list( 
              RBCconc = data$RBCconc,
              RET = data$RET,
              RBC = data$RBC,
              Hb = data$totalHb,
              vO2 = data$bloodO2,
              HbinRBC = data$HbinRBC)
  
  tmp = c(data$RBCconc, data$RET, data$totalHb, data$bloodO2, data$HbinRBC)
  return(tmp)
}

# check through a list of variables
ksyn = c(0.1e-6, 0.3e-6, 0.5e-6, 0.6e-6, 0.7e-6, 1e-6, 1.5e-6, 2e-6)

data = matrix(rep(NA, length(ksyn) * 5), nrow = length(ksyn))

for(i in 1:length(ksyn))
{
  data[i,] = synalphaimpact(ksyn[i])
}

tmp = data.frame(data)
names(tmp) = c('RBCconc', 'RET', 'Hb', 'vO2','HBinRBC')

tmp['ksyn'] = ksyn


pRBCconc <- ggplot(data = tmp, aes(x = ksyn, y = RBCconc)) + 
  geom_line(linetype = "dashed") + 
  geom_point(shape=23, fill="red", color="red", size=3) + 
  labs(y = 'RBC count (#/uL)', x = '\u03b1-globin synthesis rate (nmol.day-1.cell-1)') + theme_bw()

pHbconc <-  ggplot(data = tmp, aes(x = ksyn, y = Hb)) + 
  geom_line(linetype = "dashed") + 
  geom_point(shape=23, fill="red", color="red", size=3) + 
  labs(y = 'Hb conc in blood (g/dL)', x = '\u03b1-globin synthesis rate (nmol.day-1.cell-1)') + theme_bw()

pvO2 <-  ggplot(data = tmp, aes(x = ksyn, y = vO2)) + 
  geom_line(linetype = "dashed") + 
  geom_point(shape=23, fill="red", color="red", size=3) + 
  labs(y = 'blood O2 (ml/dL)', x = '\u03b1-globin synthesis rate (nmol.day-1.cell-1)') + theme_bw()

pHbconc2 <-  ggplot(data = tmp, aes(x = ksyn, y = HBinRBC)) + 
  geom_line(linetype = "dashed") + 
  geom_point(shape=23, fill="red", color="red", size=3) + 
  labs(y = 'Hb conc in RBC (g/L)', x = '\u03b1-globin synthesis rate (nmol.day-1.cell-1)') + theme_bw()

# grid.arrange(pRBCconc, pHbconc, pvO2, pHbconc2, ncol = 2)
grid.arrange(pRBCconc, pHbconc, pHbconc2, ncol = 3)
