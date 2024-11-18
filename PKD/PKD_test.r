rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


pkdmod <- mread("../SickleCell/mrgsolve/fullmodel2") %>% 
            # remove HSCT and adjust alpha globin synthesis rate
            param(totalCD34infused = 0, ksynalpha = 6e-7) %>%
            # change HbS parameters to HbA
            param(Kdalphabeta_sickle = 1e-3, HbS_saturation = 0.74) %>% 
            # adjust RBC lifespan
            param(tauRBCsickle = 27) %>%  
            init(LT0 = 1)

##------------------ adjust RBC lifespan only ------------------##

sim0 <- pkdmod  %>% 
  mrgsim(end = 30 * 18, delta = 1) %>% as_tibble() 

print(tail(sim0[c('totalHb', 'HbinRBC', 'RETconc', 'RBCconc','RET2RBCperc')], n=4))

##------------------ RBC lifespan scan ------------------##

pkdmod2 <- mread("../SickleCell/mrgsolve/fullmodel2") %>% 
  # remove HSCT and adjust alpha globin synthesis rate
  param(totalCD34infused = 0, ksynalpha = 6e-7) %>%
  # change HbS parameters to HbA
  param(Kdalphabeta_sickle = 1e-3, HbS_saturation = 0.74) %>% 
  # adjust RBC lifespan
#  param(tauRBCsickle = 27) %>%  
  init(LT0 = 1)

RBClifespan_scan <- function(tauRBC = 27)
{
  simtmp <- pkdmod2 %>% param(tauRBCsickle = tauRBC) %>%
    mrgsim(end = 30 * 18, delta = 1) %>% as_tibble() %>%
    select(RETconc, RBCconc, totalHb) %>% as.data.frame()
  
  tmp = tail(simtmp, n=1)
  
  tmp['tauRBC'] = tauRBC
  
  return(tmp)
  
}

RBClifespan <- c(10, 14, 18, 22, 26, 30, 34)

rbcscanresult <- RBClifespan_scan(120)

for( i in 1:length(RBClifespan) )
{
  tmp2 <- RBClifespan_scan(RBClifespan[i])
  
  rbcscanresult = rbind(rbcscanresult, tmp2)
  
  rm(tmp2)
}

rbclifespanplot1 <- ggplot(data = rbcscanresult) + 
  geom_point(aes(x = tauRBC, y = totalHb)) + 
  labs( x = 'RBC lifespan(days)', y = 'Hb (g/dL)') + theme_bw() 

rbclifespanplot2 <- ggplot(data = rbcscanresult) + 
  geom_point(aes(x = tauRBC, y = RBCconc)) + 
  labs( x = 'RBC lifespan(days)', y = 'RBC count (#/uL)') + theme_bw() 

rbclifespanplot3 <- ggplot(data = rbcscanresult) + 
  geom_point(aes(x = tauRBC, y = RETconc)) + 
  labs( x = 'RBC lifespan(days)', y = 'RET count (#/uL)') + theme_bw() 

rbclifespanplot4 <- ggplot(data = rbcscanresult) + 
  geom_point(aes(x = RBCconc, y = totalHb)) + 
  labs( x = 'RBC count (#/uL)', y = 'Hb (g/dL)') + theme_bw() 

grid.arrange(rbclifespanplot1, rbclifespanplot2, rbclifespanplot3, rbclifespanplot4, ncol = 2)

##------------------ fitting for RET lifespan ------------------##

# read in the observed RET data; 
RETcount <- read.csv('../SickleCell/data/PKD_PKac_RET.csv', header = T)

# show no obvious correlation between RET count and PK activity
plot(RETcount$PKactivity, RETcount$RETcount, ylab = 'RET count (*10^9 per L)')

# compute the median and mean RET count in PKD 
print( c( median(RETcount$RETcount), mean(RETcount$RETcount) ), digits = 4 )

# optimize RET lifespan
RETlifespan_opt <- function(RETlifespan)
{
  sim1 <- pkdmod  %>% param(tauRET = RETlifespan) %>%
    mrgsim(end = 30 * 18, delta = 1) %>% as_tibble() %>%
    select(time, RETconc, RBCconc, totalHb)
  
  tmp <- tail(sim1, n=1)
  
  diff <- abs(tmp$RETconc - 200 * 1000)/ 1000 # unit: # per uL
  
  return(diff)
}

# out <- optim(1, fn = RETlifespan_opt, lower = 0.1, upper = 3, method = "L-BFGS-B", maxit = 1000) # somehow this code does not return a value

RETlifespan_scan <- c(0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2)
RET_diff <- rep(NA, length(RETlifespan_scan))

for(i in 1: length(RETlifespan_scan))
{
  RET_diff[i] <- RETlifespan_opt(RETlifespan_scan[i])
}


sim1 <- pkdmod  %>% param(tauRET = 1) %>%
  mrgsim(end = 30 * 18, delta = 1) %>% as_tibble() 

print(tail(sim1[c('totalHb', 'RBCconc', 'RETconc', 'RET2RBCperc', 'HbinRBC')], n=1))

##------------------ sensitivity analysis, tauRET with reticulocyte & Hb ------------------##

RETlifespanscan <- function(RETlifespan)
{
  sim1 <- pkdmod  %>% param(tauRET = RETlifespan) %>%
    mrgsim(end = 30 * 18, delta = 1) %>% as_tibble() %>%
    select(RETconc, RBCconc, totalHb, HbinRBC) 
  
  tmp <- tail(sim1, n=1)
  
  return(data.frame(tmp))
}

scandata = RETlifespanscan(3)

scandata['tauRET'] = 3

RETlifespan_scan <- c(0.1, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2)

for(i in 1: length(RETlifespan_scan))
{
  tmp <- RETlifespanscan(RETlifespan_scan[i])
  
  tmp['tauRET'] = RETlifespan_scan[i]
  
  scandata = rbind(scandata, tmp)
  
  rm(tmp)
}

retrbcp <- ggplot(data = scandata) + 
  geom_point(aes(x = tauRET, y = RBCconc)) + 
  labs(x = 'RET lifespan (days)', y = 'RBC count (#/uL)') + theme_bw() 

ret_lifespan <- ggplot(data = scandata) + 
  geom_point(aes(x = tauRET, y = RETconc)) + 
  labs(y = 'RET count (#/uL)', x = 'RET lifespan(days)') + theme_bw() 

rethb1 <- ggplot(data = scandata) + 
  geom_point(aes(x = tauRET, y = totalHb)) + 
  labs(x = 'RET lifespan (days)', y = 'Hb (g/dL)') + theme_bw() 

rethb2 <- ggplot(data = scandata) + 
  geom_point(aes(x = RETconc, y = totalHb)) + 
  labs(x = 'RET count (#/uL)', y = 'Hb (g/dL)') + theme_bw() 

grid.arrange(retrbcp, ret_lifespan, rethb1, rethb2, ncol = 1)

ggplot(data = scandata) + 
  geom_point(aes(x = tauRET, y = RETconc + RBCconc)) + 
  labs(y = 'RET + RBC count (#/uL)', x = 'RET lifespan(days)') + theme_bw() 


ggplot(data = scandata) + 
  geom_point(aes(x = tauRET, y = HbinRBC)) + 
  labs(y = 'Hb in RBC (g/L)', x = 'RET lifespan(days)') + theme_bw() 


