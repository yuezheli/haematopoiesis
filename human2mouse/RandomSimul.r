rm(list = ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(mrggsave)
theme_set(theme_bw())

source('Depleted_B.r')
source('LoadObs.r')

modo <- mread("erythrocytes_Hb_lymphoid_myeloid") %>% param(depletedparam) %>% param(alpha_muN = 0.3, pN = 0.3, delta = 0.1)

N = 100
variation1 = 0.2
variation2 = 0.4
# generate random number 
x <- rnorm(N, mean=0.3, sd=0.3 * variation1)
y <- rnorm(N, mean=0.2, sd=0.2 * variation1)
z <- rnorm(N, mean=0.1, sd=0.1 * variation2)
clp <- rnorm(N, mean=1e4, sd=1e4 * variation1)

#hist(x)
#hist(y)
hist(z)
#hist(CLPtransplant)

sim_pre = modo %>% init(LT0 = 1)  %>%  mrgsim(end = 1500, delta = 0.5) %>% as_tibble() %>% tail(n=1) 
sim_all <- modo %>% init(sim_pre * 0.9) %>% init(CLP1 = 10000) %>% mrgsim(end = 150) %>% as_tibble() %>% select(time, transT)
rm(sim_pre)


for(i in 1:N )
{
  tmp_alpha = max(x[i], 0)
  tmp_pN = max(y[i], 0)
  tmp_delta = max(z[i], 0)
  tmpmod = modo %>% param(alpha_muN = tmp_alpha, pN = tmp_pN, delta = tmp_delta)
  
  sim_pre = tmpmod %>% init(LT0 = 1)  %>%  mrgsim(end = 1500, delta = 0.5) %>% as_tibble() %>% tail(n=1) 
  
  sim <- modo %>% init(sim_pre * 0.9) %>% init(CLP1 = clp[i]) %>%
    mrgsim(end = 150) %>% as_tibble() %>% select(time, transT)
  
  sim_all <- rbind(sim_all, sim)
  
  rm(tmp_alpha, tmp_pN, tmpmod, sim_pre, sim)
}



result = sim_all %>%group_by(time)    %>% 
  summarize(lo=quantile(transT, 0.1, na.rm = TRUE), hi=quantile(transT, 0.9, na.rm = TRUE), med = quantile(transT, 0.5, na.rm = TRUE) )

print(
ggplot(data = result) + 
  geom_line(aes(x = time, y = med * 100, color = 'median ')) + 
  geom_ribbon(aes(x = time,ymin = lo * 100, ymax = hi * 100), alpha = 0.2) + 
  geom_point(data = clpt_T, aes(x = time, y = cell_fraction, color = 'obs')) + 
  geom_errorbar(data=clpt_T, aes(x=time, ymin=T_lower,ymax=T_upper), alpha = 0.5) +
  labs(y = 'trans T percentage (%)', x = 'time (days)', color = '') + 
  ggtitle("CLP transplant")
)

ggsave('img/T_adjusted2.png',width = 8, height = 3, units = 'in', dpi = 300)


