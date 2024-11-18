rm(list = ls())

library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script
theme_set(theme_bw())


mod <- mread("Bcell_mus") %>%  mrgsim(end = 600, delta = 1) %>% as_tibble() %>% select(-ID)
mod2 <- mread("Bcell_mus2") %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% select(-ID)


print(ggplot() + 
        geom_line(data = mod, aes(time/4, BcellBM, color = 'B cells in bone marrow')) + 
        geom_line(data = mod, aes(time/4, Bcell_spleen, color = 'B cells in spleen')) +
        geom_line(data = mod, aes(time/4, BMrec, color = 'circulating B cells')) +
              geom_line(data = mod2, aes(time, BcellBM, color = 'B cells in bone marrow (adjusted)')) + 
              geom_line(data = mod2, aes(time, Bcell_spleen, color = 'B cells in spleen (adjusted)')) +
              geom_line(data = mod2, aes(time, BMrec, color = 'circulating B cells (adjusted)')) +          
        labs(y = 'cell count ', x = 'time (days)', color = '') +
        scale_y_continuous(trans = 'log10') + 
        theme_bw()
)

tmp = tail(mod, n=1)
tmp2 = tail(mod2, n=1)


## parameter scan, influx into preproB cells
mread("Bcell_mus2") %>%
        idata_set( expand.idata(s = c(0.5e6, 1e6, 1.2e6, 1.5e6, 2e6)) ) %>% 
        mrgsim(end = 150, delta = 1) %>% 
        plot( BMrec + Bcell_spleen + BcellBM ~ time )

## compare steady state of depleted and control case
source('DepletionParam.r')
mod3 <- mread("Bcell_mus") %>% param(depletedparam) %>%  mrgsim(end = 600, delta = 1) %>% as_tibble() %>% select(-ID)


print(ggplot() + 
              geom_line(data = mod, aes(time, BcellBM, color = 'B cells in bone marrow')) + 
              geom_line(data = mod, aes(time, Bcell_spleen, color = 'B cells in spleen')) +
              geom_line(data = mod, aes(time, BMrec, color = 'circulating B cells')) +
              geom_line(data = mod3, aes(time, BcellBM, color = 'B cells in bone marrow (depleted)')) + 
              geom_line(data = mod3, aes(time, Bcell_spleen, color = 'B cells in spleen (depleted)')) +
              geom_line(data = mod3, aes(time, BMrec, color = 'circulating B cells (depleted)')) +          
              labs(y = 'cell count ', x = 'time (days)', color = '') +
              scale_y_continuous(trans = 'log10') + 
              theme_bw()
)
