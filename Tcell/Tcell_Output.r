# Goal: test T cell dynamics in mouse

rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


# read in data from Scollay et al., 1980
bloodt <- read.csv('scollary1980.csv', header = TRUE) %>% filter(organ == "blood")

bloodt['cell_count'] = bloodt['perorgan'] /1e5 * 1e7 # scale to the total cell count in blood 

#------------ parameter scan, t cell entering and leaving periphery blood ------------#
if(FALSE)
{
        mod <- mread("Tcell5") %>% param(enter_lymph = 10)
        sim_pre <- mrgsim(mod, end = 200, delta = 1) %>% as_tibble() %>% tail(n=1)
        mod_default <- mod %>% init(sim_pre * 0.25) %>% init(cd4rec = 1, cd8rec = 1, cd4lym = 0, cd8lym = 0) 
        
        exidata <- expand.idata(exit_lymph = c(0.2, 0.6, 1, 1.2))
        sim_post <- mod_default %>% idata_set(exidata) %>% mrgsim(end = 1/3, delta = 0.01) %>% 
                    select(ID, time, cd4rec, cd8rec) %>% group_by(ID)
        
        tmp = ggplot() + 
                geom_line(data = sim_post %>% filter(ID == 1), aes(time * 24, cd4rec + cd8rec, color = 're-entering rate = 0.2 per day')) +
                geom_line(data = sim_post %>% filter(ID == 2), aes(time * 24, cd4rec + cd8rec, color = 're-entering rate = 0.6 per day')) +
                geom_line(data = sim_post %>% filter(ID == 3), aes(time * 24, cd4rec + cd8rec, color = 're-entering rate = 1 per day')) +
                geom_line(data = sim_post %>% filter(ID == 4), aes(time * 24, cd4rec + cd8rec, color = 're-entering rate = 1.2 per day')) +
                geom_point(data = bloodt, aes(time_hour, cell_count, color = 'blood T obs')) + 
                labs(y = 'cell count', x = 'time (hour)', color = '') +
                scale_y_continuous(trans = 'log10', limits = c(100, 1e6)) + theme_bw()
        
        png('FITClabeling_4.png', width = 10, height = 4, units = 'in', res = 300)
        print(tmp)
        dev.off()
} 


#------------ testing the dynamics of simulation ------------#
if(FALSE)
{
        mod <- mread("Tcell4") 

        sim_pre <- mrgsim(mod, end = 200, delta = 1) %>% as_tibble() %>% tail(n=1)
        sim_post <- mod %>% init(sim_pre * 0.25) %>%
                init(cd4rec = 1, cd8rec = 1, cd4lym = 0, cd8lym = 0) %>%
                mrgsim(end = 1/3, delta = 0.01) %>% as_tibble() 
        
        tmp = ggplot() + 
                geom_line(data = sim_post, aes(time * 24, cd4rec + cd8rec, color = 'blood T simul')) + 
                geom_point(data = bloodt, aes(time_hour, cell_count, color = 'blood T obs')) + 
                labs(y = 'cell count', x = 'time (hour)', color = '') +
                scale_y_continuous(trans = 'log10', limits = c(10, 1e6)) + theme_bw()
        
        png('FITClabeling_1.png', width = 10, height = 4, units = 'in', res = 300) 
        print(tmp)
        dev.off()
}


if(FALSE)
{
        mod <- mread("Tcell5") %>% param(enter_lymph = 10)

        sim_pre <- mrgsim(mod, end = 200, delta = 1) %>% as_tibble() %>% tail(n=1)
        sim_post <- mod %>% init(sim_pre * 0.25) %>%
                            init(cd4rec = 1, cd8rec = 1, cd4lym = 0, cd8lym = 0) %>%
                    mrgsim(end = 1/3, delta = 0.01) %>% as_tibble() 

        tmp = ggplot() + 
                geom_line(data = sim_post, aes(time * 24, cd4rec + cd8rec, color = 'blood T simul')) + 
                geom_point(data = bloodt, aes(time_hour, cell_count, color = 'blood T obs')) + 
                labs(y = 'cell count', x = 'time (hour)', color = '') +
                scale_y_continuous(trans = 'log10', limits = c(10, 1e6)) + theme_bw()
        
        #png('FITClabeling_2.png', width = 10, height = 4, units = 'in', res = 300) # without parameter change
        png('FITClabeling_2_2.png', width = 10, height = 4, units = 'in', res = 300) # with parameter adjustment
        print(tmp)
        dev.off()
}

#------------ test how different rates for cell leaving blood changes the dynamics ------------#
if(FALSE)
{
        enter_rate_scan_ratio <- function(cd4_leaving_blood_rate = 10, cd4cd8ratio = 1)
        {
                mod <- mread("Tcell6") %>% param(enter_lymph_cd4 = cd4_leaving_blood_rate, 
                                                 enter_lymph_cd8 = cd4_leaving_blood_rate/cd4cd8ratio)
                
                sim_pre <- mrgsim(mod, end = 200, delta = 1) %>% as_tibble() %>% tail(n=1)
                sim_post <- mod %>% init(sim_pre * 0.25) %>%
                                init(cd4rec = 1, cd8rec = 1, cd4lym = 0, cd8lym = 0) %>%
                                mrgsim(end = 1/3, delta = 0.01) %>% as_tibble() 
                return(sim_post)
        }
        
        cd4_leaving_blood = 10
        
        simpost1 <- enter_rate_scan_ratio(cd4_leaving_blood_rate = cd4_leaving_blood, cd4cd8ratio = 1)
        simpost2 <- enter_rate_scan_ratio(cd4_leaving_blood_rate = cd4_leaving_blood, cd4cd8ratio = 2)
        simpost5 <- enter_rate_scan_ratio(cd4_leaving_blood_rate = cd4_leaving_blood, cd4cd8ratio = 5)
        simpost10 <- enter_rate_scan_ratio(cd4_leaving_blood_rate = cd4_leaving_blood, cd4cd8ratio = 10)
        simpost_point1 <- enter_rate_scan_ratio(cd4_leaving_blood_rate = cd4_leaving_blood, cd4cd8ratio = 0.1)
        
        tmp = ggplot() + 
                geom_line(data = simpost1, aes(time * 24, cd4rec + cd8rec, color = 'ratio = 1')) +
                geom_line(data = simpost2, aes(time * 24, cd4rec + cd8rec, color = 'ratio = 2')) +
                geom_line(data = simpost5, aes(time * 24, cd4rec + cd8rec, color = 'ratio = 5')) +
                geom_line(data = simpost10, aes(time * 24, cd4rec + cd8rec, color = 'ratio = 10')) +
                geom_line(data = simpost_point1, aes(time * 24, cd4rec + cd8rec, color = 'ratio = 0.1')) +
                geom_point(data = bloodt, aes(time_hour, cell_count, color = 'blood T obs')) + 
                labs(y = 'cell count', x = 'time (hour)', color = '') +
                scale_y_continuous(limits = c(1e3, 1e5)) + theme_bw()
        
        print(tmp)
        
}

#------------ test v6 for output ------------#
if(FALSE)
{
        mod <- mread("Tcell6") 
        # pre-labeling steady state
        sim_pre <- mrgsim(mod, end = 200, delta = 1) %>% as_tibble() %>% tail(n=1)
        # labeling process; assuming efficiency = 25% as reported in Scollay et al., 1980, Table 1
        labeling_eff = 0.25
        sim_pre_prime = sim_pre * (1-labeling_eff)
        sim_pre_prime$N01 = sim_pre$N00 * labeling_eff
        sim_pre_prime$N11 = sim_pre$N10 * labeling_eff
        sim_pre_prime$N21 = sim_pre$N20 * labeling_eff
        sim_pre_prime$N31 = sim_pre$N30 * labeling_eff
        sim_pre_prime$N41 = sim_pre$N40 * labeling_eff
        sim_pre_prime$P01 = sim_pre$P00 * labeling_eff
        sim_pre_prime$P11 = sim_pre$P10 * labeling_eff
        sim_pre_prime$P21 = sim_pre$P20 * labeling_eff
        sim_pre_prime$P31 = sim_pre$P30 * labeling_eff
        sim_pre_prime$P41 = sim_pre$P40 * labeling_eff
        sim_pre_prime$P51 = sim_pre$P50 * labeling_eff
        sim_pre_prime$P61 = sim_pre$P60 * labeling_eff
        sim_pre_prime$P71 = sim_pre$P70 * labeling_eff
        sim_pre_prime$S401 = sim_pre$S400 * labeling_eff
        sim_pre_prime$S411 = sim_pre$S410 * labeling_eff
        sim_pre_prime$S421 = sim_pre$S420 * labeling_eff
        sim_pre_prime$S801 = sim_pre$S800 * labeling_eff
        sim_pre_prime$S811 = sim_pre$S810 * labeling_eff
        sim_pre_prime$S821 = sim_pre$S820 * labeling_eff
        # post-labeling simulation
        sim_post <- mod %>% init(sim_pre_prime) %>% mrgsim(end = 0.3, delta = 0.005) %>% as_tibble() %>% filter(time > 0)
        
        tmp = ggplot() + 
                geom_line(data = sim_post, aes(time * 24, cd4rec1 + cd8rec1, color = 'simul ')) +
                geom_point(data = bloodt, aes(time_hour, cell_count, color = 'blood T obs')) + 
                labs(y = 'cell count', x = 'time (hour)', color = '') +
                scale_y_continuous(limits = c(1e3, 1e5), trans = 'log10') + theme_bw()
        
        png('FITClabeling_5.png', width = 10, height = 4, units = 'in', res = 300) 
        print(tmp)
        dev.off()
        
}
