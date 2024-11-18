# Goal: test T cell dynamics in mouse

rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


##----- STEADY STATES -----##

if(FALSE)
{
    mod <- mread("Tcell4") %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
    mod <- mread("Tcell4") %>% param(exit_lymph = 1e-4) %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
    mod <- mread("Tcell5") %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
    mod <- mread("Tcell5") %>% param(enter_lymph = 10) %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
    mod <- mread("Tcell5") %>% param(enter_lymph = 10, exit_lymph = 0.2) %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
    mod <- mread("Tcell6") %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
}


##----- sequester fraction scanning -----##
if (FALSE)
{
    sequester_fraction_scan <- function(exit_lymph_rate = 0.01)
        {
            mod <- mread("Tcell4") %>% param(exit_lymph = exit_lymph_rate) %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
            return(list(totalcount = mod$tcellbloodconc, tcell_seq = mod$tcell_lymph)) # return T cell count per uL
        }
    
    # this part look for T cell count
    tmp <- matrix(rep(NA,30), ncol = 3) %>% as.data.frame()
    colnames(tmp) <- c('exit_lymph_rate', 't_count', 'tcell_seq')
    
    tmp$exit_lymph_rate = 10^(seq(from = -9, to = 0, length.out = 10))
    
    for(i in 1:length(tmp$exit_lymph_rate))
    {
        tmp2 = sequester_fraction_scan(tmp$exit_lymph_rate[i])
        tmp$t_count[i] = tmp2$totalcount
        tmp$tcell_seq[i] = tmp2$tcell_seq
        rm(tmp2)
    }
    
    png('exit_rate_scan.png', width = 10, height = 4, units = 'in', res = 300)
    par(mfrow = c(1,2))
    plot(tmp$exit_lymph_rate, tmp$t_count, xlab = 'rate of T cell re-entering blood', ylab = 'T cell count (# per uL)', log="xy", pch = 15)
    plot(tmp$exit_lymph_rate, tmp$tcell_seq, xlab = 'rate of T cell re-entering blood', ylab = 'T cell count in other organs', log="xy", pch = 15)
    dev.off()
}


##----- testing scaling the model from in Million cells to Normal -----##
if (FALSE)
{
    mod <- mread("Tcell") %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% select(-ID)
    
    mod2 <-  mread("Tcell") %>% param(sigmaN = 0.02e6) %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% select(-ID)
    
    print(ggplot() + 
            geom_line(data = mod, aes(time, DN, color = 'DN'), size=1) + 
            geom_line(data = mod, aes(time, earlyDP, color = 'early DP'), size=1) +
            geom_line(data = mod, aes(time, lateDP, color = 'late DP'), size=1) +
            geom_line(data = mod2, aes(time, DN/1e6, color = 'DN (scaled)'), linetype = "dashed", size=3) + 
            geom_line(data = mod2, aes(time, earlyDP/1e6, color = 'early DP (scaled)'), linetype = "dashed", size=3) +
            geom_line(data = mod2, aes(time, lateDP/1e6, color = 'late DP (scaled)'), linetype = "dashed", size=3) +
            labs(y = 'cell count (million)', x = 'time (days)', color = '') +
            scale_y_continuous(trans = 'log10', limits = c(1e-6, 100)) + 
            theme_bw()
    )
}

##----- testing different existing and re-entering rate of naive T cells to periphery blood -----##
if(FALSE)
{
    lymphratesscan <- function(enter_lymph = 0.1, exit_lymph = 1.2)
    {
        mod <- mread("Tcell5") %>% 
            param(enter_lymph = enter_lymph, exit_lymph = exit_lymph) %>% 
            mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
        return(list(tcellbloodconc = mod$tcellbloodconc, 
                    tcell_b_ratio = mod$tcell_ratio_blood, 
                    tcell_l_ratio = mod$tcell_lymph_ratio, 
                    tcell_lymph_count = mod$tcell_lymph)) 
    }
    
    enter_lymph = 10^seq(from = -1, to = 2, by = 0.1)
    exit_lymph = 10^seq(from = -1, to = 1, by = 0.1)
    
    tcellbloodconc = matrix(rep(NA, length(enter_lymph) * length(exit_lymph)), nrow = length(enter_lymph))
    tcell_b_ratio = matrix(rep(NA, length(enter_lymph) * length(exit_lymph)), nrow = length(enter_lymph))
    tcell_l_ratio = matrix(rep(NA, length(enter_lymph) * length(exit_lymph)), nrow = length(enter_lymph))
    tcell_lymph_count = matrix(rep(NA, length(enter_lymph) * length(exit_lymph)), nrow = length(enter_lymph))
    
    
    for(i in 1:length(enter_lymph))
    {
        for(j in 1:length(exit_lymph))
        {
            tmp = lymphratesscan(enter_lymph[i], exit_lymph[j])
            tcellbloodconc[i,j] = tmp$tcellbloodconc
            tcell_b_ratio[i,j] = tmp$tcell_b_ratio
            tcell_l_ratio[i,j] = tmp$tcell_l_ratio
            tcell_lymph_count[i,j] = tmp$tcell_lymph_count
            rm(tmp)
        }
    }
    
    png('enter_exit_rate_double_scan_v5.png', width = 14, height = 4, units = 'in', res = 300)
    par(mfrow=c(1,3)) 
    contour(enter_lymph, exit_lymph, tcellbloodconc, xlab = 'rate of leaving blood',ylab = 'rate of re-entering blood', main = 'naive T cell count per uL blood')
    contour(enter_lymph, exit_lymph, tcell_lymph_count, xlab = 'rate of leaving blood',ylab = 'rate of re-entering blood', main = 'total naive T cell in other organs')
    contour(enter_lymph, exit_lymph, tcell_b_ratio, xlab = 'rate of leaving blood',ylab = 'rate of re-entering blood', main = 'CD4+:CD8+ T cell ratio, blood')
    #contour(enter_lymph, exit_lymph, tcell_l_ratio, xlab = 'rate of leaving blood',ylab = 'rate of re-entering blood', main = 'CD4+:CD8+ T cell ratio, lymph')
    dev.off()
}


##----- testing naive T cell entering secondary lymphoid system -----##
if(FALSE)
{
    enterrate_scan <- function(enter_lymph = 0.1)
    {
        mod <- mread("Tcell5") %>% param(enter_lymph = enter_lymph) %>% mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
        return(list(tcellbloodconc = mod$tcellbloodconc, tcell_lymph_count = mod$tcell_lymph)) 
    }
    
    tcellenterrate = 10^(seq(-1, 2, by = 0.1)) # rate range consistent with Sove et al., 2020
    tblood = rep(NA, length(tcellenterrate))
    tblymph = rep(NA, length(tcellenterrate))
    
    for(i in 1:length(tcellenterrate))
    {
        tmp = enterrate_scan(tcellenterrate[i])
        tblood[i] = tmp$tcellbloodconc
        tblymph[i] = tmp$tcell_lymph_count
        rm(tmp)
    }
    
    png('enter_rate_scan.png', width = 10, height = 4, units = 'in', res = 300)
    par(mfrow = c(1,2))
    plot(tcellenterrate, tblood, xlab = 'rate of T cell leaving blood (per day)', ylab = 'T cell count (# per uL)', log="xy", pch = 15)
    abline(v = 10)
    abline(h = 24000)
    plot(tcellenterrate, tblymph, xlab = 'rate of T cell leaving blood (per day)', ylab = 'T cell count in other organs', log="xy", pch = 15)
    abline(v = 10)
    abline(h = 2e8)
    dev.off()
}


##----- testing naive T cell death rate in secondary lymphoid organs -----##
if(FALSE)
{
    deathrate_scan <- function(cd4lymphdeath = 1/31, cd8lymphdeath = 1/72)
    {
        mod <- mread("Tcell6") %>% 
            param(death_cd4l = cd4lymphdeath, death_cd8l = cd8lymphdeath) %>% 
            mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
        return(list(tcellbloodconc = mod$tcellbloodconc, 
                    tcell_b_ratio = mod$tcell_ratio_blood, 
                    tcell_l_ratio = mod$tcell_lymph_ratio, 
                    tcell_lymph_count = mod$tcell_lymph)) 
    }
    
    cd4lymphdeath = c(seq(from = 0.001, to = 0.009, by = 0.001) ,seq(from = 0.01, to = 0.15, by = 0.01))
    cd8lymphdeath = c(seq(from = 0.001, to = 0.009, by = 0.001), seq(from = 0.01, to = 0.1, by = 0.01))
    
    tcellbloodconc = matrix(rep(NA, length(cd4lymphdeath) * length(cd8lymphdeath)), nrow = length(cd4lymphdeath))
    tcell_b_ratio = matrix(rep(NA, length(cd4lymphdeath) * length(cd8lymphdeath)), nrow = length(cd4lymphdeath))
    tcell_l_ratio = matrix(rep(NA, length(cd4lymphdeath) * length(cd8lymphdeath)), nrow = length(cd4lymphdeath))
    tcell_lymph_count = matrix(rep(NA, length(cd4lymphdeath) * length(cd8lymphdeath)), nrow = length(cd4lymphdeath))
    
    
    for(i in 1:length(cd4lymphdeath))
    {
        for(j in 1:length(cd8lymphdeath))
        {
            tmp = deathrate_scan(cd4lymphdeath[i], cd8lymphdeath[j])
            tcellbloodconc[i,j] = tmp$tcellbloodconc
            tcell_b_ratio[i,j] = tmp$tcell_b_ratio
            tcell_l_ratio[i,j] = tmp$tcell_l_ratio
            tcell_lymph_count[i,j] = tmp$tcell_lymph_count
            rm(tmp)
        }
    }
    
    png('lymph_death_scan.png', width = 10, height = 4, units = 'in', res = 300)
    par(mfrow = c(1,3))
    contour(cd4lymphdeath, cd8lymphdeath, tcellbloodconc, xlab = 'CD4+ death rate in other organ',ylab = 'CD8+ death rate in other organ', main = 'naive T cell count per uL blood')
    contour(cd4lymphdeath, cd8lymphdeath, tcell_b_ratio, xlab = 'CD4+ death rate in other organ',ylab = 'CD8+ death rate in other organ', main = 'CD4+:CD8+ T cell ratio, blood', nlevels = 25)
    #contour(cd4lymphdeath, cd8lymphdeath, tcell_l_ratio, xlab = 'CD4+ death rate in other organ',ylab = 'CD8+ death rate in other organ', main = 'CD4+:CD8+ T cell ratio, lymph',  nlevels = 25)
    contour(cd4lymphdeath, cd8lymphdeath, tcell_lymph_count, xlab = 'CD4+ death rate in other organ',ylab = 'CD8+ death rate in other organ', main = 'total naive T cell in other organs')
    dev.off()
}

##----- testing different rates for naive T cell to enter other organs -----##
if(FALSE)
{
    leaving_blood_scan <- function(enter_lymph_cd4 = 10, enter_lymph_cd8 = 10)
    {
        mod <- mread("Tcell6") %>% 
            param(enter_lymph_cd4 = enter_lymph_cd4, enter_lymph_cd8 = enter_lymph_cd4) %>% 
            #param(death_cd4l = 0, death_cd8l = 0) %>% 
            mrgsim(end = 150, delta = 1) %>% as_tibble() %>% tail(n=1)
        return(list(tcellbloodconc = mod$tcellbloodconc, 
                    tcell_b_ratio = mod$tcell_ratio_blood, 
                    tcell_l_ratio = mod$tcell_lymph_ratio, 
                    tcell_lymph_count = mod$tcell_lymph)) 
    }
    
    enter_lymph_cd4 = 10^seq(from = -1, to = 1.4, by = 0.2)
    enter_lymph_cd8 = 10^seq(from = -1, to = 1.2, by = 0.2)
    
    tcellbloodconc = matrix(rep(NA, length(enter_lymph_cd4) * length(enter_lymph_cd8)), nrow = length(enter_lymph_cd4))
    tcell_b_ratio = matrix(rep(NA, length(enter_lymph_cd4) * length(enter_lymph_cd8)), nrow = length(enter_lymph_cd4))
    tcell_l_ratio = matrix(rep(NA, length(enter_lymph_cd4) * length(enter_lymph_cd8)), nrow = length(enter_lymph_cd4))
    tcell_lymph_count = matrix(rep(NA, length(enter_lymph_cd4) * length(enter_lymph_cd8)), nrow = length(enter_lymph_cd4))
    
    
    for(i in 1:length(enter_lymph_cd4))
    {
        for(j in 1:length(enter_lymph_cd8))
        {
            tmp = leaving_blood_scan(enter_lymph_cd4[i], enter_lymph_cd8[j])
            tcellbloodconc[i,j] = tmp$tcellbloodconc
            tcell_b_ratio[i,j] = tmp$tcell_b_ratio
            tcell_l_ratio[i,j] = tmp$tcell_l_ratio
            tcell_lymph_count[i,j] = tmp$tcell_lymph_count
            rm(tmp)
        }
    }
    
    png('differential_cd4cd8_leavebrate_scan_nodeath.png', width = 14, height = 4, units = 'in', res = 300)
    #png('differential_cd4cd8_leavebrate_scan_withdeath.png', width = 10, height = 10, units = 'in', res = 300)
    par(mfrow = c(1,3))
    contour(enter_lymph_cd4, enter_lymph_cd8, tcellbloodconc, xlab = 'CD4+ leaving blood',ylab = 'CD8+ leaving blood', main = 'naive T cell count per uL blood')
    contour(enter_lymph_cd4, enter_lymph_cd8, tcell_b_ratio, xlab = 'CD4+ leaving blood',ylab = 'CD8+ leaving blood', main = 'CD4+:CD8+ T cell ratio, blood')
    #contour(enter_lymph_cd4, enter_lymph_cd8, tcell_l_ratio, xlab = 'CD4+ leaving blood',ylab = 'CD8+ leaving blood', main = 'CD4+:CD8+ T cell ratio, lymph')
    contour(enter_lymph_cd4, enter_lymph_cd8, tcell_lymph_count, xlab = 'CD4+ leaving blood',ylab = 'CD8+ leaving blood', main = 'total naive T cell in other organs')
    dev.off()
}

##----- testing different rates for naive T cell to re-enter blood stream -----##
if(FALSE)
{
    reenter_blood_scan <- function(exit_lymph_rate = 0.01)
    {
        mod <- mread("Tcell6") %>% param(exit_lymph = exit_lymph_rate) %>% mrgsim(end = 200, delta = 1) %>% as_tibble() %>% tail(n=1)
        return(list(bloodt = mod$tcellbloodconc, bloodtratio = mod$tcell_ratio_blood,lympht = mod$tcell_lymph)) # return T cell count per uL
    }
    
    # this part look for T cell count
    tmp <- matrix(rep(NA,48), ncol = 4) %>% as.data.frame()
    colnames(tmp) <- c('exit_lymph_rate', 'bloodt', 'bloodtratio', 'lympht')
    
    tmp$exit_lymph_rate = seq(from = 0.1, to = 1.2, length.out = 12)
    
    for(i in 1:length(tmp$exit_lymph_rate))
    {
        tmp2 = reenter_blood_scan(tmp$exit_lymph_rate[i])
        tmp$bloodt[i] = tmp2$bloodt
        tmp$lympht[i] = tmp2$lympht
        tmp$bloodtratio[i] = tmp2$bloodtratio
        rm(tmp2)
    }
    
    png('exit_rate_scan_v6.png', width = 10, height = 4, units = 'in', res = 300)
    par(mfrow = c(1,3))
    plot(tmp$exit_lymph_rate, tmp$bloodt, xlab = 'rate of T cell re-entering blood', ylab = 'T cell count (# per uL)', log="y", pch = 15)
    plot(tmp$exit_lymph_rate, tmp$bloodtratio, xlab = 'rate of T cell re-entering blood', ylab = 'T cell count in other organs', log="y", pch = 15)
    plot(tmp$exit_lymph_rate, tmp$lympht, xlab = 'rate of T cell re-entering blood', ylab = 'T cell count in other organs', log="y", pch = 15)
    dev.off()
}