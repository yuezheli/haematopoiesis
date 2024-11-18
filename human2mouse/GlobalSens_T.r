# Global sensitivity analysis of B cell dynamics after HSCT/ MPP/ CLP

rm(list = ls())
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(sensitivity)
library(PKPDmisc)
library(mrgsim.parallel)

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


##---------------------------------- prepare system before transplant ----------------------------------##
modo <- mread("erythrocytes_Hb_lymphoid_myeloid") %>% param(alphaCLP2proB = 3.1)

source('transplant_functions.r')

mod_mppt = global_mppt()


##---------------------------------- create sampling method for global sensitivity analysis ----------------------------------##

gen_samples <- function(n, l, which = names(l), 
                        factor = c(0.01,100), N = NULL) { # NEW ARGUMENT: N, the absolute sampling number per parameter
  
  vars <- tidyselect::vars_select(names(l), !!(enquo(which)))
  
  l <- as.list(l)[vars]
  
  l <- lapply(l, function(x) x*factor )
  
  if(is.numeric(N)) { # NEW
    n <- N*2  
  } else {
    n <- length(l)*n*2
  }
  
  df <- as.data.frame(l)
  
  len <- length(df)
  
  X <- matrix(ncol=len, nrow=n)
  
  colnames(X) <- names(df)
  
  Y <- X
  
  for(i in seq(len)){
    r <- exp(runif(n, log(df[1,i]), log(df[2,i])))
    X[,i] <- r
    r <- exp(runif(n, log(df[1,i]), log(df[2,i])))
    Y[,i] <- r
  }
  
  return(list(x1 = as.data.frame(X), x2 = as.data.frame(Y)))
}

# set up simulation parameters for sobol2007 function
simulationboot = 1000 # number of simulation time in sobol2007
sampleperparam = 200 # sampling size per parameter, for sobol2007
set.seed(88771) # seed adopted from the GitHub page mentioned at the beginning

# set up different variation range for different parameters
setmppclp <- c('alphaMPP2CLP', 'alphaCLP2DN', 'betaCLP', 'kCLP')
setT <- c('enter_lymph_cd4', 'enter_lymph_cd8', 'exit_lymph') # both are guessed parameters, thus with the most uncertainty

N <- sampleperparam * length(c(setmppclp,setT))

# set different variation in each group
sets <- list(
  list(which = setmppclp, factor = c(0.9, 1.1), N = N, l = param(mod_mppt)), 
  list(which = setT, factor = c(0.25,4), N = N, l = param(mod_mppt))
)

##---------------------------------- set up global sensitivity analysis using Sobel method (MPP transplant) ----------------------------------##

# MPP transplant AUC
batch_MPPT <- function(x) {
  x <- mutate(x, ID = row_number())
  mod_mppt %>% 
    idata_set(x) %>%
    mrgsim(obsonly = TRUE) %>% 
    na.omit() %>%
    group_by(ID) %>% 
    summarise(bexp = auc_partial(time, transT, range = c(0, 90))) %>% 
    pull(bexp)
}

# not running this chunk because too time consuming
if(TRUE) # note this section can be computational intensive; don't run unless you have a cluster
{
  # create sample space
  samps <- map(sets, do.call, what = gen_samples)
  samp <- map(c(1,2), ~ map_dfc(samps, .x))
  
  print(Sys.time())
  tvglobal = sobol2007(batch_MPPT, X1=samp[[1]], X2=samp[[2]], nboot=simulationboot) 
  plot(tvglobal)
  lines(c(0,8), c(0.05, 0.05), type="l", lty = 2)
  print(Sys.time())
  
  # save a minimal set of the Sobol analysis data
  # that is sufficient for visualization
  sobolx <- list(S = tvglobal$S , T = tvglobal$T)
  saveRDS(sobolx, file = "SobolData_T.rds")
}



if(FALSE)
{
  
  # load previously saved Sobol sensitivity analysis data (minimal dataset)
  x <- readRDS("SobolData_T.rds")
  ### plot
  globSens <- tibble(Parameter = c(setmppclp,setT),
                     main = x$S$original,
                     main_lo = x$S$'min. c.i.',
                     main_hi = x$S$'max. c.i.',
                     total = x$T$original,
                     total_lo = x$T$'min. c.i.',
                     total_hi = x$T$'max. c.i.') %>%
    gather(effect, Index, -Parameter, -main_lo, -main_hi, -total_lo, -total_hi) %>%
    mutate(lo = ifelse(effect == "main", main_lo, total_lo),
           hi = ifelse(effect == "main", main_hi, total_hi),
           Effect = factor(effect))
  
  fig_sens_glob <- ggplot(data=globSens, aes(x=Parameter, y=Index, group=Effect, col=Effect)) +
    geom_point(position = position_dodge(width=0.3)) +
    geom_errorbar(aes(ymax=hi, ymin=lo), position=position_dodge(width=0.3), width=0) +
    theme_bw() +
    geom_hline(aes(yintercept=0.05), lty=2) +
    theme(legend.title = element_text(size = 15), legend.text=element_text(size=12)) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 12))
  
  png('img/GlobalSens_MPPT_T.png', width = 800, height = 300)
  print(fig_sens_glob)
  dev.off()
}
