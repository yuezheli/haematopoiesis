## set up another set of depleted parameters for B cell model
depletedparam <- list(
  gamma = 0.38,
  delta_oe = 0.25,
  mu_i = 0.53,
  delta_i_t = 0.2,
  delta_i_re = 0.158, 
  mu_t = 0.06,
  delta_t = 0.005, 
  phi_BM = 0.95,
  mu_re = 0.04,
  phi_s = 0.019,
  epsilon_spl = 0.018
)

lapply(depletedparam, "*", 4) # scale from per 6h to per day
