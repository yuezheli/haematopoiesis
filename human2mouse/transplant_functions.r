##----------------- HSC transplant -----------------##


conditioning_proliferation_hsct <- function(progenitor_conditioning_strength = 0.1, hscprof = 1/(2.5*7), LTnumber = 200, b_conditioning_strength = 0.999, t_conditioning_strength = 0.9, gm_conditioning_strength = 0.05, mod = modo, returnSteadyState = FALSE, simendtime = 150)
{

  sim_pre <- mod %>% init(LT0 = 1)  %>%  mrgsim(end = 30 * 50, delta = 0.5) %>% 
           as_tibble() %>% tail(n=1) %>% select(-c(ID, time)) 

  if (returnSteadyState)
  {
    return(sim_pre)
  }
  else
  {
      # steady state of sickle cell patients
      mod2 <- mod %>% init(sim_pre)

      # myeloablative preconditioning
      ratioleft = 1-progenitor_conditioning_strength
      
      mod3 <- mod2 %>% init(LT0 = sim_pre$LT0 * ratioleft) %>% 
                       init(ST01 = sim_pre$ST01 * ratioleft) %>% init(ST02 = sim_pre$ST02 * ratioleft) %>%
                       init(MPP01 = sim_pre$MPP01 * ratioleft) %>% init(MPP02 = sim_pre$MPP02 * ratioleft) %>%
                       init(CMP01 = sim_pre$CMP01 * ratioleft) %>% init(CMP02 = sim_pre$CMP02 * ratioleft) %>%
                       init(BFUE01 = sim_pre$BFUE01 * ratioleft) %>% init(BFUE02 = sim_pre$BFUE02 * ratioleft) %>%
                       init(CFUE01 = sim_pre$CFUE01 * ratioleft) %>% init(CFUE02 = sim_pre$CFUE02 * ratioleft)

      if("CLP0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(CLP0 = sim_pre$CLP0 * ratioleft) 
      } 
      
      if(("Boe0" %in% colnames(sim_pre)) & ("Boe1" %in% colnames(sim_pre)))
      {
        mod3 <-mod3 %>% init(Boe0 = sim_pre$Boe0 * ratioleft) %>% init(Bi0 = sim_pre$Bi0 * ratioleft) %>% 
                        init(Bt0 = sim_pre$Bt0 * (1-b_conditioning_strength)) %>% init(BMrec0 = sim_pre$BMrec0 * (1-b_conditioning_strength)) %>% 
                        init(BMspl0 = sim_pre$BMspl0 * (1-b_conditioning_strength))
      } 

      if("N00" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(N00 = sim_pre$N00 * (1-t_conditioning_strength)) %>% init(N10 = sim_pre$N10 * (1-t_conditioning_strength)) %>% init(N20 = sim_pre$N20 * (1-t_conditioning_strength)) %>% init(N30 = sim_pre$N30 * (1-t_conditioning_strength)) %>% init(N40 = sim_pre$N40 * (1-t_conditioning_strength)) %>%
                        init(P00 = sim_pre$P00 * (1-t_conditioning_strength)) %>% init(P10 = sim_pre$P10 * (1-t_conditioning_strength)) %>% init(P20 = sim_pre$P20 * (1-t_conditioning_strength)) %>% init(P30 = sim_pre$P30 * (1-t_conditioning_strength)) %>% init(P40 = sim_pre$P40 * (1-t_conditioning_strength)) %>%
                        init(P50 = sim_pre$P50 * (1-t_conditioning_strength)) %>% init(P60 = sim_pre$P60 * (1-t_conditioning_strength)) %>% init(P70 = sim_pre$P70 * (1-t_conditioning_strength)) %>% 
                        init(S800 = sim_pre$S800 * (1-t_conditioning_strength)) %>% init(S810 = sim_pre$S810 * (1-t_conditioning_strength)) %>% init(S820 = sim_pre$S820 * (1-t_conditioning_strength)) %>% 
                        init(S400 = sim_pre$S400 * (1-t_conditioning_strength)) %>% init(S410 = sim_pre$S410 * (1-t_conditioning_strength)) %>% init(S420 = sim_pre$S420 * (1-t_conditioning_strength)) %>% 
                        init(cd4rec0 = sim_pre$cd4rec0 * (1-t_conditioning_strength)) %>% init(cd8rec0 = sim_pre$cd8rec0 * (1-t_conditioning_strength)) %>% init(cd4lym0 = sim_pre$cd4lym0 * (1-t_conditioning_strength)) %>% init(cd8lym0 = sim_pre$cd8lym0 * (1-t_conditioning_strength))
      } 

      if("GMP01" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GMP01 = sim_pre$GMP01 * ratioleft) %>% init(GMP02 = sim_pre$GMP02 * ratioleft) 
      }

      if("GM0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GM0 = sim_pre$GM0 * (1-gm_conditioning_strength)) 
      }

      sim_hsct <- mod3 %>% init(LT1 = LTnumber) %>% param(rLT = hscprof) %>% mrgsim(end = simendtime) %>% as_tibble() %>% select(-ID)


      return( sim_hsct )
  }
}


##----------------- MPP transplant -----------------##

conditioning_mppt <- function(conditioning_strength = 0.1, latempp = 0.5, MPPnumber = 1000, b_conditioning_strength = 0.999, t_conditioning_strength = 0.9, gm_conditioning_strength = 0.05, mod = modo, returnSteadyState = FALSE, simendtime = 150)
{
  
  sim_pre <- mod %>% init(LT0 = 1)  %>%  mrgsim(end = 30 * 50, delta = 0.5) %>% 
    as_tibble() %>% tail(n=1) %>% select(-c(ID, time)) 
  
  if (returnSteadyState)
  {
    return(sim_pre)
  }
  else
  {
    # steady state of sickle cell patients
    mod2 <- mod %>% init(sim_pre)
    
    # myeloablative preconditioning
    ratioleft = 1-conditioning_strength
    
    mod3 <- mod2 %>% init(LT0 = sim_pre$LT0 * ratioleft) %>% 
      init(ST01 = sim_pre$ST01 * ratioleft) %>% init(ST02 = sim_pre$ST02 * ratioleft) %>%
      init(MPP01 = sim_pre$MPP01 * ratioleft) %>% init(MPP02 = sim_pre$MPP02 * ratioleft) %>%
      init(CMP01 = sim_pre$CMP01 * ratioleft) %>% init(CMP02 = sim_pre$CMP02 * ratioleft) %>%
      init(BFUE01 = sim_pre$BFUE01 * ratioleft) %>% init(BFUE02 = sim_pre$BFUE02 * ratioleft) %>%
      init(CFUE01 = sim_pre$CFUE01 * ratioleft) %>% init(CFUE02 = sim_pre$CFUE02 * ratioleft)
    
    if("CLP0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(CLP0 = sim_pre$CLP0 * ratioleft) 
      } 
      
      if(("Boe0" %in% colnames(sim_pre)) & ("Boe1" %in% colnames(sim_pre)))
      {
        mod3 <-mod3 %>% init(Boe0 = sim_pre$Boe0 * ratioleft) %>% init(Bi0 = sim_pre$Bi0 * ratioleft) %>% 
                        init(Bt0 = sim_pre$Bt0 * (1-b_conditioning_strength)) %>% init(BMrec0 = sim_pre$BMrec0 * (1-b_conditioning_strength)) %>% 
                        init(BMspl0 = sim_pre$BMspl0 * (1-b_conditioning_strength))
      } 

      if("N00" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(N00 = sim_pre$N00 * (1-t_conditioning_strength)) %>% init(N10 = sim_pre$N10 * (1-t_conditioning_strength)) %>% init(N20 = sim_pre$N20 * (1-t_conditioning_strength)) %>% init(N30 = sim_pre$N30 * (1-t_conditioning_strength)) %>% init(N40 = sim_pre$N40 * (1-t_conditioning_strength)) %>%
                        init(P00 = sim_pre$P00 * (1-t_conditioning_strength)) %>% init(P10 = sim_pre$P10 * (1-t_conditioning_strength)) %>% init(P20 = sim_pre$P20 * (1-t_conditioning_strength)) %>% init(P30 = sim_pre$P30 * (1-t_conditioning_strength)) %>% init(P40 = sim_pre$P40 * (1-t_conditioning_strength)) %>%
                        init(P50 = sim_pre$P50 * (1-t_conditioning_strength)) %>% init(P60 = sim_pre$P60 * (1-t_conditioning_strength)) %>% init(P70 = sim_pre$P70 * (1-t_conditioning_strength)) %>% 
                        init(S800 = sim_pre$S800 * (1-t_conditioning_strength)) %>% init(S810 = sim_pre$S810 * (1-t_conditioning_strength)) %>% init(S820 = sim_pre$S820 * (1-t_conditioning_strength)) %>% 
                        init(S400 = sim_pre$S400 * (1-t_conditioning_strength)) %>% init(S410 = sim_pre$S410 * (1-t_conditioning_strength)) %>% init(S420 = sim_pre$S420 * (1-t_conditioning_strength)) %>% 
                        init(cd4rec0 = sim_pre$cd4rec0 * (1-t_conditioning_strength)) %>% init(cd8rec0 = sim_pre$cd8rec0 * (1-t_conditioning_strength)) %>% init(cd4lym0 = sim_pre$cd4lym0 * (1-t_conditioning_strength)) %>% init(cd8lym0 = sim_pre$cd8lym0 * (1-t_conditioning_strength))
      } 

      if("GMP01" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GMP01 = sim_pre$GMP01 * ratioleft) %>% init(GMP02 = sim_pre$GMP02 * ratioleft) 
      }

      if("GM0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GM0 = sim_pre$GM0 * (1-gm_conditioning_strength)) 
      }

    sim_mppt <- mod3 %>% init(MPP11 = MPPnumber * (1-latempp), MPP12 = MPPnumber * latempp ) %>% mrgsim(delta = 0.5, end = simendtime) %>% as_tibble() %>% select(-ID)
    
    return( sim_mppt )
  }
}


##----------------- CLP transplant -----------------##

conditioning_clpt <- function(conditioning_strength = 0.1, CLPnumber = 1e4, b_conditioning_strength = 0.999, t_conditioning_strength = 0.9, gm_conditioning_strength = 0.05, mod = modo, returnSteadyState = FALSE, simendtime = 150)
{

  sim_pre <- mod %>% init(LT0 = 1)  %>%  mrgsim(end = 30 * 50, delta = 0.5) %>% 
           as_tibble() %>% tail(n=1) %>% select(-c(ID, time)) 

  if (returnSteadyState)
  {
    return(sim_pre)
  }
  else
  {
      # steady state of sickle cell patients
      mod2 <- mod %>% init(sim_pre)

      # myeloablative preconditioning
      ratioleft = 1-conditioning_strength

      mod3 <- mod2 %>% init(LT0 = sim_pre$LT0 * ratioleft) %>% 
                       init(ST01 = sim_pre$ST01 * ratioleft) %>% init(ST02 = sim_pre$ST02 * ratioleft) %>%
                       init(MPP01 = sim_pre$MPP01 * ratioleft) %>% init(MPP02 = sim_pre$MPP02 * ratioleft) %>%
                       init(CMP01 = sim_pre$CMP01 * ratioleft) %>% init(CMP02 = sim_pre$CMP02 * ratioleft) %>%
                       init(BFUE01 = sim_pre$BFUE01 * ratioleft) %>% init(BFUE02 = sim_pre$BFUE02 * ratioleft) %>%
                       init(CFUE01 = sim_pre$CFUE01 * ratioleft) %>% init(CFUE02 = sim_pre$CFUE02 * ratioleft)

      if("CLP0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(CLP0 = sim_pre$CLP0 * ratioleft) 
      } 
      
      if(("Boe0" %in% colnames(sim_pre)) & ("Boe1" %in% colnames(sim_pre)))
      {
        mod3 <-mod3 %>% init(Boe0 = sim_pre$Boe0 * ratioleft) %>% init(Bi0 = sim_pre$Bi0 * ratioleft) %>% 
                        init(Bt0 = sim_pre$Bt0 * (1-b_conditioning_strength)) %>% init(BMrec0 = sim_pre$BMrec0 * (1-b_conditioning_strength)) %>% 
                        init(BMspl0 = sim_pre$BMspl0 * (1-b_conditioning_strength))
      } 

      if("N00" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(N00 = sim_pre$N00 * (1-t_conditioning_strength)) %>% init(N10 = sim_pre$N10 * (1-t_conditioning_strength)) %>% init(N20 = sim_pre$N20 * (1-t_conditioning_strength)) %>% init(N30 = sim_pre$N30 * (1-t_conditioning_strength)) %>% init(N40 = sim_pre$N40 * (1-t_conditioning_strength)) %>%
                        init(P00 = sim_pre$P00 * (1-t_conditioning_strength)) %>% init(P10 = sim_pre$P10 * (1-t_conditioning_strength)) %>% init(P20 = sim_pre$P20 * (1-t_conditioning_strength)) %>% init(P30 = sim_pre$P30 * (1-t_conditioning_strength)) %>% init(P40 = sim_pre$P40 * (1-t_conditioning_strength)) %>%
                        init(P50 = sim_pre$P50 * (1-t_conditioning_strength)) %>% init(P60 = sim_pre$P60 * (1-t_conditioning_strength)) %>% init(P70 = sim_pre$P70 * (1-t_conditioning_strength)) %>% 
                        init(S800 = sim_pre$S800 * (1-t_conditioning_strength)) %>% init(S810 = sim_pre$S810 * (1-t_conditioning_strength)) %>% init(S820 = sim_pre$S820 * (1-t_conditioning_strength)) %>% 
                        init(S400 = sim_pre$S400 * (1-t_conditioning_strength)) %>% init(S410 = sim_pre$S410 * (1-t_conditioning_strength)) %>% init(S420 = sim_pre$S420 * (1-t_conditioning_strength)) %>% 
                        init(cd4rec0 = sim_pre$cd4rec0 * (1-t_conditioning_strength)) %>% init(cd8rec0 = sim_pre$cd8rec0 * (1-t_conditioning_strength)) %>% init(cd4lym0 = sim_pre$cd4lym0 * (1-t_conditioning_strength)) %>% init(cd8lym0 = sim_pre$cd8lym0 * (1-t_conditioning_strength))
      } 

      if("GMP01" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GMP01 = sim_pre$GMP01 * ratioleft) %>% init(GMP02 = sim_pre$GMP02 * ratioleft) 
      }

      if("GM0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GM0 = sim_pre$GM0 * (1-gm_conditioning_strength)) 
      }
      
      sim_hsct <- mod3 %>% init(CLP1 = CLPnumber) %>% mrgsim(end = 150) %>% as_tibble() %>% select(-ID)


      return( sim_hsct )
  }
}

##----------------- GMP transplant -----------------##


conditioning_gmpt <- function(conditioning_strength = 0.1, lategmp = 0, GMPnumber = 5e4, t_conditioning_strength = 0.9, b_conditioning_strength = 0.999, gm_conditioning_strength = 0.05, mod = modo, returnSteadyState = FALSE, simendtime = 150)
{

  sim_pre <- mod %>% init(LT0 = 1)  %>%  mrgsim(end = 30 * 50, delta = 0.5) %>% 
           as_tibble() %>% tail(n=1) %>% select(-c(ID, time)) 

  if (returnSteadyState)
  {
    return(sim_pre)
  }
  else
  {
      # steady state of sickle cell patients
      mod2 <- mod %>% init(sim_pre)

      # myeloablative preconditioning
      ratioleft = 1-conditioning_strength

      mod3 <- mod2 %>% init(LT0 = sim_pre$LT0 * ratioleft) %>% 
                       init(ST01 = sim_pre$ST01 * ratioleft) %>% init(ST02 = sim_pre$ST02 * ratioleft) %>%
                       init(MPP01 = sim_pre$MPP01 * ratioleft) %>% init(MPP02 = sim_pre$MPP02 * ratioleft) %>%
                       init(CMP01 = sim_pre$CMP01 * ratioleft) %>% init(CMP02 = sim_pre$CMP02 * ratioleft) %>%
                       init(BFUE01 = sim_pre$BFUE01 * ratioleft) %>% init(BFUE02 = sim_pre$BFUE02 * ratioleft) %>%
                       init(CFUE01 = sim_pre$CFUE01 * ratioleft) %>% init(CFUE02 = sim_pre$CFUE02 * ratioleft)

      if("CLP0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(CLP0 = sim_pre$CLP0 * ratioleft) 
      } 
      
      if(("Boe0" %in% colnames(sim_pre)) & ("Boe1" %in% colnames(sim_pre)))
      {
        mod3 <-mod3 %>% init(Boe0 = sim_pre$Boe0 * ratioleft) %>% init(Bi0 = sim_pre$Bi0 * ratioleft) %>% 
                        init(Bt0 = sim_pre$Bt0 * (1-b_conditioning_strength)) %>% init(BMrec0 = sim_pre$BMrec0 * (1-b_conditioning_strength)) %>% 
                        init(BMspl0 = sim_pre$BMspl0 * (1-b_conditioning_strength))
      } 

      if("N00" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(N00 = sim_pre$N00 * (1-t_conditioning_strength)) %>% init(N10 = sim_pre$N10 * (1-t_conditioning_strength)) %>% init(N20 = sim_pre$N20 * (1-t_conditioning_strength)) %>% init(N30 = sim_pre$N30 * (1-t_conditioning_strength)) %>% init(N40 = sim_pre$N40 * (1-t_conditioning_strength)) %>%
                        init(P00 = sim_pre$P00 * (1-t_conditioning_strength)) %>% init(P10 = sim_pre$P10 * (1-t_conditioning_strength)) %>% init(P20 = sim_pre$P20 * (1-t_conditioning_strength)) %>% init(P30 = sim_pre$P30 * (1-t_conditioning_strength)) %>% init(P40 = sim_pre$P40 * (1-t_conditioning_strength)) %>%
                        init(P50 = sim_pre$P50 * (1-t_conditioning_strength)) %>% init(P60 = sim_pre$P60 * (1-t_conditioning_strength)) %>% init(P70 = sim_pre$P70 * (1-t_conditioning_strength)) %>% 
                        init(S800 = sim_pre$S800 * (1-t_conditioning_strength)) %>% init(S810 = sim_pre$S810 * (1-t_conditioning_strength)) %>% init(S820 = sim_pre$S820 * (1-t_conditioning_strength)) %>% 
                        init(S400 = sim_pre$S400 * (1-t_conditioning_strength)) %>% init(S410 = sim_pre$S410 * (1-t_conditioning_strength)) %>% init(S420 = sim_pre$S420 * (1-t_conditioning_strength)) %>% 
                        init(cd4rec0 = sim_pre$cd4rec0 * (1-t_conditioning_strength)) %>% init(cd8rec0 = sim_pre$cd8rec0 * (1-t_conditioning_strength)) %>% init(cd4lym0 = sim_pre$cd4lym0 * (1-t_conditioning_strength)) %>% init(cd8lym0 = sim_pre$cd8lym0 * (1-t_conditioning_strength))
      } 

      if("GMP01" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GMP01 = sim_pre$GMP01 * ratioleft) %>% init(GMP02 = sim_pre$GMP02 * ratioleft) 
      }

      if("GM0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GM0 = sim_pre$GM0 * (1-gm_conditioning_strength)) 
      }
      
      sim_hsct <- mod3 %>% init(GMP11 = GMPnumber * (1-lategmp), GMP12 = GMPnumber * lategmp) %>% mrgsim(end = simendtime) %>% as_tibble() %>% select(-ID)


      return( sim_hsct )
  }
}






##-------------------------------------------------- TRANSPLANT FUNCTIONS FOR GLOBAL SENSITIVITY ANALYSIS --------------------------------------------------##

global_mppt <- function(conditioning_strength = 0.5, latempp = 0, MPPnumber = 1000, b_conditioning_strength = 0.999, t_conditioning_strength = 0.9, gm_conditioning_strength = 0.05, mod = modo, returnSteadyState = FALSE)
{
  
  sim_pre <- mod %>% init(LT0 = 1)  %>%  mrgsim(end = 30 * 50, delta = 0.5) %>% 
    as_tibble() %>% tail(n=1) %>% select(-c(ID, time)) 
  
  if (returnSteadyState)
  {
    return(sim_pre)
  }
  else
  {
    # steady state of sickle cell patients
    mod2 <- mod %>% init(sim_pre)
    
    # myeloablative preconditioning
    ratioleft = 1-conditioning_strength
    
    mod3 <- mod2 %>% init(LT0 = sim_pre$LT0 * ratioleft) %>% 
      init(ST01 = sim_pre$ST01 * ratioleft) %>% init(ST02 = sim_pre$ST02 * ratioleft) %>%
      init(MPP01 = sim_pre$MPP01 * ratioleft) %>% init(MPP02 = sim_pre$MPP02 * ratioleft) %>%
      init(CMP01 = sim_pre$CMP01 * ratioleft) %>% init(CMP02 = sim_pre$CMP02 * ratioleft) %>%
      init(BFUE01 = sim_pre$BFUE01 * ratioleft) %>% init(BFUE02 = sim_pre$BFUE02 * ratioleft) %>%
      init(CFUE01 = sim_pre$CFUE01 * ratioleft) %>% init(CFUE02 = sim_pre$CFUE02 * ratioleft)
    
    if("CLP0" %in% colnames(sim_pre))
    {
      mod3 <-mod3 %>% init(CLP0 = sim_pre$CLP0 * ratioleft) 
    } 
    
    if(("Boe0" %in% colnames(sim_pre)) & ("Boe1" %in% colnames(sim_pre)))
    {
      mod3 <-mod3 %>% init(Boe0 = sim_pre$Boe0 * ratioleft) %>% init(Bi0 = sim_pre$Bi0 * ratioleft) %>% 
        init(Bt0 = sim_pre$Bt0 * ratioleft) %>% init(BMrec0 = sim_pre$BMrec0 * ratioleft) %>% 
        init(BMspl0 = sim_pre$BMspl0 * (1-b_conditioning_strength))
    } 
    
    if("N00" %in% colnames(sim_pre))
    {
      mod3 <-mod3 %>% init(N00 = sim_pre$N00 * (1-t_conditioning_strength)) %>% init(N10 = sim_pre$N10 * (1-t_conditioning_strength)) %>% init(N20 = sim_pre$N20 * (1-t_conditioning_strength)) %>% init(N30 = sim_pre$N30 * (1-t_conditioning_strength)) %>% init(N40 = sim_pre$N40 * (1-t_conditioning_strength)) %>%
        init(P00 = sim_pre$P00 * (1-t_conditioning_strength)) %>% init(P10 = sim_pre$P10 * (1-t_conditioning_strength)) %>% init(P20 = sim_pre$P20 * (1-t_conditioning_strength)) %>% init(P30 = sim_pre$P30 * (1-t_conditioning_strength)) %>% init(P40 = sim_pre$P40 * (1-t_conditioning_strength)) %>%
        init(P50 = sim_pre$P50 * (1-t_conditioning_strength)) %>% init(P60 = sim_pre$P60 * (1-t_conditioning_strength)) %>% init(P70 = sim_pre$P70 * (1-t_conditioning_strength)) %>% 
        init(S800 = sim_pre$S800 * (1-t_conditioning_strength)) %>% init(S810 = sim_pre$S810 * (1-t_conditioning_strength)) %>% init(S820 = sim_pre$S820 * (1-t_conditioning_strength)) %>% 
        init(S400 = sim_pre$S400 * (1-t_conditioning_strength)) %>% init(S410 = sim_pre$S410 * (1-t_conditioning_strength)) %>% init(S420 = sim_pre$S420 * (1-t_conditioning_strength)) %>% 
        init(cd4rec0 = sim_pre$cd4rec0 * (1-t_conditioning_strength)) %>% init(cd8rec0 = sim_pre$cd8rec0 * (1-t_conditioning_strength)) %>% init(cd4lym0 = sim_pre$cd4lym0 * (1-t_conditioning_strength)) %>% init(cd8lym0 = sim_pre$cd8lym0 * (1-t_conditioning_strength))
    } 
    
    if("GMP01" %in% colnames(sim_pre))
    {
      mod3 <-mod3 %>% init(GMP01 = sim_pre$GMP01 * ratioleft) %>% init(GMP02 = sim_pre$GMP02 * ratioleft) 
    }
    
    if("GM0" %in% colnames(sim_pre))
    {
      mod3 <-mod3 %>% init(GM0 = sim_pre$GM0 * (1-gm_conditioning_strength)) 
    }
    
    mod_mppt <- mod3 %>% init(MPP11 = MPPnumber * (1-latempp), MPP12 = MPPnumber * latempp ) 
    
    return( mod_mppt )
  }
}


global_clpt <- function(conditioning_strength = 0.1, CLPnumber = 1e4, b_conditioning_strength = 0.999, t_conditioning_strength = 0.1, gm_conditioning_strength = 0.05, mod = modo, returnSteadyState = FALSE)
{
  
  sim_pre <- mod %>% init(LT0 = 1)  %>%  mrgsim(end = 30 * 50, delta = 0.5) %>% 
    as_tibble() %>% tail(n=1) %>% select(-c(ID, time)) 
  
  if (returnSteadyState)
  {
    return(sim_pre)
  }
  else
  {
    # steady state of sickle cell patients
    mod2 <- mod %>% init(sim_pre)
    
    # myeloablative preconditioning
    ratioleft = 1-conditioning_strength
    
    mod3 <- mod2 %>% init(LT0 = sim_pre$LT0 * ratioleft) %>% 
      init(ST01 = sim_pre$ST01 * ratioleft) %>% init(ST02 = sim_pre$ST02 * ratioleft) %>%
      init(MPP01 = sim_pre$MPP01 * ratioleft) %>% init(MPP02 = sim_pre$MPP02 * ratioleft) %>%
      init(CMP01 = sim_pre$CMP01 * ratioleft) %>% init(CMP02 = sim_pre$CMP02 * ratioleft) %>%
      init(BFUE01 = sim_pre$BFUE01 * ratioleft) %>% init(BFUE02 = sim_pre$BFUE02 * ratioleft) %>%
      init(CFUE01 = sim_pre$CFUE01 * ratioleft) %>% init(CFUE02 = sim_pre$CFUE02 * ratioleft)
    
    if("CLP0" %in% colnames(sim_pre))
    {
      mod3 <-mod3 %>% init(CLP0 = sim_pre$CLP0 * ratioleft) 
    } 
    
    if(("Boe0" %in% colnames(sim_pre)) & ("Boe1" %in% colnames(sim_pre)))
    {
      mod3 <-mod3 %>% init(Boe0 = sim_pre$Boe0 * ratioleft) %>% init(Bi0 = sim_pre$Bi0 * ratioleft) %>% 
        init(Bt0 = sim_pre$Bt0 * (1-b_conditioning_strength)) %>% init(BMrec0 = sim_pre$BMrec0 * (1-b_conditioning_strength)) %>% 
        init(BMspl0 = sim_pre$BMspl0 * (1-b_conditioning_strength))
    } 
    
    if("N00" %in% colnames(sim_pre))
    {
      mod3 <-mod3 %>% init(N00 = sim_pre$N00 * (1-t_conditioning_strength)) %>% init(N10 = sim_pre$N10 * (1-t_conditioning_strength)) %>% init(N20 = sim_pre$N20 * (1-t_conditioning_strength)) %>% init(N30 = sim_pre$N30 * (1-t_conditioning_strength)) %>% init(N40 = sim_pre$N40 * (1-t_conditioning_strength)) %>%
        init(P00 = sim_pre$P00 * (1-t_conditioning_strength)) %>% init(P10 = sim_pre$P10 * (1-t_conditioning_strength)) %>% init(P20 = sim_pre$P20 * (1-t_conditioning_strength)) %>% init(P30 = sim_pre$P30 * (1-t_conditioning_strength)) %>% init(P40 = sim_pre$P40 * (1-t_conditioning_strength)) %>%
        init(P50 = sim_pre$P50 * (1-t_conditioning_strength)) %>% init(P60 = sim_pre$P60 * (1-t_conditioning_strength)) %>% init(P70 = sim_pre$P70 * (1-t_conditioning_strength)) %>% 
        init(S800 = sim_pre$S800 * (1-t_conditioning_strength)) %>% init(S810 = sim_pre$S810 * (1-t_conditioning_strength)) %>% init(S820 = sim_pre$S820 * (1-t_conditioning_strength)) %>% 
        init(S400 = sim_pre$S400 * (1-t_conditioning_strength)) %>% init(S410 = sim_pre$S410 * (1-t_conditioning_strength)) %>% init(S420 = sim_pre$S420 * (1-t_conditioning_strength)) %>% 
        init(cd4rec0 = sim_pre$cd4rec0 * (1-t_conditioning_strength)) %>% init(cd8rec0 = sim_pre$cd8rec0 * (1-t_conditioning_strength)) %>% init(cd4lym0 = sim_pre$cd4lym0 * (1-t_conditioning_strength)) %>% init(cd8lym0 = sim_pre$cd8lym0 * (1-t_conditioning_strength))
    } 
    
    if("GMP01" %in% colnames(sim_pre))
    {
      mod3 <-mod3 %>% init(GMP01 = sim_pre$GMP01 * ratioleft) %>% init(GMP02 = sim_pre$GMP02 * ratioleft) 
    }
    
    if("GM0" %in% colnames(sim_pre))
    {
      mod3 <-mod3 %>% init(GM0 = sim_pre$GM0 * (1-gm_conditioning_strength)) 
    }
    
    mod_clpt <- mod3 %>% init(CLP1 = CLPnumber) 
    
    
    return( mod_clpt )
  }
}

global_hsct <- function(conditioning_strength = 0.1, hscprof = 1/(2.5*7), LTnumber = 200, b_conditioning_strength = 0.999, t_conditioning_strength = 0.9, gm_conditioning_strength = 0.05, mod = modo, returnSteadyState = FALSE, simendtime = 150)
{

  sim_pre <- mod %>% init(LT0 = 1)  %>%  mrgsim(end = 30 * 50, delta = 0.5) %>% 
           as_tibble() %>% tail(n=1) %>% select(-c(ID, time)) 

  if (returnSteadyState)
  {
    return(sim_pre)
  }
  else
  {
      # steady state of sickle cell patients
      mod2 <- mod %>% init(sim_pre)

      # myeloablative preconditioning
      ratioleft = 1-conditioning_strength
      
      mod3 <- mod2 %>% init(LT0 = sim_pre$LT0 * ratioleft) %>% 
                       init(ST01 = sim_pre$ST01 * ratioleft) %>% init(ST02 = sim_pre$ST02 * ratioleft) %>%
                       init(MPP01 = sim_pre$MPP01 * ratioleft) %>% init(MPP02 = sim_pre$MPP02 * ratioleft) %>%
                       init(CMP01 = sim_pre$CMP01 * ratioleft) %>% init(CMP02 = sim_pre$CMP02 * ratioleft) %>%
                       init(BFUE01 = sim_pre$BFUE01 * ratioleft) %>% init(BFUE02 = sim_pre$BFUE02 * ratioleft) %>%
                       init(CFUE01 = sim_pre$CFUE01 * ratioleft) %>% init(CFUE02 = sim_pre$CFUE02 * ratioleft)

      if("CLP0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(CLP0 = sim_pre$CLP0 * ratioleft) 
      } 
      
      if(("Boe0" %in% colnames(sim_pre)) & ("Boe1" %in% colnames(sim_pre)))
      {
        mod3 <-mod3 %>% init(Boe0 = sim_pre$Boe0 * ratioleft) %>% init(Bi0 = sim_pre$Bi0 * ratioleft) %>% 
                        init(Bt0 = sim_pre$Bt0 * (1-b_conditioning_strength)) %>% init(BMrec0 = sim_pre$BMrec0 * (1-b_conditioning_strength)) %>% 
                        init(BMspl0 = sim_pre$BMspl0 * (1-b_conditioning_strength))
      } 

      if("N00" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(N00 = sim_pre$N00 * (1-t_conditioning_strength)) %>% init(N10 = sim_pre$N10 * (1-t_conditioning_strength)) %>% init(N20 = sim_pre$N20 * (1-t_conditioning_strength)) %>% init(N30 = sim_pre$N30 * (1-t_conditioning_strength)) %>% init(N40 = sim_pre$N40 * (1-t_conditioning_strength)) %>%
                        init(P00 = sim_pre$P00 * (1-t_conditioning_strength)) %>% init(P10 = sim_pre$P10 * (1-t_conditioning_strength)) %>% init(P20 = sim_pre$P20 * (1-t_conditioning_strength)) %>% init(P30 = sim_pre$P30 * (1-t_conditioning_strength)) %>% init(P40 = sim_pre$P40 * (1-t_conditioning_strength)) %>%
                        init(P50 = sim_pre$P50 * (1-t_conditioning_strength)) %>% init(P60 = sim_pre$P60 * (1-t_conditioning_strength)) %>% init(P70 = sim_pre$P70 * (1-t_conditioning_strength)) %>% 
                        init(S800 = sim_pre$S800 * (1-t_conditioning_strength)) %>% init(S810 = sim_pre$S810 * (1-t_conditioning_strength)) %>% init(S820 = sim_pre$S820 * (1-t_conditioning_strength)) %>% 
                        init(S400 = sim_pre$S400 * (1-t_conditioning_strength)) %>% init(S410 = sim_pre$S410 * (1-t_conditioning_strength)) %>% init(S420 = sim_pre$S420 * (1-t_conditioning_strength)) %>% 
                        init(cd4rec0 = sim_pre$cd4rec0 * (1-t_conditioning_strength)) %>% init(cd8rec0 = sim_pre$cd8rec0 * (1-t_conditioning_strength)) %>% init(cd4lym0 = sim_pre$cd4lym0 * (1-t_conditioning_strength)) %>% init(cd8lym0 = sim_pre$cd8lym0 * (1-t_conditioning_strength))
      } 

      if("GMP01" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GMP01 = sim_pre$GMP01 * ratioleft) %>% init(GMP02 = sim_pre$GMP02 * ratioleft) 
      }

      if("GM0" %in% colnames(sim_pre))
      {
        mod3 <-mod3 %>% init(GM0 = sim_pre$GM0 * (1-gm_conditioning_strength)) 
      }

      sim_hsct <- mod3 %>% init(LT1 = LTnumber) %>% param(rLT = hscprof)


      return( sim_hsct )
  }
}
