rm(list = ls())
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
library(plotly)


mod <- mread("thalassemia") 
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


synbetascan <- function(synbeta = 0.8, tauRBCsickle_default = 36, totalCD34infusion = 1e8)
{
  # pre-equilibrium
  modB <- mod %>% param(totalCD34infused = 0) %>% param(ksynalpha = 6e-7) %>% # these are for SCD parameter change
    param(HbS_saturation = 0.74)  %>% # change the blood carrying capacity of mutated HbA (in SCD, it is HbS)
    init(LT0 = 1) %>% # set the system with only 1 stem cell
    param(tauRBCsickle  = tauRBCsickle_default) %>% # change the avg RBC life span 
    param(E_synbeta_effect = synbeta) 
  
  sim0 <- modB %>% mrgsim(end = 30 * 12, delta = 0.5) %>% as_tibble() 
  data = tail(sim0, n = 1)
  Data <- data[,-(1:2),drop=FALSE]  
  # steady state of sickle cell patients
  mod2 <- modB %>% init(Data)
  # myeloablative preconditioning
  ratioleft = 1-0.9
  mod3 <- mod2 %>% init(LT0 = Data$LT0 * ratioleft) %>% 
    init(ST01 = Data$ST01 * ratioleft) %>% init(ST02 = Data$ST02 * ratioleft) %>%
    init(MPP01 = Data$MPP01 * ratioleft) %>% init(MPP02 = Data$MPP02 * ratioleft) %>%
    init(CMP01 = Data$CMP01 * ratioleft) %>% init(CMP02 = Data$CMP02 * ratioleft)
  
  # post-infusion simulation
  sim_infusion <- mod3 %>% param(totalCD34infused = totalCD34infusion) %>% 
    mrgsim(delta = 0.5, end = 30 * 18) %>% as_tibble() 
  

  return( data.frame(  list(time = sim_infusion['time'], HbA_total = sim_infusion['totalHbA']) )  )
}

sim_point8_36 = synbetascan(synbeta = 0.8, tauRBCsickle_default = 36)
sim_point4_36 = synbetascan(synbeta = 0.4, tauRBCsickle_default = 36)
sim_0_36 = synbetascan(synbeta = 0.01, tauRBCsickle_default = 36)


sim_point8_12 = synbetascan(synbeta = 0.8, tauRBCsickle_default = 12)
sim_point4_12 = synbetascan(synbeta = 0.4, tauRBCsickle_default = 12)
sim_0_12 = synbetascan(synbeta = 0.01, tauRBCsickle_default = 12)


HbA_t87q_plot = ggplot() + 
  geom_line(data = sim_point8_36, aes(x = time/30, y = totalHbA, color = '\u03B2-globin syn 80%; RBC lifespan = 36 days'), alpha = 0.5) + 
  geom_line(data = sim_point4_36, aes(x = time/30, y = totalHbA, color = '\u03B2-globin syn 40%; RBC lifespan = 36 days'), alpha = 0.5) + 
  geom_line(data = sim_0_36, aes(x = time/30, y = totalHbA, color = '\u03B2-globin syn 1%; RBC lifespan = 36 days'), alpha = 0.5) + 
  geom_line(data = sim_point8_12, aes(x = time/30, y = totalHbA, color = '\u03B2-globin syn 80%; RBC lifespan = 12 days'), alpha = 0.5) + 
  geom_line(data = sim_point4_12, aes(x = time/30, y = totalHbA, color = '\u03B2-globin syn 40%; RBC lifespan = 12 days'), alpha = 0.5) + 
  geom_line(data = sim_0_12, aes(x = time/30, y = totalHbA, color = '\u03B2-globin syn 1%; RBC lifespan = 12 days'), alpha = 0.5) + 
  labs(y = expression(HbA^T87Q~conc~(g/dL)), x = 'time (months)', color = ' ') + theme_bw() 


## ------------------ steady state before transfusion ------------------ ##

preequilibrium <- function(synbeta = 0.8, tauRBCsickle_default = 36)
{
  # pre-equilibrium
  modB <- mod %>% param(totalCD34infused = 0) %>% param(ksynalpha = 6e-7) %>% # these are for SCD parameter change
    param(HbS_saturation = 0.74)  %>% # change the blood carrying capacity of mutated HbA (in SCD, it is HbS)
    init(LT0 = 1) %>% # set the system with only 1 stem cell
    param(tauRBCsickle  = tauRBCsickle_default) %>% # change the avg RBC life span 
    param(E_synbeta_effect = synbeta) 
  
  sim0 <- modB %>% mrgsim(end = 30 * 12, delta = 0.5) %>% as_tibble() %>% select('RBCconc', 'totalHb')
  data = tail(sim0, n = 1)
  
  return(data)
}

synbeta_range = c(0.01, 0.1, 0.3, 0.4, 0.5, 0.8)
tauRBC_range = c(12, 20, 36, 43)

RBCconc = matrix(rep(NA, length(synbeta_range) * length(tauRBC_range), ncol = length(tauRBC_range)), ncol = 4)
totalHb = matrix(rep(NA, length(synbeta_range) * length(tauRBC_range), ncol = length(tauRBC_range)), ncol = 4)


for(i in 1:length(synbeta_range))
{
  for(j in 1:length(tauRBC_range))
  {
    tmp = preequilibrium(synbeta = synbeta_range[i], tauRBCsickle_default = tauRBC_range[j])
    
    RBCconc[i,j] = tmp$RBCconc
    totalHb[i,j] = tmp$totalHb
  }
}



fig_rbc <- plot_ly(x = tauRBC_range, y = synbeta_range * 100, z = RBCconc, type = 'contour', 
                   autocontour = F,
                   contours = list(
                     start = 1.5e6,
                     end = 5e6,
                     size = 5e5
                   ))%>%
  layout(title = "RBC conc(#/uL)", 
         yaxis = list(title = "\u03B2-globin syn (%)"), 
         xaxis = list(title = "RBC lifespan(days)"))



fig_hb <- plot_ly(x = tauRBC_range, y = synbeta_range * 100 , , z = totalHb, type = 'contour',
                  colorscale = 'Jet',
                  autocontour = F,
                  contours = list(
                    start = 3,
                    end = 10,
                    size = 1)    
              ) %>%
  layout(yaxis = list(title = '\u03B2-globin syn (%)'), 
         xaxis = list(title = 'RBC lifespan (days)'),  
         title = 'Hb conc (g/dL)' ) 


fig <- subplot(fig_rbc, fig_hb, nrows = 1, titleY = TRUE, titleX = TRUE, margin = 0.1 )

fig <- fig %>%layout(title = ' ',
                     plot_bgcolor='#e5ecf6', 
                     xaxis = list( 
                       zerolinecolor = '#ffff', 
                       zerolinewidth = 2, 
                       gridcolor = 'ffff'), 
                     yaxis = list( 
                       zerolinecolor = '#ffff', 
                       zerolinewidth = 2, 
                       gridcolor = 'ffff'))

# Update title
annotations = list( 
  list( 
    x = 0.2,  
    y = 1,  
    text = "RBC count (#/uL)",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),
  list( 
    x = 0.8,  
    y = 1,  
    text = "Hb concentration (g/dL)",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ))

fig <- fig %>%layout(annotations = annotations) 
fig
