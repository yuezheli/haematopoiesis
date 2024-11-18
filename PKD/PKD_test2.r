rm(list = ls())

# load required packages
library(tidyverse)
library(mrgsolve)
library(gridExtra)
library(grid)
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # set the working directory at the folder that contains the script


pkdmod <- mread("PKD") %>%
          init(LT0 = 1) # %>% init(LT1 = 1)

# parameter tuning for kdeathRETendo

RETdeathscan <- function(RETdeath = 0.5, tauRBC_endo = 27)
{
  sim0 <- pkdmod  %>% param(kdeathRETendo = RETdeath, tauRBCendo = tauRBC_endo) %>%
    mrgsim(end = 30 * 18, delta = 1) %>% as_tibble() %>%
    select(time, totalHbA, totalHbF, totalHbA2, totalHb, LTHSC, RETconc, RBCconc, reticulocyte) %>%
    tail(n = 1) %>% as.data.frame()
  
  sim0['RETdeath'] = RETdeath
  
  sim0['tauRBCendo'] = tauRBC_endo
  
  return(sim0)
}


##---------- firstly look for RET death rate----------##

RETdeathrate <- (1:30)/30

tmpsim <- RETdeathscan(0)

for(i in 1:length(RETdeathrate))
{
  tmp <- RETdeathscan(RETdeathrate[i])
  
  tmpsim <- rbind(tmpsim, tmp)
  
  rm(tmp)
}

bloodcellcount <- ggplot(data = tmpsim) +
  geom_line(aes(x = RETdeath, y = RETconc, color = 'RET count')) + 
  geom_line(aes(x = RETdeath, y = RBCconc, color = 'RBC count')) + 
  labs(y = 'cell count (#/uL)', x = 'RET death rate', color = '') + theme_bw() + theme(legend.position = "bottom")

HBconc <- ggplot(data = tmpsim) +
  geom_line(aes(x = RETdeath, y = totalHbA, color = 'HbA')) + 
  geom_line(aes(x = RETdeath, y = totalHbF, color = 'HbF')) + 
  geom_line(aes(x = RETdeath, y = totalHbA2, color = 'HbA2')) + 
  geom_line(aes(x = RETdeath, y = totalHb, color = 'Hb')) + 
  labs(y = 'Hb conc (g/dL)', x = 'RET death rate', color = '') + theme_bw() + theme(legend.position = "bottom")

grid.arrange(bloodcellcount, HBconc, ncol = 2)


##---------- double scan on RBC lifespan and RET death rate -----------##

RBClifespan <- c(20, 23, 26, 28, 30, 32, 40, 50, 60)

RETdeathrate <- (1:10)/10

RBCcountholder = matrix(rep(NA, length(RBClifespan) * length(RETdeathrate)), ncol = length(RETdeathrate))
Hbconcholder = matrix(rep(NA, length(RBClifespan) * length(RETdeathrate)), ncol = length(RETdeathrate))
reticulocyteholder = matrix(rep(NA, length(RBClifespan) * length(RETdeathrate)), ncol = length(RETdeathrate))


for(i in 1:length(RBClifespan))
{
  
  for(j in 1:length(RETdeathrate))
  {
    tmp <- RETdeathscan(RETdeathrate[j], RBClifespan[i])
    
    RBCcountholder[i,j] = tmp$RBCconc
    Hbconcholder[i,j] = tmp$totalHb
    reticulocyteholder[i,j] = tmp$reticulocyte
    
    rm(tmp)
  }

}

library(plotly)

fig_rbc <- plot_ly(x = RETdeathrate, y = RBClifespan, z = RBCcountholder, type = 'contour', 
                   colorscale = 'Viridis',
                   autocontour = F,
                   contours = list(
                     start = 2e6,
                     end = 5e6,
                     size = 4e5
                   ))%>%
  layout(title = "RBC conc(#/uL)", 
         xaxis = list(title = "RET death rate (day-1)"), 
         yaxis = list(title = "RBC lifespan(days)"))

fig_hb <- plot_ly(x = RETdeathrate, y = RBClifespan, z = Hbconcholder, type = 'contour', 
                  colorscale = 'Picnic',
                   autocontour = F,
                   contours = list(
                     start = 2,
                     end = 12,
                     size = 1
                   ))%>%
  layout(title = "Hb conc(g/dL)", 
         xaxis = list(title = "RET death rate (day-1)"), 
         yaxis = list(title = "RBC lifespan(days)"))

fig_reticulocyte <- plot_ly(x = RETdeathrate, y = RBClifespan, z = reticulocyteholder, type = 'contour', 
                  colorscale = 'Greys',
                  autocontour = T,
                  contours = list(
                    start = 2,
                    end = 12,
                    size = 1
                  ))%>%
  layout(title = "reticulocyte (%)", 
         xaxis = list(title = "RET death rate (day-1)"), 
         yaxis = list(title = "RBC lifespan(days)"))


fig <- subplot(fig_rbc, fig_hb, fig_reticulocyte, nrows = 1, titleY = TRUE, titleX = TRUE, margin = 0.05)

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
    x = 0.05,  
    y = 1,  
    text = "RBC count (#/uL)",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ),
  list( 
    x = 0.45,  
    y = 1,  
    text = "Hb concentration (g/dL)",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  ), 
  list( 
    x = 0.85,  
    y = 1,  
    text = "reticulocyte (%)",  
    xref = "paper",  
    yref = "paper",  
    xanchor = "center",  
    yanchor = "bottom",  
    showarrow = FALSE 
  )
)

fig <- fig %>%layout(annotations = annotations) 
fig
