library(ggplot)


# fit for Fig 2a. 

fig2a <- read.csv('Fig2a.csv', header = TRUE)

spleenNK = 2e6

fig2a['NK'] = spleenNK * fig2a$'BrdU'

# add_nls <- nls(NK ~ c + s * time + exp(r*time), data = fig2a, start = list(s = 1, r = 1), algorithm = "plinear")
# nonlinear regression reports error, indicating the exponential part does not really exist

lm(NK ~ time , data = fig2a, method = "qr")
# (Intercept)        time 
# 208012.80    43655.74 
# estimated differentiation: rate * NKP = 43k

# fit for Fig 3a, b; fit for NK death rate in bone marrow and spleen

fig3a <- read.csv('Fig3a.csv', header = TRUE)

lm(log(BrdU) ~ time , data = fig3a, method = "qr")
# Coefficients:
#   (Intercept)         time  
#     -0.1962      -0.0251  

# NK death rate, in bone marrow
fig3b <- read.csv('Fig3b.csv', header = TRUE)

lm(log(BrdU) ~ time , data = fig3b, method = "qr")
# Coefficients:
#  (Intercept)         time  
#   -0.440259    -0.007204  

