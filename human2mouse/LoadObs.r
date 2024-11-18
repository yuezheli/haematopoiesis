# read in observed data
## RBC 
hsct_rbc <- read.csv('MouseData/Boyer2019_1M_HSC.csv', header = TRUE) %>% filter(cell_type == "RBC")
hsct_rbc[,4] <- read.csv('MouseData/Boyer2019_1M_HSC.csv', header = TRUE) %>% filter(cell_type == "RBClower") %>% select(cell_count) 
colnames(hsct_rbc)[4] <- "RBC_lower"
hsct_rbc$RBC_upper <- hsct_rbc$cell_count + (hsct_rbc$cell_count - hsct_rbc$RBC_lower)

mppt_rbc <- read.csv('MouseData/Boyer2019_1N_MPP.csv', header = TRUE) %>% filter(cell_type == "RBC")
mppt_rbc[,4] <- read.csv('MouseData/Boyer2019_1N_MPP.csv', header = TRUE) %>% filter(cell_type == "RBClower") %>% select(cell_count) 
colnames(mppt_rbc)[4] <- "RBC_lower"
mppt_rbc$RBC_upper <- mppt_rbc$cell_count + (mppt_rbc$cell_count - mppt_rbc$RBC_lower)

## B cell 
hsct_B <- read.csv('MouseData/Boyer2019_1A_HSC.csv', header = TRUE) %>% filter(cell_type == "B")
hsct_B[,4] <- read.csv('MouseData/Boyer2019_1A_HSC.csv', header = TRUE) %>% filter(cell_type == "B_upper") %>% select(cell_fraction) 
colnames(hsct_B)[4] <- "B_upper"
hsct_B$B_lower <- hsct_B$cell_fraction - (hsct_B$B_upper - hsct_B$cell_fraction)
mppt_B <- read.csv('MouseData/Boyer2019_1B_MPP.csv', header = TRUE) %>% filter(cell_type == "B")
mppt_B[,4] <- read.csv('MouseData/Boyer2019_1B_MPP.csv', header = TRUE) %>% filter(cell_type == "B_upper") %>% select(cell_fraction) 
colnames(mppt_B)[4] <- "B_upper"
mppt_B$B_lower <- mppt_B$cell_fraction - (mppt_B$B_upper - mppt_B$cell_fraction)
clpt_B <- read.csv('MouseData/Boyer2019_1G_CLP.csv', header = TRUE) %>% filter(cell_type == "B")
clpt_B[,4] <- read.csv('MouseData/Boyer2019_1G_CLP.csv', header = TRUE) %>% filter(cell_type == "B_upper") %>% select(cell_fraction) 
colnames(clpt_B)[4] <- "B_upper"
clpt_B$B_lower <- clpt_B$cell_fraction - (clpt_B$B_upper - clpt_B$cell_fraction)

## T cell
hsct_T <- read.csv('MouseData/Boyer2019_1A_HSC.csv', header = TRUE) %>% filter(cell_type == "T")
hsct_T[,4] <- read.csv('MouseData/Boyer2019_1A_HSC.csv', header = TRUE) %>% filter(cell_type == "T_upper") %>% select(cell_fraction) 
colnames(hsct_T)[4] <- "T_upper"
hsct_T$T_lower <- hsct_T$cell_fraction - (hsct_T$T_upper - hsct_T$cell_fraction)
mppt_T <- read.csv('MouseData/Boyer2019_1B_MPP.csv', header = TRUE) %>% filter(cell_type == "T")
mppt_T[,4] <- read.csv('MouseData/Boyer2019_1B_MPP.csv', header = TRUE) %>% filter(cell_type == "T_upper") %>% select(cell_fraction) 
colnames(mppt_T)[4] <- "T_upper"
mppt_T$T_lower <- mppt_T$cell_fraction - (mppt_T$T_upper - mppt_T$cell_fraction)
clpt_T <- read.csv('MouseData/Boyer2019_Supp1D_CLP.csv', header = TRUE)
