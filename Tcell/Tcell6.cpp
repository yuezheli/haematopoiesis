[PROB]

T cell dynamics in thymus, adapted of Thomas-Vaslin et al., 2008
https://pubmed.ncbi.nlm.nih.gov/18250431/

All cells are in real cell count

Additional compartment: 
CD4rec -- CD4+ cells in the blood stream
CD8rec -- CD8+ cells in the blood stream



[SET]
delta = 1

[GLOBAL]  
#define min(a,b)  ( a<b ? a : b)


[CMT] 

// double negative thymocyte
N00
N10
N20
N30
N40
N01
N11
N21
N31
N41

// double positive thymocytes
P00
P10
P20
P30
P40
P50
P60
P70
P01
P11
P21
P31
P41
P51
P61
P71

// single positive CD4 thymocyte
S400
S410
S420
S401
S411
S421

// single positive CD8 thymocyte
S800
S810
S820
S801
S811
S821

// T cells in blood
cd4rec0 
cd8rec0
cd4rec1 
cd8rec1

// T cells in the secondary lymphoid organs
cd4lym0
cd8lym0
cd4lym1
cd8lym1

[PARAM]

// input rate for thymic precursor
sigmaN = 2e4

// proliferation 
pN = 0.23   // DN cells
pP = 4.5    // DP cells
pS = 0.23   // SP cells

// selection DP into SP4/ SP8 (fraction)
alpha4 = 0.06    // DP -> SP4
alpha8 = 0.01    // DP -> SP8


// differentiation parameters
alpha_muN = 0.29 // DN 
alpha_muP = 0.2029  // DP
alpha_e = 0.994   // SP4
alpha_r = 0.48   // SP8

// exponent of differentiation function
n = 127 

// death rate for DN, DP, SP
delta = 0 

// removal rate last stage DP (LP)
muLP = 0.37

// T cell entering and leaving secondary lymphoid organs
enter_lymph_cd4 = 10
enter_lymph_cd8 = 5
exit_lymph = 1.2

// death rate of T cells in blood
death_cd4b = 1/31
death_cd8b = 1/72

// death rate of T cells in secondary lymphoid organs
death_cd4l = 0.002
death_cd8l = 0.001


[ODE]

// dynamics of DN cells
dxdt_N00 = sigmaN - (pN + delta) * N00; 
dxdt_N10 =  2 * pN * N00 - ( pN + delta +  min( pow(alpha_muN * 1, n), 100) ) * N10; 
dxdt_N20 =  2 * pN * N10 - ( pN + delta +  min( pow(alpha_muN * 2, n), 100) ) * N20; 
dxdt_N30 =  2 * pN * N20 - ( pN + delta +  min( pow(alpha_muN * 3, n), 100) ) * N30; 
dxdt_N40 =  2 * pN * N30 - ( pN + delta +  min( pow(alpha_muN * 4, n), 100) ) * N40; 

dxdt_N01 = - (pN + delta) * N01; // transduced branch, cell number can only shrink
dxdt_N11 =  2 * pN * N01 - ( pN + delta +  min( pow(alpha_muN * 1, n), 100) ) * N11; 
dxdt_N21 =  2 * pN * N11 - ( pN + delta +  min( pow(alpha_muN * 2, n), 100) ) * N21; 
dxdt_N31 =  2 * pN * N21 - ( pN + delta +  min( pow(alpha_muN * 3, n), 100) ) * N31; 
dxdt_N41 =  2 * pN * N31 - ( pN + delta +  min( pow(alpha_muN * 4, n), 100) ) * N41; 


// dynamics of DP cells
double sum_mu_N0 = min( pow(alpha_muN * 1, n), 100) * N10 + min( pow(alpha_muN * 2, n), 100) * N20 + min( pow(alpha_muN * 3, n), 100) * N30 + min( pow(alpha_muN * 4, n), 100) * N40; 
double sum_mu_N1 = min( pow(alpha_muN * 1, n), 100) * N11 + min( pow(alpha_muN * 2, n), 100) * N21 + min( pow(alpha_muN * 3, n), 100) * N31 + min( pow(alpha_muN * 4, n), 100) * N41; 

dxdt_P00 = sum_mu_N0 + 2 * pN * N40 - (pP + delta) * P00; 
dxdt_P10 = 2 * pP * P00 - ( pP + delta + min( pow(alpha_muP * 1, n), 100) ) * P10; 
dxdt_P20 = 2 * pP * P10 - ( pP + delta + min( pow(alpha_muP * 2, n), 100) ) * P20; 
dxdt_P30 = 2 * pP * P20 - ( pP + delta + min( pow(alpha_muP * 3, n), 100) ) * P30; 
dxdt_P40 = 2 * pP * P30 - ( pP + delta + min( pow(alpha_muP * 4, n), 100) ) * P40; 
dxdt_P50 = 2 * pP * P40 - ( pP + delta + min( pow(alpha_muP * 5, n), 100) ) * P50; 
dxdt_P60 = 2 * pP * P50 - ( pP + delta + min( pow(alpha_muP * 6, n), 100) ) * P60; 

dxdt_P01 = sum_mu_N1 + 2 * pN * N41 - (pP + delta) * P01; 
dxdt_P11 = 2 * pP * P01 - ( pP + delta + min( pow(alpha_muP * 1, n), 100) ) * P11; 
dxdt_P21 = 2 * pP * P11 - ( pP + delta + min( pow(alpha_muP * 2, n), 100) ) * P21; 
dxdt_P31 = 2 * pP * P21 - ( pP + delta + min( pow(alpha_muP * 3, n), 100) ) * P31; 
dxdt_P41 = 2 * pP * P31 - ( pP + delta + min( pow(alpha_muP * 4, n), 100) ) * P41; 
dxdt_P51 = 2 * pP * P41 - ( pP + delta + min( pow(alpha_muP * 5, n), 100) ) * P51; 
dxdt_P61 = 2 * pP * P51 - ( pP + delta + min( pow(alpha_muP * 6, n), 100) ) * P61; 

double sum_mu_P0 = min( pow(alpha_muP * 1, n), 100) * P10 +  min( pow(alpha_muP * 2, n), 100) * P20 + min( pow(alpha_muP * 3, n), 100) * P30 +
  min( pow(alpha_muP * 4, n), 100) * P40 +  min( pow(alpha_muP * 5, n), 100) * P50 + min( pow(alpha_muP * 6, n), 100) * P60; 

double sum_mu_P1 = min( pow(alpha_muP * 1, n), 100) * P11 +  min( pow(alpha_muP * 2, n), 100) * P21 + min( pow(alpha_muP * 3, n), 100) * P31 +
  min( pow(alpha_muP * 4, n), 100) * P41 +  min( pow(alpha_muP * 5, n), 100) * P51 + min( pow(alpha_muP * 6, n), 100) * P61;

dxdt_P70 = sum_mu_P0 + 2 * pP * P60 - muLP * P70;
dxdt_P71 = sum_mu_P1 + 2 * pP * P61 - muLP * P71;  

// dynamics of CD4+ cells
dxdt_S400 = alpha4 * muLP * P70 - (pS + delta) * S400; 
dxdt_S410 = 2 * pS * S400 - (pS + delta + min( pow(alpha_e * 1, n), 100) ) * S410; 
dxdt_S420 = 2 * pS * S410 - (delta + min( pow(alpha_e * 2, n), 100) ) * S420; 

dxdt_S401 = alpha4 * muLP * P71 - (pS + delta) * S401; 
dxdt_S411 = 2 * pS * S401 - (pS + delta + min( pow(alpha_e * 1, n), 100) ) * S411; 
dxdt_S421 = 2 * pS * S411 - (delta + min( pow(alpha_e * 2, n), 100) ) * S421; 

double sum_SP40 =  min( pow(alpha_e * 1, n), 100) * S410 +  min( pow(alpha_e * 2, n), 100) * S420; 
double sum_SP41 =  min( pow(alpha_e * 1, n), 100) * S411 +  min( pow(alpha_e * 2, n), 100) * S421; 


// dynamics of CD8+ cells

dxdt_S800 = alpha8 * muLP * P70 - (pS + delta) * S800; 
dxdt_S810 = 2 * pS * S800 - (pS + delta +  min( pow(alpha_e * 1, n), 100) ) * S810; 
dxdt_S820 = 2 * pS * S810 - (delta + min( pow(alpha_e * 2, n), 100) ) * S820; 

dxdt_S801 = alpha8 * muLP * P71 - (pS + delta) * S801; 
dxdt_S811 = 2 * pS * S801 - (pS + delta +  min( pow(alpha_e * 1, n), 100) ) * S811; 
dxdt_S821 = 2 * pS * S811 - (delta + min( pow(alpha_e * 2, n), 100) ) * S821; 

double sum_SP80 =  min( pow(alpha_e * 1, n), 100) * S810 +  min( pow(alpha_e * 2, n), 100) * S820; 
double sum_SP81 =  min( pow(alpha_e * 1, n), 100) * S811 +  min( pow(alpha_e * 2, n), 100) * S821; 


// dynamics of T cells in blood
dxdt_cd4rec0 = sum_SP40 - death_cd4b * cd4rec0 + exit_lymph * cd4lym0 - enter_lymph_cd4 * cd4rec0;
dxdt_cd8rec0 = sum_SP80 - death_cd8b * cd8rec0 + exit_lymph * cd8lym0 - enter_lymph_cd8 * cd8rec0;

dxdt_cd4rec1 = sum_SP41 - death_cd4b * cd4rec1 + exit_lymph * cd4lym1 - enter_lymph_cd4 * cd4rec1;
dxdt_cd8rec1 = sum_SP81 - death_cd8b * cd8rec1 + exit_lymph * cd8lym1 - enter_lymph_cd8 * cd8rec1;

// dynamics of T cells in secondary lymphoid system
dxdt_cd4lym0 = enter_lymph_cd4 * cd4rec0 - (exit_lymph + death_cd4l) * cd4lym0;
dxdt_cd8lym0 = enter_lymph_cd8 * cd8rec0 - (exit_lymph + death_cd8l) * cd8lym0; 

dxdt_cd4lym1 = enter_lymph_cd4 * cd4rec1 - (exit_lymph + death_cd4l) * cd4lym1;
dxdt_cd8lym1 = enter_lymph_cd8 * cd8rec1 - (exit_lymph + death_cd8l) * cd8lym1; 

[TABLE]

capture tcellbloodconc = (cd4rec0 + cd8rec0 + cd4rec1 + cd8rec1)/ 2e3;
capture tcell_ratio_blood = (cd4rec0 + cd4rec1)/ (cd8rec0 + cd8rec1);
capture tcell_lymph = cd4lym0 + cd8lym0 + cd4lym1 + cd8lym1;
