[PROB]

T cell dynamics, adapted of Thomas-Vaslin et al., 2008
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
N0
N1
N2
N3
N4

// double positive thymocytes
P0
P1
P2
P3
P4
P5
P6
P7

// single positive CD4 thymocyte
S40
S41
S42

// single positive CD8 thymocyte
S80
S81
S82

// T cells in blood
cd4rec 
cd8rec

// T cells in the secondary lymphoid organs
cd4lym
cd8lym

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
enter_lymph = 0.1
exit_lymph = 1.2

// death rate of T cells in blood
death_cd4b = 1/31
death_cd8b = 1/72


[ODE]

// dynamics of DN cells
dxdt_N0 = sigmaN - (pN + delta) * N0; 
dxdt_N1 =  2 * pN * N0 - ( pN + delta +  min( pow(alpha_muN * 1, n), 100) ) * N1; 
dxdt_N2 =  2 * pN * N1 - ( pN + delta +  min( pow(alpha_muN * 2, n), 100) ) * N2; 
dxdt_N3 =  2 * pN * N2 - ( pN + delta +  min( pow(alpha_muN * 3, n), 100) ) * N3; 
dxdt_N4 =  2 * pN * N3 - ( pN + delta +  min( pow(alpha_muN * 4, n), 100) ) * N4; 


// dynamics of DP cells
double sum_mu_N = min( pow(alpha_muN * 1, n), 100) * N1 + min( pow(alpha_muN * 2, n), 100) * N2 + min( pow(alpha_muN * 3, n), 100) * N3 + min( pow(alpha_muN * 4, n), 100) * N4; 

dxdt_P0 = sum_mu_N + 2 * pN * N4 - (pP + delta) * P0; 

dxdt_P1 = 2 * pP * P0 - ( pP + delta + min( pow(alpha_muP * 1, n), 100) ) * P1; 
dxdt_P2 = 2 * pP * P1 - ( pP + delta + min( pow(alpha_muP * 2, n), 100) ) * P2; 
dxdt_P3 = 2 * pP * P2 - ( pP + delta + min( pow(alpha_muP * 3, n), 100) ) * P3; 
dxdt_P4 = 2 * pP * P3 - ( pP + delta + min( pow(alpha_muP * 4, n), 100) ) * P4; 
dxdt_P5 = 2 * pP * P4 - ( pP + delta + min( pow(alpha_muP * 5, n), 100) ) * P5; 
dxdt_P6 = 2 * pP * P5 - ( pP + delta + min( pow(alpha_muP * 6, n), 100) ) * P6; 

double sum_mu_P = min( pow(alpha_muP * 1, n), 100) * P1 +  min( pow(alpha_muP * 2, n), 100) * P2 + min( pow(alpha_muP * 3, n), 100) * P3 +
                  min( pow(alpha_muP * 4, n), 100) * P4 +  min( pow(alpha_muP * 5, n), 100) * P5 + min( pow(alpha_muP * 6, n), 100) * P6; 

dxdt_P7 = sum_mu_P + 2 * pP * P6 - muLP * P7; 

// dynamics of CD4+ cells
dxdt_S40 = alpha4 * muLP * P7 - (pS + delta) * S40; 

dxdt_S41 = 2 * pS * S40 - (pS + delta + min( pow(alpha_e * 1, n), 100) ) * S41; 

dxdt_S42 = 2 * pS * S41 - (delta + min( pow(alpha_e * 2, n), 100) ) * S42; 

double sum_SP4 =  min( pow(alpha_e * 1, n), 100) * S41 +  min( pow(alpha_e * 2, n), 100) * S42; 

// dynamics of CD8+ cells

dxdt_S80 = alpha8 * muLP * P7 - (pS + delta) * S80; 

dxdt_S81 = 2 * pS * S80 - (pS + delta +  min( pow(alpha_e * 1, n), 100) ) * S81; 

dxdt_S82 = 2 * pS * S81 - (delta + min( pow(alpha_e * 2, n), 100) ) * S82; 

double sum_SP8 =  min( pow(alpha_e * 1, n), 100) * S81 +  min( pow(alpha_e * 2, n), 100) * S82; 

// dynamics of T cells in blood

dxdt_cd4rec = sum_SP4 - death_cd4b * cd4rec + exit_lymph * cd4lym - enter_lymph * cd4rec;
dxdt_cd8rec = sum_SP8 - death_cd8b * cd8rec + exit_lymph * cd8lym - enter_lymph * cd8rec;

// dynamics of T cells in secondary lymphoid system

dxdt_cd4lym = enter_lymph * cd4rec - exit_lymph * cd4lym;
dxdt_cd8lym = enter_lymph * cd8rec - exit_lymph * cd8lym; 

[TABLE]

capture DPdeath = muLP * P7 * (1-alpha4 - alpha8);

capture tcellbloodconc = (cd4rec + cd8rec)/ 2e3;
capture tcell_ratio_blood = cd4rec/ cd8rec;
capture tcell_lymph_ratio = cd4lym/ cd8lym; 
capture tcell_lymph = cd4lym + cd8lym;

[CAPTURE]
sum_SP4, sum_SP8