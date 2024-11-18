[PROB]

T cell dynamics, implementation of Thomas-Vaslin et al., 2008
https://pubmed.ncbi.nlm.nih.gov/18250431/

Note all the cells are modeled in Millions of Cells

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

[PARAM]

// input rate for thymic precursor
sigmaN = 0.02

// proliferation 
pN = 0.23   // DN cells
pP = 4.5    // DP cells
pS = 0.23   // SP cells

// selection DP into SP4/ SP8 (fraction)
alpha4 = 0.06    // DP -> SP4
alpha8 = 0.01    // DP -> SP8

// fraction of exported SP cells going to spleen
f4 = 0.36  // SP4 -> spleen
f8 = 0.74  // SP8 -> spleen

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

[TABLE]

capture DN = N0 + N1 + N2 + N3 + N4; 
capture earlyDP = P0 + P1 + P2 + P3 + P4 + P5 + P6;
capture lateDP = P7; 
capture DP = P0 + P1 + P2 + P3 + P4 + P5 + P6 + P7; 
capture SP4 = S40 + S41 + S42; 
capture SP8 = S80 + S81 + S82; 
capture thymocyte = DN + DP + SP4 + SP8; 

capture DP2SP4_flux = alpha4 * muLP * P7;
capture DP2SP8_flux = alpha8 * muLP * P7;
capture DN2DP_flux = sum_mu_N + pN * N4;

capture P5toP7_flux =  min( pow(alpha_muP * 5, n), 100) * P5 ; 
capture P6toP7_flux =  min( pow(alpha_muP * 6, n), 100) * P6 ; 

capture DPdeath = muLP * P7 * (1-alpha4 - alpha8);

capture SP4_blood_flux = sum_SP4; 
capture SP8_blood_flux = sum_SP8; 
capture SP4_spleen_flux = sum_SP4 * f4; 
capture SP8_spleen_flux = sum_SP8 * f8; 