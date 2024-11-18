[PROB]

This file is copied from prior work that was presented at ACoP14 (National Habor, MD; Nov 2023)
https://www.metrumrg.com/wp-content/uploads/2023/11/poster-YL.pdf

This file removed transplant arm associated with HSCT

[SET]
delta = 0.5, end = 70

[GLOBAL]  
#define min(a,b)  ( a<b ? a : b) 

[CMT] 

// blood cells
LT0
ST01
ST02
MPP01
MPP02
CMP01
CMP02
BFUE01
BFUE02
CFUE01
CFUE02
RET0
RBC0

// hemoglobins
alpha_RET0
beta_RET0
gamma_RET0
delta_RET0
alphabeta_RET0
alphagamma_RET0
alphadelta_RET0
HbA_RET0
HbF_RET0
HbA2_RET0
alpha_RBC0
beta_RBC0
gamma_RBC0
delta_RBC0
alphabeta_RBC0
alphagamma_RBC0
alphadelta_RBC0
HbA_RBC0
HbF_RBC0
HbA2_RBC0

// Lymphocytes progenitors
CLP0

// B cells in bone marrow
Boe0 // propreB cells
Bi0  // immature B cells

// spleen B cells
Bt0  // transitional B cells
BMspl0 // mature B cells in spleen

// mature recirculating cells
BMrec0 

// double negative thymocyte (4 stages)
N00
N10
N20
N30
N40

// double positive thymocytes (7 stages)
P00
P10
P20
P30
P40
P50
P60
P70

// single positive CD4 thymocyte
S400
S410
S420

// single positive CD8 thymocyte
S800
S810
S820

// T cells in blood
cd4rec0 
cd8rec0

// T cells in the lymph nodes and peripheral organs
cd4lym0
cd8lym0

cd4peripheral0
cd8peripheral0

// myeloid 
GMP01
GMP02
GM0

[PARAM]

// total mean residence time (days);
// mean residence time in each subcompartment will be computed later;
tauLT = 100
tauST = 20
tauMPP = 2
tauCMP = 4
tauBFUE = 7
tauCFUE = 7
tauRET = 3
tauRBCendo = 120
tauRBCtrans = 120 

ssLT = 1275 // total capacity for LT-HSC

// amplification parameters (A.U.);
aST = 1000
aMPP = 1000
aCMP = 36
aBFUE = 32

// replication rate for LT-HSC
rLT = (1/(2.5*7));

// hemoglobin related parameters
VRET = 0.09e-12 // reticulocyte volume; L 
VRBC = 0.09e-12 // erythrocyte volume; L
MWHh = 64500 // hemoglobin molecular weight, g/mol or unit in ng/nmol

// change MWHh, following suggestions in Pittman, 2011
// https://www.ncbi.nlm.nih.gov/books/NBK54103/
// MWHh = 64400

ksynalpha = 6e-7 // alpha globin synthesis rate; nmol/day/cell
ratiosynbeta = 0.5 // beta globin/ alpha globin synthesis rate ratio; unitless; ratio assumed based on gene copy number
ratiosyngamma = 0.03 // gamma globin/ alpha globin synthesis rate ratio
ratiosyndelta = 0.04 // delta globin/ alpha globin synthesis rate ratio
thalfmonomer = 0.25 // free monomer half life; days
kon = 1e-5*(60*60*24) // bimolecular binding on rate constant; unit nmol-1.day-1 (value adjusted for unit)
Kdalphabeta = 1e-3 // αβ dimer dissociation constant, nM; 
Kdalphagamma = 1e-5 // αγ dimer dissociation constant, nM
Kdalphadelta = 1e-2  // αδ dimer dissociation constant, nM
KdHbA = 100 // HbA (α₂β₂) tetramer dissociation constant, nM
KdHbF = 100 // HbF (α₂γ₂) tetramer dissociation constant, nM
KdHbA2 = 100 // HbA2 (α₂δ₂) tetramer dissociation rate, nM

HbA_saturation = 0.74
HbF_saturation = 0.88
HbA2_saturation = 0

// GMP and GM related dynamics
kCMP2GMP = 2.5 // differentiation rate, CMP -> GMP; unit day-1
betaGM = 3 // GM death rate; unit day-1
aGMP = 128 // GMP proliferation date, day-1
tauGMP = 0.12 //mean residence time, day-1

// CLP parameters
alphaMPP2CLP = 0.022 // differentiation rate, MPP -> CLP; unit day-1
betaCLP = 3 // CLP proliferation rate; unit day-1
alphaCLP2proB = 3 // differentiation rate, CLP -> proB; unit day-1
kCLP = 0.015 // CLP death rate, day-1 

CLP2proB_amp = 4 // dummy parameter; introduced to match the flow
CLP2DN_amp = 512 // dummy parameter; introduced to match the cell flow
alphaCLP2DN = 2.5e-4 // assumed, based on  Zlotoff and Bhandoola, 2012; https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3076003/

//--- B cell parameter ---//

K = 8.4e9  // propreB and immature B cell maximum capacity in bone marrow

// preproB cell rates (per day)
gamma = 2.2 // propreB cell self-renew
delta_oe = 0.5 * 4 // death rate of propreB cells
delta_r = 0.05 * 4 // immature B cells -> propreB cells

// immature B cell rates (per day)
mu_i = 0.1 * 4 // immature B cell death rate
delta_i_t = 0.6 * 4 // immature B cells -> transitional B cells
delta_i_re = 0.19 * 4  // immature B cells -> mature recirculating B cells 

// mature recirculating B cells (per day)
phi_s = 0.03 * 4 // mature B cell in spleen -> mature recirculting B cells
mu_re = 0.008 * 4 // mature B cell death
phi_BM = 0.94 * 4 // mature recirculating B cells -> mature B cells in spleen

// transitional B cells in spleen (per day)
mu_t = 0.03 * 4 // death rate of transitional B cells
delta_t = 0.03 * 4 // transitional B cells -> mature B cells in spleen

// mature B cells death in spleen (per day)
epsilon_spl = 0.032

//--- T cell parameters ---//

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
// delta = 0 
delta_dn = 0
delta_dp = 1e-6
delta_sp = 0

// removal rate last stage DP (LP)
muLP = 0.37

// death rate of T cells
death_naiveT_cd8 = 0.002
death_naiveT_cd4 = 0.002

// T cell leaving and re-entering the blood (all unit scaled to day-1)
enter_lymph = 0.07
enter_peripheral = 4.8
exit_lymph = 1.13
exit_peripheral_cd4= 0.55 
exit_peripheral_cd8= 1

// T cell proliferation in central and peripheral tissue
env_naiveT = 1e9 // Naive T cell Density for Half-Maximal Peripheral Proliferation
k_nT_pro = 3.2e8 // unit. cells.day-1
bloodvolume = 5  // unit in L
peripheralvolume = 60 // unit L

[MAIN]

// LT-HSC dynamics
double kLT2ST = 1/tauLT;
double deltaLT = (rLT - kLT2ST)/ssLT; 

// ST-HSC dynamics
double tau1ST = (log2(aST) -1)/log2(aST) * tauST;
double tau2ST = tauST - tau1ST;
double krep1ST = (aST/2-1)/ tau1ST;
double krep2ST = 1/ tau2ST;
double k12ST = aST/(2*tau1ST);
double kST2MPP = 2/tau2ST;

// MPP dynamics
double tau1MPP = (log2(aMPP)-1)/log2(aMPP) * tauMPP;
double tau2MPP = tauMPP - tau1MPP;
double krep1MPP = (aMPP/2-1)/ tau1MPP;
double krep2MPP = 1/ tau2MPP;
double k12MPP = aMPP/(2*tau1MPP);
double kMPP2CMP = 2/tau2MPP;

// CMP dynamics
double tau1CMP = (log2(aCMP)-1)/log2(aCMP) * tauCMP;
double tau2CMP = tauCMP - tau1CMP;
double krep1CMP = (aCMP/2-1)/ tau1CMP;
double krep2CMP = 1/ tau2CMP;
double k12CMP = aCMP/(2*tau1CMP);
double kCMP2BFUE = 2/tau2CMP;

// BFU-E dynamics
double tau1BFUE = (log2(aBFUE)-1)/log2(aBFUE) * tauBFUE;
double tau2BFUE = tauBFUE - tau1BFUE;
double krep1BFUE = (aBFUE/2-1)/ tau1BFUE;
double krep2BFUE = 1/ tau2BFUE;
double k12BFUE = aBFUE/(2*tau1BFUE);
double kBFUE2CFUE = 2/tau2BFUE;


// RET and RBC dynamics
double kRET2RBC = 1/tauRET;
double kdeathRBCendo = 1/tauRBCendo;
double kdeathRBCtrans = 1/tauRBCtrans;


/// monomer synthesis rate
double ksynbeta = ratiosynbeta * ksynalpha;
double ksyngamma = ratiosyngamma * ksynalpha;
double ksyndelta = ratiosyndelta * ksynalpha;

// monomer degredation rate
double kdeg = log(2)/ thalfmonomer; 

// dimer dissociation rate
double koffalphagamma = Kdalphagamma * kon; 
double koffalphadelta = Kdalphadelta * kon; 
double koffalphabeta = Kdalphabeta * kon; 

// tetramer dissociation rate
double koffHbA = KdHbA * kon; 
double koffHbF = KdHbF * kon; 
double koffHbA2 = KdHbA2 * kon; 

// total RET & RBC volume 
double VRET0 = VRET * RET0 ; 
double VRBC0 = VRBC * RBC0 ;

// GMP dynamics
double tau1GMP = ( log2(aGMP)-1 )/ log2(aGMP) * tauGMP;
double tau2GMP = tauGMP - tau1GMP;
double krep1GMP = ( aGMP/2-1 )/ tau1GMP;
double krep2GMP = 1/tau2GMP; 
double k12GMP = aGMP/(2*tau1GMP);
double kGMP2GM = 2/tau2GMP;

[ODE]

// hemoglobin concentration

double HbA = (HbA_RET0 + HbA_RBC0) * MWHh*1e-9/ (10*bloodvolume);  // unit in g/dL 
double HbF = (HbF_RET0 + HbF_RBC0) * MWHh*1e-9/ (10*bloodvolume); // unit in g/dL  
double HbA2 = (HbA2_RET0 + HbA2_RBC0) * MWHh*1e-9/ (10*bloodvolume);  // unit in g/dL 


// calculate O2 in the blood
double vO2 = ( HbA_saturation * HbA + HbF_saturation * HbF + HbA2_saturation * HbA2 ) * 1.34;

// calculate CFU-E amplification time
double aCFUE = 550 * exp( -0.23 * vO2 );

// CFU-E dynamics
double tau1CFUE = (log2(aCFUE)-1)/log2(aCFUE) * tauCFUE;
double tau2CFUE = tauCFUE - tau1CFUE;
double krep1CFUE = (aCFUE/2-1)/ tau1CFUE;
double krep2CFUE = 1/ tau2CFUE;
double k12CFUE = aCFUE/(2*tau1CFUE);
double kCFUE2RET = 2/tau2CFUE;


// LT-HSC compartment
dxdt_LT0 = rLT * LT0 - deltaLT * LT0 * LT0 - kLT2ST * LT0; // endogenous branch

// ST-HSC compartment
dxdt_ST01 = kLT2ST * LT0 + (krep1ST - k12ST) * ST01; // endogenous branch
dxdt_ST02 = k12ST * ST01 + (krep2ST - kST2MPP) * ST02; // endogenous branch

// MPP compartment
dxdt_MPP01 = kST2MPP * ST02 + (krep1MPP - k12MPP) * MPP01; // endogenous branch
dxdt_MPP02 = k12MPP * MPP01 + (krep2MPP - kMPP2CMP - alphaMPP2CLP) * MPP02; // endogenous branch

// CMP compartment
dxdt_CMP01 = kMPP2CMP * MPP02 + (krep1CMP - k12CMP) * CMP01;  // endogenous branch
dxdt_CMP02 = k12CMP * CMP01 + (krep2CMP - kCMP2BFUE - kCMP2GMP) * CMP02; // endogenous branch

// GMP compartment
dxdt_GMP01 = kCMP2GMP * CMP02 + (krep1GMP-k12GMP) * GMP01;
dxdt_GMP02 = k12GMP * GMP01 + (krep2GMP - kGMP2GM) * GMP02;

// GM compartment
dxdt_GM0 = kGMP2GM * GMP02 - betaGM * GM0; 

// BFU-E compartment
dxdt_BFUE01 = kCMP2BFUE * CMP02 + (krep1BFUE - k12BFUE) * BFUE01;   // endogenous branch
dxdt_BFUE02 = k12BFUE * BFUE01 + (krep2BFUE - kBFUE2CFUE) * BFUE02; // endogenous branch

// CFU-E compartment
dxdt_CFUE01 = kBFUE2CFUE * BFUE02 + (krep1CFUE - k12CFUE) * CFUE01; // endogenous branch
dxdt_CFUE02 = k12CFUE * CFUE01 + (krep2CFUE - kCFUE2RET) * CFUE02;  // endogenous branch

// RET compartment
dxdt_RET0 = kCFUE2RET * CFUE02 - kRET2RBC * RET0;  // endogenous branch

// RBC
dxdt_RBC0 = kRET2RBC * RET0 - kdeathRBCendo * RBC0;  // endogenous branch

// monomer in endogenous reticulocyte; 
dxdt_alpha_RET0 = ksynalpha * RET0 - (kdeg + kRET2RBC)*alpha_RET0 - (kon*alpha_RET0*beta_RET0 + kon*alpha_RET0*gamma_RET0 + kon*alpha_RET0*delta_RET0) + (koffalphabeta*alphabeta_RET0 + koffalphagamma*alphagamma_RET0 + koffalphadelta*alphadelta_RET0)*VRET0;
dxdt_beta_RET0  = ksynbeta  * RET0 - (kdeg + kRET2RBC)*beta_RET0  - kon * alpha_RET0 * beta_RET0  + koffalphabeta *VRET0 * alphabeta_RET0  ; 
dxdt_gamma_RET0 = ksyngamma * RET0 - (kdeg + kRET2RBC)*gamma_RET0 - kon * alpha_RET0 * gamma_RET0 + koffalphagamma      *VRET0 * alphagamma_RET0 ;
dxdt_delta_RET0 = ksyndelta * RET0 - (kdeg + kRET2RBC)*delta_RET0 - kon * alpha_RET0 * delta_RET0 + koffalphadelta      *VRET0 * alphadelta_RET0 ; 


// dimers in endogenous reticulocyte; 
dxdt_alphabeta_RET0  = kon * alpha_RET0 * beta_RET0   - kRET2RBC * alphabeta_RET0  - koffalphabeta * alphabeta_RET0  * VRET0 - 2 * kon * alphabeta_RET0  * alphabeta_RET0  + 2 * koffHbA  *VRET0 * HbA_RET0 ;
dxdt_alphagamma_RET0 = kon * alpha_RET0 * gamma_RET0  - kRET2RBC * alphagamma_RET0 - koffalphagamma      * alphagamma_RET0 * VRET0 - 2 * kon * alphagamma_RET0 * alphagamma_RET0 + 2 * koffHbF  *VRET0 * HbF_RET0 ; 
dxdt_alphadelta_RET0 = kon * alpha_RET0 * delta_RET0  - kRET2RBC * alphadelta_RET0 - koffalphadelta      * alphadelta_RET0 * VRET0 - 2 * kon * alphadelta_RET0 * alphadelta_RET0 + 2 * koffHbA2 *VRET0 * HbA2_RET0; 

// tetramer in endogenous reticulocyte; 
dxdt_HbA_RET0  = kon * alphabeta_RET0  * alphabeta_RET0  - (koffHbA  *VRET0 + kRET2RBC) * HbA_RET0;
dxdt_HbF_RET0  = kon * alphagamma_RET0 * alphagamma_RET0 - (koffHbF  *VRET0 + kRET2RBC) * HbF_RET0;
dxdt_HbA2_RET0 = kon * alphadelta_RET0 * alphadelta_RET0 - (koffHbA2 *VRET0 + kRET2RBC) * HbA2_RET0;

// monomer in endogenous RBC; 
dxdt_alpha_RBC0 = kRET2RBC * alpha_RET0 - (kdeg + kdeathRBCendo) * alpha_RBC0 - (kon*alpha_RBC0*beta_RBC0 + kon*alpha_RBC0*gamma_RBC0 + kon*alpha_RBC0*delta_RBC0) + (koffalphabeta * alphabeta_RBC0 + koffalphagamma * alphagamma_RBC0 + koffalphadelta * alphadelta_RBC0) * VRBC0; 
dxdt_beta_RBC0 =  kRET2RBC * beta_RET0  - (kdeg + kdeathRBCendo) * beta_RBC0  - kon * alpha_RBC0 * beta_RBC0  + koffalphabeta *VRBC0 * alphabeta_RBC0  ;
dxdt_gamma_RBC0 = kRET2RBC * gamma_RET0 - (kdeg + kdeathRBCendo) * gamma_RBC0 - kon * alpha_RBC0 * gamma_RBC0 + koffalphagamma      *VRBC0 * alphagamma_RBC0 ;
dxdt_delta_RBC0 = kRET2RBC * delta_RET0 - (kdeg + kdeathRBCendo) * delta_RBC0 - kon * alpha_RBC0 * delta_RBC0 + koffalphadelta      *VRBC0 * alphadelta_RBC0 ;

// dimers in endogenous RBC
dxdt_alphabeta_RBC0  = kRET2RBC * alphabeta_RET0  - kdeathRBCendo * alphabeta_RBC0  + kon * alpha_RBC0 * beta_RBC0  - koffalphabeta *VRBC0 * alphabeta_RBC0  - 2*kon* alphabeta_RBC0 *alphabeta_RBC0  + 2*koffHbA *VRBC0 * HbA_RBC0;
dxdt_alphagamma_RBC0 = kRET2RBC * alphagamma_RET0 - kdeathRBCendo * alphagamma_RBC0 + kon * alpha_RBC0 * gamma_RBC0 - koffalphagamma      *VRBC0 * alphagamma_RBC0 - 2*kon* alphagamma_RBC0*alphagamma_RBC0 + 2*koffHbF *VRBC0 * HbF_RBC0;
dxdt_alphadelta_RBC0 = kRET2RBC * alphadelta_RET0 - kdeathRBCendo * alphadelta_RBC0 + kon * alpha_RBC0 * delta_RBC0 - koffalphadelta      *VRBC0 * alphadelta_RBC0 - 2*kon* alphadelta_RBC0*alphadelta_RBC0 + 2*koffHbA2*VRBC0 * HbA2_RBC0;

// tetramers in endogenous RBC
dxdt_HbA_RBC0  = kRET2RBC * HbA_RET0  - kdeathRBCendo * HbA_RBC0  + kon*alphabeta_RBC0 *alphabeta_RBC0  - koffHbA * VRBC0 * HbA_RBC0;
dxdt_HbF_RBC0  = kRET2RBC * HbF_RET0  - kdeathRBCendo * HbF_RBC0  + kon*alphagamma_RBC0*alphagamma_RBC0 - koffHbF * VRBC0 * HbF_RBC0; 
dxdt_HbA2_RBC0 = kRET2RBC * HbA2_RET0 - kdeathRBCendo * HbA2_RBC0 + kon*alphadelta_RBC0*alphadelta_RBC0 - koffHbA2* VRBC0 * HbA2_RBC0; 

// CLP compartment
dxdt_CLP0 = alphaMPP2CLP * MPP02 - kCLP * CLP0 + betaCLP * CLP0 - alphaCLP2proB * CLP0 - alphaCLP2DN * CLP0; 

// B cell models; assume CLP -> propreB has another round of division
dxdt_Boe0 = CLP2proB_amp * alphaCLP2proB * CLP0 + ( gamma * ( 1 -  (Boe0 + BMrec0)/K ) - delta_oe ) * Boe0 + delta_r * Bi0; 

dxdt_Bi0 = delta_oe * Boe0 - (mu_i + delta_i_t + delta_r + delta_i_re) * Bi0;

dxdt_BMrec0 = delta_i_re * Bi0 + phi_s * BMspl0 - (mu_re + phi_BM) * BMrec0; 

dxdt_Bt0 = delta_i_t * Bi0 - (mu_t + delta_t) * Bt0;

dxdt_BMspl0 = delta_t * Bt0 + phi_BM * BMrec0 - (phi_s + epsilon_spl) * BMspl0; 

// dynamics of DN cells
dxdt_N00 = CLP2DN_amp * alphaCLP2DN * CLP0 - (pN + delta_dn) * N00; 
dxdt_N10 =  2 * pN * N00 - ( pN + delta_dn +  min( pow(alpha_muN * 1, n), 100) ) * N10; 
dxdt_N20 =  2 * pN * N10 - ( pN + delta_dn +  min( pow(alpha_muN * 2, n), 100) ) * N20; 
dxdt_N30 =  2 * pN * N20 - ( pN + delta_dn +  min( pow(alpha_muN * 3, n), 100) ) * N30; 
dxdt_N40 =  2 * pN * N30 - ( pN + delta_dn +  min( pow(alpha_muN * 4, n), 100) ) * N40; 

// dynamics of DP cells
double sum_mu_N0 = min( pow(alpha_muN * 1, n), 100) * N10 + min( pow(alpha_muN * 2, n), 100) * N20 + min( pow(alpha_muN * 3, n), 100) * N30 + min( pow(alpha_muN * 4, n), 100) * N40; 

dxdt_P00 = sum_mu_N0 + 2 * pN * N40 - (pP + delta_dp) * P00; 
dxdt_P10 = 2 * pP * P00 - ( pP + delta_dp + min( pow(alpha_muP * 1, n), 100) ) * P10; 
dxdt_P20 = 2 * pP * P10 - ( pP + delta_dp + min( pow(alpha_muP * 2, n), 100) ) * P20; 
dxdt_P30 = 2 * pP * P20 - ( pP + delta_dp + min( pow(alpha_muP * 3, n), 100) ) * P30; 
dxdt_P40 = 2 * pP * P30 - ( pP + delta_dp + min( pow(alpha_muP * 4, n), 100) ) * P40; 
dxdt_P50 = 2 * pP * P40 - ( pP + delta_dp + min( pow(alpha_muP * 5, n), 100) ) * P50; 
dxdt_P60 = 2 * pP * P50 - ( pP + delta_dp + min( pow(alpha_muP * 6, n), 100) ) * P60; 

double sum_mu_P0 = min( pow(alpha_muP * 1, n), 100) * P10 +  min( pow(alpha_muP * 2, n), 100) * P20 + min( pow(alpha_muP * 3, n), 100) * P30 +
                   min( pow(alpha_muP * 4, n), 100) * P40 +  min( pow(alpha_muP * 5, n), 100) * P50 + min( pow(alpha_muP * 6, n), 100) * P60; 

dxdt_P70 = sum_mu_P0 + 2 * pP * P60 - muLP * P70;

// dynamics of CD4+ cells
dxdt_S400 = alpha4 * muLP * P70 - (pS + delta_sp) * S400; 
dxdt_S410 = 2 * pS * S400 - (pS + delta_sp + min( pow(alpha_e * 1, n), 100) ) * S410; 
dxdt_S420 = 2 * pS * S410 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S420; 

double sum_SP40 =  min( pow(alpha_e * 1, n), 100) * S410 +  min( pow(alpha_e * 2, n), 100) * S420; 

// dynamics of CD8+ cells

dxdt_S800 = alpha8 * muLP * P70 - (pS + delta_sp) * S800; 
dxdt_S810 = 2 * pS * S800 - (pS + delta_sp +  min( pow(alpha_e * 1, n), 100) ) * S810; 
dxdt_S820 = 2 * pS * S810 - (delta_sp + min( pow(alpha_e * 2, n), 100) ) * S820; 

double sum_SP80 =  min( pow(alpha_e * 1, n), 100) * S810 +  min( pow(alpha_e * 2, n), 100) * S820; 

// dynamics of T cells in blood
dxdt_cd4rec0 = sum_SP40 - death_naiveT_cd4 * cd4rec0 + exit_lymph * cd4lym0 - enter_lymph * cd4rec0 + exit_peripheral_cd4 * cd4peripheral0 - enter_peripheral * cd4rec0 + k_nT_pro/bloodvolume * (cd4rec0/bloodvolume)/(cd4rec0/bloodvolume + env_naiveT);
dxdt_cd8rec0 = sum_SP80 - death_naiveT_cd8 * cd8rec0 + exit_lymph * cd8lym0 - enter_lymph * cd8rec0 + exit_peripheral_cd8 * cd8peripheral0 - enter_peripheral * cd8rec0 + k_nT_pro/bloodvolume * (cd8rec0/bloodvolume)/(cd8rec0/bloodvolume + env_naiveT);

// dynamics of T cells in peripheral tissues
dxdt_cd4peripheral0 = k_nT_pro/peripheralvolume * (cd4peripheral0/peripheralvolume)/(cd4peripheral0/peripheralvolume + env_naiveT) - death_naiveT_cd4 * cd4peripheral0 - exit_peripheral_cd4 * cd4peripheral0 + enter_peripheral * cd4rec0;
dxdt_cd8peripheral0 = k_nT_pro/peripheralvolume * (cd8peripheral0/peripheralvolume)/(cd8peripheral0/peripheralvolume + env_naiveT) - death_naiveT_cd8 * cd8peripheral0 - exit_peripheral_cd8 * cd8peripheral0 + enter_peripheral * cd8rec0;

// dynamics of T cells in secondary lymphoid system
dxdt_cd4lym0 = enter_lymph * cd4rec0 - exit_lymph * cd4lym0;
dxdt_cd8lym0 = enter_lymph * cd8rec0 - exit_lymph * cd8lym0; 

[TABLE]
// progenitors in thymus
capture DN = N00 + N10 + N20 + N30 + N40; 
capture DP = P00 + P10 + P20 + P30 + P40 + P50 + P60 + P70; 
capture SP4thymus = S400 + S410 + S420;
capture SP8thymus = S800 + S810 + S820;

// B cells in spleen
capture splenicB = Bt0 + BMspl0;

// T cells in other ograns
capture nT_lymphnodes = cd4lym0 + cd8lym0;
capture nT_peripheral = cd4peripheral0 + cd8peripheral0 ;

// cells & proteins in peripheral blood
capture RBCconc = (RBC0)/ (bloodvolume * 1e6); // RBC per uL
capture RETconc = (RET0)/ (bloodvolume * 1e6); // RET per uL
capture Bconc = (BMrec0)/(bloodvolume * 1e6); // naive B per uL
capture Tconc = (cd4rec0 +  cd8rec0)/(bloodvolume * 1e6); // naive T per uL
capture GMconc = (GM0)/(bloodvolume * 1e6); // granulocytes per uL
capture WBCconc = Bconc + Tconc + GMconc; // white blood cell count per uL blood

capture totalHb = HbA + HbF + HbA2; // unit: g/dL
capture HbinRBC = (HbA_RBC0 + HbF_RBC0 + HbA2_RBC0) * MWHh*1e-9 / (VRBC0); // unit: g/L


// fluxes
capture CLP2proB0 = CLP2proB_amp * alphaCLP2proB * CLP0; // CLP -> propreB cells; 
capture CLPexport2thymus = alphaCLP2DN * CLP0; 

// fluxes related to thymus
capture thymic_output = sum_SP40 + sum_SP80;
capture naiveTprof = k_nT_pro/peripheralvolume * (cd4peripheral0/peripheralvolume)/(cd4peripheral0/peripheralvolume + env_naiveT) + 
                     k_nT_pro/peripheralvolume * (cd8peripheral0/peripheralvolume)/(cd8peripheral0/peripheralvolume + env_naiveT) + 
                     k_nT_pro/bloodvolume * (cd4rec0/bloodvolume)/(cd4rec0/bloodvolume + env_naiveT) + 
                     k_nT_pro/bloodvolume * (cd8rec0/bloodvolume)/(cd8rec0/bloodvolume + env_naiveT);
                     
[CAPTURE]
aCFUE, sum_SP40, sum_SP80
