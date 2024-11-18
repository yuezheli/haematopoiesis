[PROB]

This file is copied from human_ery_lymp_mk_neutrophil.cpp

This file incorporates the ADC-induced toxicity on myelocytes (basd on AGS5-mcMMAF) and CFU-MK (based on T-DM1). 

PK backbone was based on T-DM1 (from Scheuher et al., 2023). 

[SET]
delta = 0.5, end = 70

[GLOBAL]  
#define min(a,b)  ( a<b ? a : b) 
#define max(a,b)  ( a<b ? b : a) 

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
GMP
myeloblast
promyelocyte
myelocyte01
myelocyte02
metamyelocyte
bandcell
neutrophil_blood
neutrophil_pool

// megakaryocyte 
BFUMK01
BFUMK02
CFUMK01
CFUMK02
MK0
platelet0
TPO

// ADC 
adc_central
adc_peripheral
pl_central
pl_peripheral
c_sTAA
c_sTAA_adc
p_sTAA
p_sTAA_adc

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
betaGM = 2.5 // GM death rate; unit day-1; estimated from half life of 6.7h; https://pubmed.ncbi.nlm.nih.gov/28303306/
tauGMP = 0.12 //mean residence time, day-1; https://www.nature.com/articles/nature14242
a_myelocyte = 32  // back calculate from myelocyte -> metamyelocyte transition time; https://www.nature.com/articles/nature14242
tau_myeloblast = 2.08 // mean residence time of proliferating neutrophil progenitors, [day], https://pubmed.ncbi.nlm.nih.gov/28303306/
tau_promyelocyte = 2.08 // mean residence time of proliferating neutrophil progenitors, [day], https://pubmed.ncbi.nlm.nih.gov/28303306/
tau_myelocyte = 2.08 // mean residence time of proliferating neutrophil progenitors, [day], https://pubmed.ncbi.nlm.nih.gov/28303306/
tau_metamyelocyte = 0.5 // mean residence time for nonproliferating progenitor cells for neutrophil (metamyelocyte & band cell), [day],
k_pool_neutrophil = 10 // rate for neutrophil to bind to/ unbind from blood vessel, [day-1], assumed, https://pubmed.ncbi.nlm.nih.gov/28303306/

// MK and platelet related dynamics 
kCMP2BFUMK = 0.3  // estimated from https://ghe.metrumrg.com/yuezhel/Thrombocytopenia-ADC
aBFUMK = 8
aCFUMK0 = 8
tauBFUMK = 9.0
tauCFUMK = 12.0
tauMK = 0.1
tauPlatelet = 7.0
thalfTPO = 1.25
kMKplatelet = 480.
alpha1 = -4.34
beta1 = 0.009
beta2 = 4.72E-3

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

// ADC
IC50_adc_cfumk = 50 // [nM]
n_adc = 1.13
k_in_adc_cfumk = 0.37 // ADC in-rate into CFU-MK; [day-1]
IC50_adc_myelocyte = 14 // [nM]
k_in_adc_myelocyte = 3 // [day-1]

// PK; from https://link.springer.com/article/10.1007/s10928-023-09884-6 unless specified
Vcent = 3 // [L]
Vperi = 13 // [L]
Q = 0.196 // [L/day]
Q_sTAA = 0.46  // [L/day]
thalf_adc = 11.6 // [day]
kdec = 0.00144 // T-Dxd deconjugation; [day-1]
kon_TAA = 8.64
kd_TAA = 0.3 // [nM]
thalf_sTAA = 0.21 // [day]
thalf_sTAA_adc = 11.6 // assumed to be the same as ADC; [day]
thalf_pl = 0.85/24 // Dxd half life; tuned [day]
base_sTAA = 1.16 // tuned; [nM]
KpBM = 0.07 // partition coefficient, plasma:bone marrow
MW_sTAA = 1e5
MW_pl = 500
DAR = 8
// t_dist_12_pl = 0.1/24  // distribution half life of Dxd, [day]
P_dist_12_pl = 83.02  // partition coef of Dxd
// k12_pl = log(2) / (0.1/24) * 83.02 / ( 83.02 + 3/13 )
// k21_pl = log(2) / (0.1/24) * 3/13 / ( 83.02 + 3/13 )
Q_pl = 166 // computed from distribution half life [L/day]
IC50_pl_neut = 9.54 // [nM]
emax_pl_neut = 0.201/24 // [day]
k_in_pl = 46.08/24 // [day]; https://www.metrumrg.com/wp-content/uploads/2023/10/poster-icsb-2023-YL-compressed.pdf
k_out_pl = 32.32/24 // [day]; https://www.metrumrg.com/wp-content/uploads/2023/10/poster-icsb-2023-YL-compressed.pdf

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

// BFU-MK dynamics 
double tau1BFUMK = (log2(aBFUMK)-1)/log2(aBFUMK) * tauBFUMK;
double tau2BFUMK = tauBFUMK - tau1BFUMK;
double krep1BFUMK = (aBFUMK/2-1)/ tau1BFUMK;
double krep2BFUMK = 1/ tau2BFUMK;
double k12BFUMK = aBFUMK/(2*tau1BFUMK);
double kBFUMK2CFUMK = 2/tau2BFUMK;

// CFU-MK dynamics
double aCFUMK = aCFUMK0 * log(2.73 + beta2 * TPO) / log(1.38) * ( pow(IC50_adc_cfumk, n_adc)/( pow(IC50_adc_cfumk, n_adc) + pow(max(0, adc_central * KpBM), n_adc) ) ); 
double krep1CFUMK = 0; 
double krep2CFUMK = 0;
double k12CFUMK = 0;
double kCFUMK2MK = 0;
if (aCFUMK > 2.1){
  double tau1CFUMK = (log2(aCFUMK)-1)/log2(aCFUMK) * tauCFUMK;
  double tau2CFUMK = tauCFUMK - tau1CFUMK;
  krep1CFUMK = (aCFUMK/2-1)/ tau1CFUMK;
  krep2CFUMK = 1/ tau2CFUMK;
  k12CFUMK = aCFUMK/(2*tau1CFUMK);
  kCFUMK2MK = 2/tau2CFUMK;
} 
  
  
// MK and platelet dynamics 
double kdeathMK = 1/tauMK;
double kdeathPlatelet = 1/tauPlatelet;

// TPO dynamics
double ksynTPO = exp(alpha1 - beta1 * platelet0/bloodvolume * 1E-9);  // [nmol]
double kdegTPO = log(2)/thalfTPO;

// myelocyte dynamics 
double amp_myelocyte = a_myelocyte * pow(IC50_adc_myelocyte, n_adc) / ( pow(IC50_adc_myelocyte, n_adc) + pow(max(0, adc_central * KpBM), n_adc) ); 
double krep1_myelocyte = 0; 
double krep2_myelocyte = 0; 
double k12_myelocyte = 0; 
double k_myelocyte2metamyelocyte = 0; 
if (amp_myelocyte > 2.1){
  double tau1_myelocyte = (log2(amp_myelocyte)-1)/log2(amp_myelocyte) * tau_myelocyte;
  double tau2_myelocyte = tau_myelocyte - tau1_myelocyte;
  krep1_myelocyte = (amp_myelocyte/2-1)/ tau1_myelocyte;
  krep2_myelocyte = 1/ tau2_myelocyte;
  k12_myelocyte = amp_myelocyte/(2*tau1_myelocyte);
  k_myelocyte2metamyelocyte = 2/tau2_myelocyte;
}

// ADC PK 
double kelim = log(2)/ thalf_adc; 
double kelim_pl = log(2)/ thalf_pl; 
double koff_taa = kon_TAA * kd_TAA;
double kdeg_sTAA = log(2)/thalf_sTAA;
double kdeg_sTAA_adc = log(2)/thalf_sTAA_adc;
double ksyn_sTAA = kdeg_sTAA * base_sTAA;

// initial condition 
c_sTAA_0 = 8E3/MW_sTAA; 


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
dxdt_CMP02 = k12CMP * CMP01 + (krep2CMP - kCMP2BFUE - kCMP2GMP - kCMP2BFUMK) * CMP02; // endogenous branch

// MK dynamics 
dxdt_BFUMK01 = kCMP2BFUMK * CMP02 + (krep1BFUMK - k12BFUMK) * BFUMK01;  
dxdt_BFUMK02 = k12BFUMK * BFUMK01 + (krep2BFUMK - kBFUMK2CFUMK) * BFUMK02; 
dxdt_CFUMK01 = kBFUMK2CFUMK * BFUMK02 + (krep1CFUMK - k12CFUMK) * CFUMK01; 
dxdt_CFUMK02 = k12CFUMK * CFUMK01 + (krep2CFUMK - kCFUMK2MK) * CFUMK02;
dxdt_MK0   = kCFUMK2MK * CFUMK02 - kdeathMK * MK0;  
dxdt_platelet0   = kMKplatelet * MK0 - kdeathPlatelet * platelet0; 
dxdt_TPO = ksynTPO - kdegTPO * TPO;

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

// neutrophil compartment
dxdt_GMP = kCMP2GMP * CMP02 - 1/tauGMP * GMP;

dxdt_myeloblast = 2 * 1/tauGMP * GMP - 1/tau_myeloblast * myeloblast - emax_pl_neut*pl_central/(IC50_pl_neut + pl_central)*myeloblast;

dxdt_promyelocyte = 2 * 1/tau_myeloblast * myeloblast - 1/tau_promyelocyte * promyelocyte - emax_pl_neut*pl_central/(IC50_pl_neut + pl_central)*promyelocyte; 

dxdt_myelocyte01 = 2 * 1/tau_promyelocyte * promyelocyte + (krep1_myelocyte - k12_myelocyte) * myelocyte01 - emax_pl_neut*pl_central/(IC50_pl_neut + pl_central)*myelocyte01; 
dxdt_myelocyte02 = k12_myelocyte * myelocyte01 + (krep2_myelocyte - k_myelocyte2metamyelocyte) * myelocyte02 - emax_pl_neut*pl_central/(IC50_pl_neut + pl_central)*myelocyte02; 

dxdt_metamyelocyte = k_myelocyte2metamyelocyte * myelocyte02 - 1/tau_metamyelocyte * metamyelocyte;
dxdt_bandcell = 1/tau_metamyelocyte * metamyelocyte - 1/tau_metamyelocyte * bandcell;
dxdt_neutrophil_blood = 1/tau_metamyelocyte * bandcell - betaGM * neutrophil_blood - k_pool_neutrophil * neutrophil_blood + k_pool_neutrophil * neutrophil_pool - emax_pl_neut*pl_central/(IC50_pl_neut + pl_central)*neutrophil_blood;
dxdt_neutrophil_pool = k_pool_neutrophil * neutrophil_blood - k_pool_neutrophil * neutrophil_pool - emax_pl_neut*pl_central/(IC50_pl_neut + pl_central)*neutrophil_pool;

// ADC PK 
dxdt_adc_central = -Q/Vcent*(adc_central - adc_peripheral) - kelim*adc_central - kdec*adc_central - kon_TAA*adc_central*c_sTAA + koff_taa*c_sTAA_adc ;
dxdt_adc_peripheral = Q/Vperi*(adc_central - adc_peripheral) - kelim*adc_peripheral - kdec*adc_peripheral - kon_TAA*adc_peripheral*p_sTAA + koff_taa*p_sTAA_adc;
dxdt_c_sTAA = ksyn_sTAA - kdeg_sTAA*c_sTAA - kon_TAA*adc_central*c_sTAA + koff_taa*c_sTAA_adc -Q_sTAA/Vcent*(c_sTAA - p_sTAA);
dxdt_c_sTAA_adc = kon_TAA*adc_central*c_sTAA - koff_taa*c_sTAA_adc - kdeg_sTAA_adc*c_sTAA_adc -Q/Vcent*(c_sTAA_adc - p_sTAA_adc);
dxdt_p_sTAA = ksyn_sTAA - kdeg_sTAA*p_sTAA - kon_TAA*adc_peripheral*p_sTAA + koff_taa*p_sTAA_adc + Q_sTAA/Vperi*(c_sTAA - p_sTAA);
dxdt_p_sTAA_adc = kon_TAA*adc_peripheral*p_sTAA - koff_taa*p_sTAA_adc - kdeg_sTAA_adc*p_sTAA_adc + Q/Vperi*(c_sTAA_adc - p_sTAA_adc);
dxdt_pl_central = DAR*(kelim*adc_central + kdec*adc_central + kdeg_sTAA_adc*c_sTAA_adc) - Q_pl/Vcent*(pl_central - pl_peripheral/P_dist_12_pl) - kelim_pl*pl_central; 
dxdt_pl_peripheral = DAR*(kelim*adc_peripheral + kdec*adc_peripheral + kdeg_sTAA_adc*p_sTAA_adc) + Q_pl/Vperi*(pl_central - pl_peripheral/P_dist_12_pl) * Vcent/Vperi - kelim_pl*pl_peripheral;

[TABLE]
// cells & proteins in peripheral blood
capture RBCconc = (RBC0)/ (bloodvolume * 1e6); // RBC per uL
capture RETconc = (RET0)/ (bloodvolume * 1e6); // RET per uL
capture Bconc = (BMrec0)/(bloodvolume * 1e6); // naive B per uL
capture Tconc = (cd4rec0 + cd8rec0)/(bloodvolume * 1e6); // naive T per uL
capture Neuconc = (neutrophil_blood)/(bloodvolume * 1e6); // neutrophil per uL

capture totalHb = HbA + HbF + HbA2; // unit: g/dL
capture HbinRBC = (HbA_RBC0 + HbF_RBC0 + HbA2_RBC0) * MWHh*1e-9 / (VRBC0); // unit: g/L

capture PLTconc = platelet0 /(bloodvolume * 1e6); // platelet T per uL
capture TPO_nM = TPO / bloodvolume; // thrombopoietin concentration in blood [nM]

// fluxes related to thymus
capture thymic_output = sum_SP40 + sum_SP80;
capture naiveTprof = k_nT_pro/peripheralvolume * (cd4peripheral0/peripheralvolume)/(cd4peripheral0/peripheralvolume + env_naiveT) + 
                     k_nT_pro/peripheralvolume * (cd8peripheral0/peripheralvolume)/(cd8peripheral0/peripheralvolume + env_naiveT) + 
                     k_nT_pro/bloodvolume * (cd4rec0/bloodvolume)/(cd4rec0/bloodvolume + env_naiveT) + 
                     k_nT_pro/bloodvolume * (cd8rec0/bloodvolume)/(cd8rec0/bloodvolume + env_naiveT);
                     
[CAPTURE]
k_myelocyte2metamyelocyte
