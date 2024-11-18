[PROB] 

This model is from Shahaf et al., 2016. 
https://www.frontiersin.org/articles/10.3389/fimmu.2016.00077/full

This models BrdU labeling

[SET]
delta = 0.1

[CMT]

// B cells in bone marrow
UBoe 
UBi  
LBoe
LBi

// spleen B cells
UBt  
UBMspl 
LBt  
LBMspl 

// mature recirculating cells
UBMrec 
LBMrec

[PARAM]

s = 3e5 // cell influx into B cells, count of cells per 6h
K = 6e6  // propreB and immature B cell maximum capacity in bone marrow

// preproB cell rates (per 6h)
gamma = 0.3  // propreB cell self-renew
delta_oe = 0.5  // death rate of propreB cells
delta_r = 0.05  // immature B cells -> propreB cells

// immature B cell rates (per 6h)
mu_i = 0.1  // immature B cell death rate
delta_i_t = 0.6  // immature B cells -> trasitional B cells
delta_i_re = 0.19  // immature B cells -> mature recirculating B cells 

// mature recirculating B cells (per 6h)
phi_s = 0.03 // mature B cell in spleen -> mature recirculting B cells
mu_re = 0.008  // mature B cell death
phi_BM = 0.94  // mature B cell in spleen -> mature recirculating B cells

// transitional B cells in spleen (per 6h)
mu_t = 0.03 // death rate of transitional B cells
delta_t = 0.03 // transitional B cells -> mature B cells in spleen

// mature B cells in spleen (per 6h)
epsilon_spl = 0.008

// add additional dummy parameter for output
dvtype = 1 // this parameter needs to be changed

[ODE]

// total number of B cells

double Boe = UBoe + LBoe; 
double Bi = UBi + LBi; 
double Bt = UBt + LBt; 
double BMrec = UBMrec + LBMrec; 
double BMspl = UBMspl + LBMspl;

// ODE for unlabeled cells

dxdt_UBoe = s - ( gamma * ( 1 -  (Boe + BMrec)/K ) + delta_oe ) * UBoe + delta_r * UBi; 

dxdt_UBi = delta_oe * UBoe - (mu_i + delta_i_t + delta_r + delta_i_re) * UBi;

dxdt_UBMrec = delta_i_re * UBi + phi_s * UBMspl - (mu_re + phi_BM) * UBMrec; 

dxdt_UBt = delta_i_t * UBi - (mu_t + delta_t) * UBt;

dxdt_UBMspl = delta_t * UBt + phi_BM * UBMrec - (phi_s + epsilon_spl) * UBMspl; 

// ODE for labeled cells

dxdt_LBoe = ( 2*UBoe + LBoe ) * gamma * ( 1 -  (Boe + BMrec)/K )  - delta_oe * LBoe + delta_r * LBi; 

dxdt_LBi = delta_oe * LBoe - (mu_i + delta_i_t + delta_r + delta_i_re) * LBi;

dxdt_LBMrec = delta_i_re * LBi + phi_s * LBMspl - (mu_re + phi_BM) * LBMrec; 

dxdt_LBt = delta_i_t * LBi - (mu_t + delta_t) * LBt;

dxdt_LBMspl = delta_t * LBt + phi_BM * LBMrec - (phi_s + epsilon_spl) * LBMspl; 

[TABLE]
capture BrdUimmature = LBi/Bi; 
capture BrdUmaturerec = LBMrec/BMrec; 
capture BrdUtransit =  LBt/ Bt;
capture BrdUsplenicmature = LBMspl/ BMspl;

double DV;
if(dvtype == 1){ DV =  LBi/Bi; }
else if (dvtype == 2){ DV =  LBMrec/BMrec; }
else if (dvtype == 3){ DV =  LBt/ Bt; }
else if (dvtype == 4){ DV =  LBMspl/ BMspl; }

[CAPTURE]

DV