[PROB] 

Humanized version of B cell model. 

changes: 
1. K: 6 millions -> 8.4 billions
2. gamma: 1.2 -> 2.2
3. s: 1.2M -> 22M (based on the scaled up model)

Rational see mrgsolve/mouse2human readme.md

[SET]
delta = 1


[CMT]

// B cells in bone marrow
Boe // propreB cells
Bi  // immature B cells

// spleen B cells
Bt  // transitional B cells
BMspl // mature B cells in spleen

// mature recirculating cells
BMrec 

[PARAM]

s = 22e6 // cell influx into B cells, count of cells per day
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

// mature B cells in spleen (per day)
epsilon_spl = 0.008 * 4

[ODE]

dxdt_Boe = s + ( gamma * ( 1 -  (Boe + BMrec)/K ) - delta_oe ) * Boe + delta_r * Bi; 

dxdt_Bi = delta_oe * Boe - (mu_i + delta_i_t + delta_r + delta_i_re) * Bi;

dxdt_BMrec = delta_i_re * Bi + phi_s * BMspl - (mu_re + phi_BM) * BMrec; 

dxdt_Bt = delta_i_t * Bi - (mu_t + delta_t) * Bt;

dxdt_BMspl = delta_t * Bt + phi_BM * BMrec - (phi_s + epsilon_spl) * BMspl; 

[TABLE]

capture Bcell_spleen = Bt + BMspl;
capture BcellBM = Boe + Bi;