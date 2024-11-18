# date: May 2, 2024
# author: Yuezhe Li 
# purpose of this code: to model platelet development in human (add TPO synthesis)
# model based on ABM poster (see doc folder)
# model assumes synthesis rate of TPO is exponentially negatively correlated with platelet number (Engle et al., 1999; https://pubmed.ncbi.nlm.nih.gov/10354155/)
# model also assume that TPO logrithmic promote CFU-MK proliferation (Broudy et al., 1996; https://pubmed.ncbi.nlm.nih.gov/8822921/)

const V_plasma = 3.126 # [L]
const V_blood = 5 # [L]

using DifferentialEquations, ComponentArrays
using Parameters: @unpack

function mk_dev!(du, u, p, t)
    @unpack LT0, ST01, ST02, MPP01, MPP02, CMP01, CMP02, BFUMK01, BFUMK02, CFUMK01, CFUMK02, MK0, platelet0, TPO = u; 
    @unpack tauLT, tauST, tauMPP, tauCMP, tauBFUMK, tauCFUMK, tauMK, tauPlatelet, thalfTPO, ssLT, rLT, aST, aMPP, aCMP, aBFUMK, aCFUMK0, kMKplatelet, alpha1, beta1, beta2 = p

    # LT-HSC dynamics
    kLT2ST = 1/tauLT;
    deltaLT = (rLT - kLT2ST)/ssLT; 
    # ST-HSC dynamics
    tau1ST = (log2(aST) -1)/log2(aST) * tauST;
    tau2ST = tauST - tau1ST;
    krep1ST = (aST/2-1)/ tau1ST;
    krep2ST = 1/ tau2ST;
    k12ST = aST/(2*tau1ST);
    kST2MPP = 2/tau2ST;
    # MPP dynamics
    tau1MPP = (log2(aMPP)-1)/log2(aMPP) * tauMPP;
    tau2MPP = tauMPP - tau1MPP;
    krep1MPP = (aMPP/2-1)/ tau1MPP;
    krep2MPP = 1/ tau2MPP;
    k12MPP = aMPP/(2*tau1MPP);
    kMPP2CMP = 2/tau2MPP;
    # CMP dynamics
    tau1CMP = (log2(aCMP)-1)/log2(aCMP) * tauCMP;
    tau2CMP = tauCMP - tau1CMP;
    krep1CMP = (aCMP/2-1)/ tau1CMP;
    krep2CMP = 1/ tau2CMP;
    k12CMP = aCMP/(2*tau1CMP);
    kCMP2BFUMK = 2/tau2CMP;
    # BFU-MK dynamics
    tau1BFUMK = (log2(aBFUMK)-1)/log2(aBFUMK) * tauBFUMK;
    tau2BFUMK = tauBFUMK - tau1BFUMK;
    krep1BFUMK = (aBFUMK/2-1)/ tau1BFUMK;
    krep2BFUMK = 1/ tau2BFUMK;
    k12BFUMK = aBFUMK/(2*tau1BFUMK);
    kBFUMK2CFUMK = 2/tau2BFUMK;
    # CFU-MK dynamics 
    aCFUMK = aCFUMK0 * log(2.73 + beta2 * TPO) / log(1.38); 
    if aCFUMK > 2.1 
        tau1CFUMK = (log2(aCFUMK)-1)/log2(aCFUMK) * tauCFUMK;
        tau2CFUMK = tauCFUMK - tau1CFUMK;
        krep1CFUMK = (aCFUMK/2-1)/ tau1CFUMK;
        krep2CFUMK = 1/ tau2CFUMK;
        k12CFUMK = aCFUMK/(2*tau1CFUMK);
        kCFUMK2MK = 2/tau2CFUMK;
    else
        krep1CFUMK = 0. 
        krep2CFUMK = 0. 
        k12CFUMK = 0. 
        kCFUMK2MK = 0.
    end
    # MK and platelet dynamics 
    kdeathMK = 1/tauMK
    kdeathPlatelet = 1/tauPlatelet
    # TPO dynamics 
    ksynTPO = exp(alpha1 - beta1 * platelet0/V_blood * 1E-9)
    kdegTPO = log(2)/thalfTPO
    # ODEs
    du.LT0 = rLT * LT0 - deltaLT * LT0 * LT0 - kLT2ST * LT0;
    du.ST01   = kLT2ST * LT0 + (krep1ST - k12ST) * ST01;
    du.ST02   = k12ST * ST01 + (krep2ST - kST2MPP) * ST02; 
    du.MPP01  = kST2MPP * ST02 + (krep1MPP - k12MPP) * MPP01; 
    du.MPP02  = k12MPP * MPP01 + (krep2MPP - kMPP2CMP) * MPP02; 
    du.CMP01  = kMPP2CMP * MPP02 + (krep1CMP - k12CMP) * CMP01; 
    du.CMP02  = k12CMP * CMP01 + (krep2CMP - kCMP2BFUMK) * CMP02;
    du.BFUMK01 = kCMP2BFUMK * CMP02 + (krep1BFUMK - k12BFUMK) * BFUMK01;  
    du.BFUMK02 = k12BFUMK * BFUMK01 + (krep2BFUMK - kBFUMK2CFUMK) * BFUMK02; 
    du.CFUMK01 = kBFUMK2CFUMK * BFUMK02 + (krep1CFUMK - k12CFUMK) * CFUMK01; 
    du.CFUMK02 = k12CFUMK * CFUMK01 + (krep2CFUMK - kCFUMK2MK) * CFUMK02;
    du.MK0   = kCFUMK2MK * CFUMK02 - kdeathMK * MK0;  
    du.platelet0   = kMKplatelet * MK0 - kdeathPlatelet * platelet0; 
    du.TPO = ksynTPO - kdegTPO * TPO;
end
