# date: May 10, 2024
# author: Yuezhe Li 
# purpose of this code: combine T-DM1 PK model with a simplified platelet development model

const V_plasma = 3.126 # [L]
const V_blood = 5 # [L]

using DifferentialEquations, ComponentArrays
using Parameters: @unpack

function mk_tdm1_simp!(du, u, p, t)
    @unpack CFUMK01, CFUMK02, MK0, platelet0, TPO, 
            cent_tdm1, peri_tdm1, c_sher2, c_sher2_adc, p_sher2, p_sher2_adc = u; 
    @unpack tauCFUMK, tauMK, tauPlatelet, thalfTPO, aCFUMK0, kMKplatelet, alpha1, beta1, beta2, IC50_tdm1, n_tdm1, 
            Vcent, Vperi, Q, Q_sher2, thalf_adc, kdec, kon_her2, kd_her2, thalf_sher2, thalf_sher2_adc, base_sher2, KpBM, infusion = p
    # define BFU-MK -> CFU-MK steady state influx (obtained from simulation on the full model)
    INFLUX = 1.55E8
    # CFU-MK dynamics 
    aCFUMK = aCFUMK0 * log(2.73 + beta2 * TPO)/log(1.38) * ( IC50_tdm1^n_tdm1./(IC50_tdm1^n_tdm1 .+ max(cent_tdm1*KpBM, 0.).^n_tdm1) )
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
    # dynamics of T-DM1 
    kelim = log(2)/ thalf_adc
    # dynamics of sHER2 
    koff_her2 = kon_her2 * kd_her2
    kdeg_sher2 = log(2)/thalf_sher2
    kdeg_sher2_adc = log(2)/thalf_sher2_adc
    ksyn_sher2 = kdeg_sher2 * base_sher2
    # MK dynamics
    du.CFUMK01 = INFLUX + (krep1CFUMK - k12CFUMK) * CFUMK01; 
    du.CFUMK02 = k12CFUMK * CFUMK01 + (krep2CFUMK - kCFUMK2MK) * CFUMK02;
    du.MK0   = kCFUMK2MK * CFUMK02 - kdeathMK * MK0;  
    du.platelet0   = kMKplatelet * MK0 - kdeathPlatelet * platelet0; 
    du.TPO = ksynTPO - kdegTPO * TPO;
    # T-DM1 PK
    du.cent_tdm1 = -Q/Vcent*(cent_tdm1 - peri_tdm1) - kelim*cent_tdm1 - kdec*cent_tdm1 - kon_her2*cent_tdm1*c_sher2 + koff_her2*c_sher2_adc + infusion
    du.peri_tdm1 = Q/Vperi*(cent_tdm1 - peri_tdm1) - kelim*peri_tdm1 - kdec*peri_tdm1 - kon_her2*peri_tdm1*p_sher2 + koff_her2*p_sher2_adc
    du.c_sher2_adc = kon_her2*cent_tdm1*c_sher2 - koff_her2*c_sher2_adc - kdeg_sher2_adc*c_sher2_adc -Q/Vcent*(c_sher2_adc - p_sher2_adc)
    du.c_sher2 = ksyn_sher2 - kdeg_sher2*c_sher2 - kon_her2*cent_tdm1*c_sher2 + koff_her2*c_sher2_adc -Q_sher2/Vcent*(c_sher2 - p_sher2)
    du.p_sher2_adc = kon_her2*peri_tdm1*p_sher2 - koff_her2*p_sher2_adc - kdeg_sher2_adc*p_sher2_adc + Q/Vperi*(c_sher2_adc - p_sher2_adc)
    du.p_sher2 = ksyn_sher2 - kdeg_sher2*p_sher2 - kon_her2*peri_tdm1*p_sher2 + koff_her2*p_sher2_adc + Q_sher2/Vperi*(c_sher2 - p_sher2)
end

