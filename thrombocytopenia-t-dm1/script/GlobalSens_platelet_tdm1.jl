# date: May 15, 2024
# author: Yuezhe Li 
# purpose of this code: global sensitivity analysis of simplified combined model

const MW_ADC = 15E4; # g/mol
const BW = 70; #kg

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames, DataFramesMeta
using Parameters: @unpack
using GlobalSensitivity, QuasiMonteCarlo, CairoMakie, Statistics

# read in simplified model 
include("model/platelet_dynamics_simp_tdm1.jl");

p_default =  ComponentArray(tauCFUMK = 12.0, tauMK = 0.1, tauPlatelet = 7.0, thalfTPO = 1.25, aCFUMK0 = 8, 
    kMKplatelet = 480., alpha1 = -2.86, beta1 = 0.009, beta2 = 4.72E-3, IC50_tdm1 = 50., n_tdm1 = 1.13, 
    Vcent = 3., Vperi = 13., Q = 0.196, Q_sher2 = 0.46, thalf_adc = 11.6, kdec = 0.073, kon_her2 = 8.64, kd_her2 = 0.3, 
    thalf_sher2 = 0.21, thalf_sher2_adc = 8.16, base_sher2 = 1.16, KpBM = 0.07, infusion = 0.); 

dose_mgkg = 3.6

u0_default = ComponentArray(CFUMK01=1.45E9, CFUMK02=5E9, MK0=3.86E8, platelet0=1.3E12, TPO=0.0099, 
    cent_tdm1 = dose_mgkg*1E6*BW/MW_ADC/V_plasma, peri_tdm1 = 0., c_sher2 = p_default.base_sher2, c_sher2_adc = 0., p_sher2 = p_default.base_sher2, p_sher2_adc = 0.);

tspan = (0.0, 21.0);  
prob = ODEProblem(mk_tdm1_simp!, u0_default, tspan, p_default);

f1 = function (p)
    p_tmp = deepcopy(p_default); 
    p_tmp.tauPlatelet = p[1];
    p_tmp.kMKplatelet = p[2];
    p_tmp.thalfTPO = p[3];
    p_tmp.IC50_tdm1 = p[4];
    sol_pre = solve(remake(prob; p = p_tmp, tspan = (0., 100)), Tsit5(); saveat = 1.); 
    u02 = sol_pre.u[end];
    prob1 = remake(prob; p = p_tmp, u0 = u02); 
    sol = solve(prob1, Tsit5(); saveat = 1.)
    return [minimum(sol[4, :]), mean(sol[4, :]), maximum(sol[4, :])]
end

bounds = [[3, 12], [240, 1E5], [0.83, 1.66], [2, 50]]; 

sobol_sens = gsa(f1, Sobol(), bounds, samples = 500)

fig = Figure(resolution = (600, 400));
barplot(fig[1,1], [1,2,3,4], sobol_sens.S1[1, :], axis = (xlabel = "", xticks = ([1,2,3,4], ["PLT residence time", "PLT/MK/day", "t1/2, TPO", "IC50, T-DM1"]), ylabel = "First order"));
barplot(fig[2,1], [1,2,3,4], sobol_sens.ST[1, :], axis = (xlabel = "", xticks = ([1,2,3,4], ["PLT residence time", "PLT/MK/day", "t1/2, TPO", "IC50, T-DM1"]), ylabel = "Total order"));
fig

save("deliv/figure/global_sens.png", fig);
