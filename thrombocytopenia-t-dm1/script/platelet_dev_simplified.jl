# author: Yuezhe Li 
# date: May 6, 2024
# purpose of this script: to code for a simplified platelet development model 
# BFU-MK -> CFU-MK was replaced by a steady state flow 

const V_plasma = 3.126 # [L]
const V_blood = 5 # [L]

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames
using Parameters: @unpack
using Plots

function mk_dev_simp!(du, u, p, t)
    @unpack CFUMK01, CFUMK02, MK0, platelet0, TPO = u; 
    @unpack tauCFUMK, tauMK, tauPlatelet, thalfTPO, aCFUMK0, kMKplatelet, alpha1, beta1, beta2 = p
    # define BFU-MK -> CFU-MK steady state influx (obtained from simulation on the full model)
    INFLUX = 1.55E8
    # CFU-MK dynamics 
    aCFUMK = aCFUMK0 * log(2.73 + beta2 * TPO) / log(1.38); 
    tau1CFUMK = (log2(aCFUMK)-1)/log2(aCFUMK) * tauCFUMK;
    tau2CFUMK = tauCFUMK - tau1CFUMK;
    krep1CFUMK = (aCFUMK/2-1)/ tau1CFUMK;
    krep2CFUMK = 1/ tau2CFUMK;
    k12CFUMK = aCFUMK/(2*tau1CFUMK);
    kCFUMK2MK = 2/tau2CFUMK;
    # MK and platelet dynamics 
    kdeathMK = 1/tauMK
    kdeathPlatelet = 1/tauPlatelet
    # TPO dynamics 
    ksynTPO = exp(alpha1 - beta1 * platelet0/V_blood * 1E-9)
    kdegTPO = log(2)/thalfTPO
    # ODEs
    du.CFUMK01 = INFLUX + (krep1CFUMK - k12CFUMK) * CFUMK01; 
    du.CFUMK02 = k12CFUMK * CFUMK01 + (krep2CFUMK - kCFUMK2MK) * CFUMK02;
    du.MK0   = kCFUMK2MK * CFUMK02 - kdeathMK * MK0;  
    du.platelet0   = kMKplatelet * MK0 - kdeathPlatelet * platelet0; 
    du.TPO = ksynTPO - kdegTPO * TPO;
end

p_mk_simp =  ComponentArray(tauCFUMK = 12.0, tauMK = 0.1, tauPlatelet = 7.0, thalfTPO = 1.25, aCFUMK0 = 8, 
    kMKplatelet = 480., alpha1 = -2.86, beta1 = 0.009, beta2 = 4.72E-3
); 

init_mk_simp = ComponentArray(CFUMK01=1.45E9, CFUMK02=5E9, MK0=3.86E8, platelet0=1.3E12, TPO=0.0099);

sol_simp = solve(ODEProblem(mk_dev_simp!, init_mk_simp, (0.0, 100.0), p_mk_simp), alg = AutoTsit5(Rosenbrock23()), saveat=1.0);
sdf_simp = DataFrame(sol_simp); rename!(sdf_simp, Symbol.(names(sdf_simp)[2:end]) .=> collect(keys(sol_simp.u[end])));

p_steadystate = plot(xlabel = "Time (Day)", ylabel = "Cell count (#)", legend = :outerright, ylims = [1E8, 1E13],yaxis = :log10);
plot!(sdf_simp.timestamp, sdf_simp.CFUMK01 .+ sdf_simp.CFUMK02, label = "CFU-MK");
plot!(sdf_simp.timestamp, sdf_simp.MK0, label = "MK");
plot!(sdf_simp.timestamp, sdf_simp.platelet0, label = "Platelet");
display(p_steadystate);
