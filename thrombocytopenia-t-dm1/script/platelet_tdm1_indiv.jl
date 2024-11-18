# date: May 10, 2024
# author: Yuezhe Li 
# purpose: to test combined T-DM1 and platelet dynamics (simplified)

const MW_sHER2 = 1E5; # g/mol
const MW_ADC = 15E4; # g/mol
const BW = 70; #kg

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames, DataFramesMeta
using Parameters: @unpack
using Plots

# observed data; obtained from Krop et al., 2010, Fig S1;  # https://ascopubs.org/doi/10.1200/JCO.2009.26.2071
# a patient with low platelet count (navy)
obs1 = DataFrame(
    time_day = [2.7,6.6,15.3,21.9,28.5,42.2,48.8,56.9,62.9,71.1,77.7,84.7,93.9,104.5,115.2,118.7,127.3,133.9,141.6,148.6,154.7,161.3,168.9,175.5,197.3,206.4,225.2,267.3], 
    platelet_k_per_uL = [139.4,57.7,109.4,56.6,131.1,73.4,123.9,113.2,63.9,108.4,103.6,62.8,100,61.6,89.3,94.1,71.3,96.6,105,55.7,78.6,83.4,72.6,81.1,87.1,69.1,100.4,97]
);
# a patient with high platelet count (purple)
obs2 = DataFrame(
    time_day = [6.7,18.1,30.1,35.8,41.4,50.9,57.1,63.2,71.7,78.4,86,93,107.8],
    platelet_k_per_uL = [143.1,349.8,145.5,343.9,363.1,154,330.7,343.9,190.1,300.7,304.3,191.4,304.4]
);
# a patient with medium platelet count (deepskyblue)
obs3 = DataFrame(
    time_day = [5.2,22,44.4,51.4,66.7,74.2,86.9,93.9,101.1,115.3,123.4,129.5,136.1,143.6,157.9,164.9,179.2,186.2,199.9,206.5,214.2,228.3,242,249.6,261.3,271.4,284.7,298.3,306.4,312.9,326.2,333.2], 
    platelet_k_per_uL = [116.6,244.1,274.2,162.4,255,186.5,225,137.3,184.2,193.9,159,204.7,198.7,144.7,214.4,151.9,194,146,158.1,141.2,208.6,126.9,190.6,179.8,195.5,172.7,280.9,165.6,176.4,135.6,187.3,139.2]
);

# read in simplified model 
include("model/platelet_dynamics_simp_tdm1.jl");

p_default =  ComponentArray(tauCFUMK = 12.0, tauMK = 0.1, tauPlatelet = 7.0, thalfTPO = 1.25, aCFUMK0 = 8, 
    kMKplatelet = 480., alpha1 = -2.86, beta1 = 0.009, beta2 = 4.72E-3, IC50_tdm1 = 50., n_tdm1 = 1.13, 
    Vcent = 3., Vperi = 13., Q = 0.196, Q_sher2 = 0.46, thalf_adc = 11.6, kdec = 0.073, kon_her2 = 8.64, kd_her2 = 0.3, 
    thalf_sher2 = 0.21, thalf_sher2_adc = 8.16, base_sher2 = 1.16, KpBM = 0.07, infusion = 0.); 

u0_default = ComponentArray(CFUMK01=1.45E9, CFUMK02=5E9, MK0=3.86E8, platelet0=1.3E12, TPO=0.0099, 
    cent_tdm1 = 0., peri_tdm1 = 0., c_sher2 = 8E3/MW_sHER2, c_sher2_adc = 0., p_sher2 = 8E3/MW_sHER2, p_sher2_adc = 0.);

function RepeatedDosing(dose_mgkg, infusion_time, AddDose = [0., 21.])
    global cbs0 = [];
    if length(AddDose) >= 1
        for i in 1:length(AddDose)
            function affect_infusion_on!(integrator)
                integrator.p.infusion = dose_mgkg*1E6*BW/MW_ADC/V_plasma/infusion_time; # [nM]
            end
            function affect_infusion_off!(integrator)
                integrator.p.infusion = 0.;
            end
            cb01 = PresetTimeCallback(AddDose[i],affect_infusion_on!);
            cb02 = PresetTimeCallback( (AddDose[i]+infusion_time),affect_infusion_off!);
            global cbs0 = push!(cbs0, cb01);
            global cbs0 = push!(cbs0, cb02);
        end
    end
    cbset0 = CallbackSet(cbs0...);
    return cbset0
end

function ModelCalibration(p_mk_simp = p_default, simsdays = 105., u0_default = u0_default, Dose_mgkg = 3.6, infusion_time = 0.125, )
    # simulation to a steady state before dosing 
    sol_pre = solve(ODEProblem(mk_tdm1_simp!, u0_default, (0., simsdays), p_mk_simp), alg = AutoTsit5(Rosenbrock23()), saveat=1.0);
    u0_ss = sol_pre.u[end];
    # set up repeated dosing as suggested in Krop et al., 2010; # https://ascopubs.org/doi/10.1200/JCO.2009.26.2071
    repeated_doses = floor((simsdays-1)/21);
    AddDose = collect(0:repeated_doses) * 21. 
    cbset0 = RepeatedDosing(Dose_mgkg, infusion_time, AddDose);
    # simulation after dosing 
    sol_simp = solve(ODEProblem(mk_tdm1_simp!, u0_ss, (-0.001, simsdays), p_mk_simp), alg = AutoTsit5(Rosenbrock23()), saveat=1.0, callback = cbset0);
    sdf_simp = DataFrame(sol_simp); rename!(sdf_simp, Symbol.(names(sdf_simp)[2:end]) .=> collect(keys(sol_simp.u[end])));
    @rsubset!(sdf_simp, :timestamp >= 0.);
    return sdf_simp;
end

# individual tunning: platelet residence time, platelet/mk/day, 
p_low = deepcopy(p_default); p_low.tauPlatelet = 3.; p_low.kMKplatelet = 500.; p_low.IC50_tdm1 = 10.
sdf1 = ModelCalibration(p_low, 200.);
p1 = plot(xlabel = "Time (Day)", ylabel = "Platelet count (1000 #/uL)", ylim = [0, 400], xlim = [0., 105.]);
plot!(sdf1.timestamp, sdf1.platelet0/1E3/(V_blood*1E6), label = "pred");
scatter!(obs1.time_day, obs1.platelet_k_per_uL, label = "obs");
plot!(xticks = ([0, 21, 42, 63, 84, 105], ["0", "21", "42", "63", "84", "105"]));

p_high = deepcopy(p_default); p_high.tauPlatelet = 7; p_high.kMKplatelet = 1200.; p_high.IC50_tdm1 = 2.
sdf2 = ModelCalibration(p_high, 200.);
p2 = plot(xlabel = "Time (Day)", ylabel = "Platelet count (1000 #/uL)", ylim = [0, 600], xlim = [0., 105.]);
plot!(sdf2.timestamp, sdf2.platelet0/1E3/(V_blood*1E6), label = "pred");
scatter!(obs2.time_day, obs2.platelet_k_per_uL, label = "obs");
plot!(xticks = ([0, 21, 42, 63, 84, 105], ["0", "21", "42", "63", "84", "105"]));

p_med = deepcopy(p_default); p_med.tauPlatelet = 5.; p_med.kMKplatelet = 620.; p_med.IC50_tdm1 = 30.
sdf3 = ModelCalibration(p_med, 334.);
p3 = plot(xlabel = "Time (Day)", ylabel = "Platelet count (1000 #/uL)", ylim = [0, 600], xlim = [0., 105.]);
plot!(sdf3.timestamp, sdf3.platelet0/1E3/(V_blood*1E6), label = "pred");
scatter!(obs3.time_day, obs3.platelet_k_per_uL, label = "obs");
plot!(xticks = ([0, 21, 42, 63, 84, 105], ["0", "21", "42", "63", "84", "105"]));

savefig(plot(p2, p3, ncol = 2), "deliv/figure/indiv_comparison.png");
