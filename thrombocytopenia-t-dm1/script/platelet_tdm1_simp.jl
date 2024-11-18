# date: May 10, 2024
# author: Yuezhe Li 
# purpose: to test combined T-DM1 and platelet dynamics (simplified)

const MW_sHER2 = 1E5; # g/mol
const MW_ADC = 15E4; # g/mol
const BW = 70; #kg

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames
using Parameters: @unpack
using Plots

include("model/platelet_dynamics_simp_tdm1.jl");

p_default =  ComponentArray(tauCFUMK = 12.0, tauMK = 0.1, tauPlatelet = 7.0, thalfTPO = 1.25, aCFUMK0 = 8, 
    kMKplatelet = 480., alpha1 = -2.86, beta1 = 0.009, beta2 = 4.72E-3, IC50_tdm1 = 20., n_tdm1 = 1.13, 
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

cbset0 = RepeatedDosing(3.6, 0.125, [0., 21., 42., 63., 84.]);

sol_simp = solve(ODEProblem(mk_tdm1_simp!, u0_default, (-0.001, 105.), p_default), alg = AutoTsit5(Rosenbrock23()), saveat=1.0, callback = cbset0);
sdf_simp = DataFrame(sol_simp); rename!(sdf_simp, Symbol.(names(sdf_simp)[2:end]) .=> collect(keys(sol_simp.u[end])));

p_plt = plot(xlabel = "Time (Day)", ylabel = "Platelet count (1000 #/uL)", ylim = [150, 300], xlim = [0., 105.]);
plot!(sdf_simp.timestamp, sdf_simp.platelet0/1E3/(V_blood*1E6), label = false);
plot!(xticks = ([0, 21, 42, 63, 84, 105], ["0", "21", "42", "63", "84", "105"]));

p_pk = plot(xlabel = "Time (Day)", ylabel = "Plasma T-DM1 (nM)", ylim = [0, 600], xlim = [0., 105.]);
plot!(sdf_simp.timestamp, sdf_simp.cent_tdm1, label = false);
plot!(xticks = ([0, 21, 42, 63, 84, 105], ["0", "21", "42", "63", "84", "105"]));

display(plot(p_plt, p_pk, ncol = 2));

savefig(plot(p_plt, p_pk, ncol = 2), "deliv/figure/simplifed_combined.png");
