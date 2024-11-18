# date: May 15, 2024
# author: Yuezhe Li 
# purpose of this code: to test for alternative dose in Girish et al., 2012

const MW_ADC = 15E4; # g/mol
const BW = 70; #kg

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames, DataFramesMeta
using Parameters: @unpack
using Plots

# read in simplified model 
include("model/platelet_dynamics_simp_tdm1.jl");

p_default =  ComponentArray(tauCFUMK = 12.0, tauMK = 0.1, tauPlatelet = 7.0, thalfTPO = 1.25, aCFUMK0 = 8, 
    kMKplatelet = 480., alpha1 = -2.86, beta1 = 0.009, beta2 = 4.72E-3, IC50_tdm1 = 20., n_tdm1 = 1.13, 
    Vcent = 3., Vperi = 13., Q = 0.196, Q_sher2 = 0.46, thalf_adc = 11.6, kdec = 0.073, kon_her2 = 8.64, kd_her2 = 0.3, 
    thalf_sher2 = 0.21, thalf_sher2_adc = 8.16, base_sher2 = 1.16, KpBM = 0.07, infusion = 0.); 

u0_default = ComponentArray(CFUMK01=1.45E9, CFUMK02=5E9, MK0=3.86E8, platelet0=1.3E12, TPO=0.0099, 
    cent_tdm1 = 0., peri_tdm1 = 0., c_sher2 = p_default.base_sher2, c_sher2_adc = 0., p_sher2 = p_default.base_sher2, p_sher2_adc = 0.);

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

function ModelSims(p_mk_simp = p_default, Dose_mgkg = 3.6, simsdays = 105., u0_default = u0_default, infusion_time = 0.125, )
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

# test parameters that makes platelet count in the range of literature range 
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4914488/
p1 = deepcopy(p_default); p1.tauPlatelet = 4.; sdf1 = ModelSims(p1, 0.); sdf1.platelet0[end]/1E3/(V_blood*1E6)  # 148.48
p2 = deepcopy(p_default); p2.tauPlatelet = 11.; sdf2 = ModelSims(p2, 0.); sdf2.platelet0[end]/1E3/(V_blood*1E6)  # 408.3

# simulation on multiple doses
sdf1_4point8 = ModelSims(p1, 4.8);
sdf1_3point6 = ModelSims(p1, 3.6);
sdf1_2point4 = ModelSims(p1, 2.4);
sdf1_1point2 = ModelSims(p1, 1.2);

sdf2_4point8 = ModelSims(p2, 4.8);
sdf2_3point6 = ModelSims(p2, 3.6);
sdf2_2point4 = ModelSims(p2, 2.4);
sdf2_1point2 = ModelSims(p2, 1.2);

plt1 = plot(xlabel = "Time (Day)", ylabel = "Platelet count (1000 #/uL)", ylim = [0, 410], xlim = [0., 105.], legend = :outerright);
plot!(sdf1_4point8.timestamp, sdf1_4point8.platelet0/1E3/(V_blood*1E6), label = "low base PLT, 4.8mg/kg");
plot!(sdf1_3point6.timestamp, sdf1_3point6.platelet0/1E3/(V_blood*1E6), label = "low base PLT, 3.6mg/kg");
plot!(sdf1_2point4.timestamp, sdf1_2point4.platelet0/1E3/(V_blood*1E6), label = "low base PLT, 2.4mg/kg");
plot!(sdf1_1point2.timestamp, sdf1_1point2.platelet0/1E3/(V_blood*1E6), label = "low base PLT, 1.2mg/kg");
plot!(sdf2_4point8.timestamp, sdf2_4point8.platelet0/1E3/(V_blood*1E6), label = "high base PLT, 4.8mg/kg", linestyle = :dashdot);
plot!(sdf2_3point6.timestamp, sdf2_3point6.platelet0/1E3/(V_blood*1E6), label = "high base PLT, 3.6mg/kg", linestyle = :dashdot);
plot!(sdf2_2point4.timestamp, sdf2_2point4.platelet0/1E3/(V_blood*1E6), label = "high base PLT, 2.4mg/kg", linestyle = :dashdot);
plot!(sdf2_1point2.timestamp, sdf2_1point2.platelet0/1E3/(V_blood*1E6), label = "high base PLT, 1.2mg/kg", linestyle = :dashdot);
hline!([25., 45., 75., 90], linestyle = :dash, color = :black, alpha = 0.3, label = false);
annotate!(105, 20, ("Grade IV TCP", 8, :left, :black));
annotate!(105, 40, ("Grade III TCP", 8, :left, :black));
annotate!(105, 70, ("Grade II TCP", 8, :left, :black));
annotate!(105, 85, ("Grade I TCP", 8, :left, :black));
plot!(xticks = ([0, 21, 42, 63, 84, 105], ["0", "21", "42", "63", "84", "105"]));
display(plt1);

savefig(plt1, "deliv/figure/alt_dose_tcp.png");
