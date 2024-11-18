# date: 5/9/24
# author: Yuezhe Li 
# purpose of this code: to verify T-DM1 PK 
# parameters were obtained from Scheuher et al., 2023; # https://link.springer.com/article/10.1007/s10928-023-09884-6
# tumor compartment was dropped to simplified the model
# in addition, HER2-expression cells were also removed to simplify the model 

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames, DataFramesMeta
using Parameters: @unpack
using Plots
using Optimization, OptimizationBBO

const MW_sHER2 = 1E5; # g/mol
const MW_ADC = 15E4; # g/mol
const h2s = 3600;
const d2h = 24.;
const d2s = d2h*h2s;
const BW = 70; #kg

# observed data obtained from Girish et al., 2012; # https://link.springer.com/article/10.1007/s00280-011-1817-3
obs = DataFrame(
    time_day = [0.1,0.25,1,2,3,0.1,0.25,1,2,3,7,0.1,0.25,1,2,3,7,0.25,1,2,3,7,10,14,17,21,0.1,0.25,1,2,3,7,10,17,21,0.1,0.25,1,2,3,7,14,17,21], 
    T_DM1_ugperml = [9.66134,7.1674,4.22606,2.54967,1.11527,13.02307,10.35053,7.33393,5.69654,1.37138,0.52264,20.14835,20.14835,11.08888,7.33393,3.208,0.49918,74.61359,55.35313,45.01597,33.39572,17.15592,10.35053,3.68201,2.02644,0.40596,74.61359,66.51853,48.22717,34.96559,27.15909,13.95206,9.44196,2.27305,0.77229,126.54482,118.11885,89.66407,78.12105,65.00812,29.09648,12.43836,7.1674,3.3588],
    Dose_mgkg = [0.3,0.3,0.3,0.3,0.3,0.6,0.6,0.6,0.6,0.6,0.6,1.2,1.2,1.2,1.2,1.2,1.2,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,2.4,3.6,3.6,3.6,3.6,3.6,3.6,3.6,3.6,3.6,4.8,4.8,4.8,4.8,4.8,4.8,4.8,4.8,4.8]
);

# PK model 
function tdm1_pk!(du, u, p, t)
    @unpack cent_tdm1, peri_tdm1, c_sher2, c_sher2_adc, p_sher2, p_sher2_adc  = u;
    @unpack Vcent, Vperi, Pdist_adc, t_dist_adc, P_dist_sher2, t_dist_sher2, thalf_adc, kdec, kon_her2, kd_her2, thalf_sher2, thalf_sher2_adc, base_sher2 = p
    # dynamics of T-DM1 
    k12_adc = log(2)/t_dist_adc * Pdist_adc / (Pdist_adc + Vcent/Vperi)
    k21_adc = log(2)/t_dist_adc * Vcent/Vperi / (Pdist_adc + Vcent/Vperi)
    kelim = log(2)/ thalf_adc
    # dynamics of sHER2 
    k12_sher2 = log(2)/t_dist_sher2 * P_dist_sher2 / (P_dist_sher2 + Vcent/Vperi)
    k21_sher2 = log(2)/t_dist_sher2 * Vcent/Vperi / (P_dist_sher2 + Vcent/Vperi)
    koff_her2 = kon_her2 * kd_her2
    kdeg_sher2 = log(2)/thalf_sher2
    kdeg_sher2_adc = log(2)/thalf_sher2_adc
    ksyn_sher2 = kdeg_sher2 * base_sher2
    # ODEs
    du.cent_tdm1 = -k12_adc*cent_tdm1 + k21_adc*peri_tdm1 - kelim*cent_tdm1 - kdec*cent_tdm1 - kon_her2*cent_tdm1*c_sher2 + koff_her2*c_sher2_adc
    du.peri_tdm1 = k12_adc*cent_tdm1 - k21_adc*peri_tdm1 - kelim*peri_tdm1 - kdec*peri_tdm1 - kon_her2*peri_tdm1*p_sher2 + koff_her2*p_sher2_adc 
    du.c_sher2_adc = kon_her2*cent_tdm1*c_sher2 - koff_her2*c_sher2_adc - kdeg_sher2_adc*c_sher2_adc - k12_adc*c_sher2_adc + k21_adc*p_sher2_adc
    du.c_sher2 = ksyn_sher2 - kdeg_sher2*c_sher2 - kon_her2*cent_tdm1*c_sher2 + koff_her2*c_sher2_adc - k12_sher2*c_sher2 + k21_sher2*p_sher2
    du.p_sher2_adc = kon_her2*peri_tdm1*p_sher2 - koff_her2*p_sher2_adc - kdeg_sher2_adc*p_sher2_adc + k12_adc*c_sher2_adc - k21_adc*p_sher2_adc
    du.p_sher2 = ksyn_sher2*P_dist_sher2 - kdeg_sher2*p_sher2 - kon_her2*peri_tdm1*p_sher2 + koff_her2*p_sher2_adc + k12_sher2*c_sher2 - k21_sher2*p_sher2
end

p_tmd1 =  ComponentArray(Vcent = 3., Vperi = 13., Pdist_adc = 0.187, t_dist_adc = 10/d2h, thalf_adc = 11.6, P_dist_sher2 = 1., t_dist_sher2 = 8/d2h, 
                         kdec = 8.5E-7*d2s, kon_her2 = 1E-4*d2s, kd_her2 = 0.3, thalf_sher2 = 5/d2h, thalf_sher2_adc = 11.6, 
                         base_sher2 = 8E3/MW_sHER2); # units: nM, day
        
u0 = ComponentArray(cent_tdm1 = 0., peri_tdm1 = 0., c_sher2 = 8E3/MW_sHER2, c_sher2_adc = 0., p_sher2 = 8E3/MW_sHER2, p_sher2_adc = 0.); 

function sims_dose(dose_mgkg, sims_days = 21)
    u0_tmp = deepcopy(u0); 
    u0_tmp.cent_tdm1 = BW * dose_mgkg/MW_ADC*1E6/p_tmd1.Vcent; # [nM]
    sol = solve(ODEProblem(tdm1_pk!, u0_tmp, (0.0, sims_days), p_tmd1), alg = AutoTsit5(Rosenbrock23()), saveat=1.);
    sdf = DataFrame(sol); rename!(sdf, Symbol.(names(sdf)[2:end]) .=> collect(keys(sol.u[end])));
    sdf.time_day .= sdf.timestamp
    sdf.ADC_ugml .= sdf.cent_tdm1 * MW_ADC * 1E-6
    sdf.Dose_mgkg .= dose_mgkg
    return sdf
end

sdf = vcat(sims_dose(0.3), sims_dose(0.6), sims_dose(1.2), sims_dose(2.4), sims_dose(3.6), sims_dose(4.8));

p_pk = plot(xlabel = "Time (Day)", ylabel = "ADC conc (nM)", legend = :outerright, legendtitle = "Dose Group (mg/kg)", palette = :Set2_6); 
plot!(sdf.time_day, sdf.cent_tdm1, group = sdf.Dose_mgkg, linewidth = 2); 
scatter!(obs.time_day, obs.T_DM1_ugperml*1E6/MW_ADC, group = obs.Dose_mgkg, markerstrokewidth=0, markersize=6, alpha = 0.8); 
plot!(yaxis = :log, xticks = ([0, 7, 14, 21], ["0", "7", "14", "21"]));
display(p_pk);

# simply PK 
function tdm1_simp_pk!(du, u, p, t)
    @unpack cent_tdm1, peri_tdm1, c_sher2, c_sher2_adc, p_sher2, p_sher2_adc  = u;
    @unpack Vcent, Vperi, Q, Q_sher2, thalf_adc, kdec, kon_her2, kd_her2, thalf_sher2, thalf_sher2_adc, base_sher2 = p
    # dynamics of T-DM1 
    kelim = log(2)/ thalf_adc
    # dynamics of sHER2 
    koff_her2 = kon_her2 * kd_her2
    kdeg_sher2 = log(2)/thalf_sher2
    kdeg_sher2_adc = log(2)/thalf_sher2_adc
    ksyn_sher2 = kdeg_sher2 * base_sher2
    # ODEs
    du.cent_tdm1 = -Q/Vcent*(cent_tdm1 - peri_tdm1) - kelim*cent_tdm1 - kdec*cent_tdm1 - kon_her2*cent_tdm1*c_sher2 + koff_her2*c_sher2_adc
    du.peri_tdm1 = Q/Vperi*(cent_tdm1 - peri_tdm1) - kelim*peri_tdm1 - kdec*peri_tdm1 - kon_her2*peri_tdm1*p_sher2 + koff_her2*p_sher2_adc
    du.c_sher2_adc = kon_her2*cent_tdm1*c_sher2 - koff_her2*c_sher2_adc - kdeg_sher2_adc*c_sher2_adc -Q/Vcent*(c_sher2_adc - p_sher2_adc)
    du.c_sher2 = ksyn_sher2 - kdeg_sher2*c_sher2 - kon_her2*cent_tdm1*c_sher2 + koff_her2*c_sher2_adc -Q_sher2/Vcent*(c_sher2 - p_sher2)
    du.p_sher2_adc = kon_her2*peri_tdm1*p_sher2 - koff_her2*p_sher2_adc - kdeg_sher2_adc*p_sher2_adc + Q/Vperi*(c_sher2_adc - p_sher2_adc)
    du.p_sher2 = ksyn_sher2 - kdeg_sher2*p_sher2 - kon_her2*peri_tdm1*p_sher2 + koff_her2*p_sher2_adc + Q_sher2/Vperi*(c_sher2 - p_sher2)
end

function sims_simp_dose(dose_mgkg, obs_time = 1., sims_days = 21, p_s_tmd1 = p_s_tmd1)
    u0_tmp = deepcopy(u0); 
    u0_tmp.c_sher2 = p_s_tmd1.base_sher2; 
    u0_tmp.p_sher2 = p_s_tmd1.base_sher2; 
    u0_tmp.cent_tdm1 = BW * dose_mgkg/MW_ADC*1E6/p_s_tmd1.Vcent; # [nM]
    sol = solve(ODEProblem(tdm1_simp_pk!, u0_tmp, (0.0, sims_days), p_s_tmd1), alg = AutoTsit5(Rosenbrock23()), saveat=obs_time);
    sdf = DataFrame(sol); rename!(sdf, Symbol.(names(sdf)[2:end]) .=> collect(keys(sol.u[end])));
    sdf.time_day .= sdf.timestamp
    sdf.ADC_ugml .= sdf.cent_tdm1 * MW_ADC * 1E-6
    sdf.Dose_mgkg .= dose_mgkg
    return sdf
end

p_s_tmd1 =  ComponentArray(Vcent = 3., Vperi = 13., Q = 0.196, Q_sher2 = 1.4, thalf_adc = 11.6, 
                         kdec = 0.073, kon_her2 = 8.64, kd_her2 = 0.3, thalf_sher2 = 0.21, thalf_sher2_adc = 11.6, 
                         base_sher2 = 4.); # units: nM, day

# simulate for steady state for partition coefficient 
function auc(t, conc)
    auc0 = 0; 
    for i in 2:length(t)
        auc_tmp = (conc[i] + conc[i-1]) * (t[i] - t[i-1])/2
        auc0 += auc_tmp
    end
    return auc0
end

function ss_exp(dose_mgkg, Q, obs_time = 0.1, sims_days = 50., p_s_tmd1 = p_s_tmd1)
    p_tmp = deepcopy(p_s_tmd1);
    p_tmp.Q = Q;
    p_tmp.base_sher2 = 0;
    p_tmp.kdec = 0;
    u0_tmp = deepcopy(u0); 
    u0_tmp.c_sher2 = 0; 
    u0_tmp.p_sher2 = 0; 
    u0_tmp.cent_tdm1 = BW * dose_mgkg/MW_ADC*1E6/p_s_tmd1.Vcent; # [nM]
    sol = solve(ODEProblem(tdm1_simp_pk!, u0_tmp, (0.0, sims_days), p_tmp), alg = AutoTsit5(Rosenbrock23()), saveat=obs_time);
    sdf = DataFrame(sol); rename!(sdf, Symbol.(names(sdf)[2:end]) .=> collect(keys(sol.u[end])));
    cent_exp = auc(sdf.timestamp, sdf.cent_tdm1);
    peri_exp = auc(sdf.timestamp, sdf.peri_tdm1);
    return peri_exp/cent_exp;
end

# ss_exp(2.4, 0.196); # 0.187

function loss1(p, obs, opt = true)
    p_tmp = deepcopy(p_s_tmd1); p_s_tmd1.Q_sher2 = p[1]; p_s_tmd1.thalf_sher2_adc = p[2]; p_s_tmd1.base_sher2 = p[3];
    pk_point3 = @rsubset(obs, :Dose_mgkg .== 0.3); 
    pk_point6 = @rsubset(obs, :Dose_mgkg .== 0.6); 
    pk_1point2 = @rsubset(obs, :Dose_mgkg .== 1.2); 
    pk_2point4 = @rsubset(obs, :Dose_mgkg .== 2.4); 
    pk_3point6 = @rsubset(obs, :Dose_mgkg .== 3.6); 
    pk_4point8 = @rsubset(obs, :Dose_mgkg .== 4.8); 
    sdf_point3 = sims_simp_dose(0.3, pk_point3.time_day, 21., p_tmp);
    sdf_point6 = sims_simp_dose(0.6, pk_point6.time_day, 21., p_tmp);
    sdf_1point2 = sims_simp_dose(1.2, pk_1point2.time_day, 21., p_tmp);
    sdf_2point4 = sims_simp_dose(2.4, pk_2point4.time_day, 21., p_tmp);
    sdf_3point6 = sims_simp_dose(3.6, pk_3point6.time_day, 21., p_tmp);
    sdf_4point8 = sims_simp_dose(4.8, pk_4point8.time_day, 21., p_tmp);
    diff = sum((sdf_point3.ADC_ugml .- pk_point3.T_DM1_ugperml).^2) + sum((sdf_point6.ADC_ugml .- pk_point6.T_DM1_ugperml).^2) + 
           sum((sdf_1point2.ADC_ugml .- pk_1point2.T_DM1_ugperml).^2) + sum((sdf_2point4.ADC_ugml .- pk_2point4.T_DM1_ugperml).^2) +
           sum((sdf_3point6.ADC_ugml .- pk_3point6.T_DM1_ugperml).^2) + sum((sdf_4point8.ADC_ugml .- pk_4point8.T_DM1_ugperml).^2) 
    if opt 
        return diff
    else
        return vcat(sdf_point3, sdf_point6, sdf_1point2, sdf_2point4, sdf_3point6, sdf_4point8)
    end
end

#sol2 = solve(OptimizationProblem(OptimizationFunction(loss1, Optimization.AutoForwardDiff()), [2.1, 11.6, 0.08], obs, ub = [4., 11.6, 5.], lb = [0., 0., 0.]), BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 2E5, maxtime = 1E3, progress = true ); 
# sol2 = [0.46, 8.16, 1.16]
p_s_tmd1.Q_sher2 = 0.46
p_s_tmd1.thalf_sher2_adc = 8.16
p_s_tmd1.base_sher2 = 1.16

sdf2 = vcat(sims_simp_dose(0.3), sims_simp_dose(0.6), sims_simp_dose(1.2), sims_simp_dose(2.4), sims_simp_dose(3.6), sims_simp_dose(4.8));

p2_pk = plot(xlabel = "Time (Day)", ylabel = "ADC conc (nM)", legend = :bottomleft, legendtitle = "Dose Group (mg/kg)", legendtitlefontsize = 8, palette = :Set2_6); 
plot!(sdf2.time_day, sdf2.cent_tdm1, group = sdf2.Dose_mgkg, linewidth = 2); 
scatter!(obs.time_day, obs.T_DM1_ugperml*1E6/MW_ADC, group = obs.Dose_mgkg, markerstrokewidth=0, markersize=6, alpha = 0.8); 
plot!(yaxis = :log, ylim = [1E-4, 2E3],  xticks = ([0, 7, 14, 21], ["0", "7", "14", "21"]));
display(p2_pk);

# secondary comparison with data from Krop et al., 2010; # https://pubmed.ncbi.nlm.nih.gov/20421541/

tdm1_3point6mgkg = DataFrame(
    time_hr = [3,24,48,72,96,168,240,338,408,504],
    T_DM1_ngperml = [63865, 44989.7, 32692.7, 26940.1, 20194.1, 12811.5, 8387.4, 3102.3, 2030.6, 662.3]
);

sdf3 = sims_simp_dose(3.6); 

p3_pk = plot(xlabel = "Time (Day)", ylabel = "ADC conc (ng/mL)", legend = :bottomleft); 
plot!(sdf3.time_day, sdf3.ADC_ugml*1E3, glinewidth = 2, label = "pred", color = "deepskyblue"); 
scatter!(tdm1_3point6mgkg.time_hr/d2h, tdm1_3point6mgkg.T_DM1_ngperml, markerstrokewidth=0, markersize=6, alpha = 0.6, label = "Krop et al., 2010", color = "deepskyblue"); 
plot!(yaxis = :log, ylim = [1E3, 1E5], xticks = ([0, 7, 14, 21], ["0", "7", "14", "21"]));
display(p3_pk);

savefig(plot(p2_pk, p3_pk, ncol = 2), "deliv/figure/T-DM1-pk-fitting.png");
