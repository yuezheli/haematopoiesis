# date: May 7, 2024
# author: Yuezhe Li 
# purpose of this script: to develop in vitro model on how T-DM1 impede MK development

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames
using Parameters: @unpack
using Plots
using Optimization, OptimizationBBO

# add observed data; digitized from Uppal et al., 2015, Fig 3C (timeline scaled to 0); https://pubmed.ncbi.nlm.nih.gov/25370470/
control = DataFrame(
    time_day = [0., 3., 6., 9., 14.],
    normalized_mk = [1, 1, 1.22, 2.07, 2.37]
); 
tdm1_treated = DataFrame(
    time_day = [0., 3., 6., 9., 14.],
    normalized_mk = [1, 0.71, 0.64, 0.37, 0.3]
); 

# define simplified CFU-MK, MK, and T-DM1 
function tdm1_mk_dev_simp!(du, u, p, t)
    @unpack CFUMK01, CFUMK02, MK0, cyto_tdm1, medium_tdm1 = u; 
    @unpack tauCFUMK, tauMK, aCFUMK0, k_in_tdm1, IC50_tdm1, n_tdm1 = p
    # CFU-MK dynamics (simplified)
    aCFUMK = aCFUMK0 * ( IC50_tdm1^n_tdm1./(IC50_tdm1^n_tdm1 .+ max(cyto_tdm1, 0.).^n_tdm1) )
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
    # MK dynamics 
    kdeathMK = 1/tauMK
    # ODEs
    du.CFUMK01 = (krep1CFUMK - k12CFUMK) * CFUMK01; 
    du.CFUMK02 = k12CFUMK * CFUMK01 + (krep2CFUMK - kCFUMK2MK) * CFUMK02;
    du.MK0   = kCFUMK2MK * CFUMK02 - kdeathMK * MK0;  
    du.cyto_tdm1 = k_in_tdm1 * medium_tdm1 - kCFUMK2MK * cyto_tdm1
    du.medium_tdm1 = 0.  # assuming medium T-DM1 is overly abundant and uptake from cells would not impact the medium conc
end

p0 =  ComponentArray(tauCFUMK = 12.0, tauMK = 0.1, aCFUMK0 = 11.0, k_in_tdm1 = 0., IC50_tdm1 = 20, n_tdm1 = 1.); # IC50_tdm1 in nM

init_0 = ComponentArray(CFUMK01=5.2E5, CFUMK02=0., MK0=2.3E5, cyto_tdm1 = 0., medium_tdm1 = 0.); # medium T-DM1 conc in nM

# optimization of doubling time and doubling time (to suit for the culture condition)
function loss_fit_0(p, obs = control, opt = true, p0 = p0, init_0 = init_0)
    p_tmp = deepcopy(p0); p_tmp.tauCFUMK = p[1]; p_tmp.aCFUMK0 = p[2]; p_tmp.tauMK = p[3];
    init_tmp = deepcopy(init_0); init_tmp.CFUMK01 = p[4]; init_tmp.CFUMK02 = 5.2E5 - p[4]; 
    sol0 = solve(ODEProblem(tdm1_mk_dev_simp!, init_tmp, (0.0, 14.0), p_tmp), alg = AutoTsit5(Rosenbrock23()), saveat=obs.time_day);
    sdf0 = DataFrame(sol0); rename!(sdf0, Symbol.(names(sdf0)[2:end]) .=> collect(keys(sol0.u[end])));
    sdf0.norm_mk .= sdf0.MK0 ./ sdf0.MK0[1]
    if opt
        return sum((sdf0.norm_mk .- obs.normalized_mk).^2)
    else
        return sdf0
    end
end

# sol0 = solve(OptimizationProblem(OptimizationFunction(loss_fit_0, Optimization.AutoForwardDiff()), [6., 11., 0.1, 1.2E5], control, ub = [24., 1000., 2, 5.2E5], lb = [0., 0., 0., 0.]), BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 2E5, maxtime = 1E3, progress = true ); 

pred0 = loss_fit_0([24., 15.21, 2.1, 5.2E5], control, false);

p_noTDM1 = plot(xlabel = "Time (Days)", ylabel = "Normalized MK cell count", legend = :outerright, ylim = [0, 3]);
plot!(pred0.timestamp, pred0.norm_mk, label = "predicted, no T-DM1");
plot!(control.time_day, control.normalized_mk, seriestype=:scatter, label="Uppal et al., 2015");
display(p_noTDM1);

# optimization for T-DM1 uptake and drug effect 
p1 = deepcopy(p0);
p1.tauCFUMK = 24.
p1.tauMK = 2.1
p1.aCFUMK0 = 15.21
init_1 = ComponentArray(CFUMK01=5.2E5, CFUMK02=0., MK0=2.3E5, cyto_tdm1 = 0., medium_tdm1 = 41.6); # medium T-DM1 conc in nM

function loss_fit_1(p, obs = tdm1_treated, opt = true, p1 = p1, init_1 = init_1)
    p_tmp = deepcopy(p1); p_tmp.k_in_tdm1 = p[1]; p_tmp.n_tdm1 = p[2]; 
    sol0 = solve(ODEProblem(tdm1_mk_dev_simp!, init_1, (0.0, 14.0), p_tmp), alg = AutoTsit5(Rosenbrock23(autodiff = false)), saveat=obs.time_day);
    sdf0 = DataFrame(sol0); rename!(sdf0, Symbol.(names(sdf0)[2:end]) .=> collect(keys(sol0.u[end])));
    sdf0.norm_mk .= sdf0.MK0 ./ sdf0.MK0[1]
    if opt
        return sum((sdf0.norm_mk .- obs.normalized_mk).^2)
    else
        return sdf0
    end
end

# sol1 = solve(OptimizationProblem(OptimizationFunction(loss_fit_1, Optimization.AutoForwardDiff()), [0.1, 5], tdm1_treated, ub = [0.6, 10], lb = [0., 0.02]), BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 2E5, maxtime = 1E3, progress = true ); 
# sol1 = [0.37, 1.13]
pred1 = loss_fit_1([0.37, 1.13], tdm1_treated, false)

p_TDM1 = plot(xlabel = "Time (Days)", ylabel = "Normalized MK cell count", legend = :outerright, ylim = [0, 3]);
plot!(pred0.timestamp, pred0.norm_mk, label = "predicted, no T-DM1", color = :blue);
plot!(control.time_day, control.normalized_mk, seriestype=:scatter, label="Uppal et al., 2015 (control)", color = :blue);
plot!(pred1.timestamp, pred1.norm_mk, label = "predicted, T-DM1", color = :red);
plot!(tdm1_treated.time_day, tdm1_treated.normalized_mk, seriestype=:scatter, label="Uppal et al., 2015 (T-DM1)", color = :red);
display(p_TDM1);

savefig(p_TDM1, "deliv/figure/T-DM1-CFU-MK-fitting.png");

