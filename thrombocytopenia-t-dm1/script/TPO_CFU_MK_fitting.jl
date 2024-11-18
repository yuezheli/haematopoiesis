# date: May 3, 2024
# author: Yuezhe Li 
# purpose of this code: to fit for TPO's effect on CFU-MK amplification 

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames
using Parameters: @unpack
using Plots

# model also assume that TPO logrithmic promote CFU-MK proliferation (Broudy et al., 1996; https://pubmed.ncbi.nlm.nih.gov/8822921/)
# data from Broudy et al., 1996; obs = DataFrame(TPO_U_mL = [15, 60, 350, 1000, 2000], CFU_MK_colonies = [3.55, 13.75, 24.9, 34.96, 40.85]);
# assuming 129.4U TPO is 1mg TPO (https://www.researchgate.net/publication/15494676_Defining_units_for_thrombopoietin)
# and TPO molecular weight to be 38kDa 

dat = DataFrame(TPO_uM = [15, 60, 350, 1000, 2000]/129.4/(38E3)*1E6,  # [uM]
                CFU_MK_colonies = [3.55, 13.75, 24.9, 34.96, 40.85]); 

function cfumk_prof!(du, u, p, t)
    @unpack CFUMK01, CFUMK02 = u;
    @unpack tauCFUMK, aCFUMK0, alpha2, beta2, gamma2, TPO = p;
    aCFUMK = aCFUMK0 * log(2.73 + alpha2 + beta2 * TPO)/log(gamma2); 
    tau1CFUMK = (log2(aCFUMK)-1)/log2(aCFUMK) * tauCFUMK;
    tau2CFUMK = tauCFUMK - tau1CFUMK;
    krep1CFUMK = (aCFUMK/2-1)/ tau1CFUMK;
    krep2CFUMK = 1/ tau2CFUMK;
    k12CFUMK = aCFUMK/(2*tau1CFUMK);
    kCFUMK2MK = 2/tau2CFUMK;
    du.CFUMK01 = (krep1CFUMK - k12CFUMK) * CFUMK01; 
    du.CFUMK02 = k12CFUMK * CFUMK01 + (krep2CFUMK - kCFUMK2MK) * CFUMK02;
end

p0 = ComponentArray(tauCFUMK = 12.0,  aCFUMK0 = 11., alpha2 = 0., beta2 = 0., gamma2 = 2.73, TPO = 0.);
u0 = ComponentArray(CFUMK01=2E5, CFUMK02=0);

function murine_cfu_mk_culture(p0 = p0, u0 = u0, sims_days = 5.)
    sol = solve(ODEProblem(cfumk_prof!, u0, (0.0, sims_days), p0), alg = AutoTsit5(Rosenbrock23()), saveat=1.0);
    # sdf = DataFrame(sol); rename!(sdf, Symbol.(names(sdf)[2:end]) .=> collect(keys(sol.u[end])));
    return (sol.u[end][:CFUMK01] + sol.u[end][:CFUMK02])/sum(u0);
end

function loss_fit(p, obs, opt = true, u0 = u0, p0 = p0, sims_days = 5.)
    cfu_amp = [];
    for i in 1:length(obs.TPO_uM)
        p_tmp = deepcopy(p0); 
        p_tmp.tauCFUMK = p[1]
        p_tmp.alpha2 = p[2]
        p_tmp.beta2 = p[3]
        p_tmp.aCFUMK0 = p[4]
        p_tmp.gamma2 = p[5]
        p_tmp.TPO = obs.TPO_uM[i]
        append!(cfu_amp, murine_cfu_mk_culture(p_tmp, u0, sims_days));
    end
    if opt
        return sum( (cfu_amp .- obs.CFU_MK_colonies).^2 )
    else
        return cfu_amp;
    end
end

# optimization 
using Optimization, OptimizationBBO

f_loss = OptimizationFunction(loss_fit, Optimization.AutoForwardDiff());
sol = solve(OptimizationProblem(f_loss, [1., 0., 0., 11., 2.73], dat, ub = [12., 10, 1E20, 100., 100], lb = [0., 0., 0., 1., 1.1]), BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 2E5, maxtime = 1E3, progress = true ); 
# sol = [12., 2.36E-11, 0.24, 66.6, 1.38]

sol_out = loss_fit(sol, dat, false)

p = plot(xlabel = "TPO (uM)", ylabel = "#CFU-Meg colonies/200000 cells");
plot!(dat.TPO_uM, sol_out, label = "fitted");
scatter!(dat.TPO_uM, dat.CFU_MK_colonies, label = "obs");
display(p);

savefig(p, "deliv/figure/TPO-CFUMK-prof-fit.png");

# fitting between platelet number and TPO level 
# data obtained from Engel et al., 1999, Fig 3; # https://pubmed.ncbi.nlm.nih.gov/10354155/
dat2 = DataFrame(platelet_e9_L = [78.61, 301.7], TPO_pg_mL = [240.3, 31.17]); # this is 2 dots from the linear regression line 
beta1 = ( log(dat2.TPO_pg_mL[2]/38E3) - log(dat2.TPO_pg_mL[1]/38E3) ) / (dat2.platelet_e9_L[2] - dat2.platelet_e9_L[1]); # TPO in nM
# beta1 = -0.009
alpha1 = log(dat2.TPO_pg_mL[2]/38E3) - beta1 * dat2.platelet_e9_L[2]
# alpha1 = -4.34
# note this alpha1 is not the same as the one in `model/platelet_dynamics_TPO1.jl`. That one needs to adjust for kdegTPO