# date: May 2, 2024
# author: Yuezhe Li 
# purpose of this code: to troubleshoot platelet development model
# TPO half life obtained from https://www.ncbi.nlm.nih.gov/books/NBK12518/
# normal TPO and platelet number: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4339937/
# normal TPO level in serum: 3.18E-12 mol/L => 3.18E-3 nM (ranges in 2.13E-3 and 6.25E-3 nM)
# normal platelet level in blood: 193E9 L-1 (ranges in 130E9/L and 332E9/L)
# beta1 was black calculated based on data from Engel et al., 1999 # https://pubmed.ncbi.nlm.nih.gov/10354155/

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames
using Plots
using Optimization, OptimizationBBO

include("model/platelet_dynamics_TPO1.jl")

p_mk =  ComponentArray(tauLT = 100.0, tauST = 20.0, tauMPP = 2.0, tauCMP = 4.0, tauBFUMK = 9.0, tauCFUMK = 12.0, tauMK = 0.1, tauPlatelet = 7.0, thalfTPO = 1.25, 
    ssLT = 1275.0, rLT = 1/(2.5*7), aST = 1000.0, aMPP = 380.0, aCMP = 4.0, aBFUMK = 8., aCFUMK0 = 11., 
    kMKplatelet = 450., alpha1 = -3.5, beta1 = 0.009, beta2 = 4.72E-3
); 

init_mk = ComponentArray(LT0 = 1E3, ST01 = 0, ST02=0, MPP01=0, MPP02=0, CMP01=0, CMP02=0, BFUMK01=0, BFUMK02=0, CFUMK01=0, CFUMK02=0, MK0=0, platelet0=0, TPO=0);

function loss_fit(p, obs, opt = true, p_mk = p_mk, init_mk = init_mk)
    p_tmp = deepcopy(p_mk); 
    p_tmp.alpha1 = p[1]
    p_tmp.aCFUMK0 = p[2]
    p_tmp.kMKplatelet = p[3]
    # simulate to steady state
    sol = solve(ODEProblem(mk_dev!, init_mk, (0.0, 600.0), p_tmp), alg = AutoTsit5(Rosenbrock23()), saveat=1.0);
    diff1 = (sol.u[end].platelet0/V_blood - obs[1])/obs[1]
    diff2 = (sol.u[end].TPO/V_plasma - obs[2]) / obs[2]  # normalized TPO diff (TPO in nM)
    diff3 = (sol.u[end].MK0 - obs[3]) / obs[3] 
    if opt 
        return diff1*diff1 + diff2*diff2 + diff3*diff3
    else
        return sol.u[end].platelet0/V_blood, sol.u[end].TPO/V_plasma, sol.u[end].MK0, sol.u[end]
    end
end

# sol = solve(OptimizationProblem(OptimizationFunction(loss_fit, Optimization.AutoForwardDiff()), [-3.5, 11., 480.], [193E9, 3.18E-3, 1.2E8], ub = [0, 12., 600.], lb = [-10, 8., 480.]), BBO_adaptive_de_rand_1_bin_radiuslimited(), maxiters = 2E5, maxtime = 1E3, progress = true ); 
# sol = [-2.86, 8.0, 480.];

platelet_perL, TPO_nM, MK_base, BFUMK2_base = loss_fit(sol, [193E9, 3.18E-3, 1.2E8], false)
#(2.59e11, 0.00317, 3.86e8, (LT0 = 1274.9999954564246, ST01 = 229.4124441777085, ST02 = 12793.774433450413, MPP01 = 22524.453214759884, MPP02 = 565353.7307239509, CMP01 = 9.689999645645952e6, CMP02 = 1.937999921549127e7, BFUMK01 = 1.1627999336736813e8, BFUMK02 = 2.3255998448042998e8, CFUMK01 = 1.4595757028927734e9, CFUMK02 = 5.000386348496645e9, MK0 = 3.867557677623322e8, platelet0 = 1.2994992845111394e12, TPO = 0.009939413585182959)

#=
tau1BFUMK = (log2(p_mk.aBFUMK)-1)/log2(p_mk.aBFUMK) * p_mk.tauBFUMK;
tau2BFUMK = p_mk.tauBFUMK - tau1BFUMK;
kBFUMK2CFUMK = 2/tau2BFUMK;
kBFUMK2CFUMK * BFUMK2_base # 1.55E8
=#