# date: May 1, 2024
# author: Yuezhe Li 
# purpose of this code: to troubleshoot platelet development model

using Pkg; Pkg.activate("");

using DifferentialEquations, ComponentArrays, DataFrames
using Plots

include("model/platelet_dynamics_noTPO.jl")

p_mk =  ComponentArray(tauLT = 100.0, tauST = 20.0, tauMPP = 2.0, tauCMP = 4.0, tauBFUMK = 9.0, tauCFUMK = 12.0, tauMK = 0.1, tauPlatelet = 4.0,
    ssLT = 1275.0, rLT = 1/(2.5*7), aST = 1000.0, aMPP = 380.0, aCMP = 4.0, aBFUMK = 8., aCFUMK = 11., 
    kMKplatelet = 450.
); 

init_mk = ComponentArray(LT0 = 1E3, ST01 = 0, ST02=0, MPP01=0, MPP02=0, CMP01=0, CMP02=0, BFUMK01=0, BFUMK02=0, CFUMK01=0, CFUMK02=0, MK0=0, platelet0=0);

prob = ODEProblem(mk_dev!, init_mk, (0.0, 600.0), p_mk);
sol = solve(prob, alg = AutoTsit5(Rosenbrock23()), saveat=1.0);
sdf = DataFrame(sol); rename!(sdf, Symbol.(names(sdf)[2:end]) .=> collect(keys(sol.u[end])));

p_steadystate = plot(xlabel = "Time (Day)", ylabel = "Cell count (#)", legend = :outerright, ylims = [1E2, 1E12],yaxis = :log10);
plot!(sdf.timestamp, sdf.LT0, label = "Long-term HSC");
plot!(sdf.timestamp, sdf.ST01 .+ sdf.ST02, label = "Short-term HSC");
plot!(sdf.timestamp, sdf.MPP01 .+ sdf.MPP02, label = "MPP");
plot!(sdf.timestamp, sdf.CMP01 .+ sdf.CMP02, label = "CMP");
plot!(sdf.timestamp, sdf.BFUMK01 .+ sdf.BFUMK02, label = "BFU-MK");
plot!(sdf.timestamp, sdf.CFUMK01 .+ sdf.CFUMK02, label = "CFU-MK");
plot!(sdf.timestamp, sdf.MK0, label = "MK");
plot!(sdf.timestamp, sdf.platelet0, label = "Platelet");
display(p_steadystate);

savefig(p_steadystate, "deliv/figure/v1-steadystatetest.png");

# check steady state platelet count 
sdf.platelet0[end]/5E3  # platelet count per mL blood
# 61E6

#=
sdf[end, 2:end]
DataFrameRow
 Row │ LT0      ST01     ST02     MPP01    MPP02      CMP01    CMP02    BFUMK01   BFUMK02   CFUMK01    CFUMK02   MK0        platelet0  
     │ Float64  Float64  Float64  Float64  Float64    Float64  Float64  Float64   Float64   Float64    Float64   Float64    Float64    
─────┼─────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
 601 │  1275.0  229.412  12793.8  22524.5  5.65354e5   9.69e6  1.938e7  1.1628e8  2.3256e8  1.32268e9  2.9579e9  1.70544e8  3.06979e11
=#

# estimation on MK and platelet number 
# assuming MK cell is 1 in 1E4 per bone marrow cell (https://pubmed.ncbi.nlm.nih.gov/1060175/)
# bone marrow cellularity = 0.4E6 mm^3 (https://academic.oup.com/ajcp/article-abstract/61/2/199/1770391?redirectedFrom=PDF)
# bone marrow volume = 3L (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5738992/)
# this translate into total MK cell number 1.2E8

# simulation on tauPlatelet = 7. (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7589436/)
p_mk2 = deepcopy(p_mk); p_mk2.tauPlatelet = 7.
sol = solve(ODEProblem(mk_dev!, init_mk, (0.0, 600.0), p_mk2), alg = AutoTsit5(Rosenbrock23()), saveat=1.0);
sol.u[end][:platelet0]/5E3 # 1.07e8

p_mk3 = deepcopy(p_mk); p_mk3.tauPlatelet = 10.
sol = solve(ODEProblem(mk_dev!, init_mk, (0.0, 600.0), p_mk3), alg = AutoTsit5(Rosenbrock23()), saveat=1.0);
sol.u[end][:platelet0]/5E3 # 1.5e8

p_mk4 = deepcopy(p_mk); p_mk4.tauPlatelet = 10.; p_mk4.kMKplatelet = 900.
sol = solve(ODEProblem(mk_dev!, init_mk, (0.0, 600.0), p_mk4), alg = AutoTsit5(Rosenbrock23()), saveat=1.0);
sol.u[end][:platelet0]/5E3 # 3e8

