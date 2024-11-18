# Import packages required to run the model
# note: these packages need to be imported only ONCE in each Julia session.
using DifferentialEquations
using RecursiveArrayTools
using Plots
using DataFrames
using CSV

include("scdODEs.jl");

tspan = (0.0, 24.0); # unit: hr


u0 = vcat([1], zeros(188));


prob = ODEProblem(scdODEs!, u0, tspan); 

sol = solve(prob, Tsit5(), reltol = 1e-2);

LT_HSC = [sol.u[i][1] for i in 1:length(sol.t)];
ST_HSC = [(sol.u[i][2] + sol.u[i][3]) for i in 1:length(sol.t)];

# read in simulation result from mrgsolve for comparison
df = CSV.read("../data/sickle_simul.csv", DataFrame)

plot_lthsc = plot(sol.t, LT_HSC, linewidth = 1, xlabel = "time (h)", legend=:bottomright, 
                 ylabel = "LT-HSC cell number", label = "Julia simulation")
plot!(df[:,1], df[:,2], label = "mrgsolve solution", color = 2, linestyle = :dot, linewidth = 4)


plot_sthsc = plot(sol.t, ST_HSC, linewidth = 1, xlabel = "time (h)", legend=:bottomright, 
                 ylabel = "ST-HSC cell number", label = "Julia simulation")
plot!(df[:,1], df[:,3] + df[:,4], label = "mrgsolve solution", color = 2, linestyle = :dot, linewidth = 4)

plotd = plot(plot_lthsc, plot_sthsc, layout = @layout [a; b])

savefig(plotd,"../img/julia_HSC.png")
