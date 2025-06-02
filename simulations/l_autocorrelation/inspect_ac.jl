using DelimitedFiles
using Statistics
using Plots
using Dates
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")

Tc = 2/log(1+sqrt(2)) # = 2.269
beta_values = 1 ./ generate_T_intervals(4.0,Tc,20)
copies = 20
p = plot()
L = 95
l = 3
beta_index = 19
p = plot(legend = false)
q = plot(legend = false)

for i in 1:copies
    data = readdlm("L=$(L)_l=$(l)_copy$(i)_autocorrelation_unconnected.csv", ',')
    row = data[beta_index,:]
    row = row .* 2 .- 1
    plot!(p, row)


    data = readdlm("L=$(L)_l=$(l)_copy$(i)_autocorrelation_E.csv", ',')
    row = data[beta_index,:]
    row = row ./ row[1]
    plot!(q, row)


end
title!(q, "T = $(round(1 / beta_values[beta_index], digits = 2)), L=$(L), l = $(l)")
title!(p, "T = $(round(1 / beta_values[beta_index], digits = 2)), L=$(L), l = $(l)")
xlabel!(p, "Autocorrelation time, MC steps")
xlabel!(q, "Autocorrelation time, MC steps")
ylabel!(p, "Autocorrelation")
ylabel!(q, "Autocorrelation")
savefig(p, "inspect_unc.png")
savefig(q, "inspect_E.png")
