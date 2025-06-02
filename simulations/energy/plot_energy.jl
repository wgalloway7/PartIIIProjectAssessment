using DelimitedFiles
using Plots
using Statistics
using SpecialFunctions
using LaTeXStrings


include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")
#include("../../experiments/singleflip/singleflip_fit.jl")

beta_values = 1 ./ generate_T_intervals(4.0, 1.0, 100)
global p_values = [16]
l_values = [2*i - 1 for i in 1:50]
q = plot()
plot2 = plot()



#colors = (cgrad(:RdBu, length(p_values)))
mid_red    = RGB(0.9, 0.4, 0.4)    # warm, mid-dark red
mid_blue   = RGB(0.2, 0.2, 0.8)    # cool, mid-dark blue
gradient = cgrad([mid_blue, mid_red])
N = length(l_values)
colors = [gradient[i] for i in range(0, 1, length=N)]

q = plot()
plot2 = plot()
plateaux_energies = []
for (j,l) in enumerate(l_values)
    for (i,p) in enumerate(p_values)
        d = readdlm("l=$(l)_p=$(p)_energy_p_$(l)_.csv")
        d = reshape(d, length(d))
        if l in [1,l_values[end]]
            plot!(q, beta_values, d, label = "l = $l", color = colors[j])
            plot!(plot2, beta_values, d - onsager_energy.(beta_values, 100), label = "l = $l", color = colors[j])
        else
            plot!(q, beta_values, d, label = "", color = colors[j])
            plot!(plot2, beta_values, d - onsager_energy.(beta_values, 100), label = "", color = colors[j])
        end



        
    end
    d = readdlm("l=$(l)_p=16_energy_p_$(l)_.csv")
    E_plateau = d[end]
    push!(plateaux_energies, E_plateau)
end
plot!(q, beta_values, onsager_energy.(beta_values, 100), color = :black, label  = "")
vline!(q,[log(1+sqrt(2))/2], color = :black, linestyle = :dash, label = L"T_c")
xlabel!(q, L"\beta = 1/T")
ylabel!(q, L"E")
plot!(q, colorbar=true, cbar_ticks=1:length(colors), cbar_label="Color Bar")
savefig(q, "energy_vs_beta.png")

vline!(plot2, [log(1+sqrt(2))/2], color = :black, linestyle = :dash, label = L"T_c")
xlabel!(plot2, L"\beta = 1/T")
ylabel!(plot2, L"E_{Onsager} - E")
plot!(plot2, colorbar=true, cbar_ticks=1:length(colors), cbar_label="Color Bar")
savefig(plot2, "energy_vs_beta_log.png")

writedlm("plateaux_energies.csv", plateaux_energies, ',')
plot3 = plot()
plot!(plot3, l_values, plateaux_energies, label = "Plateaux Energies", color = :black)
savefig(plot3, "plateaux_energies.png")