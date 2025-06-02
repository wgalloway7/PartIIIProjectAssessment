using DelimitedFiles
using Plots
using Statistics
using SpecialFunctions
using LaTeXStrings


include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")
#include("../../experiments/singleflip/singleflip_fit.jl")


beta_values = 1 ./ generate_T_intervals(10.0,0.65,100)

global p_values = [100]

move = "l chain flip"
copies = 20
println(p_values)
L_values = [i for i in 95:105]
L_l_pairs = Dict{Int64, Vector{Int64}}() 
for L in L_values
    for l in 1:L
        if valid_l(L,l) && l % 2 == 1 && l != 1
            if haskey(L_l_pairs, L)
                push!(L_l_pairs[L], l)
            else
                L_l_pairs[L] = [l]
            end
        end
    end
end
L_l_pairs = sort(collect(L_l_pairs), by = x -> x[1])

for (L,l_values) in L_l_pairs
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
            d = readdlm("l=$(l)_p=$(p)_energy_$(L)_$(l)_.csv")
            
            d = reshape(d, length(d))
            if l in [1,l_values[end]]
                plot!(q, beta_values, d, label = "l = $l", color = colors[j])
                plot!(plot2, beta_values, d - onsager_energy.(beta_values, 100), label = "l = $l", color = colors[j])
            else
                plot!(q, beta_values, d, label = "", color = colors[j])
                plot!(plot2, beta_values, d - onsager_energy.(beta_values, 100), label = "", color = colors[j])
            end



            
        end
        d = readdlm("l=$(l)_p=100_energy_$(L)_$(l)_.csv")
        E_plateau = d[end]
        push!(plateaux_energies, E_plateau)
    end
    plot!(q, beta_values, onsager_energy.(beta_values, 100), color = :black, label  = "")
    vline!(q,[log(1+sqrt(2))/2], color = :black, linestyle = :dash, label = L"T_c")
    xlabel!(q, L"\beta = 1/T")
    ylabel!(q, L"E")
    plot!(q, colorbar=true, cbar_ticks=1:length(colors), cbar_label="Color Bar")
    savefig(q, "L=$(L)_energy_vs_beta.png")

    vline!(plot2, [log(1+sqrt(2))/2], color = :black, linestyle = :dash, label = L"T_c")
    xlabel!(plot2, L"\beta = 1/T")
    ylabel!(plot2, L"E_{Onsager} - E")
    plot!(plot2, colorbar=true, cbar_ticks=1:length(colors), cbar_label="Color Bar")
    savefig(plot2, "L=$(L)_energy_vs_beta_log.png")

    writedlm("L=$(L)_plateaux_energies.csv", plateaux_energies, ',')
    plot3 = plot()
    plot!(plot3, l_values, plateaux_energies, label = "Plateaux Energies", color = :black)
    savefig(plot3, "L=$(L)_plateaux_energies.png")
end