using DelimitedFiles
using Plots
using LaTeXStrings

include("../../core/lattice.jl")
L_values = [i for i in 95:105]
L_l_pairs = Dict{Int64, Vector{Int64}}() 
for L in L_values
    #L_l_pairs[L] = [1]
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
for (L, l_values) in L_l_pairs
    for l in l_values
        data,header = readdlm("l=$(l)_saddles_simulation_$(L).csv",'\t', header = true)
        E_row = data[:,1]
        saddle_row = data[:,2]
        p = plot(legend = false)
        scatter!(p, E_row, saddle_row, markersize = 1,color =RGB(0.2, 0.2, 0.8), markerstrokewidth = 0, alpha= 0.8)

        xlabel!(p, latexstring("Energy, \$E\$"))
        ylabel!(p, latexstring("Saddle index, \$K\$"))

        savefig(p, "energy_saddles_$(L)_$(l).png")

    

    end
end


L = 95
l_values = [1,9,27,63]

mid_red    = RGB(0.9, 0.4, 0.4)    # warm, mid-dark red
mid_blue   = RGB(0.2, 0.2, 0.8)    # cool, mid-dark blue
gradient = cgrad([mid_blue, mid_red])
N = length(l_values)
colors = [gradient[i] for i in range(0, 1, length=N)]

p = plot(legend = :topleft)
for (i,l) in enumerate(l_values)
    data,header = readdlm("l=$(l)_saddles_simulation_$(L).csv",'\t', header = true)
    E_row = data[:,1]
    saddle_row = data[:,2]

    scatter!(p, E_row, log10.(saddle_row), markersize = 1, markerstrokewidth = 0, alpha= 0.8, label = "l = $l", color = colors[i])

    xlabel!(p, latexstring("Energy, \$E\$"))
    ylabel!(p, latexstring("log_{10}\\,(\\mathrm{Saddle\\ index}),\\ \\log_{10} K"))


    savefig(p, "energy_saddles_95_selected_l.png")

end
