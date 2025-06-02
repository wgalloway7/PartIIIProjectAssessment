using DelimitedFiles
using Plots
using LaTeXStrings

include("../../core/lattice.jl")
#onset energies
#sigmoid fit energies
#plateau energies (already got for L = 95=105)

#maybe VFT energies from VFT fit?


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

all_sig_energies = []
all_plateau_energies = []
all_l_ratios = []
all_k2_energies = []
all_k0_energies = []
all_new_onset_energies = []
for (L, l_values) in L_l_pairs
    sigmoid_fit_data,_ = readdlm("sigmoid_fit_k2_L=$(L).csv", '\t',header=true)
    sigmoid_energies = sigmoid_fit_data[:,3]
    plateau_data = readdlm("L=$(L)_plateaux_energies.csv")
    
    onset_data,_ = readdlm("zOnset_energies_$L.csv", ',', header=true)
    onset_min_energies = onset_data[:,2]
    onset_saddles_energies = onset_data[:,3]
    println(onset_min_energies)

    new_onset_data = readdlm("z_new_onset_E_$(L).csv", ',')
    new_onset_energies = new_onset_data[:,1]
    for (i,l) in enumerate(l_values)
        plateau_energy = plateau_data[i]
        sigmoid_energy = sigmoid_energies[i+1]
        println("L = $L, l = $l, plateau energy = $plateau_energy, sigmoid energy = $sigmoid_energy")
        push!(all_sig_energies, sigmoid_energy)
        push!(all_plateau_energies, plateau_energy)
        push!(all_l_ratios, l / L)
        push!(all_k0_energies, onset_min_energies[i])
        push!(all_k2_energies, onset_saddles_energies[i])
        push!(all_new_onset_energies, new_onset_energies[i])

    end
end
p = plot()

#p = scatter!(p, all_l_ratios, all_k0_energies, label=latexstring("E_{-}"), marker = :pentagon)
p = scatter!(p, all_l_ratios, all_new_onset_energies, label=latexstring("E_{onset}"), marker = :diamond)
p = scatter(p, all_l_ratios, all_sig_energies, label=latexstring("E_{sig}"), marker = :circle)

p = scatter!(p, all_l_ratios, all_plateau_energies, label=latexstring("E*"), marker = :square)

#p = scatter!(p, all_l_ratios, all_k2_energies, label=latexstring("E_{+}"), marker = :utriangle)
xlabel!(p, L"l/L")
ylabel!(p, latexstring("Energy, \$E\$"))
savefig(p, "topo.png")
