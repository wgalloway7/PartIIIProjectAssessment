using DelimitedFiles
using Plots
using LaTeXStrings
include("../../core/lattice.jl")
L_values = [i for i in 95:105]
L_l_pairs = Dict{Int64, Vector{Int64}}() 
Tc_inf = 2 / log(1 + sqrt(2))  # Critical temperature for the 2D Ising model


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
    plateau_data = readdlm("L=$(L)_plateaux_energies.csv")
    new_onset_data = readdlm("z_new_onset_E_$(L).csv", ',')

    writedlm("zOnset_energies_$L.csv", [["l", "E_K0", "E_K2p"]], ',')
    for (i,l) in enumerate(l_values)

        new_onset_E = new_onset_data[i, 1]

        data = readdlm("proportions_l=$(l)_saddles_simulation_$(L).csv", '\t'; skipstart = 1)
        energy_set = data[:, 1]
        pK0    = data[:, 2]
        pK1    = data[:, 3]
        pK2p   = data[:, 4]

        sorted_indices = sortperm(energy_set)

        function first_sustained_energy(energies::Vector{Float64}, values::Vector{Float64}; direction::Symbol)

            sorted_pairs = sort(collect(zip(energies, values)), by = x -> x[1])
            if direction == :down
                sorted_pairs = reverse(sorted_pairs)
            end
        
            streak = 0
            for i in 1:length(sorted_pairs)
                _, val = sorted_pairs[i]
                if val > 0
                    streak += 1
                    if streak == 3
                        return sorted_pairs[i - 2][1]
                    end
                else
                    streak = 0
                end
            end
            return nothing
        end
        
        E_K2p_first = first_sustained_energy(energy_set, pK2p; direction = :up)
        E_K0_first  = first_sustained_energy(energy_set, pK0; direction = :down)


        nbins = 200
        emin, emax = minimum(energy_set), maximum(energy_set)
        edges = range(emin, emax; length = nbins + 1) |> collect
        bin_centers = (edges[1:end-1] .+ edges[2:end]) ./ 2


        sum_pK0  = zeros(nbins)
        sum_pK1  = zeros(nbins)
        sum_pK2p = zeros(nbins)
        count_in_bin = zeros(Int, nbins)


        for i in eachindex(energy_set)
            E = energy_set[i]
            # find bin index: edges[b] ≤ E < edges[b+1], except E == emax → bin = nbins
            bin = findfirst(b -> E ≥ edges[b] && E < edges[b+1], 1:nbins)
            if bin === nothing
                if E == emax
                    bin = nbins
                else
                    continue
                end
            end
            sum_pK0[bin]  += pK0[i]
            sum_pK1[bin]  += pK1[i]
            sum_pK2p[bin] += pK2p[i]
            count_in_bin[bin] += 1
        end


        binned_centers = Float64[]
        avg_pK0  = Float64[]
        avg_pK1  = Float64[]
        avg_pK2p = Float64[]

        for bin in 1:nbins
            N = count_in_bin[bin]
            if N > 0
                push!(binned_centers, bin_centers[bin])
                push!(avg_pK0,  sum_pK0[bin]  / N)
                push!(avg_pK1,  sum_pK1[bin]  / N)
                push!(avg_pK2p, sum_pK2p[bin] / N)
            end
        end


        scatter(
            binned_centers, avg_pK2p;
            label       = L"K ≥ 2",
            markershape = :diamond,
            markersize = 3,
        )
        
        scatter!(
            binned_centers, avg_pK1;
            label       = L"K = 1",
            markershape = :circle,
            markersize = 3,
        )
        if l < L/2
            plot!(legend = :right)
        else
            plot!(legend = :left)
        end
        scatter!(
            binned_centers, avg_pK0;
            label       = L"K = 0",
            markershape = :square,
            markersize = 3,
        )
        
        vline!([plateau_data[i]], label = "Plateau energy")

        xlabel!(L"Energy, $E$")
        ylabel!(latexstring("Saddle index proportion, \$p_K(E)\$"))
        title!(latexstring("\$L = $(L),\\ l = $(l)\$"))

        #E_VFT = onsager_energy(1 / VFT_temps[i], L)
        #vline!([E_VFT], label = "VFT energy", color = :red, linestyle = :dash)


        #if E_K0_first !== nothing
        #    vline!([E_K0_first];  label = L"E_{on}",  linestyle = :dash,
        #                          color = :black)
        #    annotate!([E_K0_first+0.075], [0.6], text(L"E_{on}", 8, :black))
        #end                         
        #if E_K2p_first !== nothing
        #    vline!([E_K2p_first]; label = L"E_{min}", linestyle = :dash,
        #                          color = :black)    
        #    annotate!([E_K2p_first-0.075], [0.6], text(L"E_{min}", 8, :black))
        #end
        
        vline!([new_onset_E]; label = L"E_{on}", color = :black, linestyle = :dash)
        #annotate!([new_onset_E+0.075], [0.6], text(L"E_{on}", 8, :black))
        vline!([plateau_data[i]]; label = L"E^{*}", color = :red, linestyle = :dash)
        E_onsa = -sqrt(2) - 0.622/L

        vline!([E_onsa]; label = L"E_{mag}", color = :blue, linestyle = :dash)


        savefig("proportions_histogram_$(L)_$(l).png")

        open("zOnset_energies_$L.csv", "a") do io
            writedlm(io, [[l,E_K0_first, E_K2p_first]], ',')
        end
    end
end

