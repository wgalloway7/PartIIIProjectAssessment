using DelimitedFiles
using Plots
using LaTeXStrings
include("../../core/lattice.jl")
L_values = vcat(10:50,95:105)
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
l_choice = 5
L_choice = [12,27,48,97]

mid_red    = RGB(0.9, 0.4, 0.4)
mid_blue   = RGB(0.2, 0.2, 0.8)   
gradient = cgrad([mid_blue, mid_red])
colors = [gradient[i] for i in range(0, 1, length(L_choice))]
#colors = ["red","blue", "green"]

p = plot(legend = :right)
global color_count = 1
for (L, l_values) in L_l_pairs
    plateau_data = readdlm("L=$(L)_plateaux_energies.csv")
    

    for (i,l) in enumerate(l_values)
        if l == l_choice && L in L_choice
            sigmoid_fit_data,header = readdlm("sigmoid_fit_k2_L=$(L).csv", '\t',header=true)
            x0_data = sigmoid_fit_data[:, 3]
            alpha_data = sigmoid_fit_data[:, 2]
            x0 = x0_data[i]

            data = readdlm("proportions_l=$(l)_saddles_simulation_$(L).csv", '\t'; skipstart = 1)
            energy_set = data[:, 1]
            energy_set .-= x0  
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
                        if streak == 2
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


            nbins = 50
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

        
            scatter!(p,
                binned_centers, avg_pK2p;
                label       = "L = $(L)",
                markershape = :diamond,
                markersize = 3, color = colors[color_count]
            )
            vline!([plateau_data[i]], label = "Plateau energy")

            xlabel!(L"Rescaled energy, $E - E_{sig}$")
            ylabel!(latexstring("Saddle index proportion, \$p_{K>1}(E)\$"))
            title!(latexstring("\$l = $(l)\$"))
            xlims!(p, (-1.0, 1.0))
            xdata = values = range(-1.0, stop=1.0, length=100)
    
            sigmoid(x, p) = @. 1 / (1 + exp(- p[1] * (x - p[2])))
            plot!(p, xdata, sigmoid(xdata, [alpha_data[i], 0.0]), label = "", color = colors[color_count])
            println(color_count)
            global color_count += 1
            
        end
            


        

            

            

    end
end

savefig(p,"zproportions_histogram_zero_center.png")
