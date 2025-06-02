using DelimitedFiles
using Plots
using LaTeXStrings
using LsqFit
include("../../core/lattice.jl")
L_values = [10,20,23,25,46,47,50,52]
L_values = 10:50
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


lchoice = 5
plot_obj = plot(legend = :right)
for (L, l_values) in L_l_pairs
    
    sigmoid_fit_data,header = readdlm("sigmoid_fit_k2_L=$(L).csv", '\t',header=true)
    writedlm("zOnset_energies_$L.csv", [["l", "E_K0", "E_K2p"]], ',')
    for (i,l) in enumerate(l_values)

        alpha = sigmoid_fit_data[i, 2]
        x0 = sigmoid_fit_data[i, 3]
        n_exponent = sigmoid_fit_data[i, 4]

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


        function fss(x)
            return alpha * (x -x0)
        end
        
        if l == lchoice
            scatter!(plot_obj,
                fss.(binned_centers), avg_pK2p;
                label       = "L = $L",
                markershape = :diamond,
                markersize = 3,
                markerstrokewidth = 0.5
            )
            #vline!([plateau_data[i]], label = "Plateau energy")
            xlims!(plot_obj, (fss(minimum(binned_centers)), 10))
            xlabel!(plot_obj, latexstring("\$\\theta(L)(E - E_{sig})\$"))
            ylabel!(plot_obj, latexstring("Saddle index proportion \$ p_{K > 1}(E)\$"))


            #E_VFT = onsager_energy(1 / VFT_temps[i], L)
            #vline!([E_VFT], label = "VFT energy", color = :red, linestyle = :dash)
        end

    

        
    end
    
end
title!(plot_obj, latexstring("\$l = $(lchoice)\$"))
savefig(plot_obj,"zfss_$lchoice.png")


q = plot()
params_dict = Dict{Int64, Tuple{Vector{Float64}, Vector{Float64}}}()
for L in L_values
    try
        sigmoid_fit_data,header = readdlm("sigmoid_fit_k2_L=$(L).csv", '\t',header=true)
        l_values = sigmoid_fit_data[:, 1]
        alpha_values = sigmoid_fit_data[:, 2]
        x0_values = sigmoid_fit_data[:, 3]
        for (i,l) in enumerate(l_values)
            alpha = alpha_values[i]
            x0 = x0_values[i]
            if !(l  in keys(params_dict))
                params_dict[l] = (Float64[], Float64[])
            end
            push!(params_dict[l][1], alpha)
            push!(params_dict[l][2], float(L))
        end
    catch e
        println("Error reading data for L = $L: $e")
        continue
    end
end

for l in [3, 5, 7, 11]
    alphas, Ls = params_dict[l]          # same as before

    idx           = sortperm(Ls)         # sort by L
    Ls_sorted     = Ls[idx]
    alphas_sorted = alphas[idx]

    scatter!(q, Ls_sorted, alphas_sorted;
             label="ℓ = $l",
             markershape=:diamond,
             markersize=3)

    #model(x, p) = p[1] .* (x .- p[2])    # p = [a, b]

    #p0  = [0.1, 10.0]                    # initial guess
    #fit = curve_fit(model, Ls_sorted, alphas_sorted, p0)

    #x_fit = range(0, last(Ls_sorted); length=100)
    #y_fit = model(x_fit, fit.param)      # predict with fitted params

    #plot!(q, x_fit, y_fit; label="Fit ℓ = $l", linestyle=:dash)
end
xlabel!(q, L"System size, $L$")
ylabel!(q, L"Sigmoid slope, $\theta$")
savefig(q, "zalphas.png")


