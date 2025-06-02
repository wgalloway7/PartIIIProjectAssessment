using DelimitedFiles
using Plots
using LaTeXStrings
include("../../core/lattice.jl")
L_values = [10,20,23,25,46,47,50,52]
L_values = 10:50
L_values = vcat(10:50, 95:105)
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
# (L, l_values) in L_l_pairs
#    println("L = $L, l values = $l_values")
#end

L_l_pairs = sort(collect(L_l_pairs), by = x -> x[1])
p = plot()
for l_choice in [2i + 1 for i in 1:5]
    L_to_plot = Float64[]
    spread_to_plot = Float64[]
    for (L,l_values) in L_l_pairs
        onset_data,_ = readdlm("zOnset_energies_$L.csv", ',', header = true)
        EK0_data = onset_data[:, 2]
        EK2p_data = onset_data[:, 3]
        for (i, l) in enumerate(l_values)
            if l == l_choice
                E_K0 = EK0_data[i]
                E_K2p = EK2p_data[i]
                spread = E_K0 - E_K2p
                push!(spread_to_plot, spread)
                push!(L_to_plot, L)
            end
        end
    end

    scatter!(1 ./ L_to_plot, spread_to_plot, label = "l = $l_choice", xlabel = L"1/L", ylabel = latexstring("Crossover spread, \$E_{on} - E_{min}\$"),
     markersize = 3, markerstrokewidth  = 0.5)

    #function fit_func(x, a, b)
    #    return a * x .+ b
    #end
    #fit_params = curve_fit(fit_func, 1 ./ L_to_plot, spread_to_plot)
    #fit_x = range(0, stop = 1 / minimum(L_to_plot), length = 100)
    #fit_y = fit_func(fit_x, fit_params.param[1], fit_params.param[2])
    #plot!(fit_x, fit_y, label = "Fit for l = $l_choice", linestyle = :dash, color = :black)

end
savefig(p, "zfss_energy_spread_choice.png")



L = 46
mid_red    = RGB(0.9, 0.4, 0.4)
mid_blue   = RGB(0.2, 0.2, 0.8)
gradient = cgrad([mid_blue, mid_red])
colors = [gradient[i] for i in range(0, 1, length(keys(L_l_pairs)))]
q = plot(legend = false)
global count_color = 1
for (L, l_values) in L_l_pairs
    onset_data,_ = readdlm("zOnset_energies_$L.csv", ',', header = true)
    
    plot!(q, onset_data[:,1] ./ L, abs.(onset_data[:,2] - onset_data[:,3]), marker = :circle, label = "L = $L",color = colors[count_color])
    global count_color += 1
end
xlabel!(q, "l/L")
ylabel!(q, L"Energy spread, \$E_{on} - E_{min}\$")
savefig(q, "zfss_energy_spread_L.png")
    