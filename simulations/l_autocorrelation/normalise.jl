using DelimitedFiles
using Statistics
using Plots
using Dates
using LaTeXStrings

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")
Tc = 2/log(1+sqrt(2)) # = 2.269
println(Tc)
beta_values = 1 ./ generate_T_intervals(4.0,Tc,20)
copies = 20
L_values = [i for i in 95:98]
L_l_pairs = Dict{Int64, Vector{Int64}}() 
for L in L_values
    L_l_pairs[L] = [1]
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
        if 1 == 1
            # Load all autocorrelation data first
            E_all = Dict{Int, Vector{Vector{Float64}}}()
            M_all = Dict{Int, Vector{Vector{Float64}}}()

            for (i, beta) in enumerate(beta_values)
                E_all[i] = Vector{Vector{Float64}}()
                M_all[i] = Vector{Vector{Float64}}()
            end
            println("averaging")
            for j in 1:copies
                E_ac = readdlm("L=$(L)_l=$(l)_copy$(j)_autocorrelation_E.csv", ',')
                M_ac = readdlm("L=$(L)_l=$(l)_copy$(j)_autocorrelation_m.csv", ',')
                #normalise the autocorrelation functions
                for (i, beta) in enumerate(beta_values)
                    E_row = E_ac[i, :]
                    M_row = M_ac[i, :]

                    E_row = E_row ./ E_row[1]
                    M_row = M_row ./ M_row[1]

                    push!(E_all[i], E_row)
                    push!(M_all[i], M_row)
                end
            end

            E0 = zeros(length(beta_values), length(E_all[1][1]))
            M0 = zeros(length(beta_values), length(M_all[1][1]))
            # compute average, binning the NaN values
            for i in 1:length(beta_values)
                for t in 1:size(E0, 2)
                    E_values = [row[t] for row in E_all[i] if !isnan(row[t])]
                    M_values = [row[t] for row in M_all[i] if !isnan(row[t])]

                    E0[i, t] = isempty(E_values) ? NaN : mean(E_values)
                    M0[i, t] = isempty(M_values) ? NaN : mean(M_values)
                end
            end

            tau_E_values = zeros(length(beta_values))
            tau_M_values = zeros(length(beta_values))
            for (i,beta) in enumerate(beta_values)
                E_row = E0[i, :]
                M_row = M0[i, :]

                e_index_E = findfirst(x -> x < 1 / exp(1), E_row)
                e_index_M = findfirst(x -> x < 1 / exp(1), M_row)

                tau_E = isnothing(e_index_E) ? (length(E_row)) : e_index_E
                tau_M = isnothing(e_index_M) ? (length(M_row)) : e_index_M

                tau_E_values[i] = tau_E
                tau_M_values[i] = tau_M
            end
            q3 = plot()
            plot!(q3, beta_values, tau_E_values, label = "tau_E", color = "blue")
            plot!(q3, beta_values, tau_M_values, label = "tau_M", color = "red")
            savefig(q3, "a_L=$(L)_l=$(l)_tau_exp-1.png")
            println("Fitting")

            E_tau_fit_values = zeros(length(beta_values))
            M_tau_fit_values = zeros(length(beta_values))
            E_alpha_fit_values  = zeros(length(beta_values))
            M_alpha_fit_values  = zeros(length(beta_values))
            for (i, beta) in enumerate(beta_values)
                E_row = E0[i, :]
                M_row = M0[i, :]
                
                E_fit_params = fit_stretched_exponential_from_autocorrelation_function(E_row)
                M_fit_params = fit_stretched_exponential_from_autocorrelation_function(M_row)
                E_tau_fit_values[i] = E_fit_params[1]
                M_tau_fit_values[i] = M_fit_params[1]
                E_alpha_fit_values[i] = E_fit_params[2]
                M_alpha_fit_values[i] = M_fit_params[2]
            
            end


            output_file = "a_L=$(L)_l=$(l)_fit_results.csv"
            open(output_file, "w") do io
                writedlm(io, [["beta" "tau_E" "alpha_E" "tau_M" "alpha_M"]], ',')
                    for i in 1:length(beta_values)
                        writedlm(io, [[beta_values[i], E_tau_fit_values[i], E_alpha_fit_values[i], M_tau_fit_values[i], M_alpha_fit_values[i]]], ',')
                    end
            end

            q4 = plot()
            mask  = beta_values .< 0.5
            cropped_beta_values = beta_values[mask]
            E_tau_fit_cropped =E_tau_fit_values[mask]
            M_tau_fit_cropped =M_tau_fit_values[mask]
            plot!(q4, 1 ./ cropped_beta_values, log10.(E_tau_fit_cropped), label = "tau_E", color = "blue")
            plot!(q4, 1 ./ cropped_beta_values, log10.(M_tau_fit_cropped), label = "tau_M", color = "red")
            xlabel!(q4, "T")
            ylabel!(q4, "log10(tau)")
            savefig(q4, "a_L=$(L)_l=$(l)_tau_fit.png")

            


            q5 = plot()
            plot!(q5, beta_values, E_alpha_fit_values, label = "alpha_E", color = "blue")
            plot!(q5, beta_values, M_alpha_fit_values, label = "alpha_M", color = "red")
            savefig(q5, "a_L=$(L)_l=$(l)_alpha_fit.png")



            mid_red    = RGB(0.9, 0.4, 0.4)    # warm, mid-dark red
            mid_blue   = RGB(0.2, 0.2, 0.8)    # cool, mid-dark blue
            gradient = cgrad([mid_blue, mid_red])
            N = length(beta_values)
            colors = [gradient[i] for i in range(0, 1, length=N)]

            function stretched_exponential(x, tau, alpha)
                return exp(-((x / tau) ^ alpha))
            end
            q1 = plot(legend = :right, size = (900, 600))
            xlabel!(q1, "Time, MC steps")
            ylabel!(q1, "Autocorrelation")

            for i in 1:length(beta_values)
                if 1 == 1
                    if (i == 1) || (i == length(beta_values))
                        plot!(q1, E0[i, :], label = latexstring("T = $(round(1 / beta_values[i], digits=2))"), color = colors[i])
                    else
                        plot!(q1, E0[i, :], label = "", color = colors[i])
                    end
                    try
                        plot!(q1, [i for i in 0:(length(E0[i, :])-1)], stretched_exponential.([i for i in 0:(length(E0[i, :])-1)], E_tau_fit_values[i], E_alpha_fit_values[i]), label = "", color = colors[i], linestyle = :dash)
                    catch e
                        println("Error in plotting fit for beta = $(beta_values[i]): $e")
                        println("tau = $(E_tau_fit_values[i]), alpha = $(E_alpha_fit_values[i])")
                    end
                end 
            end

            bb       = bbox(0.7, 0.0, 0.3, 0.3)
            xvals    = 1.0 ./ beta_values
            τvals    = E_alpha_fit_values

            plot!(q1, xvals, τvals,
              inset   = (1, bb),   
              subplot = 2,         
              legend  = false,
              xlabel  = "T",
              ylabel  = L"\alpha",
              ylims = (0,1.0),
              
              ms      = 3)

            savefig(q1, "a_L=$(L)_l=$(l)_average_autocorrelation_E.png")
            q2 = plot(legend = false)
            for i in 1:length(beta_values)
                if 1==1
                    if (i == 1) || (i == length(beta_values))
                        plot!(q2, M0[i, :], label = latexstring("T = $(round(1 / beta_values[i], digits=2))"), color = colors[i])
                    else
                        plot!(q2, M0[i, :], label = "", color = colors[i])
                    end

                    try
                        plot!(q2, [i for i in 0:(length(M0[i, :])-1)], stretched_exponential.([i for i in 0:(length(M0[i, :])-1)], M_tau_fit_values[i], M_alpha_fit_values[i]), label = "fit", color = colors[i], linestyle = :dash)
                    catch e
                        println("Error in plotting fit for beta = $(beta_values[i]): $e")
                        println("tau = $(M_tau_fit_values[i]), alpha = $(M_alpha_fit_values[i])")
                    end
                end
            end
            savefig(q2, "a_L=$(L)_l=$(l)_average_autocorrelation_m.png")
        end
    end
end