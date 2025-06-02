using DelimitedFiles
using Statistics
using Plots
using Dates
using SpecialFunctions
using LaTeXStrings

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

Tc = 2 / log(1 + sqrt(2)) # = 2.269
println(Tc)


beta_values = 1 ./ generate_T_intervals(4.0, Tc, 20)   
copies       = 20                                      
L_values     = 95:97                             


L_l_pairs = Dict{Int, Vector{Int}}()
for L in L_values
    L_l_pairs[L] = [1]
    for l in 1:L
        if valid_l(L, l) && isodd(l) && l != 1
            push!(L_l_pairs[L], l)
        end
    end
end
L_l_pairs = sort(collect(L_l_pairs), by = x -> x[1])


onsager_mag(beta) = (1 - sinh(2β)^(-4))^(1/8)

for (L, l_values) in L_l_pairs



    for l in l_values
        println("L = $(L), l = $(l)")

        unc_all = Dict(i => Vector{Vector{Float64}}() for (i, _) in enumerate(beta_values))


        for j in 1:copies
            ac = readdlm("L=$(L)_l=$(l)_copy$(j)_autocorrelation_unconnected.csv", ',')
            for (i, _) in enumerate(beta_values)
                push!(unc_all[i], ac[i, :])
            end
        end


        for i in eachindex(beta_values)
            for v in unc_all[i]
                v .= (v .- 0.5) .* 2
            end
        end


        unc0 = zeros(length(beta_values), length(unc_all[1][1]))
        for (i, rows) in unc_all
            for t in 1:size(unc0, 2)
                vals = filter(!isnan, map(r -> r[t], rows))
                unc0[i, t] = isempty(vals) ? NaN : mean(vals)
            end
        end


        unc_tau_fit_values   = zeros(length(beta_values))
        unc_alpha_fit_values = zeros(length(beta_values))

        for (i, _) in enumerate(beta_values)
            unc_row            = unc0[i, :]
            τ, α               = fit_stretched_exponential_from_autocorrelation_function(unc_row)
            unc_tau_fit_values[i]   = τ
            unc_alpha_fit_values[i] = α
        end

        output_file = "z_L=$(L)_l=$(l)_unc_tau_fit.csv"
        open(output_file, "w") do io
            writedlm(io, [["beta" "tau_unc" "alpha_unc"]], ',')
            for i in 1:length(beta_values)
                writedlm(io, [[beta_values[i], unc_tau_fit_values[i], unc_alpha_fit_values[i]]], ',')
            end
        end

        q4 = plot(1 ./ beta_values, log10.(unc_tau_fit_values),
                  label = "", color = :blue)
        xlabel!(q4, L"T")
        ylabel!(q4, L"\log_{10} \\tau")
        annotate!(q4, (Tc, 0), text(L"T_c", 12))
        savefig(q4, "z_L=$(L)_l=$(l)_unc_tau_fit.png")

        q5 = plot(beta_values, unc_alpha_fit_values, label = "α_unc", color = :blue)
        savefig(q5, "z_L=$(L)_l=$(l)_unc_alpha_fit.png")


        mid_red   = RGB(0.8, 0.2, 0.2)
        mid_blue  = RGB(0.2, 0.2, 0.8)
        gradient  = cgrad([mid_blue, mid_red])
        N         = length(beta_values)
        colors    = [gradient[i] for i in range(0, 1, length = N)]


        stretched_exponential(t, τ, α) = exp(-((t / τ)^α))


        q1 = plot(legend = :right, size = (900, 600))
        xlabel!(q1, "Time, MC steps")
        ylabel!(q1, "Autocorrelation")


        for i in 1:length(beta_values)
            plot!(q1, unc0[i, :],
                  label = (i == 1 || i == length(beta_values)) ?
                          latexstring("T = $(round(1 / beta_values[i], digits = 2))") : "",
                  color = colors[i])

            try
                plot!(q1,
                      stretched_exponential.(1:length(unc0[i, :]),
                                              unc_tau_fit_values[i],
                                              unc_alpha_fit_values[i]),
                      linestyle = :dash, color = colors[i], label = "")
            catch e
                @warn "Error plotting fit for β=$(beta_values[i]): $e"
            end
        end


        bb       = bbox(0.7, 0.0, 0.3, 0.3)  # **top‑right**
        xvals    = 1.0 ./ beta_values
        τvals    = unc_alpha_fit_values


        plot!(q1, xvals, τvals,
              inset   = (1, bb),   
              subplot = 2,         
              xlabel  = "T",
              ylabel  = L"\alpha",
              ylims = (0,1.0),
              
              ms      = 3)

        savefig(q1, "z_L=$(L)_l=$(l)_unc_autocorrelation.png")
    end
end
