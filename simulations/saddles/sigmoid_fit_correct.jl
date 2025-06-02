using DelimitedFiles
using Plots
using LsqFit

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

function fit_sigmoid(filename::String="sigmoid_fit",l::Int64=1, L::Int=95)
    data = readdlm("proportions_l=$(l)_saddles_simulation_$(L).csv")
    E    = Float64.(data[2:end, 1])
    k2   = Float64.(data[2:end, 4])

    sigmoid(x,p) = @. 1  - (1 + exp(p[1] * (x - p[2])))^(-1)

    idx = sortperm(E)
    E = E[idx]
    k2 = k2[idx]

    if (all(k2 .== 0.0)) || (all(k2 .== 1.0))
        return NaN, NaN
    end
    Δ = 0.0
    Eshift = E .+ Δ
    p0 = [5.0, -2.0 + Δ, 2.0]
    lower = [1e-8, -2.0 + Δ, 1e-8]
    upper = [1000.0, 0.0 + Δ, 10.0]

    fit = curve_fit(sigmoid, Eshift, k2, p0; lower=lower, upper=upper)
    #fit = curve_fit(sigmoid, Eshift, k2, p0)
    alpha, x0, n = fit.param
    x0 -= Δ  # shift back

    p = plot(E, k2, label="data")
    plot!(p, E, sigmoid(E, [alpha, x0, n]), label="fit", color=:red)
    savefig(p, "$(filename)_$(L)_$(l).png")
    return alpha, x0, n
end




L_values = 95:105
L_l_pairs = Dict{Int64, Vector{Int64}}()
for L in L_values
    L_l_pairs[L] = [1]
    for l in 1:L
        if valid_l(L, l) && l % 2 == 1 && l != 1
            push!(get!(L_l_pairs, L, Int64[]), l)
        end
    end
end

for (L, l_values) in sort(collect(L_l_pairs), by = x -> x[1])
    println("Processing L = $L")

    writedlm("sigmoid_fit_k2_L=$(L).csv", [["l", "alpha", "x0", "n"]])
    for l in l_values
        alpha, x0, n = fit_sigmoid("sigmoid_fit", l, L)
        println("L = $L, l = $l, alpha = $alpha, x0 = $x0, n = $n")

        open("sigmoid_fit_k2_L=$(L).csv", "a") do file
            writedlm(file, [[l, alpha, x0, n]])
        end
    end
end

