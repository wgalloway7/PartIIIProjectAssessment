using DelimitedFiles
using Plots
using LsqFit

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")


function fit_sigmoid(l::Int; rev::Bool=false, filename::String="sigmoid_fit", L::Int=95)
    data = readdlm("proportions_l=$(l)_saddles_simulation_$(L).csv")
    E    = Float64.(data[2:end, 1])
    k2   = Float64.(data[2:end, rev ? 2 : 4])
    fwd  = rev ? -1.0 : 1.0

    sigmoid(x, p) = @. 1 / (1 + exp(-fwd * p[1] * (x - p[2])))

    idx = sortperm(E)
    E = E[idx]
    k2 = k2[idx]


    if (rev && all(k2 .== 0.0)) || (!rev && all(k2 .== 1.0))
        return NaN, NaN
    end

    Δ = 5.0
    Eshift = E .+ Δ
    p0 = [2.0, -2.0 + Δ]
    lower = [1.0, -2.0 + Δ]
    upper = [50.0, 0.0 + Δ]

    fit = curve_fit(sigmoid, Eshift, k2, p0; lower=lower, upper=upper)
    k, x0 = fit.param
    x0 -= Δ  


    p = plot(E, k2, label="data")
    plot!(p, E, sigmoid(E, [k, x0]), label="fit", color=:red)
    savefig(p, "$(filename)_$l.png")

    return k, x0
end


L_values = [10,20,23,25,46,47,50,52]
L_l_pairs = Dict{Int64, Vector{Int64}}()
for L in L_values
    for l in 1:L
        if valid_l(L, l) && l % 2 == 1 && l != 1
            push!(get!(L_l_pairs, L, Int64[]), l)
        end
    end
end


for (L, l_values) in sort(collect(L_l_pairs), by = x -> x[1])
    println("Processing L = $L")


    fwd_file = "sigmoid_fit_k2_$L.csv"
    rev_file = "sigmoid_fit_k0_$L.csv"
    writedlm(fwd_file, [["l" "x0" "k"]], ',')
    writedlm(rev_file, [["l" "x0" "k"]], ',')


    for l in l_values
        println("  Fitting l = $l")


        k_rev, x0_rev = fit_sigmoid(l; rev=true, filename="sigmoid_fit_k0_$L", L=L)
        open(rev_file, "a") do io
            writedlm(io, [l x0_rev k_rev], ',')
        end

s
        k_fwd, x0_fwd = fit_sigmoid(l; rev=false, filename="sigmoid_fit_k2_$L", L=L)
        open(fwd_file, "a") do io
            writedlm(io, [l x0_fwd k_fwd], ',')
        end
    end
end
