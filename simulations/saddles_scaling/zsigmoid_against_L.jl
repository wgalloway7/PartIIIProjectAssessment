using DelimitedFiles
using Plots
using Statistics
using LsqFit

include("../../core/lattice.jl")
using LaTeXStrings

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
L_l_pairs = sort(collect(L_l_pairs), by = x -> x[1])

l_choices = [3,5,7,9,11,13,15,19]

mid_red    = RGB(0.9, 0.4, 0.4)
mid_blue   = RGB(0.2, 0.2, 0.8)
gradient = cgrad([mid_blue, mid_red])
colors = [gradient[i] for i in range(0, 1, length(l_choices))]


R2_linear = Float64[]
R2_nonlinear = Float64[]

p = plot(legend = :topleft, size = (900, 600))
for (i,l_choice) in enumerate(l_choices)
    L_output = []
    x0_output = []
    for (L, l_values) in L_l_pairs
        for (j,l) in enumerate(l_values)
            if l == l_choice
                sigmoid_fit_data,header = readdlm("sigmoid_fit_k2_L=$(L).csv", '\t', header=true)
                x0_data = sigmoid_fit_data[:, 3]
                x0 = x0_data[j]
                push!(L_output, L)
                push!(x0_output, x0)
            end
        end
    end

    plot!(p, 1 ./ L_output, x0_output, label = "l = $l_choice", color = colors[i], marker = :circle)

    inv_L = 1.0 ./ L_output

    function linear(x,p)
        return p[1] .* x .+ p[2]
    end

    function nonlinear(x,p)
        return p[1] .+ p[2] .* x .^ p[3]
    end

    initial_guess = [1.0, 0.0]
    lin_fit = curve_fit(linear, inv_L, x0_output, initial_guess)

    non_lin_initial_guess = [-2.0, 1.0, 1.0]
    lower_bound = [-10.0, 0.0, 0.0]
    upper_bound = [0.0, 5.0, 5.0]
    non_lin_fit = curve_fit(nonlinear, inv_L, x0_output, non_lin_initial_guess, lower=lower_bound)

    y = x0_output
    ŷ_lin = linear(inv_L, lin_fit.param)
    ŷ_non = nonlinear(inv_L, non_lin_fit.param)
    sst = sum((y .- mean(y)).^2)
    sse_lin = sum((y .- ŷ_lin).^2)
    sse_non = sum((y .- ŷ_non).^2)
    push!(R2_linear, 1 - sse_lin / sst)
    push!(R2_nonlinear, 1 - sse_non / sst)

    x_fit = range(0, 0.12, length = 100)
    plot!(p, x_fit, linear(x_fit, lin_fit.param), label = "", color = colors[i]*0.5, linestyle = :dot, linewidth = 1.0)
    plot!(p, x_fit, nonlinear(x_fit, non_lin_fit.param), label = "", color = colors[i], linestyle = :dash, linewidth = 1.0)
end

xlims!(p, 0, 0.12)
ylims!(p, -2.2, 0.0)
hline!(p, [-2.0], label = L"E_0", color = :red, linewidth = 2.0, linestyle = :dash)
xlabel!(p, "1/L")
ylabel!(p, L"E_{sig}")
savefig(p, "zsigmoid_against_L_non_lin.png")


println("\nR² values for each l (order matches l_choices = $l_choices):")
println("Linear fit     R² = ", R2_linear)
println("Non-linear fit R² = ", R2_nonlinear)

bb       = bbox(0.7, 0.1, 0.3, 0.3)
xvals    = l_choices = [3,5,7,9,11,13,15,19]
yvals    = R2_linear

scatter!(p, [xvals,xvals], [R2_linear, R2_nonlinear],
    label = ["Linear" "Power law"],
    inset   = (1, bb),   
    subplot = 2,        
    legend  = true,
    xlabel  = "l",
    ylabel  = "R²",
    ylims = (0.95,1.05),
    markerstrokewidth = 0.5,
    xticks=false)

savefig(p, "zzinset.png")