using DelimitedFiles
using Plots
using LaTeXStrings
using LsqFit

include("../../core/lattice.jl")

L = 95
l_values = [1, 3, 7, 9, 21, 27]
Tc = 2 / log(1 + sqrt(2))  # Critical temperature for the 2D Ising model

mid_red  = RGB(0.9, 0.4, 0.4)    # warm, mid-dark red
mid_blue = RGB(0.2, 0.2, 0.8)    # cool, mid-dark blue
gradient = cgrad([mid_blue, mid_red])
N = length(l_values)
colors = [gradient[i] for i in range(0, 1, length=N)]

data_0, _ = readdlm("z_L=95_l=1_unc_tau_fit.csv", ',', header=true)
tau_0 = data_0[:, 2]

Tf_vals = []
VFT_fit_vals = [Float64[], Float64[], Float64[], Int64[]]  # C, B, Tf
CSD_fit_vals = [Float64[], Float64[], Int64[]]             # C, v
VFT_fit_vals[4] = l_values
CSD_fit_vals[3] = l_values
p = plot(legend=:right)
for (i, l) in enumerate(l_values)
    data, _ = readdlm("z_L=95_l=$(l)_unc_tau_fit.csv", ',', header=true)
    beta = data[:, 1]
    tau = data[:, 2]
    alpha = data[:, 3]

    mask = tau .> 1.0
    tau_mask = tau[mask]
    logtau = log10.(tau_mask)
    T_values = 1 ./ beta[mask]

    # VFT function
    function VFT(T, p)
        C, B, Tf = p
        return C .+ B ./ (T .- Tf)
    end

    guess = [1.0, 1.0, 1e-7]
    lower_bounds = [-5, 0.0, 1e-8]
    upper_bounds = [2.0, 5.0, Tc]

    fit = curve_fit(VFT, T_values, logtau, guess, lower=lower_bounds, upper=upper_bounds)
    println("l=$l VFT fit params: ", fit.param)

    # Critical slowing down function (only fit for T > Tc)
    function crit_slowing(T, p)
        C, v= p
        return C .* ((T .- Tc) .^ (-v))
    end

    valid = T_values .> Tc + 1e-6  # exclude values ≤ Tc for this fit
    T_fit = T_values[valid]
    tau_fit = tau_mask[valid]
    println(tau_mask)
    println(T_fit)

    CSD_guess = [1.0, 1.0, 0.1]
    CSD_lower_bounds = [1e-8, 1e-8, 0.0]
    CSD_upper_bounds = [1e8, 20.0,2.0]

    CSD_fit = curve_fit(crit_slowing, T_fit, tau_fit, CSD_guess, lower=CSD_lower_bounds, upper=CSD_upper_bounds)
    println("l=$l CSD fit params: ", CSD_fit.param)

    # Plotting
    if l != 1
        plot!(p, 1 ./ beta, log10.(tau), label=latexstring("l = $l"), color=colors[i])
        plot!(p, T_values, VFT(T_values, fit.param), color=colors[i], linestyle=:dash, label="")
    else
        # Plot only data points for l=1, no fit curve
        plot!(p, 1 ./ beta, log10.(tau), label=latexstring("l = $l"), color=:red, linestyle=:dot, linewidth=2)
    end

    # Collect fit params
    push!(Tf_vals, fit.param[3])
    push!(VFT_fit_vals[1], fit.param[1])
    push!(VFT_fit_vals[2], fit.param[2])
    push!(VFT_fit_vals[3], fit.param[3])
    push!(CSD_fit_vals[1], CSD_fit.param[1])
    push!(CSD_fit_vals[2], CSD_fit.param[2])
end

open("VFT_fit_L=$(L).csv", "w") do io
    println(io, "C,B,Tf,l")
    writedlm(io, hcat(VFT_fit_vals...), ',')
end

open("CSD_fit_L=$(L).csv", "w") do io
    println(io, "C,v,l")
    writedlm(io, hcat(CSD_fit_vals...), ',')
end


xlabel!(p, L"T")
ylabel!(p, L"\log_{10}(\mathrm{Autocorrelation\ time}\ ), \log_{10}(\tau)")
vline!(p, [Tc], linestyle=:dash, color=:black, label="")
annotate!(p, Tc, -0.5, text(L"T_c", :black, 10, :center))

savefig(p, "tau_vs_T_L=$L.png")

using Plots.PlotMeasures

# Create inset plot of Tf vs l (excluding l=1)
k = plot()
scatter!(k, l_values[1:end], Tf_vals[1:end], label="", color=:blue)
xlabel!(k, L"l")
ylabel!(k, L"T_f")
xlims!(k, minimum(l_values[1:end])-1, maximum(l_values[2:end]) + 1)
ylims!(k, minimum(Tf_vals[1:end]) - 0.1, maximum(Tf_vals[1:end]) + 0.1)
savefig(k, "Tf_vs_l_L=$L.png")

bb = bbox(0.4, 0.05, 0.35, 0.35)


scatter!(
    p,
    l_values[l_values .!= 1],
    Tf_vals[l_values .!= 1],
    inset   = (1, bb),
    subplot = 2,
    legend  = false,
    xlabel  = "l",
    ylabel  = L"T_f",
    ms      = 4,
    xlims   = (minimum(l_values)-1, maximum(l_values) + 1),
    ylims   = (minimum(Tf_vals) - 0.1, maximum(Tf_vals) + 0.1),
    marker  = :diamond,
    color   = :blue,
)

# Then, plot the point where l = 1 in red
scatter!(
    p,
    l_values[l_values .== 1],
    Tf_vals[l_values .== 1],
    subplot = 2,
    marker  = :diamond,
    ms      = 4,
    color   = :red,
)
savefig(p, "tau_vs_T_L=$(L)_with_inset.png")
