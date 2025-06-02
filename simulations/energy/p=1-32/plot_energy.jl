using DelimitedFiles
using Plots
using Statistics
using SpecialFunctions
using LaTeXStrings


include("../../../core/lattice.jl")
include("../../../core/montecarlo.jl")
include("../../../core/main.jl")
#include("../../experiments/singleflip/singleflip_fit.jl")

beta_values = 1 ./ generate_T_intervals(4.0, 1.0, 100)
global p_values = [i for i in 16:32]

q = plot()
plot2 = plot()

l = 1
mid_red    = RGB(0.9, 0.4, 0.4)    # warm, mid-dark red
mid_blue   = RGB(0.2, 0.2, 0.8)    # cool, mid-dark blue
gradient = cgrad([mid_red, mid_blue])
N = length(p_values)
colors = [gradient[i] for i in range(0, 1, length=N)]

q = plot()
plot2 = plot()
for (i,p) in enumerate(p_values)
    d = readdlm("k=$(l)_p=$(p)_energy_p.csv")
    d = reshape(d, length(d))
    if p in [16,32]
        plot!(q, beta_values, d, label = "p = $p", color = colors[i])
        plot!(plot2, beta_values, d - onsager_energy.(beta_values, 100), label = "p = $p", color = colors[i])
    else
        plot!(q, beta_values, d, label = "", color = colors[i])
        plot!(plot2, beta_values, d - onsager_energy.(beta_values, 100), label = "", color = colors[i])
    end
end
#plot!(q, beta_values, onsager_energy.(beta_values, 100), color = :springgreen, label  = "Analytical Result", linestyle = :dot, linewidth = 3)
vline!(q,[log(1+sqrt(2))/2], color = :black, linestyle = :dash, label = L"T_c")
plot!(q)
xlabel!(q, L"\beta = 1/T")
ylabel!(q, "E")
plot!(q, colorbar=true, cbar_ticks=1:length(colors), cbar_label="Color Bar")
savefig(q, "E_beta.png")

vline!(plot2, [log(1+sqrt(2))/2], color = :blue, linestyle = :dash, label = "Tc")
plot!(plot2, background_color = :gray)
xlabel!(plot2, "Beta")
ylabel!(plot2, "log(onsager_energy - energy)")
plot!(plot2, colorbar=true, cbar_ticks=1:length(colors), cbar_label="Color Bar")
savefig(plot2, "E_vs_beta_log.png")

# Create a standalone colorbar figure
colorbar_plot = plot(legend=false, ticks=nothing, border=:none, size=(200, 400))

for (i, p) in enumerate(p_values)
    scatter!(colorbar_plot, [1], [i],
        markersize=15,
        color=colors[i],
        label="p = $p",
        markerstrokecolor=:black)
end

yticks!(colorbar_plot, 1:length(p_values), ["p = $p" for p in p_values])
xlims!(colorbar_plot, 0.5, 1.5)
ylims!(colorbar_plot, 0, length(p_values) + 1)
xaxis!(colorbar_plot, false)
title!(colorbar_plot, "Color Key")

savefig(colorbar_plot, "colorbar_legend.png")
