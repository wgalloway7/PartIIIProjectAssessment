using Plots
using DelimitedFiles
using LaTeXStrings


include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")
#include("../../experiments/singleflip/singleflip_fit.jl")


data = readdlm("saddle_critical_energies.csv")
println(size(data))
E_lower_raw = data[2:end, 2]
E_upper_raw = data[2:end, 3]
ls = data[2:end, 1]

E_lower = Float64[]
for val in E_lower_raw
    if val == "nothing"
        push!(E_lower, NaN)
    elseif isa(val, Float64)
        push!(E_lower, val)
    else
        push!(E_lower, parse(Float64, val))
    end
end

E_upper = Float64[]
for val in E_upper_raw
    if val == "nothing"
        push!(E_upper, NaN)
    elseif isa(val, Float64)
        push!(E_upper, val)
    else
        push!(E_upper, parse(Float64, val))
    end
end




beta_values = 1 ./ generate_T_intervals(4.0, 1.0, 100)
global p_values = [16]
l_values = [2*i - 1 for i in 1:50]
q = plot()
plot2 = plot()



#colors = (cgrad(:RdBu, length(p_values)))
mid_red    = RGB(0.9, 0.4, 0.4)    # warm, mid-dark red
mid_blue   = RGB(0.2, 0.2, 0.8)    # cool, mid-dark blue
gradient = cgrad([mid_blue, mid_red])
N = length(l_values)
colors = [gradient[i] for i in range(0, 1, length=N)]

q = plot()
plot2 = plot()
plateaux_energies = []
plat_tc_energies = []

for (j,l) in enumerate(l_values)
    if l in [17,67]
        for (i,p) in enumerate(p_values)
            d = readdlm("l=$(l)_p=$(p)_energy_p_$(l)_.csv")
            d = reshape(d, length(d))
            if l in [1,l_values[end]]
                plot!(q, beta_values, d, label = "l = $l", color = colors[j])
                plot!(plot2, beta_values, d - onsager_energy.(beta_values, 100), label = "l = $l", color = colors[j])
            else
                plot!(q, beta_values, d, label = "", color = colors[j])
                plot!(plot2, beta_values, d - onsager_energy.(beta_values, 100), label = "", color = colors[j])
            end
        end
        E_l = E_lower[l]
        E_u = E_upper[l]
        println(ls[l], E_l, E_u)
        hline!(q, [E_l], label = "E_lower", color = colors[j], linestyle = :dash)
        hline!(q, [E_u], label = "E_upper", color = colors[j], linestyle = :dot)

    end

end

plot!(q, beta_values, onsager_energy.(beta_values, 100), color = :black, label  = "")
vline!(q,[log(1+sqrt(2))/2], color = :black, linestyle = :dash, label = L"T_c")
xlabel!(q, L"\beta = 1/T")
ylabel!(q, L"E")
plot!(q, colorbar=true, cbar_ticks=1:length(colors), cbar_label="Color Bar")
savefig(q, "energy_vs_beta_topo.png")

