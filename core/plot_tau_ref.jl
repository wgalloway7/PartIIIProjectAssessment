using DelimitedFiles
using Plots
using LaTeXStrings
data = readdlm("average_tau_wide1.csv", ',')

beta_values = data[1, :]
tau_values = data[2, :]
p = plot()
plot!(p, 1 ./ beta_values, tau_values,label = "", color = :blue)
vline!(p, [2 / log(1 + sqrt(2))], label = L"$T_c$", color = :red, linestyle = :dash)
xlabel!(p, "Temperature")
ylabel!(p, L"Autocorrelation time $\tau$ (MC steps)")
savefig(p, "tau_SF.png")
