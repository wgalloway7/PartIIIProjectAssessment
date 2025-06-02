using DelimitedFiles
using Plots
using LaTeXStrings
using LsqFit
using Statistics

Tc = 2 / log(1 + sqrt(2))  # Critical temperature for the 2D Ising model
function VFT(T, p)
    C, B, Tf = p
    return C .+ B ./ (T .- Tf)
end
l = 1
data_1,_= readdlm("z_L=95_l=$(l)_unc_tau_fit.csv", ',', header=true)
plot3 = plot()
T_values = 1 ./ data_1[:, 1]
tau_values = data_1[:, 2]
inv_red_t = 1 ./ (T_values .- Tc)

mask = tau_values .> 1.0
tau_values = tau_values[mask]
inv_red_t = inv_red_t[mask]
T_values = T_values[mask]

logt = log10.(inv_red_t[1:end-1])
logtau = log10.(tau_values[1:end-1])
println(size(logt), size(logtau))
p = plot()
function model(x,p)
    a,b = p 
    return a .+ b .* x
end
initial_params = [0.0, 1.0]
fit = curve_fit(model, logt, logtau, initial_params)

best_fit_params = coef(fit)
println("Goodness of fit (R²): ", 1 - sum((logtau .- model(logt, best_fit_params)).^2) / sum((logtau .- mean(logtau)).^2))

plot!(p, logt, model(logt, best_fit_params), label="Best fit line")


scatter!(p, logt, logtau)
savefig(p, "CSD.png")




function VFT(T, p)
    C, B, Tf = p
    return C .+ B ./ (T .- Tf)
end

initial_params = [0.0, 1.0, 1.0]
logtau_values = log10.(tau_values)
fit = curve_fit(VFT, T_values, logtau_values, initial_params)
println(fit.param)

q = plot()
T_array = range(minimum(T_values), stop=maximum(T_values), length=100)
plot!(q, T_array, VFT(T_array, fit.param), label="VFT Fit", color=:red)
scatter!(q, T_values, logtau_values, label="Data", color=:blue)
savefig(q, "VFT.png")

# Compute goodness of fit (R²) manually for the VFT fit
predicted_logtau = VFT(T_values, fit.param)
residuals = logtau_values .- predicted_logtau
total_variation = logtau_values .- mean(logtau_values)
r_squared = 1 - sum(residuals.^2) / sum(total_variation.^2)

println("Goodness of fit (R²) for VFT fit: ", r_squared)
