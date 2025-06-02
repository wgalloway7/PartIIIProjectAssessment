using DelimitedFiles
using Plots
using LaTeXStrings
using LsqFit
using Statistics

Tc_inf = 2 / log(1 + sqrt(2))  # Critical temperature for the 2D Ising model
Tc = Tc_inf  + Tc_inf * (0.3603 / 95)
function VFT(T, p)
    C, B, Tf = p
    return C .+ B ./ (T .- Tf)
end
VFT_r = []
CSD_r = []
l_values = [1, 3, 7, 9, 21, 27]
for l in l_values
    data_1,_= readdlm("z_L=95_l=$(l)_unc_tau_fit.csv", ',', header=true)
    plot3 = plot()
    T_values = 1 ./ data_1[1:end-1, 1]
    
    tau_values = data_1[1:end-1, 2]

    mask = tau_values .> 1.0
    tau_values = tau_values[mask]
    T_values = T_values[mask]
    println("minmax $(minimum(T_values)), $(maximum(T_values))")

    logtau = log10.(tau_values)
    p = plot()
    function model(x,p)
        a,b = p 
        return a .+ b .* x
    end

    function CSD(T,p)
        C,v,Tc_fit = p
        return C .* ((T .- Tcrit) .^ (-2))
    end
    initial_params = [0.0, 1.0,2.2646]
    upper_bounds = [1e8, 10.0, Tcrit]
    lower_bounds = [-1e8, 1e-8, 2.2646]
    #fit = curve_fit(model, logt, logtau, initial_params)
    fit = curve_fit(CSD, T_values, tau_values, initial_params, upper = upper_bounds, lower = lower_bounds)
    fit_params = fit.param
    println("CSD fit parameters for l=$(l): ", fit_params)
    

    predicted_tau = CSD(T_values, fit_params)
    residuals = tau_values .- predicted_tau
    total_variation = tau_values .- mean(tau_values)
    r_squared = 1 - sum(residuals.^2) / sum(total_variation.^2)
    push!(CSD_r, r_squared)
    
    T_array = range(minimum(T_values), stop=maximum(T_values), length=100)
    plot!(p, T_array, CSD(T_array, fit_params), label="Best fit line")
    xlabel!("Temperature (T)")
    ylabel!("Tau")


    scatter!(p, T_values, tau_values)
    savefig(p, "CSD_$(l).png")




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
    xlabel!(q, "Temperature (T)")
    ylabel!(q, "log10(tau)")
    scatter!(q, T_values, logtau_values, label="Data", color=:blue)
    savefig(q, "VFT_$(l).png")


    predicted_logtau = VFT(T_values, fit.param)
    residuals = logtau_values .- predicted_logtau
    total_variation = logtau_values .- mean(logtau_values)
    r_squared = 1 - sum(residuals.^2) / sum(total_variation.^2)
    push!(VFT_r, r_squared)

    println("Goodness of fit (R²) for VFT fit: ", r_squared)
end

println("CSD R² values: ", CSD_r)
println("VFT R² values: ", VFT_r)

q = plot()
scatter!(q, l_values, CSD_r, label="CSD R²", color=:blue)
scatter!(q, l_values, VFT_r, label="VFT R²", color=:red)
xlabel!(q, "l")
ylabel!(q, "R²")
savefig(q, "R_squared_vs_l.png")