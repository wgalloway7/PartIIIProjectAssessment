using DelimitedFiles
using Plots
using LsqFit
l_values = [1,3,7,9,21,27]
Tc = 2 / (log(1 + sqrt(2)))
q = plot()
T_VFT = []
for l in l_values
    data_tuple = readdlm("z_L=95_l=$(l)_unc_tau_fit.csv", ',', header = true)
    data = data_tuple[1]

    beta_values = data[:,1]
    tau_E_values = data[:,2]
    bc = log(1 + sqrt(2)) / 2

    mask = tau_E_values .> 1.0
    beta_values = beta_values[mask]
    tau_E_values = tau_E_values[mask]

    logtauE = log.(tau_E_values)
    T_values = 1 ./ beta_values

    function VFT(T, p)
        C, B, Tf = p
        return C .+ B ./ (T .- Tf)
    end

    guess = [1.0, 1.0, 1e-7]
    lower_bounds = [-5, 0.0, 1e-8]


    fit = curve_fit(VFT, T_values, logtauE, guess, lower = lower_bounds)

    params = fit.param
    Tf = params[3]
    println(params)
    push!(T_VFT, Tf)

    p = plot()
    plot!(p, T_values, logtauE, label="data")
    plot!(p, T_values, VFT(T_values, params), label="Fit $l")
    vline!(p, [Tf], label = "$l Tf")
    vline!(p, [Tc], label  ="Tc")
    savefig(p, "zVFT_connected_$(l).png")

    plot!(q, T_values, logtauE, label="$l")
    plot!(q, T_values, VFT(T_values, params), label="$l fit")
    vline!(q, [Tf], label = "$l Tf")
end
vline!(q, [Tc], label = "Tc")
savefig(q, "zVFT_all_l.png")


writedlm("VFT_temps.csv", T_VFT, ',')
