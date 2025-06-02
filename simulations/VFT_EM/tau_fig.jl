using DelimitedFiles
using Plots
using LaTeXStrings
using LsqFit

include("../../core/lattice.jl")

L = 95
l_values = [1,3,7,9,21,27]
Tc_inf = 2 / log(1 + sqrt(2))  # Critical temperature for the 2D Ising model
Tc = Tc_inf  + Tc_inf * (0.3603 / 95)

mid_red    = RGB(0.9, 0.4, 0.4)    
mid_blue   = RGB(0.2, 0.2, 0.8)    
gradient = cgrad([mid_blue, mid_red])
N = length(l_values)
colors = [gradient[i] for i in range(0, 1, length=N)]


data_0,_ = readdlm("1.csv", ',', header = true)
tau_0 = data_0[:, 2]

Tf_vals = []
p = plot()
for (i,l) in enumerate(l_values)
    data,_ = readdlm("$(l).csv", ',', header = true)
    beta = data[:, 1]
    tau = data[:, 2]
    alpha = data[:, 3]

    
    mask = tau .> 1.0

    logtau = log10.(tau[mask])
    T_values = 1 ./ beta[mask]

    function VFT(T, p)
        C, B, Tf = p
        return C .+ B ./ (T .- Tf)
    end

    guess = [1.0, 1.0, 1e-7]
    lower_bounds = [-5, 0.0, 1e-8]
    upper_bounds = [5.0, 5.0, 10.0]


    fit = curve_fit(VFT, T_values, logtau, guess, lower = lower_bounds, upper = upper_bounds)
    println(fit.param)

    if l != 1
        plot!(p, 1 ./ beta, log10.(tau), label = latexstring("l = $l"), color = colors[i])
        plot!(p, T_values, VFT(T_values, fit.param), color = colors[i], linestyle = :dash, label = "")
    else
        plot!(p, 1 ./ beta, log10.(tau), label = latexstring("l = $l"), color = :red, linestyle = :dot, linewidth = 2)  
    end


    
    if l != 1
        push!(Tf_vals, fit.param[3])
    end
end


xlabel!(p, L"T")
ylabel!(p, L"\log_{10}(\mathrm{Autocorrelation\ time}\ ), \log_{10}(\tau)")
vline!(p, [Tc], linestyle = :dash, color = :black, label = "")
annotate!(p, Tc, -0.5, text(L"T_c", :black, 10, :center))

savefig(p, "tau_vs_T_L=$(L)_connected.png")


q = plot()
scatter!(q, l_values[2:end], Tf_vals)
xlabel!(q, L"l")
ylabel!(q, L"T_f")
savefig(q, "Tf_vs_l_L=$(L)_connected.png")

