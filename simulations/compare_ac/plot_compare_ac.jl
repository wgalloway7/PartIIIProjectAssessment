using DelimitedFiles
using Plots
using LaTeXStrings
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")


L = 100
res = 100
beta_values = 1 ./ generate_T_intervals(4.0, 0.3, res)
E_matrix = readdlm("E_1_compare_ac_new.csv", ',')
M_matrix = readdlm("M_1_compare_ac_new.csv", ',')
S_matrix = readdlm("S_1_compare_ac_new.csv", ',')

E0 = zeros(size(E_matrix))
M0 = zeros(size(M_matrix))
S0 = zeros(size(S_matrix))


copies_list = [i for i in 1:12]
for i in copies_list
    #println(i)
    #println(size(readdlm("E_$(i)_compare_ac_new.csv", ',')))
    #println(size(readdlm("M_$(i)_compare_ac_new.csv", ',')))
    #println(size(readdlm("S_$(i)_compare_ac_new.csv", ',')))
    E0 .+= readdlm("E_$(i)_compare_ac_new.csv", ',')
    M0 .+= readdlm("M_$(i)_compare_ac_new.csv", ',')
    S0 .+= readdlm("S_$(i)_compare_ac_new.csv", ',')
end
E0 = E0 ./ length(copies_list)
M0 = M0 ./ length(copies_list)
S0 = S0 ./ length(copies_list)

p1 = plot()
p2 = plot()
p3 = plot()

E_tau_values = []
M_tau_values = []
S_tau_values = []
for (i,b) in enumerate(beta_values)

    plot!(p1, E0[i,:], label = round(b, digits = 2))
    plot!(p2, M0[i,:], label = round(b, digits = 2))
    plot!(p3, S0[i,:], label = round(b, digits = 2))
    Mrow = M0[i,:]
    Srow = S0[i,:]
    Erow = E0[i,:]

    tau = length(Mrow)
    for (j,a) in enumerate(Mrow)
        if a < exp(-1)
            tau = j
            break
        end
    end
    push!(M_tau_values, tau)

    tau = length(Srow)
    for (j,a) in enumerate(Srow)
        if a < exp(-1)
            tau = j
            break
        end
    end
    push!(S_tau_values, tau)
    tau = length(Erow)

    for (j,a) in enumerate(Erow)
        if a < exp(-1)
            tau = j
            break
        end
    end
    push!(E_tau_values, tau)
    
end
savefig(p1, "E_ac.png")
savefig(p2, "M_ac.png")
savefig(p3, "S_ac.png")

q = plot()
println(E_tau_values)
println(M_tau_values)
println(S_tau_values)
plot!(q, 1 ./ beta_values, log10.(E_tau_values), label = L"\tau_E", linestyle = :solid)
plot!(q, 1 ./ beta_values, log10.(M_tau_values), label = L"\tau_m", linestyle  =:dash, color = :red)
plot!(q, 1 ./ beta_values, log10.(S_tau_values), label = L"\tau_{\sigma}", linestyle = :dot, color = :green)
vline!(q, [2/(log(1 + sqrt(2)))], label = L"T_c")

beta_ac = log(L^2) / 8
vline!(q, [1 / beta_ac], label = L"T_{\textrm{ordered}}", color = :black)

xlabel!(q, L"T")
ylabel!(q, L"log(Autocorrelation time) $\log(\tau)$ (MC steps)")

savefig(q, "plot_compare_ac.png")