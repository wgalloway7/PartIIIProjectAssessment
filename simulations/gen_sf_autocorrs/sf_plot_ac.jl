using DelimitedFiles
using Plots
using Statistics
using LsqFit
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")   
init_d = readdlm("L=96_l=1_copy1_autocorrelation_E.csv", ',')
beta_values = 1 ./ generate_T_intervals(10.0,0.65,100)

L_values = [i for i in 95:105]
#L_values = [95,96]
rows, cols = size(init_d)
copies = 8
for L in L_values
    d0 = zeros(rows, cols)
    for j in 1:copies
        d = readdlm("L=$(L)_l=1_copy$(j)_autocorrelation_E.csv", ',')
        for i in 1:rows
            d[i, :] = d[i,:]  ./ d[i,1]
        end
        d0 = d0 .+ d
        if NaN in d
            println("NaN found in L=$(L)_l=1_copy$(j)_autocorrelation_E.csv")
        end
    end
    d0 = d0 ./ copies

    p = plot()
    for i in 1:rows
        plot!(p, d0[i, :])
    end
    savefig(p, "zL=$(L)_l=1_autocorrelation_E.png")

    function stretched_exponential(x,p)
        tau,beta = p
        return exp( - (x/tau)^(beta))
    end

    tau_values = []

    for i in 1:rows
        fit = fit_stretched_exponential_from_autocorrelation_function(d0[i, :])
        push!(tau_values, fit[1])
    end
        

            
  
    #println(tau_values)
    q = plot()
    plot!(q, beta_values, tau_values, label = "L=$(L)_l=1")
    xlabel!("beta")
    ylabel!("tau")
    title!("L=$(L)_l=1")
    savefig(q, "zL=$(L)_l=1_autocorrelation_E_tau.png")

    writedlm("zL=$(L)_l=1_autocorrelation_E_tau.csv", (beta_values, tau_values), ',')


end