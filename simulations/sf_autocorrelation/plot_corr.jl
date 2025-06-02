using Random
using Plots
using Statistics
using Dates
using DelimitedFiles
using Base.Threads

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")
beta_values= 1 ./ generate_T_intervals(4.0,1.0,100)
d = readdlm("l=1_copy1_autocorrelation.csv", ',')
data_sum = zeros(size(d))
for j in 1:12
    data = readdlm("l=1_copy$(j)_autocorrelation.csv", ',')
    data_sum .+= data
end
data = data_sum ./ 12

params = readdlm("l=1_fit_autocorrelation.csv", ',', skipstart = 1)
function stretched_exponential(t, p)
    tau, alpha, epsilon = p
    return exp.(-abs.(t ./ tau).^alpha) .+ epsilon
end


beta_indexes_to_plot = [1,10,20,30,50,70,100]
for i in 1:size(data)[1]
    if i in beta_indexes_to_plot
        p = plot()
        row = data[i, :]
        x_data = length(row)
        plot!(p, row)
        plot!(p, stretched_exponential(0:x_data-1, params[i, :]), color = :red, linestyle = :dash)
        xlims!(0,50)
        xlabel!("Time")
        ylabel!("Correlation")
        savefig(p, "correlation_$(beta_values[i]).png")
    end
end

tau = params[:,1]
q = plot()
plot!(q, beta_values, tau)
vline!([log(1+sqrt(2))/2], color = :red, linestyle = :dash, label = "Tc")
xlabel!("Beta")
ylabel!("Energy autocorrelation tau (MC steps)")

savefig(q, "tau_vs_beta.png")

