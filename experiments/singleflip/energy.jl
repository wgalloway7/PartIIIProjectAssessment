using DelimitedFiles

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")
data = readdlm("single_flips_fixed.csv", ',')
beta_values = 1 ./ generate_T_intervals(4.0, 0.25, 100)

rows, columns = size(data)
println(rows)
p = plot()
for i in 1:rows
    plot!(p,1 ./ beta_values, data[i, :], label="Run $i", xlabel="Temperature (T)", ylabel="Energy (E)")
end
vline!(p, [2 / log(1 + sqrt(2))], label="Critical temperature", color=:red, linestyle=:dash)
savefig(p,"single_flips_fixed.png")