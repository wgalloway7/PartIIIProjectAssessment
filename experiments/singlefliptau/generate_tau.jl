using Random
using Plots
using Statistics
using Dates
using DelimitedFiles
using Base.Threads
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")


function run_parallel_simulation(lattice::Lattice, beta_values::Vector{Float64}, filename::String, equilib_steps::Int64, measurement_steps::Int64, max_lag::Int64, copies::Int64, l::Int64 = 1, move::String = "single flip")
    tau_runs = []


    lock_obj = ReentrantLock()
    @threads for i in 1:copies
        println("Starting copy $(i)")
        time = now()


        thread_lattice = Lattice(lattice.L)


        tau_values = generate_tau(thread_lattice, beta_values, true, l = l, move = move, filename = "$(l)_$(i)_$(filename)", max_lag = max_lag, equilib_steps = equilib_steps, measurement_steps = measurement_steps)


        lock(lock_obj)
        println("Finished copy $(i) in $(now() - time)")
        try
            push!(tau_runs, tau_values)
        finally
            unlock(lock_obj)
        end
    end

    return tau_runs
end
L = 100
res = 100
lattice = Lattice(L)
beta_values= 1 ./ generate_T_intervals(4.0,1.0,res)
filename = "autocorrelation_wide.csv"
equilib_steps = 5000
measurement_steps = 5000
max_lag = 1000
copies = 12
l_values = [1]
move = "l chain flip"

p = plot()
vline!(p, [2 / log(1 + sqrt(2))], label = "Tc")
xlabel!(p,"Temperature")
ylabel!(p,"Tau (MC steps)")


q = plot()
vline!(q, [log(1 + sqrt(2)) / 2], label = "Tc")
xlabel!(q,"Beta")
ylabel!(q,"Tau (MC steps)")



l = 1
time = now()
println("Starting l = $(l)")
tau_runs = run_parallel_simulation(lattice, beta_values, filename, equilib_steps, measurement_steps, max_lag, copies, l, move)
average_tau = mean(hcat(tau_runs...), dims = 2)[:]
println(now() - time)
writedlm("average_tau_wide$(l).csv", (beta_values, average_tau), ',')

    
plot!(p, 1 ./ beta_values, average_tau, label = "l = $(l)")
plot!(q, beta_values, average_tau, label = "l = $(l)")
    


savefig(p, "tau_wide.png")
savefig(q, "tau_wide_beta.png")