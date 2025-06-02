using Random
using Plots
using Dates
using DelimitedFiles
using Base.Threads

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

function run_parallel_simulation(lattice::Lattice, beta_values::Vector{Float64}, filename::String, equilib_steps::Int64, measurement_steps::Int64, max_lag::Int64, copies::Int64)
    open(filename * ".log", "w") do io
        println(io, "Starting simulation at $(now())")
        println(io, "Generating tau values for spin, mangetisation and energy")
        println(io, "Lattice size: ", lattice.L)
        println(io, "Beta values: ", beta_values)
        println(io, "Copies: $copies")
        println(io, "Equilibriating for $(equilib_steps) steps")
        println(io, "Measuring for $(measurement_steps) steps")
        println(io, "Max lag: $(max_lag)")
        println(io, "----------------------------------")
    end
    
    
    E_tau_sum =zeros(length(beta_values))
    M_tau_sum =zeros(length(beta_values))
    S_tau_sum =zeros(length(beta_values))

    lock_obj = ReentrantLock()
    @threads for i in 1:copies
        println("Starting copy $(i)")
        open(filename * ".log", "a") do io
            println(io, "Starting copy: $i")
        end
        start_time = now()


        thread_lattice = Lattice(lattice.L)  


        E_tau_values, M_tau_values, S_tau_values = generate_all_tau(thread_lattice, beta_values, true, filename = "$(i)_$(filename)", max_lag = max_lag, equilib_steps = equilib_steps, measurement_steps = measurement_steps)

        lock(lock_obj)
        println("Finished copy $(i) in $(now() - start_time)")
        open(filename * ".log", "a") do io
            println(io, "Finished AC generation for copy $i at $(now()) in $(now() - start_time)")
        end
        try
            E_tau_sum .+= E_tau_values
            M_tau_sum .+= M_tau_values
            S_tau_sum .+= S_tau_values
        finally
            unlock(lock_obj)
        end
    end
    return E_tau_sum ./ copies, M_tau_sum ./ copies, S_tau_sum ./ copies
end

L = 100
res = 100
lattice = Lattice(L)
beta_values = 1 ./ generate_T_intervals(4.0, 0.3, res)

filename = "compare_ac_new.csv"
equilib_steps = 2000
measurement_steps = 2000
max_lag = 250
copies = 12


E_tau_runs, M_tau_runs, S_tau_runs = run_parallel_simulation(lattice, beta_values, filename, equilib_steps, measurement_steps, max_lag, copies)

p = plot()
plot!(p, beta_values, E_tau_runs, label = "E_tau", title = "E_tau vs Beta", xlabel = "Beta", ylabel = "E_tau")


plot!(p, beta_values, M_tau_runs, label = "M_tau", title = "M_tau vs Beta", xlabel = "Beta", ylabel = "M_tau")


plot!(p, beta_values, S_tau_runs, label = "S_tau", title = "S_tau vs Beta", xlabel = "Beta", ylabel = "S_tau")

savefig(p, "compare_ac.png")
