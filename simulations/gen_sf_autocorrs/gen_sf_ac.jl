using DelimitedFiles
using Base.Threads
using Plots
using Statistics
using Dates
#finds value (l,L) pairs for which configuration space matches single spin flip
# then runs simulation for each pair, using single spin flips to equilibrate

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

function run_parallel_decorrelation(lattice::Lattice, beta_values::Vector{Float64}, filename::String, equilib_steps::Int64, measurement_steps::Int64, max_lag::Int64, copies::Int64, l::Int64 = 1, move::String = "single flip")

    open("$(filename)_$(l)_$(lattice.L)" * ".log", "w") do io
        println(io, "Starting simulation at $(now())")
        println(io, "Lattice size: ", lattice.L)
        println(io, "Beta values: ", beta_values)
        println(io, "Equilibration steps: $equilib_steps")
        println(io, "Measurement steps: $measurement_steps")
        println(io, "Max lag: $max_lag")
        println(io, "Copies: $copies")
        println(io, "Move type: $move")
        println(io, "l = $l")
        println(io, "L = $(lattice.L)")
        println(io, "----------------------------------")
    end

    E_autocorr_sum = [zeros(max_lag+1) for _ in 1:length(beta_values)]
    M_autocorr_sum = [zeros(max_lag+1) for _ in 1:length(beta_values)]
    unconnected_autocorr_sum = [zeros(max_lag+1) for _ in 1:length(beta_values)]
    E_autocorr_count = zeros(Int, length(beta_values))
    M_autocorr_count = zeros(Int, length(beta_values))
    lock_obj = ReentrantLock()

    @threads for i in 1:copies
        println("Starting copy $(i)")
        open("$(filename)_$(l)_$(lattice.L)" * ".log", "a") do io
            println(io, "Current copy: $i")
        end

        time = now()
        thread_lattice = Lattice(lattice.L)
        # Here we should use generate_E_m_autocorrelations_with_sf as annealing should be done with single spin flips
        # but give the same move = "single flip" this is okay
        _, E_autocorr_beta, M_autocorr_beta, unconnected_autocorr_beta = generate_E_m_autocorrelations(thread_lattice, beta_values, l = l, move = move, max_lag = max_lag, measurement_steps = measurement_steps, equilib_steps = equilib_steps, normalise = false)
        lock(lock_obj)
        try
            open("L=$(lattice.L)_l=$(l)_copy$(i)_$(filename)_E.csv", "w") do io
                for entry in E_autocorr_beta
                    writedlm(io, [entry], ',')
                end
            end
            open("L=$(lattice.L)_l=$(l)_copy$(i)_$(filename)_m.csv", "w") do io
                for entry in M_autocorr_beta
                    writedlm(io, [entry], ',')
                end
            end
            open("L=$(lattice.L)_l=$(l)_copy$(i)_$(filename)_unconnected.csv", "w") do io
                for entry in unconnected_autocorr_beta
                    writedlm(io, [entry], ',')
                end
            end
        finally
            unlock(lock_obj)
        end        
    end
end
Tc = log(1+sqrt(2))/2 # = 2.269
beta_values = 1 ./ generate_T_intervals(10.0,0.65,100)
println(beta_values)
equilib_steps = 1000
measurement_steps = 2500
max_lag = 500
copies = 8


move = "single flip"
filename = "autocorrelation"
 


for L in L_values
    lattice = Lattice(L)
    run_parallel_decorrelation(lattice, beta_values, filename, equilib_steps, measurement_steps, max_lag, copies, 1, move)
end
