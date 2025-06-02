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

    @threads :dynamic for i in 1:copies
        println("Starting copy $(i)")
        open("$(filename)_$(l)_$(lattice.L)" * ".log", "a") do io
            println(io, "Current copy: $i")
        end

        time = now()
        thread_lattice = Lattice(lattice.L)
        _, E_autocorr_beta, M_autocorr_beta, unconnected_autocorr_beta, E_avg_beta, M_avg_beta = generate_E_m_autocorrelations_with_sf(thread_lattice, beta_values, l = l, move = move, max_lag = max_lag, measurement_steps = measurement_steps, equilib_steps = equilib_steps, sim_copy = i, normalise = false, show_configs = false)
        
        
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
        open("L=$(lattice.L)_l=$(l)_copy$(i)_$(filename)_E_avg.csv", "w") do io 
            for entry in E_avg_beta
                writedlm(io, [entry], ',')
            end
        end
        open("L=$(lattice.L)_l=$(l)_copy$(i)_$(filename)_M_avg.csv", "w") do io 
            for entry in M_avg_beta
                writedlm(io, [entry], ',')
        end

        open("$(filename)_$(l)_$(lattice.L)" * ".log", "a") do io
            println(io, "Finished at time $(now())")
        end

        
    end

               
    end
end
Tc = 2/log(1+sqrt(2)) # = 2.269
println(Tc)
beta_values = 1 ./ generate_T_intervals(4.0,Tc,20)
println(beta_values)
equilib_steps = 1000
measurement_steps = 5001
max_lag = 1000
copies = 20


move = "l chain flip"
filename = "autocorrelation"
L_values = [i for i in 95:105]
L_l_pairs = Dict{Int64, Vector{Int64}}() 
for L in L_values
    L_l_pairs[L]=[1]
    for l in 1:L
        if valid_l(L,l) && l % 2 == 1 && l != 1
            if haskey(L_l_pairs, L)
                push!(L_l_pairs[L], l)
            else
                L_l_pairs[L] = [l]
            end
        end
    end
end
L_l_pairs = sort(collect(L_l_pairs), by = x -> x[1])
for (L, l_values) in L_l_pairs
    println("L = $L, l = $l_values")
end

completed_pairs = [[95,1],[95,3],[95,7],[95,9],[95,21],[95,27], [95,63],[96,1],[97,1],[97,3],[97,5]]
for (L, l_values) in L_l_pairs
    lattice = Lattice(L)
    for l in l_values
        if !([L,l] in completed_pairs)
            run_parallel_decorrelation(lattice, beta_values, filename, equilib_steps, measurement_steps, max_lag, copies, l, move)
        end
    end
end
        
