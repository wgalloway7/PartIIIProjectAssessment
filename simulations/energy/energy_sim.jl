using DelimitedFiles
using Plots
using Statistics
using Dates
using Base.Threads

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")




function process_copy(i, lattice, beta_values, l, local_p_values, move, tau_vector_file, energy_sum, lock_obj, filename)
    println("Starting copy $(i)")
    open("progress_$(filename).log", "w") do io
        println(io, "Current copy: $i")
    end
    time = now()
    thread_lattice = Lattice(lattice.L, tau_vector_file)
    for (j, p) in enumerate(local_p_values)
        energies = generate_energies(thread_lattice, beta_values, l, lattice.tau_values * p, move)
        lock(lock_obj) do
            for idx in eachindex(energy_sum[j])
                energy_sum[j][idx] += energies[idx]
            end
        end
    end
    println("Finished generating energies copy $(i) in $(now() - time)")
    open(filename * ".log", "a") do io
        println(io, "Finished energy generation for copy $i at $(now()) in $(now() - time)")
    end
end

function run_parallel_energy(lattice::Lattice, beta_values::Vector{Float64}, l::Int64,
                             p_values::Vector{Int64}, filename::String, copies::Int64,
                             move::String = "single flip", do_plot::Bool = false)
    println("started")
    open(filename * ".log", "w") do io
        println(io, "Starting simulation at $(now())")
        println(io, "Lattice size: ", lattice.L)
        println(io, "Beta values: ", beta_values)
        println(io, "Copies: $copies")
        println(io, "Move type: $move")
        println(io, "l = ", l)
        println(io, "Tau values: ", lattice.tau_values)
        println(io, "Equilibriating for p*tau steps, where p = ", p_values)
        println(io, "----------------------------------")
    end


    energy_sum = [zeros(length(beta_values)) for i in 1:length(p_values)]
    lock_obj = ReentrantLock()


    local_p_values = p_values
    local_tau_vector_file = tau_vector_file

    Threads.@threads for i in 1:copies
        process_copy(i, lattice, beta_values, l, local_p_values, move, local_tau_vector_file, energy_sum, lock_obj, filename)
    end

    energy_avg = [energy_sum[i] ./ copies for i in 1:length(p_values)]


    if do_plot
        p = plot()
        colors = (cgrad(:RdBu, length(p_values)))
        xlabel!(p, "Beta", xlabelcolor=:white)
        ylabel!(p, "Energy", ylabelcolor=:white)
        title!(p, "Energy vs Beta for various p, L = $(lattice.L), copies = $copies, $move", titlecolor=:white)
    end

    for (i, p_val) in enumerate(p_values)
        open("l=$(l)_p=$(p_val)_$(filename).csv", "w") do io
            writedlm(io, energy_avg[i], ',')
        end

        if do_plot
            plot!(p, beta_values, energy_avg[i], label = "p = $p_val", color = colors[i])
        end
    end

    if do_plot
        savefig(p, filename * ".png")
    end
end




tau_vector_file = "../../core/average_tau_wide1.csv"
lattice = Lattice(100, tau_vector_file)

beta_values = 1 ./ generate_T_intervals(4.0, 1.0, 100)

global p_values = [16]

move = "l chain flip"
copies = 20
println(p_values)
l_values = [2i - 1 for i in 1:50]

for l in l_values
    run_parallel_energy(lattice, beta_values, l, p_values, "energy_p_$(l)_", copies, move, true)
end

        


        




