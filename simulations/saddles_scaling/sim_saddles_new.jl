using DelimitedFiles
using Plots
using Base.Threads
using Dates
using StatsBase

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

function saddles_simulation(lattice::Lattice, beta_values::Vector{Float64}, l_values::Vector{Int64}, copies::Int64, filename::String, move::String, p::Int64 = 30, do_plot::Bool = false)
    println("started")
    L = lattice.L
    open(filename * ".log", "w") do io
        println(io, "Starting saddles simulation at $(now())")
        println(io, "Lattice size: ", lattice.L)
        println(io, "Beta values: ", beta_values)
        println(io, "Copies: $copies")
        println(io, "Move type: $move")
        println(io, "----------------------------------")
    end

    energy_list        = [Vector{Vector{Float64}}() for _ in 1:length(l_values)]
    saddle_index_list  = [Vector{Vector{Float64}}() for _ in 1:length(l_values)]
    beta_list          = [Vector{Vector{Float64}}() for _ in 1:length(l_values)]

    lock_obj = ReentrantLock()
    tau_to_use = ceil.(Int, lattice.tau_values)
    println("tau_to_use: ", tau_to_use)
    @threads for i in 1:copies
        #println("Starting copy $(i)")
        thread_lattice = Lattice(lattice.L)
        K_vals, E_vals, b_vals = generate_saddles_parallel(thread_lattice, beta_values, l_values, tau_to_use, 100, move, i)

        lock(lock_obj) do
            for (j,l) in enumerate(l_values)
                push!(energy_list[j], E_vals[j,:])
                push!(saddle_index_list[j], K_vals[j,:])
                push!(beta_list[j], b_vals[j,:])
            end

            open(filename * ".log", "a") do io
                println(io, "Finished copy $(i) at $(now())")
            end

        end
    end
    for (j,l) in enumerate(l_values)
        open("l=$(l)_$(filename)", "w") do io
            writedlm(io, ["Energy" "SaddleIndex" "Beta"])
            for (E, K, B) in zip(vcat((energy_list[j])...), vcat((saddle_index_list[j])...), vcat((beta_list[j])...))
                writedlm(io, [E K B])
            end
        end
    end
    E_lower_crit_vals = []
    E_upper_crit_vals = []
    for (j,l) in enumerate(l_values)
        
        all_energies = vcat((energy_list[j])...)
        all_saddles = vcat((saddle_index_list[j])...)
        all_betas = vcat((beta_list[j])...)
        println("First few  l = $(l) energies: ", all_energies[1:5])
        println("First few l = $(l) saddles: ", all_saddles[1:5])
        println("First few l = $(l) betas: ", all_betas[1:5])


        energy_groups = Dict{Float64, Vector{Float64}}()
        for (E, K) in zip(all_energies, all_saddles)
            if haskey(energy_groups, E)
                push!(energy_groups[E], K)
            else
                energy_groups[E] = [K]
            end
        end

        proportions = Dict{Float64, Tuple{Float64, Float64, Float64}}()
        k_mean = Dict{Float64, Float64}()

        for (E, K_values) in energy_groups
            total = length(K_values)
            p_K0 = count(x -> x == 0, K_values) / total
            p_K1 = count(x -> x == 1, K_values) / total
            p_K2 = count(x -> x >= 2, K_values) / total
            proportions[E] = (p_K0, p_K1, p_K2)

            k_mean[E] = mean(K_values)
        end

        open("proportions_l=$(l)_$(filename)", "w") do io
            writedlm(io, ["Energy" "p(K=0)" "p(K=1)" "p(K>=2)"])
            data_rows = [[E, p_K0, p_K1, p_K2] for (E, (p_K0, p_K1, p_K2)) in proportions]
            writedlm(io, data_rows)

        end

        if do_plot

            p = plot()
            xlabel!(p, "Energy", xlabelcolor=:white)
            ylabel!(p, "Saddle index", ylabelcolor=:white)
            title!(p, "Saddle index proportion vs Energy for L = $(lattice.L), copies = $copies, $move", titlecolor=:white)
            first_K0 = true
            first_K1 = true
            first_K2 = true
            for (E, (p_K0, p_K1, p_K2)) in proportions
                scatter!(p, [E], [p_K0], label=first_K0 ? "K = 0" : "", color="blue", markersize=1, markerstrokewidth=0.1)
                scatter!(p, [E], [p_K1], label=first_K1 ? "K = 1" : "", color="red", markersize=1, markerstrokewidth=0.1)
                scatter!(p, [E], [p_K2], label=first_K2 ? "K >= 2" : "", color="green", markersize=1, markerstrokewidth=0.1)
                first_K0 = false
                first_K1 = false
                first_K2 = false
            end
        

            Es = sort(collect(keys(energy_groups)))
            

            p_K0_vals = [proportions[E][1] for E in Es]
            p_K1_vals = [proportions[E][2] for E in Es]
            p_K2_vals = [proportions[E][3] for E in Es]

            E_lower_critical_index = findfirst(x -> x > 0.5, p_K2_vals)
            E_upper_critical_index = findlast(x -> x < 0.5, p_K2_vals)
            #println("Lower critical index: ", E_lower_critical_index)
            #println("Upper critical index: ", E_upper_critical_index)
            E_lower_critical = E_lower_critical_index !== nothing ? Es[E_lower_critical_index] : nothing
            E_upper_critical = E_upper_critical_index !== nothing ? Es[E_upper_critical_index] : nothing

        #println("Lower critical energy: ", E_lower_critical)
        #println("Upper critical energy: ", E_upper_critical)


        #p2 = plot()
        #plot!(p2, Es, p_K2_vals, label="K = 0", color="blue", markersize=1, markerstrokewidth=0.1)
        #if E_lower_critical !== nothing
        #    vline!(p2, [E_lower_critical], label="Lower critical energy", color="black", linestyle=:dash)
        #end
        #if E_upper_critical !== nothing
        #    vline!(p2, [E_upper_critical], label="Upper critical energy", color="black", linestyle=:dash)
        #end
        #savefig(p2, "saddle_proportion_$(k)_line.png")


        
            if E_lower_critical !== nothing
                vline!(p, [E_lower_critical], label="Lower critical energy", color="black", linestyle=:dash)
            end
            if E_upper_critical !== nothing
                vline!(p, [E_upper_critical], label="Upper critical energy", color="black", linestyle=:dash)
            end
            
            
            ylabel!(p, "Saddle index proportion")
            xlabel!(p, "Energy")
            savefig(p, "saddle_proportion_$(L)_$(l).png")


       

            q = plot(legend = false)
            for (E, K_values) in energy_groups
                scatter!(q, [E for _ in 1:length(K_values)], K_values, label="E = $E", color="blue", markersize=1, markerstrokewidth=0)
            end
            xlabel!(q, "Energy", xlabelcolor=:white)
            ylabel!(q, "Saddle index", ylabelcolor=:white)
            savefig(q,"energy_saddles_$(L)_$(l).png")


            r = plot()
            xlabel!(r, "Beta")
            ylabel!(r, "Energy")
            scatter!(r, all_betas, all_energies, xlabel="Beta", ylabel="Energy", title="Energy vs Beta for L = $(lattice.L), l = $l", color="blue", markersize=1, markerstrokewidth=0)
            plot!(r, sort(all_betas), onsager_energy.(sort(all_betas), lattice.L), label="Onsager energy", color="red", linestyle=:dash)
            savefig(r, "energy_vs_beta_$(L)_$(l)_" * filename * ".png")
        end
        push!(E_lower_crit_vals, E_lower_critical)
        push!(E_upper_crit_vals, E_upper_critical)
    end
    return E_lower_crit_vals, E_upper_crit_vals
end

# main execution

#beta_values = 1 ./ generate_T_intervals(4.0, 1.0, 100)
beta_values = 1 ./ generate_T_intervals(10.0,0.65,100)

move = "l chain flip"
copies = 600

E_lower_crit_vals_vals = []
E_upper_crit_vals_vals = []

L_values = 10:50
L_l_pairs = Dict{Int64, Vector{Int64}}() 
for L in L_values
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





for (L,l_values) in L_l_pairs
    filename = "saddles_simulation_$(L).csv"

    open("saddle_critical_energies$(L).csv", "w") do io
        writedlm(io, [["l_values", "E_lower_critical", "E_upper_critical"]])
    end
    

    
    lattice = Lattice(L)


    
    println("Running simulation for L = $(L), l = $(l_values)")
    E_lower_vals, E_upper_vals = saddles_simulation(lattice, beta_values, l_values, copies, filename, move, 100, true)
    push!(E_lower_crit_vals_vals, E_lower_vals)
    push!(E_upper_crit_vals_vals, E_upper_vals)

    open("saddle_critical_energies_$(L).csv", "a") do io
        writedlm(io, [[l_values, E_lower_vals, E_upper_vals]])
    end
end


