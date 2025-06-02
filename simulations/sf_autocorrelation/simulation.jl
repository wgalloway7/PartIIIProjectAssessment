using DelimitedFiles
using Base.Threads
using Plots
using Statistics
using Dates

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

        
        




# runs simulation to generate autocorrelation function from time series average on parrallel threads
# writes autocorrelation function run to file, computes average over copies
function run_parallel_decorrelation(lattice::Lattice, beta_values::Vector{Float64}, filename::String, equilib_steps::Int64, measurement_steps::Int64, max_lag::Int64, copies::Int64, l::Int64 = 1, move::String = "single flip")
    
    open(filename * ".log", "w") do io
        println(io, "Starting simulation at $(now())")
        println(io, "Lattice size: ", lattice.L)
        println(io, "Beta values: ", beta_values)
        println(io, "Equilibration steps: $equilib_steps")
        println(io, "Measurement steps: $measurement_steps")
        println(io, "Max lag: $max_lag")
        println(io, "Copies: $copies")
        println(io, "Move type: $move")
        println(io, "----------------------------------")
    end
    
    

    
    autocorr_sum = [zeros(max_lag+1) for i in 1:length(beta_values)]
    lock_obj = ReentrantLock()
    @threads for i in 1:copies

        println("Starting copy $(i)")
        open("progress_$(filename).log", "w") do io
            println(io, "Current copy: $i")
        end


        time = now()
        thread_lattice = Lattice(lattice.L)
        autocorr_beta = generate_autocorrelation(thread_lattice, beta_values, l = l, move = move, max_lag = max_lag, measurement_steps = measurement_steps, equilib_steps = equilib_steps)
        open("l=$(l)_copy$(i)_$(filename)", "w") do io
            for entry in autocorr_beta
                autocorr = entry[2]
                writedlm(io, [autocorr], ',')
            end
        end
        lock(lock_obj)
        println("Finished generating AC copy $(i) in $(now() - time)")
        open(filename * ".log", "w") do io
            println(io, "Finished AC generation for copy $i at $(now()) in $(now() - time)")
        end
        
        try
            autocorr_sum .+= [autocorr_beta[i][2] for i in 1:length(beta_values)]
        finally
            unlock(lock_obj)
        end
    end
    autocorr_avg = autocorr_sum ./ copies
    time = now()
    println("Fitting averaged autocorrelation function")
    stretched_fit_params = generate_stretched_fit(autocorr_avg, beta_values, plot_fit = false)
    #VFT_fit_params = fit_VFT_from_stretched_exponential_fit(beta_values, stretched_fit_params)
    open("l=$(l)_fit_$(filename)", "w") do io

        writedlm(io, [["tau", "alpha", "epsilon"]],',')
        for i in 1:length(beta_values)
            writedlm(io, [stretched_fit_params[i]],',')
        end
    end
    println("Finished fitting averaged autocorrelation function in $(now() - time)")

end

beta_values = 1 ./ generate_T_intervals(4.0,1.0,100)
equilib_steps = 2500
measurement_steps = 2500
max_lag = 500
copies = 12
l_values = [1]
move = "l chain flip"
filename = "autocorrelation.csv"



lattice = Lattice(100)
run_parallel_decorrelation(lattice, beta_values, filename, equilib_steps, measurement_steps, max_lag, copies, 1, move)


    


