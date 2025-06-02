# singleflip.jl
using Random
using Plots
using Statistics
using DelimitedFiles
using Dates
using ForwardDiff
using QuadGK
using SpecialFunctions
using LaTeXStrings
using Printf

include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
include("../../core/main.jl")

function calculate_C_from_derivative(energy_runs::Vector{Vector{Float64}}, beta_values::Vector{Float64})

    energies = mean(hcat(energy_runs...), dims=2)[:]
    dE_dBeta = Float64.(diff(energies) ./ diff(beta_values))
    beta_mid = Float64.((beta_values[1:end-1] .+ beta_values[2:end]) ./ 2)
    C = - dE_dBeta .* beta_mid.^2
    return C
end


function calculate_C_from_fluctuations(energy_runs::Vector{Vector{Float64}}, beta_values::Vector{Float64})

    C_values = Float64[]


    for i in 1:length(beta_values)
        beta = beta_values[i]
        

        energies_at_beta = [energy_runs[j][i] for j in 1:length(energy_runs)]
        

        E_avg = mean(energies_at_beta)
        E_squared_avg = mean(energies_at_beta .^ 2)
        

        C = (E_squared_avg - E_avg^2) * beta^2
        

        push!(C_values, C)
    end


    return C_values
end


function entropy_from_heat_capacity(C::Vector{Float64}, beta_values::Vector{Float64})
    beta_values = beta_values[1:end-1]
    T_values = 1.0 ./ beta_values
    reverse!(T_values)
    reverse!(C)

    S = Float64[0.0]
    for i in 2:length(T_values)
        dT = T_values[i] - T_values[i-1]  
        avg_C_over_T = (C[i]/T_values[i] + C[i-1]/T_values[i-1]) / 2  
        push!(S, S[end] + avg_C_over_T * dT) 
    end
    return S, T_values

end




function onsager_free_energy(b, L)
    J  = 1

    k = 2 * sinh(2*b*J) / (cosh(2*b*J)^2)
    

    function integrand(theta)
        return log(1 + sqrt(1 - k^2 * sin(theta)^2))
    end
    

    I, _ = quadgk(integrand, 0, pi/2)
    I /= pi
    

    return (log(sqrt(2) * cosh(2*b*J)) + I) / (-b)
end


function onsager_energy(b, L)
    J = 1
    k = 2 * sinh(2*b*J) / (cosh(2*b*J)^2)
    K_k = ellipk(k^2)

    return -J * coth(2 * b * J) * (1+2/pi * K_k * (2 * tanh(2*b*J)^2 -1))
    
    
end


function onsager_specific_heat(b, L)
    J = 1
    k = 2 * sinh(2*b*J) / (cosh(2*b*J)^2)
    K_k = ellipk(k^2)
    E_k = ellipe(k^2)
    k_prime = 2 * tanh(2*b*J)^2 -1
    return (b*J*coth(2*b*J))^2 * 2/pi * (2 * K_k - 2 * E_k - (1 - k_prime)*(pi/2 + k_prime * K_k))


end


function onsager_entropy(b, L)
    U = onsager_energy(b, L)
    F = onsager_free_energy(b, L)
    S = b * (U - F)
    return S
end



function figure_C(energy_runs::Vector{Vector{Float64}}, beta_values::Vector{Float64}, filename::String, folder::String, L::Int64)
    C_derivative = calculate_C_from_derivative(energy_runs, beta_values)
    C_fluctuations = calculate_C_from_fluctuations(energy_runs, beta_values)
    
    p = plot()
    plot!(p, beta_values[1:end-1], abs.(C_derivative), label = "Derivative", lw = 2)
    plot!(p, beta_values, onsager_specific_heat.(beta_values,L), label = "Onsager specific heat capacity", lw = 2)
    plot!(p, beta_values, abs.(C_fluctuations * L^2), label = "Fluctuations", lw = 2)
    
    xlabel!(p, "Beta", xlabelcolor = :white)
    ylabel!(p, "C", ylabelcolor = :white)
    title!(p, "Heat capacity from single flip energy runs, L = $L", titlecolor = :white)
    vline!(p, [log(1+sqrt(2))/2], label="beta_c = ln(1+sqrt2)/2", color=:red, linestyle=:dash)
    savefig(p,joinpath(folder, filename))

    p1 = plot(legend = :topleft)
    plot!(p1, 1 ./ beta_values[1:end-1], abs.(C_derivative), label = "Metropolis-Hastings Simulation", lw = 2, colour = :blue)
    plot!(p1, 1 ./ beta_values, onsager_specific_heat.(beta_values,L), label = "Analytical result", lw = 2, colour = :red)
    #plot!(p, beta_values, abs.(C_fluctuations * L^2), label = "Fluctuations", lw = 2)
    
    xlabel!(p1, L"$ \mathrm{T}$", xlabelcolor = :white)
    ylabel!(p1, L"Average Heat Capacity per site $\langle \mathrm{C(T)} \rangle$", ylabelcolor = :white)
    #title!(p1, "Heat capacity from single flip energy runs, L = $L", titlecolor = :white)
    vline!(p1, [2/log(1+sqrt(2))], label=L"$T_c = 2/ \ln (1+\sqrt{2})$", color=:purple, linestyle=:dash)
    savefig(p1,joinpath(folder, "beta_$filename"))
end

L = 100
beta_values = 1 ./ generate_T_intervals(4.0, 0.25, 100)
#E_runs = readdlm("experiments\\singleflip\\single_flips.csv", ',', Float64)
E_runs = readdlm("single_flips_fixed.csv", ',', Float64)
E_runs_vector = [collect(row) for row in eachrow(E_runs)]
figure_C(E_runs_vector, beta_values, "C_singleflip_many.png", "", L)







C = calculate_C_from_derivative(E_runs_vector, beta_values)
S_values,T_vals = entropy_from_heat_capacity(C, beta_values)
energies = mean(hcat(E_runs_vector...), dims=2)[:]


#plot2 = plot()
#plot!(plot2, T_vals, S_values)  # Plot entropy vs. temperature
#xlabel!(plot2, "T", xlabelcolor = :white)
#ylabel!(plot2, "S(T)", ylabelcolor = :white)
#title!(plot2, "Entropy from heat capacity, L = $L", titlecolor = :white)
#savefig(plot2,"experiments\\singleflip\\entropy.png")


plot3 = plot()
#plot!(plot3, 1 ./beta_values, energies, label = "Single flip energy", lw = 2)
plot!(plot3, beta_values, onsager_free_energy.(beta_values,L), label = "Onsager free energy", lw = 2)
xlabel!(plot3, "Beta", xlabelcolor = :white)
ylabel!(plot3, "F(T)", ylabelcolor = :white)
title!(plot3, "onsager free energy, L = $L", titlecolor = :white)
savefig(plot3,"onsager free energy_many.png")


plot4 = plot()
plot!(plot4, beta_values, onsager_energy.(beta_values,L), label = "Onsager energy", lw = 2)
xlabel!(plot4, "Beta", xlabelcolor = :white)
ylabel!(plot4, "E(T)", ylabelcolor = :white)
title!(plot4, "onsager energy, L = $L", titlecolor = :white)
savefig(plot4,"onsager energy_many.png")

#plot5 = plot()
#T_c = 2 / log(1+sqrt(2))
#println(T_c)
#plot!(plot5, beta_values, onsager_specific_heat.(beta_values,L), label = "Onsager specific heat capacity", lw = 2)
#plot!(plot5, beta_values[1:end-1], C, label = "Specific heat capacity from single flip", lw = 2)
#vline!(plot5, [1/T_c], label="beta_c = ln(1+sqrt2)/2", color=:red, linestyle=:dash)
#xlabel!(plot5, "Beta", xlabelcolor = :white)
#ylabel!(plot5, "C(T)", ylabelcolor = :white)
#title!(plot5, "onsager heat capacity, L = $L", titlecolor = :white)
#savefig(plot5,"onsager heat capacity.png")

plot6 = plot(legend = :topleft)
plot!(plot6, 1 ./beta_values, onsager_entropy.(beta_values,L), label = "Analytical Result", lw = 2, linestyle=:dash, color=:red)
plot!(plot6, T_vals, S_values, label = "Metropolis-Hastings Simulation", lw = 2, color = :blue)
hline!(plot6, [log(2)], label=L"\ln(2)", color=:green, linestyle=:dash)
vline!(plot6, [2/log(1+sqrt(2))], label=L"$T_c = 2/ \ln (1+\sqrt{2})$", color=:purple, linestyle=:dash)
xlabel!(plot6, L"\mathrm{T}", xlabelcolor = :white)
ylabel!(plot6, L"Average entropy per site $\langle \mathrm{S(T)} \rangle $", ylabelcolor = :white)
#title!(plot6, "onsager entropy, L = $L", titlecolor = :white)
xticks!(plot6, 1:4, ["1.0", "2.0", "3.0", "4.0"])
savefig(plot6,"onsager entropy_many.png")

plot7 = plot()
plot!(plot7, 1 ./beta_values, energies, label = "Metropolis-Hastings Simulation", lw = 2, color = :blue)
plot!(plot7, 1 ./beta_values, onsager_energy.(beta_values,L), label = "Analytical Result", lw = 2, linestyle=:dash, color=:red)
xlabel!(plot7, L"\mathrm{T}", xlabelcolor = :white)
ylabel!(plot7, L"Average energy per site $\langle \mathrm{E(T)} \rangle $", ylabelcolor = :white)
vline!(plot7, [2/log(1+sqrt(2))], label=L"$T_c = 2/ \ln (1+\sqrt{2})$", color=:purple, linestyle=:dash)
#title!(plot7, "Energy", titlecolor = :white)
xticks!(plot7, 1:4, ["1.0", "2.0", "3.0", "4.0"])
savefig(plot7,"single_flip_energy_many.png")
