#lattice.jl
using Random
using LinearAlgebra
using Distributions
using Combinatorics
using DelimitedFiles
using ForwardDiff
using QuadGK
using SpecialFunctions

mutable struct Lattice
    # Configuration of the lattice (2D array of binary spins)
    grid::Matrix{Int64}
    L::Int64
    tau_values::Vector{Int64}  # Vector to store the tau values read from a CSV

    # Constructor to initialize the lattice with solved configuration and tau values
    function Lattice(L::Int64, csv_file::String = "../../core/average_tau_wide1.csv")
        grid = solved_configuration(L)
        tau_values = read_tau_values(csv_file)
        new(grid, L, tau_values)
    end
end

@inline function solved_configuration(L::Int64)
    # minimal energy state
    # all spins aligned
    # acessilbe assuming ergodic
    return fill(1,L,L)
end

function read_tau_values(csv_file::String)
    # Read the CSV file
    data =readdlm(csv_file, ',')
    return Int64.(ceil.(data[2, :]))
end

function random_configuration(L::Int, m::Float64)
    # random configuration of spins
    # weighted probability of spin up or down
    if m < -1 || m > 1
        throw(ArgumentError("m out of range"))
    end

    choices = [-1, 1]
    p = 0.5 * (m + 1.0)
    probabilities = [1.0 - p, p]
    spin_distribution = Categorical(probabilities)
    return [rand(spin_distribution) == 2 ? 1 : -1 for _ in 1:L, _ in 1:L]
end


function energy(lattice::Lattice)
    # for each lattice sites
    # calculate sum of nearest neighbours
    # by shifting lattice and summing
    L = lattice.L
    grid = lattice.grid
    total_energy = 0
    # vectorised energy calculation
    right = circshift(grid, (1, 0))
    left = circshift(grid, (-1, 0))
    up = circshift(grid, (0, 1))
    down = circshift(grid, (0, -1))
    total_energy = - sum(grid .* (right + left + up + down))
    return total_energy / (2 * L^2)
end


function generate_moves(lattice::Lattice, move::String, l::Int64 = 1)
    # generate candidate site for flip moves
    # based on move type
    # returns tuple of x and y coordinates
    L  = lattice.L
    if move == "single flip"
        # randomly select a site
        x = rand(1:lattice.L)
        y = rand(1:lattice.L)
        return ([x],[y])

    elseif move == "l chain flip"
        # randomly select an initial site
        x_0 = rand(1:lattice.L)
        y_0 = rand(1:lattice.L)
        # pick subsequent sites on x axis
        return ([mod1(x_0 + i, L) for i in 0:l-1], [y_0 for i in 0:l-1])


    elseif move == "l line flip"
        # randomly select a line
        y = rand(1:lattice.L)
        # within line randomly select l sites
        return (sample(1:L, l, replace=false), fill(y, l))
      
    elseif move == "unconstrained l flip"
        # randomly select l sites
        x_sites = sample(1:L, l, replace=false)
        y_sites = sample(1:L, l, replace=false)
        return (x_sites, y_sites)
        
    elseif move == "l square flip"
        # randomly select an initial site
        x0 = rand(1:lattice.L)
        y0 = rand(1:lattice.L)

        # select whether flip has a hole (l^2 - 1 flips) or not (l^2 flips)
        has_hole = rand(Bool)
        
        #generate indices for l x l square
        indices = [(mod1(x0 + i, lattice.L), mod1(y0 + j, lattice.L)) for i in 0:l-1, j in 0:l-1]

        #remove hole if needed
        if has_hole && !isempty(indices)
            hole_index = rand(1:length(indices))
            deleteat!(indices, hole_index)
        end
        xvals = [coord[1] for coord in indices]
        yvals = [coord[2] for coord in indices]
        return (xvals, yvals)
    else
        throw(ArgumentError("Invalid move type"))
    end
end

@inline function energy_change(lattice::Lattice, flips::Tuple{Vector{Int64}, Vector{Int64}})
    L = lattice.L
    
    # unpack proposed flips
    xvals, yvals = flips

    if length(xvals) != length(yvals)
        throw(ArgumentError("xvals and yvals must have the same length"))
    end

    # mask proposed flips
    # takes into account the fact that flipping two adjacent spins
    # doesn't change their interaction energy
    # ie +1,+1 has same energy as -1,-1
    # so by masking we acknowledge that the energy change is 0
    flip_set = Set(zip(xvals,yvals))
    total_energy_change = 0
    @inbounds for (x, y) in zip(xvals, yvals)
        local sum_neighbours = 0
        for (dx,dy) in ((1,0),(-1,0),(0,1),(0,-1))
            xp = mod1(x+dx, L)
            yp = mod1(y+dy, L)
            if (xp,yp) ∉ flip_set
                sum_neighbours += lattice.grid[xp,yp]
            end
        end
        total_energy_change += lattice.grid[x, y] * 2 * sum_neighbours
    end
    return total_energy_change
end

function do_flips(lattice::Lattice, flips::Tuple{Vector{Int64}, Vector{Int64}})
    # given list of candidate sites
    # flip the spins of these sites
    # same site might be flipped multiple times?
    xvals, yvals = flips
    @inbounds for (x,y) in zip(xvals, yvals)
        lattice.grid[x,y] *= -1
    end
end

function explore_moves(lattice::Lattice, l::Int64, move::String)
    # explore all possible moves from a given configuration
    #returns the saddle index ie the number of moves that decrease the energy
    moves = []
    # Generate moves based on the type specified
    if move == "single flip"
        moves = [([i], [j]) for i in 1:lattice.L, j in 1:lattice.L]   

    elseif move == "l chain flip"
        for y in 1:lattice.L
            for x in 1:(lattice.L - l + 1)
                x_coords = [(x + i - 1) % lattice.L + 1 for i in 1:l]
                y_coords = fill(y, l)
                push!(moves, (x_coords, y_coords))
            end
        end

    elseif move == "l line flip"
        for y in 1:lattice.L
            for comb in collect(combinations(1:lattice.L, l))  # Avoid lazy evaluation issues
                push!(moves, (collect(comb), fill(y, l)))
            end
        end
    
    elseif move == "unconstrained l flip"
        indices = [(i, j) for i in 1:lattice.L for j in 1:lattice.L]  # Flatten grid
        for comb in collect(combinations(indices, l))
            x_coords, y_coords = unzip(comb)  # Efficient tuple splitting
            push!(moves, (collect(x_coords), collect(y_coords)))
        end



    else
        throw(ArgumentError("Invalid move type"))
        return 0
    end

    saddles = count(move -> energy_change(lattice, move) < 0, moves)
    return saddles
    
end


function valid_l(L::Int64, l::Int64)
    # Check ergodicity condition for l,L
    # to share configuration space of l = 1
    # note does not check if l is even
    valid = false
    if (L - ceil(l/2)) % l == 0
        valid = true
    end
    if (L- floor(l/2)) % l == 0
        valid = true
    end
    return valid
end


function onsager_free_energy(b, L)
    # F = lnZ / beta
    # where Z is the partition function
    # set k_B = 1
    J  = 1
    # Define κ = 2 sinh(2K) / cosh(2K)^2
    k = 2 * sinh(2*b*J) / (cosh(2*b*J)^2)
    
    # Define the integrand:
    # ln[(1 + sqrt(1 - κ² sin²θ)) / 2]
    function integrand(theta)
        return log(1 + sqrt(1 - k^2 * sin(theta)^2))
    end
    
    # Compute the integral from 0 to π and include the 1/(2π) factor
    I, _ = quadgk(integrand, 0, pi/2)
    I /= pi
    
    # Standard form of Onsager's free energy:
    # -b * f = ln(2 cosh(2b)) + I
    # So, f = -1/b * [ln(2 cosh(2b)) + I]
    return (log(sqrt(2) * cosh(2*b*J)) + I) / (-b)
end


function onsager_energy(b, L)
    #calculate Onsager energy per spin
    J = 1
    k = 2 * sinh(2*b*J) / (cosh(2*b*J)^2)
    K_k = ellipk(k^2)

    return -J * coth(2 * b * J) * (1+2/pi * K_k * (2 * tanh(2*b*J)^2 -1))
    
    
end

function onsager_specific_heat(b, L)  
    # Onsager specific heat per spin: C = b² * dE/db
    J = 1
    k = 2 * sinh(2*b*J) / (cosh(2*b*J)^2)
    K_k = ellipk(k^2)
    E_k = ellipe(k^2)
    k_prime = 2 * tanh(2*b*J)^2 -1
    return (b*J*coth(2*b*J))^2 * 2/pi * (2 * K_k - 2 * E_k - (1 - k_prime)*(pi/2 + k_prime * K_k))


end

function onsager_entropy(b, L)
    # Onsager entropy per spin: S = b * (E - f)
    U = onsager_energy(b, L)
    F = onsager_free_energy(b, L)
    S = b * (U - F)
    return S
end