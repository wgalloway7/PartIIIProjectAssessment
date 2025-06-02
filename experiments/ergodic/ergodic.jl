using Plots 
include("../../core/lattice.jl")
include("../../core/montecarlo.jl")
Lmax = 20
L_l = zeros(Int64, Lmax, Lmax - 1)

for L in 1:Lmax
    for l in 1:(L-1)
        if valid_l(L,l) && l %2 == 1
            L_l[L,l] = 1
        end
    end
end


using Plots

x = Int[]
y = Int[]
for i in 1:size(L_l, 1) 
    for j in 1:size(L_l, 2)  
        if L_l[i, j] != 0
            push!(x, i)
            push!(y, j) 
        end
    end
end

p = scatter(x, y, marker = :x, markersize = 2, label = "")
Lvalues = [i for i in 1:Lmax]
colours = [:red, :blue, :green, :orange, :purple, :cyan]
for x in 1:6
    plot!(Lvalues, Lvalues .* 2 /(2*x + 1) .+ 1/(2x + 1), label = "n = $x", color = colours[x])
    plot!(Lvalues, Lvalues .* 2 /(2*x + 1) .- 1/(2x + 1), label = "", color = colours[x])
end



xticks_all = minimum(Lvalues):maximum(Lvalues)
yticks_all = floor(Int, minimum(y)):ceil(Int, maximum(y))


xticks_labels = collect(minimum(Lvalues):2:maximum(Lvalues))
yticks_labels = collect(floor(Int, minimum(y)):2:ceil(Int, maximum(y)))


xticks_tuple = (xticks_all, [v in xticks_labels ? string(v) : "" for v in xticks_all])
yticks_tuple = (yticks_all, [v in yticks_labels ? string(v) : "" for v in yticks_all])


plot!(
    xlabel = "L",
    ylabel = "l",
    xticks = xticks_tuple,
    yticks = yticks_tuple,
    grid = :both 
)

savefig(p, "L_l_scatter.png")