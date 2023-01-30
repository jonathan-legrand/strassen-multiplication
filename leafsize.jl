using LinearAlgebra, Plots, StatsPlots
using Strassen

# Range of q for input matrix size
Q_MAX = 10
Q_MIN = 5

y = []
sizes = []
for q in range(Q_MIN, Q_MAX)

    # Matrices must have power of two size
    local n = 2 ^ q

    ### Randomly ggnerate matrices
    local A = rand(n, n);
    local B = rand(n, n);
    local perf
    
    # Measure perf for different leaf sizes
    _, perf = benchmarkstrass(A, B, 0, q)
    
    # Pad arrays to have consistent Q_MAX length
    for _ in range(q + 1, Q_MAX)
        push!(perf, 0)
    end

    push!(y, perf)
    push!(sizes, Symbol(n))
end

# Leaf size vector
x = [d for d in range(0, Q_MAX)]

perfmat = reduce(hcat, y)
xexpanded = reshape(x, size(x)[1], 1)

# We use logarithmic time to enhance readability
timeframe = DataFrame(hcat(xexpanded, map(log, perfmat)), :auto)
colnames = append!([Symbol("leaf")], sizes)
rename!(timeframe, colnames)
@df timeframe plot(
                   :leaf,
                   #map(log, cols(2: Q_MAX - Q_MIN + 2)),
                   cols(2: Q_MAX - Q_MIN + 2),
                   legend=:bottomright,
                   title="Strass performance",
                   xlabel="Leaf size (d = log_2(n_min))",
                   ylabel="Log execution time",
                   legend_title="Matrix size (n)",
                   markershape = :auto
                  )

# Output
savefig("leafstrass-$Q_MAX-$Q_MIN.png")
optimal_leaf = x[argmin(sum(perfmat, dims=2))]
println("Optimal leaf size = $optimal_leaf")
