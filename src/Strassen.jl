module Strassen
export strass, naive_multiplication, benchmarkcontrol, benchmarkstrass

say_hello() = println("Hello revised!")

function naive_multiplication(A, B, n)
    C = fill(0.0, (n, n))
        for i in 1 : n
                for j in 1 : n
                        for k in 1 : n
                                C[i, j] += A[i, k] * B[k, j]
                        end
                end
        end
        return C
    end

function strass(A, B, n, n_min)
        if n <= n_min
                return naive_multiplication(A, B, n)
        else
                m = div(n, 2)
                u = 1 : m
                v = m + 1 : n

                P_1 = strass(
                        A[u, u] + A[v, v], B[u, u] + B[v, v], m, n_min
                )
                P_2 = strass(A[v, u] + A[v, v], B[u, u], m, n_min)
                P_3 = strass(A[u, u], B[u, v] - B[v, v], m, n_min)
                P_4 = strass(A[v, v], B[v, u] - B[u, u], m, n_min)
                P_5 = strass(A[u, u] + A[u, v], B[v, v], m, n_min)
                P_6 = strass(A[v, u] - A[u, u], B[u, u] + B[u, v], m, n_min)
                P_7 = strass(A[u, v] - A[v, v], B[v, u] + B[v, v], m, n_min)

                C = fill(0.0, (n, n))
                C[u, u] = P_1 + P_4 - P_5 + P_7
                C[u, v] = P_3 + P_5
                C[v, u] = P_2 + P_4
                C[v, v] = P_1 + P_3 -P_2 + P_6
                return C
        end
    end

"""
Performs benchmark for every leaf size between d_min and d_max.
"""
function benchmarkstrass(A, B, d_min, d_max)::Tuple{Vector{Int}, Vector{Number}}
    leaf_perfs = []
    leaf_sizes = []
    n = size(A)[1]
    for d in range(d_min, d_max)
        n_min = 2 ^ d
        strass_results = @timed strass(A, B, n, n_min);

        # Debug printing for time and exactness
        #println("strass, leaf size = $(2^d): ", strass_results[2], "s")
        #println("isequal: ", isapprox(strass_results[1], ref, atol=1e-5))

        push!(leaf_perfs, strass_results[2])
        push!(leaf_sizes, d)
    end
    return leaf_sizes, leaf_perfs
end

"""
Performs benchmark for a given leaf size.
"""
function benchmarkstrass(A, B, d)::Float64
    n = size(A)[1]
    n_min = 2 ^ d
    strass_results = @timed strass(A, B, n, 2^d);

    # Debug printing for time and exactness
    #println("strass, leaf size = $(2^d): ", strass_results[2], "s")
    #println("isequal: ", isapprox(strass_results[1], ref, atol=1e-5))

    # return time only
    return strass_results[2]
end

function benchmarkcontrol(A, B)
    n = size(A)[1]

    ref, t_ref = @timed A * B;
    naive, t_naive = @timed naive_multiplication(A, B, n);

    ### Debug printing
    #=
    println("n = $n")
    println("built-in matmul: ", t_ref, "s")
    println()
    println("naive: ", t_naive)
    println("isequal: ", isapprox(naive, ref, atol=1e-5))
    println()
    println()
    =#
    return t_naive, t_ref
end


end #module
