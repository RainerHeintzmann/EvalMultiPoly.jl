using EvalMultiPoly
using BenchmarkTools
using Printf

function run_benchmark()
    numvars = 2
    n_order = 4

    Ns = (1, 10, 100, 1_000, 10_000, 100_000)

    cs = get_identity_multipoly_coeffs(
        Val(numvars),
        Val(n_order),
    )

    h = get_multi_poly(
        Val(numvars),
        Val(n_order),
    )

    println()
    println("EvalMultiPoly benchmark")
    println("numvars = $numvars, polynomial order = $n_order")
    println()

    println("--------------------------------------------------------------------------------")
    @printf(
        "%10s %18s %18s %12s %12s\n",
        "N",
        "total [ns]",
        "ns / element",
        "allocs",
        "bytes",
    )
    println("--------------------------------------------------------------------------------")

    for N in Ns
        cids = Tuple.(CartesianIndices((N, 1)))
        out = similar(cids, NTuple{numvars, Float32})
        rf = Ref(cs)

        # Warm-up
        broadcast!(h, out, cids, rf)

        trial = @benchmark broadcast!($h, $out, $cids, $rf)
        est = median(trial)

        @printf(
            "%10d %18.2f %18.4f %12d %12d\n",
            N,
            est.time,
            est.time / N,
            est.allocs,
            est.memory,
        )
    end

    println("--------------------------------------------------------------------------------")
end

run_benchmark()

# EvalMultiPoly benchmark
# numvars = 2, polynomial order = 4

# --------------------------------------------------------------------------------
#          N         total [ns]       ns / element       allocs        bytes
# --------------------------------------------------------------------------------
#          1              23.59            23.5944            0            0
#         10             147.43            14.7431            0            0
#        100            1340.00            13.4000            0            0
#       1000           13300.00            13.3000            0            0
#      10000          136800.00            13.6800            0            0
#     100000         1406750.00            14.0675            0            0
# --------------------------------------------------------------------------------