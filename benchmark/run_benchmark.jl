using EvalMultiPoly
using BenchmarkTools
using Printf


# Benchmark configuration
numvars = 2
n_order = 4
Ns = (1, 10, 100, 1_000, 10_000, 100_000)

# Construct polynomial coefficients and evaluator
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

    # Warm up this exact specialization before benchmarking.
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