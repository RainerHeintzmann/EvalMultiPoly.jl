using EvalMultiPoly
using BenchmarkTools

function main()
    numvars = 2
    n_order = 4

    cs = get_identity_multipoly_coeffs(
        Val(numvars),
        Val(n_order),
    )

    h = get_multi_poly(
        Val(numvars),
        Val(n_order),
    )

    cids = Tuple.(CartesianIndices((100_000, 1)))
    out = similar(cids, NTuple{numvars, Float32})
    rf = Ref(cs)

    # Warm-up
    broadcast!(h, out, cids, rf)

    @btime broadcast!($h, $out, $cids, $rf)

    return nothing
end

function test_allocations()
    cs = get_identity_multipoly_coeffs(Val(2), Val(4))
    h = get_multi_poly(Val(2), Val(4))

    cids = Tuple.(CartesianIndices((100_000, 1)))
    out = similar(cids, NTuple{2, Float32})
    rf = Ref(cs)

    broadcast!(h, out, cids, rf)  # compile

    return @allocated broadcast!(h, out, cids, rf)
end

println(test_allocations())

main()

# 0
#   1.884 ms (0 allocations: 0 bytes)