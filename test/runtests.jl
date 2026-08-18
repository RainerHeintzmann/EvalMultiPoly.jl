using Test
using EvalMultiPoly

function multipoly_allocations(N=1_000)
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

    cids = Tuple.(CartesianIndices((N, 1)))
    out = similar(cids, NTuple{numvars, Float32})
    rf = Ref(cs)

    # Compile before measuring.
    broadcast!(h, out, cids, rf)

    alloc = @allocated broadcast!(h, out, cids, rf)

    return alloc, out
end

@testset "Allocation-free multivariate polynomial evaluation" begin
    alloc, out = multipoly_allocations()

    @test alloc == 0
    @test out[1] == (1.0f0, 1.0f0)
    @test out[end] == (1000.0f0, 1.0f0)
end