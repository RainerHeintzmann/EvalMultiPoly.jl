export get_identity_poly_coeffs
export get_identity_multipoly_coeffs

"""
    get_poly_exponents(::Val{numvars}, ::Val{N})

Return the exponent tuple corresponding to every coefficient of an
N-th order polynomial with `numvars` variables.

For example:

    get_poly_exponents(Val(2), Val(3))

returns

    [(0,0),
     (1,0),
     (2,0),
     (3,0),
     (0,1),
     (1,1),
     (2,1),
     (0,2),
     (1,2),
     (0,3)]

corresponding to

    1, x, x², x³, y, xy, x²y, y², xy², y³
"""
function get_poly_exponents(::Val{numvars}, ::Val{N}) where {numvars, N}
    _poly_exponents(numvars, N, numvars)
end


function _poly_exponents(numvars, N, activevars)
    z = ntuple(_ -> 0, numvars)

    N == 0 && return [z]

    exps = typeof(z)[z]

    for n in 1:activevars
        for e in _poly_exponents(numvars, N - 1, n)
            push!(exps, Base.setindex(e, e[n] + 1, n))
        end
    end

    return exps
end

"""
    get_identity_poly_coeffs(::Val{numvars}, ::Val{N}, variable; T=Float64)

Return coefficients for the polynomial

    f(x₁, ..., xₙ) = x_variable
"""
function get_identity_poly_coeffs(
    ::Val{numvars},
    ::Val{N},
    variable::Integer;
    T=Float32,
) where {numvars, N}

    N >= 1 || throw(ArgumentError(
        "Identity transform requires polynomial order >= 1"
    ))

    1 <= variable <= numvars || throw(ArgumentError(
        "variable must be between 1 and $numvars"
    ))

    exps = get_poly_exponents(Val(numvars), Val(N))

    target = ntuple(i -> i == variable ? 1 : 0, numvars)

    return Tuple(
        e == target ? one(T) : zero(T)
        for e in exps
    )
end

"""
    get_identity_multipoly_coeffs(::Val{numvars}, ::Val{N}; T=Float64)

Return coefficients for the identity transformation

    (x₁, x₂, ..., xₙ) -> (x₁, x₂, ..., xₙ)
"""
function get_identity_multipoly_coeffs(
    ::Val{numvars},
    ::Val{N};
    T=Float32,
) where {numvars, N}

    return Tuple(
        c
        for variable in 1:numvars
        for c in get_identity_poly_coeffs(
            Val(numvars),
            Val(N),
            variable;
            T=T,
        )
    )
end