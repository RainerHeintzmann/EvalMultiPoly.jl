"""
    polynomial(::Val{0}, ::T1,  c::T2, ::Val{numvars}=Val(1)) where {NV, TS, T1 <: NTuple{NV}, T2 <: NTuple{TS}, numvars}

Create a polynomial of order 0 with numvars variables.

returned is a function that takes a tuple of variables and a tuple of coefficients and returns the value of the polynomial
    and the number of coefficients required (here 1).
"""
function _polynomial(::Val{0}, ::T1,  c::T2, ::Val{numvars}=Val(1)) where {NV, TS, T1 <: NTuple{NV}, T2 <: NTuple{TS}, numvars}
    # println("c: $(c) $(length(c))");
    return c[1], Base.tail(c)
end  #, (t,c) -> ntuple(n->c[1], Val(numvars))

"""
    polynomial(::Val{N}, t::T1, c::T2, ::Val{numvars}=Val(length(t))) where {N, NV, TS, T1 <: NTuple{NV}, T2 <: NTuple{TS}, numvars}

Represents a polynomial of order N with numvars variables (also implicitely defined via the length of the NTuple `t`).
Note that `numvars` is needed for the internal workings of the polynomial generator, but notmally not by the user.

E.g. to represent a polynomial of order 1 with 2 variables, the coefficients are ordered as follows:
c = (c0, c1, c2) where the polynomial is: c0 + c1*x + c2*y
or for a polynomial of order 2 with 2 variables:
c = (c0, c1, c2, c3, c4, c5) where the polynomial is c0 + c1*x + c2*x^2 + c3*y + c4*x*y + c5*y^2
Note that the coefficients are ordered not by the multiples in which they appear in the polynomial, but by the order of the variables.

Example:
```jldoctest
>julia polynomial(Val(2), (2, 20),  (1f0, 1f0, 1f0, 1f0, 1f0, 1f0))
467.0f0
>julia cs = Tuple(Float32.(collect(1:27)))
  (1.0f0, 2.0f0, 3.0f0, 4.0f0, 5.0f0, 6.0f0, 7.0f0, 8.0f0, 9.0f0, 10.0f0, 11.0f0, 12.0f0, 13.0f0, 14.0f0, 15.0f0, 16.0f0, 17.0f0, 18.0f0, 19.0f0, 20.0f0, 21.0f0, 22.0f0, 23.0f0, 24.0f0, 25.0f0, 26.0f0, 27.0f0)
>julia res = zeros(Float32, 200,200)
>julia @time res .= polynomial.(Ref(Val(2)), Tuple.(CartesianIndices((200,200))), Ref(cs)); 
  0.054017 seconds (147.45 k allocations: 10.168 MiB, 99.75% compilation time)
>julia @time res .= polynomial.(Ref(Val(2)), Tuple.(CartesianIndices((200,200))), Ref(cs)); 
  0.000089 seconds (3 allocations: 184 bytes)
```
"""
@generated function _polynomial(::Val{N}, t::T1, c::T2, ::Val{numvars}=Val(length(t))) where {N, NV, TS, T1 <: NTuple{NV}, T2 <: NTuple{TS}, numvars} 
    quote
        res = res = c[1]
        c = Base.tail(c)
        Base.Cartesian.@nexprs $numvars n -> begin
            p, c = _polynomial(Val(N-1), t, c, Val(n))
            res += t[n] * p
        end
        return res, c
    end
end

function evalmultipoly(::Val{N}, t::T1, c::T2, ::Val{numvars}=Val(length(t)))::RT where {N, NV, TS, RT, T1 <: NTuple{NV}, T2 <: NTuple{TS, RT}, numvars} 
    _polynomial(Val(N), t, c, Val(numvars))[1]
end

function get_multi_poly(::Val{numvars}, ::Val{N}; verbose=false) where {numvars, N}
    if (verbose)
        @info "Creating polynomials with $(numvars) variables of order $(N). Required constants: $(get_num_multipoly_vars(Val(numvars), Val(N)))"
    end
    # this defines one dimension of the multivariate polynomial
    # p = (t,c) -> evalmultipoly(Val(N), NTuple(t), c)

    function mpol(t, c::NTuple{M, RT})::NTuple{numvars, RT} where {M, RT}
        return ntuple(n->evalmultipoly(Val(N), Tuple(t), split_tuple(c, Val(numvars))[n]), Val(numvars))
        # return ntuple(n->p(t, split_tuple(c, Val(numvars))[n]), Val(numvars))
    end

    function mpol(t, c::NTuple{M, RT}, n)::RT where {M, RT}
        return evalmultipoly(Val(N), Tuple(t), split_tuple(c, Val(numvars))[n])
        # return p(t, split_tuple(c, Val(numvars))[n])
    end
    return mpol 
end

"""
    get_num_poly_vars(::Val{numvars}, ::Val{N}) where {numvars, N}

Get the number of coefficients required for a single polynomial of order N with numvars variables.
"""
function get_num_poly_vars(::Val{numvars}, ::Val{N}) where {numvars, N}
    return binomial(numvars+N, N) 
end

"""
    get_num_multipoly_vars(::Val{numvars}, ::Val{N}) where {numvars, N}

Get the number of coefficients required for a multipolynomial of order N with numvars variables.
"""
function get_num_multipoly_vars(::Val{numvars}, ::Val{N}) where {numvars, N}
    return numvars*get_num_poly_vars(Val(numvars), Val(N))
end
