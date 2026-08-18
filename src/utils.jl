export split_tuple

"""
    split_tuple(t::NTuple{S,T},::Val{numvars}) where {S,T,numvars}

Split a tuple into `numvars` parts packed into a tuple of tuples. The tuple `t` is assumed to have a length that is a multiple of `numvars`.

Example:
```juliadoc
julia> t = (1, 2, 3, 4, 5, 6)

julia> split_tuple(t, Val(2))
((1, 2, 3), (4, 5, 6))

julia> split_tuple(t, Val(3))
((1, 2), (3, 4), (5, 6))
```
"""
@generated function split_tuple(
    t::NTuple{S,T},
    ::Val{numvars}
) where {S,T,numvars}

    S % numvars == 0 ||
        error("Tuple length $S is not divisible by $numvars")
    L = S ÷ numvars

    chunks = [
        Expr(
            :tuple,
            [
                :(getfield(t, $i))
                for i in ((j - 1) * L + 1):(j * L)
            ]...
        )
        for j in 1:numvars
    ]
    return Expr(:tuple, chunks...)
end

# Small util function for the "tail" of a tuple
tail_new(t) = (ntuple(i -> t[i+1], Val(length(t)-1)))
