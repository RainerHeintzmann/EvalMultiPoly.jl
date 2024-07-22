# EvalMultiPoly.jl
Evaluates a multivariate polynomial of defined order and number of variables attempting minimal memory consumption.

It is basically a multivariate version of the `evalpoly()` function in `Base`. 
It is programmed using a recursive scheme of generating polynomial expressions. 
This means you probably should NOT use it for large degrees of polynomials.

The main point is that (apart from a fixed overhead depending on the order) it is allocation-free.
This yields substantial speed advantages, when used in broadcasting expressions.

The package `DataToFunctions.jl` using this package to define coordinate-transforms of higher order.
Here is an example demonstrating the equivalence to `evalpoly()`

```julia
    julia> p = evalpoly(2, (1, 2, 3))
    julia> p = evalmultipoly(Val(2), (2,), (1, 2, 3))
```

, where a polynomial of order two is created `p = 1 + 2x + 3x²`.
In `evalmultipoly` the first argument is the order of the multivariate polynomial, the second argument is the value of each variable as a tuple, and the third argument is the tuple of coefficients.
In both cases we obtain the result: `17`.

We can apply such a polynomial to an array as follows:

```julia
julia> sz = (200,200)
julia> cids = Tuple.(CartesianIndices((200,200)))

julia> cidsx = [ci[1] for ci in cids]
julia> res = zeros(Float32, sz...)
julia> res .= evalpoly.(cidsx, Ref((1, 2, 3)));

julia> res .= evalmultipoly.(Ref(Val(2)), cids, Ref((1, 2, 3, 0, 0, 0))); 
```

Note that the last line operates on the full two-dimensional coordinate array as tuples.
The corresponding twodimensional multivariate polynomial includes 6 coefficients:
`p = c1 + c2 x + c3 x² + c4 y + c5 x y + c6 y²`
Some of which were set to zero to yield the same result as the one-dimensional polynomial.
The number of required coefficients for a single such expression can be obtained by
`get_num_poly_vars(Val(2), Val(2))`.

The generator `get_multi_poly(::Val{numvars}, ::Val{N})` returns a function which can be applied to
an array of cartesian indices (or `NTuples`) and returns an array of tuples.

```julia
julia> f = get_multi_poly(Val(2), Val(2); verbose=true)
julia> nc = get_num_multipoly_vars(Val(2), Val(2))
julia> cids = Tuple.(CartesianIndices((200,200)))
julia> cs = Tuple(rand(Float32, nc))
julia> f.(cids, Ref(cs))
```

Note that here the type of the indices will determine the type of the final array of tuples.
If you want individual access to the result dimensions, you can use `f.(cids, Ref(cs), n)` with the dimension `n`.
