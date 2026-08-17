function test_poly_allocations()
    # get_num_poly_vars(Val(2), Val(3)) # 10 indices
    # p = get_polynomial(Val(2), Val(3))  
    # @time p.(Tuple.(CartesianIndices((200,200))),Ref((1.1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27)));

    evalpoly(2, (1, 2, 3))
    evalmultipoly(Val(2), (2,), (1, 2, 3))

    evalmultipoly(Val(2), (2, 20),  (1, 1, 1, 1, 1, 1)) == 467

    cids = Tuple.(CartesianIndices((200,200)))

    cidsx = [ci[1] for ci in cids]
    res = zeros(Float32, 200,200)
    @time res .= evalpoly.(cidsx, Ref((1,2,3)));

    # cfds = map((t)->Tuple(Float32.([t...])), cids) 
    cs = Tuple(Float32.(collect(1:27)))
    get_num_poly_vars(Val(2), Val(2)) # 6 indices required
    @time res .= evalmultipoly.(Ref(Val(2)), cids, Ref(cs)); # 2 orders, two variables
    # 0.000128 seconds (3 allocations: 168 bytes)

    get_num_poly_vars(Val(3), Val(2)) # 10 indices required
    @time res .= evalmultipoly.(Ref(Val(3)), cids, Ref(cs)); # 3 orders, two variables
    #  0.000163 seconds (3 allocations: 168 bytes)

    get_num_poly_vars(Val(4), Val(2)) # 15 indices required
    @time res .= evalmultipoly.(Ref(Val(4)), cids, Ref(cs)); # 2 orders, two variables
    #  0.000200 seconds (3 allocations: 168 bytes)

    get_num_poly_vars(Val(5), Val(2)) # 21 indices required
    @time res .= evalmultipoly.(Ref(Val(5)), cids, Ref(cs)); # 2 orders, two variables
    # 0.000291 seconds (3 allocations: 168 bytes)

    @time evalmultipoly.(Ref(Val(0)), cids, Ref(cs));
    # 0.000038 seconds (5 allocations: 156.461 KiB)

end

# Univariate example of the polynomials evaluation
univariate_polynomial_coeffs = get_num_poly_vars(Val(1), Val(3))
# 4
# 1 + x + x^2 + x^3

# Multivariate example of the polynomials evaluation
multivariate_polynomial_coeffs = get_num_poly_vars(Val(2), Val(3))
# 10
# [1, x, x², x³, y, xy, x²y, y², xy², y³]
# coeffs = (
#     c_const,
#     c_x,
#     c_x2,
#     c_x3,
#     c_y,
#     c_xy,
#     c_x2y,
#     c_y2,
#     c_xy2,
#     c_y3,
# )

# a test of the coeffs ordering:
x = 2.0
y = 3.0

c = (
    1.0,      # 1
    10.0,     # x
    100.0,    # x^2
    1000.0,   # x^3
    10000.0,  # y
    100000.0, # xy
    1e6,      # x^2 y
    1e7,      # y^2
    1e8,      # x y^2
    1e9       # y^3
)

r_package = evalmultipoly(Val(3), (x, y), c)

r_manual =
    c[1] +
    c[2]  * x +
    c[3]  * x^2 +
    c[4]  * x^3 +
    c[5]  * y +
    c[6]  * x*y +
    c[7]  * x^2*y +
    c[8]  * y^2 +
    c[9]  * x*y^2 +
    c[10] * y^3

@show r_package
@show r_manual
@show r_package == r_manual


# Multivariate example of the Multi-polynomials evaluation
multivariate_polynomial_coeffs = get_num_multipoly_vars(Val(2), Val(3))
# 20
# x′ = 1 + x + x^2 + x^3 + y + y^2 + y^3 + xy + xy^2 + yx^2
# y′ = 1 + x + x^2 + x^3 + y + y^2 + y^3 + xy + xy^2 + yx^2

# an example of how to work with the identity transforms:

cids = Tuple.(CartesianIndices((100, 100)))

numvars = length(cids[1])   # 2
n_order = 6

cs = get_identity_multipoly_coeffs(
    Val(numvars),
    Val(n_order)
)

h = get_multi_poly(
    Val(numvars),   # first = number of variables
    Val(n_order)    # second = polynomial order
)

# @track_allocs  h.(cids, Ref(cs));

# benchmarking
cids_res = similar(cids, NTuple{numvars, Float32});
rf = Ref(cs);

@b cids_res .= h.(cids, rf)
@track_allocs cids_res .= h.(cids, rf)


results = map((1, 10, 100, 1_000, 10_000)) do N
    c = Tuple.(CartesianIndices((N, 1)))
    rf = Ref(cs)
    out = similar(c, NTuple{2,Float32})

    N => (@b broadcast!($h, $out, $c, $rf))
end;
foreach(println, results)
# julia> foreach(println, results)
# 1 => Sample(evals=631, time=4.3106180665610143e-8)
# 10 => Sample(evals=72, time=3.8194444444444445e-7)
# 100 => Sample(evals=7, time=3.7571428571428575e-6)
# 1000 => Sample(time=3.77e-5)
# 10000 => Sample(time=0.0003818)


for N in (1, 10, 100, 1_000, 10_000, 100_000)
    c = Tuple.(CartesianIndices((N, 1)))
    rf = Ref(cs)
    out = similar(c, NTuple{2,Float32})

    alloc = @ballocated broadcast!($h, $out, $c, $rf)
    println("N = $N: $alloc bytes")
end
# N = 1: 0 bytes
# N = 10: 0 bytes
# N = 100: 0 bytes
# N = 1000: 0 bytes
# N = 10000: 0 bytes
# N = 100000: 0 bytes

