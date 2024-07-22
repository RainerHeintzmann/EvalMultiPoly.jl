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
