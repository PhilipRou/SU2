include("test_head.jl")


function delta_S_gauge_hex_test()
    # NX = N_x>>1
    test_field_1 = hexfield_SU2(N_x, N_t, true)
    test_field_2 = deepcopy(test_field_1)
    # μ = rand(1:2)
    # x = rand(1:N_x)
    # t = rand(1:N_t)
    μ, x, t = rand(chess_hex_link_coords(N_x, N_t))
    test_field_2[μ,x,t] = ran_SU2(rand())
    true_dif =  action_hex(test_field_2,β)-action_hex(test_field_1,β)
    test_dif = delta_S_gauge_hex(test_field_1, μ,x,t, test_field_1[μ,x,t], test_field_2[μ,x,t],β)
    @assert isapprox(true_dif, test_dif) "Something went wrong while testing delta_S_gauge_hex"
    # print("μ = ", μ,", t = ",t,", x = ",x,",    true_dif = ", true_dif, ",    test_dif = ", test_dif)
    return true
end

delta_S_gauge_hex_test()

#=
for i = 1:100
    push!(actions, action_hex(test_hex))
    chess_metro_hex!(test_hex,0.2)
end

plot(actions)
=#

function overrelax_hex!_test()
    test_field = hexfield_SU2(N_x,N_t,true)
    old_field = deepcopy(test_field)
    old_action = action_hex(test_field,β)
    # for i = 1:10
    #     t = rand(1:8)
    #     x = rand(1:8)
    #     μ = rand(1:3)
    #     overrelax_hex!(test_field,μ,x,t)
    # end
    coords = rand(chess_hex_link_coords(N_x, N_t), 10)
    for coord in coords
        µ, x, t = coord
        overrelax_hex!(test_field,µ,x,t)
    end
    @assert old_field != test_field "There was no update in overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action_hex(test_field,β)) "overrelax! didn't leave the action invariant"
    return true
end

overrelax_hex!_test()