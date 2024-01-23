include("test_head.jl")


#=
function loop_mat_hex_test()
    test_field = hexfield_SU2(N_x, N_t, true)
    old_mat = loop_mat_hex(U,rand(1:N_x),rand(1:N_t))
    new_mat = loop_mat_hex(temp_gauge_hex(U),rand(1:N_x),rand(1:N_t))
    x, t = rand(half_chess_coords(N_x,N_t))
    @assert isapprox(old, new_loop) "rhomb_half_loop() not yet correct at x = $x, t = $t"
    return true
end
=#


function RT_loop_hex_test()
    test_field = hexfield_SU2(N_x, N_t, true)
    x, t = rand(half_chess_coords(N_x,N_t))
    R = rand(1:N_x)
    T = rand(1:N_t)
    old_loop = tr(RT_loop_hex(test_field,R,T,x,t))
    new_loop = tr(RT_loop_hex(temp_gauge(test_field),R,T,x,t))
    @assert isapprox(old_loop, new_loop) "rhomb_half_loop() not yet correct at x = $x, t = $t"
    return true
end

RT_loop_hex_test()

function rhomb_half_loop_test()
    test_field = hexfield_SU2(N_x, N_t, true)
    x, t = rand(half_chess_coords(N_x,N_t))
    old_loop = tr(rhomb_half_loop(test_field,x,t))
    new_loop = tr(rhomb_half_loop(temp_gauge(test_field),x,t))
    @assert isapprox(old_loop, new_loop) "rhomb_half_loop() not yet correct at x = $x, t = $t"
    return true
end

rhomb_half_loop_test()

function edge_loop_hex_test()
    test_field = hexfield_SU2(N_x, N_t, true)
    x, t = rand(half_chess_coords(N_x,N_t))
    old_loop = tr(edge_loop_hex(test_field,x,t))
    new_loop = tr(edge_loop_hex(temp_gauge(test_field),x,t))
    @assert isapprox(old_loop, new_loop) "edge_loop_hex() not yet correct at x = $x, t = $t"
    return true
end

edge_loop_hex_test()

function rhomb_loop_test()
    test_field = hexfield_SU2(N_x, N_t, true)
    x, t = rand(half_chess_coords(N_x,N_t))
    old_loop = tr(rhomb_loop(test_field,x,t))
    new_loop = tr(rhomb_loop(temp_gauge(test_field),x,t))
    @assert isapprox(old_loop, new_loop) "rhomb_loop() not yet correct at x = $x, t = $t"
    return true
end

rhomb_loop_test()

function L_loop_hex_test()
    test_field = hexfield_SU2(N_x, N_t, true)
    x, t = rand(half_chess_coords(N_x,N_t))
    old_loop = tr(L_loop_hex(test_field,x,t))
    new_loop = tr(L_loop_hex(temp_gauge(test_field),x,t))
    @assert isapprox(old_loop, new_loop) "L_loop_hex() not yet correct at x = $x, t = $t"
    return true
end

L_loop_hex_test()


