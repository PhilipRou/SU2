include("test_head.jl")


#=
# Not in use anymore, after some patch it became inefficient for large N_t
function mod1_homebrew(t,N_t)
    return (t-1)%N_t +1
end

function mod1_homebrew_test()
    N_t = 4
    for i = 1:12
        @assert mod1(i,N_t) == mod1_homebrew(i,N_t) "Homebrew mod1 does not seem to work"
    end
    return true
end

mod1_homebrew_test()
=#

# function test_mywrite()
#     rm("D:\\Physik Uni\\julia_projects\\SU2\\test_data\\test_mywrite.txt")
#     for i = 1:10
#         mywrite("D:\\Physik Uni\\julia_projects\\SU2\\test_data\\test_mywrite.txt",i)
#     end
#     bla = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\test_data\\test_mywrite.txt")
#     @assert Int.(bla[:,1]) == collect(1:10)
#     return true
# end

# test_mywrite()

#=
function mywrite_last_conf_test()
    V = gaugefield_SU2(N_t,N_x,true)
    mywrite("D:\\Physik Uni\\julia_projects\\SU2\\data\\test_data\\mywrite_config_test.txt", V)
    W = read_last_conf("D:\\Physik Uni\\julia_projects\\SU2\\data\\test_data\\mywrite_config_test.txt")
    for μ = 1:2
        for t = 1:N_t
            for x = 1:N_x
                @assert isapprox(V[μ,t,x], W[μ,t,x])
            end
        end
    end
    return true
end


mywrite_last_conf_test()
=#

function plaq_test()
    test_field = gaugefield_SU2(N_x,N_t,true)
    plaquette = coeffs2grp(plaq(test_field, rand(1:N_x), rand(1:N_t)))
    @assert is_SU2(plaquette) "Something went wrong while testing plaq"
    return true
end

plaq_test()

#=
function gen_rect_test()
    U = gaugefield_SU2(N_t,N_x,true)
    @assert isapprox(plaq(U,1,1) , gen_rect(U,1,1,1,1)) "Something wrong with gen_rect"
    return true
end

gen_rect_test()
=#

#=
function loop_mat_test()
    l_t = rand(1:N_t)
    l_x = rand(1:N_x)
    U = gaugefield_SU2(N_t,N_x,true)
    to_be_tested = loop_mat(U,l_t,l_x)
    mat = Matrix{coeffs_SU2}(undef,N_t,N_x)
    for i = 1:8
        for j = 1:8
            mat[i,j] = gen_rect(U,i,j,l_t,l_x)
        end
    end
    for i = 1:8
        for j = 1:8
            @assert isapprox(mat[i,j], to_be_tested[i,j]) "Something went wrong while testing loop_mat"
        end
    end
    return true
end

loop_mat_test()



function measure_RT_loops_test()
    N_t = N_x = 32
    test_field = gaugefield_SU2(N_t, N_x, true)
    T1 = rand(1:N_t)
    T2 = rand(1:N_t)
    R1 = rand(1:N_t)
    R2 = rand(1:N_t)
    loops = [[T1,R1], [T2,R2]]
    @assert isapprox(measure_RT_loops_corrs(test_field, loops, 0, 0.0)[2], measure_RT_loops(test_field, loops, 0, 0.0))
    return true
end

measure_RT_loops_test()

=#

function loop_mat_test()
    test_field = gaugefield_SU2(N_x, N_t, true)
    R = rand(1:N_x)
    T = rand(1:N_t)
    x = rand(1:N_x)
    t = rand(1:N_t)
    old_loop = tr(loop_mat(test_field,R,T,)[x,t])
    new_loop = tr(loop_mat(temp_gauge(test_field),R,T)[x,t])
    @assert isapprox(old_loop, new_loop) "loop_mat() not yet correct at x = $x, t = $t"
    return true
end

loop_mat_test()

function loop_2x3_square_test()
    test_field = gaugefield_SU2(N_x, N_t, true)
    x, t = [rand(1:N_x), rand(1:N_t)]
    @assert isapprox(tr(loop_2x3_square(test_field,x,t)), tr(loop_2x3_square(temp_gauge(test_field),x,t)))
end

loop_2x3_square_test()

function RT_loop_test()
    test_field = gaugefield_SU2(N_x, N_t, true)
    x, t = [rand(1:N_x), rand(1:N_t)]
    R, T = [rand(1:N_x), rand(1:N_t)]
    old_loop = tr(RT_loop(test_field,R,T,x,t))
    new_loop = tr(RT_loop(temp_gauge(test_field),R,T,x,t))
    @assert isapprox(old_loop, new_loop) "RT_loop_test() needs a little attention"
    
    @assert isapprox(RT_loop(test_field, 2, 3, x, t), loop_2x3_square(test_field, x, t))
    return true
end

RT_loop_test()

function edge_loop_1_test()
    test_field = gaugefield_SU2(N_x, N_t, true)
    x, t = [rand(1:N_x), rand(1:N_t)]
    old_loop = tr(edge_loop_1(test_field,x,t))
    new_loop = tr(edge_loop_1(temp_gauge(test_field),x,t))
    @assert isapprox(old_loop, new_loop) "edge_loop_1_test() needs a little attention"
    return true
end

edge_loop_1_test()

function L_loop_1_test()
    test_field = gaugefield_SU2(N_x, N_t, true)
    x, t = [rand(1:N_x), rand(1:N_t)]
    old_loop = tr(L_loop_1(test_field,x,t))
    new_loop = tr(L_loop_1(temp_gauge(test_field),x,t))
    @assert isapprox(old_loop, new_loop) "L_loop_1_test() needs a little attention"
    return true
end

L_loop_1_test()

function rhomb_half_loop_test()
    test_field = gaugefield_SU2(N_x, N_t, true)
    x, t = [rand(1:N_x), rand(1:N_t)]
    old_loop = tr(rhomb_half_loop(test_field,x,t))
    new_loop = tr(rhomb_half_loop(temp_gauge(test_field),x,t))
    @assert isapprox(old_loop, new_loop) "rhomb_half_loop_test() needs a little attention"
    return true
end

rhomb_half_loop_test()

function rhomb_loop_test()
    test_field = gaugefield_SU2(N_x, N_t, true)
    x, t = [rand(1:N_x), rand(1:N_t)]
    old_loop = tr(rhomb_loop(test_field,x,t))
    new_loop = tr(rhomb_loop(temp_gauge(test_field),x,t))
    @assert isapprox(old_loop, new_loop) "rhomb_loop_test() needs a little attention"
    return true
end

rhomb_loop_test()

# Old version of loop_mat_mike(...) that has already been "tested
# (i.e. already obtained results indicate that it works)
function loop_mat_mike_old(U, R, T, β)
    NX = size(U,2)
    NT = size(U,3)
    # staples = [staple(U, μ, t, x) for μ = 1:2, t = 1:NT, x = 1:NX]    #⭕ staple_dag???
    avg_U = Array{coeffs_SU2}(undef, 2, NX, NT)
    res = Matrix{coeffs_SU2}(undef, NX, NT)
    for t = 1:NT
        for x = 1:NX
            res[x,t] = coeffs_Id_SU2()
            for μ = 1:2
                stap = staple(U,μ,x,t)
                d = sqrt(det(stap))     #⭕ sqrt(abs(det(stap))) ???
                avg_U[μ,x,t] = (besseli(2,β*d) / (besseli(1,β*d) * d)) * stap
                # avg_U[μ,t,x] = U[μ,t,x] # for debugging
            end
        end
    end

    x_arr = collect(1:NX)
    t_arr = collect(1:NT)

    if R > 2
        res = U[1,x_arr,t_arr]
        for i = 2:R-1
            circshift!(x_arr,-1)    # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
            res = res .* avg_U[1,x_arr,t_arr]
        end
        circshift!(x_arr,-1)
        res = res .* U[1,x_arr,t_arr]
        circshift!(x_arr,-1)
    else
        for i = 1:R
            res = res .* U[1,x_arr,t_arr]
            circshift!(x_arr,-1)
        end
    end

    if T > 2
        res = res .* U[2,x_arr,t_arr]
        for i = 2:T-1
            circshift!(t_arr,-1)   
            res = res .* avg_U[2,x_arr,t_arr]
        end
        circshift!(t_arr,-1)   
        res = res .* U[2,x_arr,t_arr]
        circshift!(t_arr,-1)   
    else
        for i = 1:T
            res = res .* U[2,x_arr,t_arr]
            circshift!(t_arr,-1)   
        end
    end
    
    if R > 2
        circshift!(x_arr,1)
        res = res .* adjoint.(U[1,x_arr,t_arr])
        for i = 2:R-1
            circshift!(x_arr,1)
            res = res .* adjoint.(avg_U[1,x_arr,t_arr])
        end
        circshift!(x_arr,1)
        res = res .* adjoint.(U[1,x_arr,t_arr])
    else
        for i = 1:R
            circshift!(x_arr,1)
            res = res .* adjoint.(U[1,x_arr,t_arr])
        end
    end
    
    if T > 2
        circshift!(t_arr,1)
        res = res .* adjoint.(U[2,x_arr,t_arr])
        for i = 2:T-1
            circshift!(t_arr,1)
            res = res .* adjoint.(avg_U[2,x_arr,t_arr])
        end
        circshift!(t_arr,1)
        res = res .* adjoint.(U[2,x_arr,t_arr])
    else
        for i = 1:T
            circshift!(t_arr,1)
            res = res .* adjoint.(U[2,x_arr,t_arr])
        end
    end

    return res
end

function loop_mat_mike_test()
    test_field = gaugefield_SU2(N_x, N_t, true)
    x, t = [rand(1:N_x), rand(1:N_t)]
    R, T = [rand(1:N_x), rand(1:N_t)]
    old_mat = loop_mat_mike_old(test_field,R,T,β)
    new_mat = loop_mat_mike(test_field,R,T,β)
    @assert isapprox(old_mat[x,t], new_mat[x,t]) "loop_mat_mike_test() needs a little attention"
    return true
end

loop_mat_mike_test()

# auto_corr, auto_corr_norm, auto_corr_time