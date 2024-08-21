# We want:
# staple, local metropolis, chess metropolis
include("D:\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")

function staple_dag(U, μ, t, x)
    NT = size(U,2)
    NX = size(U,3)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    t_p = mod1(t+1, NT) # t%NT +1                 
    x_p = mod1(x+1, NX) # x%NX +1                 
    t_m = mod1(t-1, NT) # (t + NT -2)%NT +1   
    x_m = mod1(x-1, NX) # (x + NX -2)%NX +1   

    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    if μ == 1
        a = U[2,t_p,x] * adjoint(U[1,t,x_p]) * adjoint(U[2,t,x])
        b = adjoint(U[2,t_p,x_m]) * adjoint(U[1,t,x_m]) * U[2,t,x_m]
    else #if μ == 2
        a = U[1,t,x_p] * adjoint(U[2,t_p,x]) * adjoint(U[1,t,x])
        b = adjoint(U[1,t_m,x_p]) * adjoint(U[2,t_m,x]) * U[1,t_m,x]
    end
    return a + b 
end

# ❌ Not in use anymore! Only for testing purposes ❌
function delta_S_gauge(U, μ, t, x, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2, β)
    return β*0.5*(tr((old_coeffs + (-1) * new_coeffs) * staple_dag(U,μ,t,x)))
end

function metro!(U, μ, t, x, step, β, acc)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,t,x]
    staple_d = staple_dag(U,μ,t,x)
    S_old = β*0.5*tr(U[μ,t,x] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function overrelax!(U, μ, t, x)
    v = proj2man!(staple_dag(U,μ,t,x))
    U[μ,t,x] = adjoint(v *  U[μ,t,x] * v)
    return nothing
end

function lexico_metro!(U, step, β, acc)
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:2
        for t = 1:NT
            for x = 1:NX
                metro!(U,μ,t,x,step, β, acc)
            end
        end
    end
    return nothing
end

function chess_metro!(U, step, β, acc)
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:2
        for trip = 1:2
            for x = 1:NX
                for t = (1+mod(x+trip,2)):2:NT
                    metro!(U,μ,t,x,step, β, acc)
                end
            end
        end
    end
    return nothing
end

#
function lexico_overrelax!(U, acc)
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:2
        for t = 1:NT
            for x = 1:NX
                overrelax!(U, μ, t, x)
            end
        end
    end
    acc[1] += 2*NT*NX
    return nothing
end

#
function chess_overrelax!(U, acc)
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    overrelax!(U,μ,t,x)
                end
            end
        end
    end
    acc[1] += 2*NT*NX
    return nothing
end





####    Hexagonal shenanigans    ####     🚧👷 under construction 👷🚧





function staple_dag_hex(U, μ, t, x)
    # NX = N_x>>1
    NT = size(U,2)
    NX = size(U,3)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    t_p = mod1(t+1, NT) # t%NT +1                 
    x_p = mod1(x+1, NX) # x%NX +1                 
    t_m = mod1(t-1, NT) # (t + NT -2)%NT +1   
    x_m = mod1(x-1, NX) # (x + NX -2)%NX +1  
    t_pp = mod1(t+2, NT) # (t+1)%NT +1
    t_mm = mod1(t-2, NT) # (t + NT - 3)%NT +1

    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    if μ == 1
        # n a a a n
        # 23123
        # t, t_pp, t_pp, t_p, t_p
        # x, x_p, x_p, x, x
        a = U[2,t,x] * adjoint(U[3,t_pp,x_p]) * adjoint(U[1,t_pp,x_p]) * adjoint(U[2,t_p,x]) * U[3,t_p,x]
        # n a a a n
        # 32132
        # t, t_mm, t_mm, t_m, t_m
        # x, x_m, x_m, x_m, x_m
        b = U[3,t,x] * adjoint(U[2,t_mm,x_m]) * adjoint(U[1,t_mm,x_m]) * adjoint(U[3,t_m,x_m]) * U[2,t_m,x_m]
    elseif μ == 2
        # n n a a a
        # 13213
        # t_p, t_p, t_m, t_m ,t
        # x_p, x_p, x, x, x
        a = U[1,t_p,x_p] * U[3,t_p,x_p] * adjoint(U[2,t_m,x]) * adjoint(U[1,t_m,x]) * adjoint(U[3,t,x]) 
        # a a a n n
        # 31231
        # t_pp, t_pp, t_p, t_p, t
        # x_p, x_p, x, x, x
        b = adjoint(U[3,t_pp,x_p]) * adjoint(U[1,t_pp,x_p]) * adjoint(U[2,t_p,x]) * U[3,t_p,x] * U[1,t,x]
    elseif μ == 3
        # n n a a a 
        # 12312
        # t_m, t_m, t_p, t_p, t
        # x, x, x_p, x_p, x
        a = U[1,t_m,x] * U[2,t_m,x] * adjoint(U[3,t_p,x_p]) * adjoint(U[1,t_p,x_p]) * adjoint(U[2,t,x])
        # a a a n n
        # 21321
        # t_mm, t_mm, t_m, t_m, t
        # x_m, x_m, x_m, x_m, x
        b = adjoint(U[2,t_mm,x_m]) * adjoint(U[1,t_mm,x_m]) * adjoint(U[3,t_m,x_m]) * U[2,t_m,x_m] * U[1,t,x]
    end

    return a+b
end

# Just for testing purposes, not to be actually used
function delta_S_gauge_hex(U, μ, t, x, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2, β)
    return β/2*tr((old_coeffs + (-1) * new_coeffs) * staple_dag_hex(U,μ,t,x))
end

function metro_hex!(U, μ, t, x, step, β, acc)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,t,x]
    staple_d = staple_dag_hex(U,μ,t,x)
    S_old = β*0.5*tr(U[μ,t,x] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function chess_metro_hex!(U, step, β, acc)
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:3
        for trip = 1:2
            for x = 1:NX
                for t = (1+mod(x+trip,2)):2:NT
                    metro_hex!(U,μ,t,x,step, β, acc)
                end
            end
        end
    end
    return nothing
end

function lexico_metro_hex!(U, step, β, acc)
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:3
        for x = 1:NX
            for t = 1:NT
                metro_hex!(U,μ,t,x,step, β, acc)
            end
        end
    end
    return nothing
end

function ran_metro_hex!(U, step, β, acc)
    NT = size(U,2)
    NX = size(U,3)
    for i = 1:3*NT*NX
        μ = rand(1:3)
        t = rand(1:NT)
        x = rand(1:NX)
        metro_hex!(U,μ,t,x,step, β, acc)
    end    
    return nothing
end

#=
for L in [16,32,64,96]
    N_t = N_x = L
    acc = [0]
    β = 6.0
    ϵ = 0.2
    n_meas = 50000
    counter = 0
    # for i = 1:100
    #     chess_metro_hex!(test_hex, 0.2)
    # end
    test_lex = hexfield_SU2(N_t,N_x,false)
    test_chess = deepcopy(test_lex)
    test_ran = deepcopy(test_lex)
    actions_lex = []
    actions_chess = []
    actions_ran = []
    for i = 1:n_meas
        if i%(Int(n_meas/100)) == 1
            println(" ")
            println("We're already ", counter, "% deep in the HEXAGONAL simulation with N_t = $N_t, N_x = $N_x, β = $β and ϵ = $ϵ !")
            # mywrite(last_conf_path, U)
            # bla = open(last_conf_path, "a")
            # write(bla, "Progress in simulation: $counter %")
            # close(bla)
            counter += 1
        end
        lexico_metro_hex!(test_lex, ϵ)
        chess_metro_hex!(test_chess, ϵ)
        ran_metro_hex!(test_ran, ϵ)
        push!(actions_lex, action_hex(test_lex))
        push!(actions_chess, action_hex(test_chess))
        push!(actions_ran, action_hex(test_ran))
    end
    image = plot(actions_chess[1:200], label = "chess")
    image = plot!(actions_lex[1:200], label = "lexicog.")
    image = plot!(actions_ran[1:200], label = "random")
    image = plot!(title = "Action Time Series, Different Metrop. Updates, Hex.
    N_t = N_x = $N_t,  β = $β,  ϵ = $ϵ",
    xlabel = "MC time")
    display(image)
    savefig("D:\\Physik Uni\\julia_projects\\SU2\\hex_data\\auto_corr\\actions_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.pdf")
    
    # auto_chess = []
    # auto_lex = []
    # auto_ran = []
    # counter = 0
    # for start = 50:100:n_meas-50
    #     if start%(Int(n_meas/50)) == 1
    #         # println(" ")
    #         println("Progress for autos: $counter%")
    #         # mywrite(last_conf_path, U)
    #         # bla = open(last_conf_path, "a")
    #         # write(bla, "Progress in simulation: $counter %")
    #         # close(bla)
    #         counter += 2
    #     end
    #     push!(auto_chess, auto_corr_time(actions_chess[start:n_meas]))
    #     push!(auto_lex, auto_corr_time(actions_lex[start:n_meas]))
    #     push!(auto_ran, auto_corr_time(actions_ran[start:n_meas]))
    # end
    # image = plot(50:100:n_meas-50, auto_chess, label = "chess")
    # image = plot!(50:100:n_meas-50, auto_lex, label = "lexicog.")
    # image = plot!(50:100:n_meas-50, auto_ran, label = "ran.")
    # image = plot!(title = "Int. Autoc. Time, Different Metrop. Updates, Hex.
    # N_t = N_x = $N_t,  β = $β,  ϵ = $ϵ",
    # xlabel = "MC time at which summation began")
    # display(image)
    # savefig("D:\\Physik Uni\\julia_projects\\SU2\\hex_data\\auto_corr\\auto_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.pdf")
    # println("chess: ", mean(auto_chess), ", ", std(auto_chess)/sqrt(length(auto_chess)))
    # println("lexicog.: ", mean(auto_lex),  ", ", std(auto_lex)/sqrt(length(auto_lex)))
    # println("random: ", mean(auto_ran),  ", ", std(auto_ran)/sqrt(length(auto_ran)))
    # println(" ")

    writedlm("D:\\Physik Uni\\julia_projects\\SU2\\hex_data\\auto_corr\\actions_2_lex_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.txt", actions_lex)
    writedlm("D:\\Physik Uni\\julia_projects\\SU2\\hex_data\\auto_corr\\actions_2_chess_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.txt", actions_chess)
    writedlm("D:\\Physik Uni\\julia_projects\\SU2\\hex_data\\auto_corr\\actions_2_ran_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.txt", actions_ran)

    # writedlm("D:\\Physik Uni\\julia_projects\\SU2\\hex_data\\auto_corr\\auto_lex_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.txt", auto_lex)
    # writedlm("D:\\Physik Uni\\julia_projects\\SU2\\hex_data\\auto_corr\\auto_chess_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.txt", auto_chess)
    # writedlm("D:\\Physik Uni\\julia_projects\\SU2\\hex_data\\auto_corr\\auto_ran_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.txt", auto_ran)
end
=#


# test_field = gaugefield_SU2(32,32,true)
# test_chess = deepcopy(test_field)
# actions_lex = []
# actions_chess = []
# for i = 1:300
#     lexico_metro!(test_field, 0.2)
#     chess_metro!(test_chess, 0.2)
#     push!(actions_lex, action(test_field))
#     push!(actions_chess, action(test_chess))
# end
# auto_corr_time(actions_chess[1:300])
# auto_corr_time(actions_lex[1:300])
# plot(actions_chess)
# plot!(actions_lex)


function overrelax_hex!(U, μ, t, x)
    v = proj2man!(staple_dag_hex(U,μ,t,x))
    U[μ,t,x] = adjoint(v *  U[μ,t,x] * v)
    return nothing
end

function chess_overrelax_hex!(U, acc)
    # NX = N_x>>1
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:3
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    overrelax_hex!(U,μ,t,x)
                end
            end
        end
    end
    acc[1] += 3*NT*NX
    return nothing
end





####    D-dimensional shenanigans   #### 🚧👷 under construction 👷🚧 





# 
function staple_dag_3d(U, μ, t, x, y)
    NT = size(U,2)
    NX = size(U,3)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    c = coeffs_SU2(0.0,0.0,0.0,0.0)
    d = coeffs_SU2(0.0,0.0,0.0,0.0)
    t_p = mod1(t+1, NT) # t%NT +1                 
    x_p = mod1(x+1, NX) # x%NX +1          
    y_p = mod1(y+1, NX)       
    t_m = mod1(t-1, NT) # (t + NT -2)%NT +1   
    x_m = mod1(x-1, NX) # (x + NX -2)%NX +1   
    y_m = mod1(y-1, NX)

    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    if μ == 1
        a = U[2,t_p,x,y] * adjoint(U[1,t,x_p,y]) * adjoint(U[2,t,x,y])
        b = adjoint(U[2,t_p,x_m,y]) * adjoint(U[1,t,x_m,y]) * U[2,t,x_m,y]
        c = U[3,t_p,x,y] * adjoint(U[1,t,x,y_p]) * adjoint(U[3,t,x,y])
        d = adjoint(U[3,t_p,x,y_m]) * adjoint(U[1,t,x,y_m]) * U[3,t,x,y_m]
    elseif μ == 2
        a = U[1,t,x_p,y] * adjoint(U[2,t_p,x,y]) * adjoint(U[1,t,x,y])
        b = adjoint(U[1,t_m,x_p,y]) * adjoint(U[2,t_m,x,y]) * U[1,t_m,x,y]
        c = U[3,t,x_p,y] * adjoint(U[2,t,x,y_p]) * adjoint(U[3,t,x,y])
        d = adjoint(U[3,t,x_p,y_m]) * adjoint(U[2,t,x,y_m]) * U[3,t,x,y_m]
    elseif μ == 3
        # n a a / a a n 
        # 1 3 1 / 1 3 1 
        # t tp t / tm tm tm 
        # x: const
        # yp y y / yp y y
        a = U[1,t,x,y_p] * adjoint(U[3,t_p,x,y]) * adjoint(U[1,t,x,y])
        b = adjoint(U[1,t_m,x,y_p]) * adjoint(U[3,t_m,x,y]) * U[1,t_m,x,y]
        # n a a / a a n 
        # 2 3 2 / 2 3 2
        # t: const
        # x xp x / xm xm xm 
        # yp y y / yp y y
        c = U[2,t,x,y_p] * adjoint(U[3,t,x_p,y]) * adjoint(U[2,t,x,y])
        d = adjoint(U[2,t,x_m,y_p]) * adjoint(U[3,t,x_m,y]) * U[2,t,x_m,y]
    end
    return a + b + c + d
end


# ❌ Not in use! Only for testing purposes ❌
function delta_S_gauge_3d(U, μ, t, x, y, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2, β)
    return β*0.5*(tr((old_coeffs + (-1) * new_coeffs) * staple_dag_3d(U,μ,t,x,y)))
end

# N_t = N_x = 16
# β = 1.0
# μ = rand(1:3)
# t = rand(1:N_t)
# x = rand(1:N_x)
# y = rand(1:N_x)
# old_field = gaugefield_SU2_3d(N_t, N_x, true)
# new_field = deepcopy(old_field)
# new_field[μ,t,x,y] = ran_SU2(rand())
# delta_S_gauge_3d(old_field, μ, t, x, y, old_field[μ,t,x,y], new_field[μ,t,x,y])
# action_3d(old_field) - action_3d(new_field)


# metro! but for D-dim. configs
function metro_3d!(U, μ, t, x, y, step, β, acc)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,t,x,y]
    staple_d = staple_dag_3d(U,μ,t,x,y)
    S_old = β*0.5*tr(U[μ,t,x,y] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x,y] = new_coeffs
        acc[1] += 1
    end
    return nothing
end


# 
function chess_metro_3d!(U, step, β, acc)
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:3
        for trip = 1:2
            for y = 1:NX
                for x = 1:NX
                    for t = (1+mod(x+y+trip,2):2:NT)
                        metro_3d!(U, μ, t, x, y, step, β, acc)
                    end
                end
            end
        end
    end
    return nothing
end

# NX = NT = 5
# for trip = 1:2
#     for y = 1:NX
#         for x = 1:NX
#             for t = (1+mod(x+y+trip,2):2:NT)
#                 println([t,x,y])
#             end
#         end
#     end
# end

# N_t = 4
# N_x = 4
# list_test = []
# for trip = 1:2
#     for y = 1:N_x
#         for x = 1:N_x
#             for t = (1+mod(x+y+trip,2):2:N_t)
#                 push!(list_test, [t,x,y])
#             end
#         end
#     end
# end

# length(list_test)
# length(list)

# list = D_dim_coords(N_t,N_x,3)
# for ind in list
#     if !(ind in list_test)
#         println(ind, " is not in list_test")
#     end
# end
