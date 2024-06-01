include("test_head.jl")


function delta_S_gauge_test()
    test_field_1 = gaugefield_SU2(N_t,N_x,true)
    test_field_2 = deepcopy(test_field_1)
    μ = rand(1:2)
    t = rand(1:N_t)
    x = rand(1:N_x)
    test_field_2[μ,t,x] = ran_SU2(rand())
    true_dif =  action(test_field_2,β)-action(test_field_1,β)
    test_dif = delta_S_gauge(test_field_1, μ,t,x, test_field_1[μ,t,x], test_field_2[μ,t,x], β)
    @assert isapprox(true_dif, test_dif) "Something went wrong while testing delta_S_gauge"
    return true
end

delta_S_gauge_test()

# function staple_dag_D_test()
#     test_field = gaugefield_SU2(64, 64, true)
#     μ = rand(1:2)
#     t = rand(1:64)
#     x = rand(1:64)
#     @assert isapprox(staple_dag(test_field, μ,t,x), staple_dag_D(test_field, μ, [t,x]))
#     return true
# end

# staple_dag_D_test()

acc = [0.0]
function local_metro!_test()
    test_field = gaugefield_SU2(N_t,N_x,true)
    test_test_field = deepcopy(test_field)
    for i = 1:50
        metro!(test_field, 1,1,1,0.1,β,acc,"SU2")
    end
    @assert test_test_field != test_field "local_metro! didn't change the config it was given"
    @assert acc[1] != 0 "local_metro! didn't increment the acceptance"
    return true
end

local_metro!_test()

# function gaugefield_SU2_acc_count_test()
#     test_field = gaugefield_SU2(8, 8, true)
#     for i = 1:50
#         local_metro!(test_field, 1, 1, 1, 0.1)
#     end
#     @assert test_field.acc_count != 0
#     return true
# end

# gaugefield_SU2_acc_count_test()


# Only the indices of chess_metro! are tested, as the only remaining thing,
# local_metro!, has already been tested above.
function chess_indices_test()
    indices = []
    for μ = 1:2
        for t = 1:N_t
            for x = 1:N_x
                push!(indices, [μ,t,x])
            end
        end
    end

    test_indices = []
    for dir = 1:2
        for pass = 1:2
            for t = 1:N_t
                for x = (1+mod(t+pass,2)):2:N_x
                    push!(test_indices, [dir,t,x])
                end
            end
        end
    end

    @assert issetequal(indices, test_indices) "The indices of chess_metro! seem weird"
    return true
end

chess_indices_test()

function overrelax!_test_SU2()
    test_field = gaugefield_SU2(N_t,N_x,true)
    old_field = deepcopy(test_field)
    old_action = action(test_field,β)
    for i = 1:100
        t = rand(1:N_t)
        x = rand(1:N_x)
        μ = rand(1:2)
        overrelax!(test_field,μ,t,x,[0.0])
    end
    @assert old_field != test_field "There was no update in overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action(test_field,β)) "overrelax! didn't leave the action invariant"
    return true
end

overrelax!_test_SU2()

function overrelax!_test_U2()
    test_field = gaugefield_U2(N_t,N_x,true)
    old_field = deepcopy(test_field)
    old_action = action(test_field,β)
    for i = 1:100
        t = rand(1:N_t)
        x = rand(1:N_x)
        μ = rand(1:2)
        overrelax!(test_field,μ,t,x,[0.0])
    end
    @assert old_field != test_field "There was no update in overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action(test_field,β)) "overrelax! didn't leave the action invariant"
    return true
end

overrelax!_test_U2()

acc = [0]
function lexico_overrelax!_test()
    test_field = gaugefield_SU2(N_t,N_x,true)
    old_field = deepcopy(test_field)
    old_action = action(test_field,β)
    lexico_overrelax!(test_field,acc)
    @assert old_field != test_field "There was no update in lexico_overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action(test_field,β)) "lexico_overrelax! didn't leave the action invariant"
    return true
end

lexico_overrelax!_test()

function chess_overrelax!_test()
    test_field = gaugefield_SU2(N_t,N_x,true)
    old_field = deepcopy(test_field)
    old_action = action(test_field,β)
    chess_overrelax!(test_field,acc)
    @assert old_field != test_field "There was no update in overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action(test_field,β)) "lexico_overrelax! didn't leave the action invariant"
    return true
end

lexico_overrelax!_test()

function insta_U2_log_test()
    for i = 1:100
        q = rand(Vector(Int(-N_x*N_t/2):Int(N_x*N_t/2)))
        @assert !(false in isapprox.(exp_u2.(insta_U2_log(N_x,N_t,q)), insta_U2(N_x,N_t,q)))
        @assert round(Int,top_charge_U2(exp_u2.(insta_U2_log(N_x,N_t,q)))) == q "We have $q VS $(round(Int,top_charge_U2(exp_u2.(insta_U2_log(N_x,N_t,q)))))"
        @assert action(exp_u2.(insta_U2_log(N_x,N_t,q)),1) - action(insta_U2(N_x, N_t, q),1) < 1.0e-12 "For q = $q there's ΔS = $(action(exp_u2.(insta_U2_log(N_x,N_t,q)),1) - action(insta_U2(N_x, N_t, q),1))"
    end
    return true
end

insta_U2_log_test()