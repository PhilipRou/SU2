include("updates_head.jl")



# 
function staple_dag_cube(U, μ, x, y, t)
    NX = size(U,2)
    NT = size(U,4)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    c = coeffs_SU2(0.0,0.0,0.0,0.0)
    d = coeffs_SU2(0.0,0.0,0.0,0.0)
    xp = mod1(x+1, NX) # x%NX +1          
    yp = mod1(y+1, NX)       
    tp = mod1(t+1, NT) # t%NT +1                 
    xm = mod1(x-1, NX) # (x + NX -2)%NX +1   
    ym = mod1(y-1, NX)
    tm = mod1(t-1, NT) # (t + NT -2)%NT +1   

    # 🐌 More efficient: only use adj_SU2 once 🐌 (but less human-readable, no?)
    if μ == 1
        # n  a  a / a  a  n 
        # 2  1  2 / 2  1  2
        # xp x  x / xp x  x
        # y  yp y / ym ym ym
        # t: const
        a = U[2,xp,y,t] * adj_SU2(U[1,x,yp,t]) * adj_SU2(U[2,x,y,t])
        b =  adj_SU2(U[2,xp,ym,t]) * adj_SU2(U[1,x,ym,t]) * U[2,x,ym,t]
        # n  a  a / a  a  n 
        # 3  1  3 / 3  1  3
        # xp x  x / xp x  x
        # y: const
        # t  tp t / tm tm tm
        c = U[3,xp,y,t] * adj_SU2(U[1,x,y,tp]) * adj_SU2(U[3,x,y,t])
        d = adj_SU2(U[3,xp,y,tm]) * adj_SU2(U[1,x,y,tm]) * U[3,x,y,tm]
    elseif μ == 2
        # n  a  a / a  a  n 
        # 1  2  1 / 1  2  1
        # x  xp x / xm xm xm
        # yp y  y / yp y  y
        # t: const.
        a = U[1,x,yp,t] * adj_SU2(U[2,xp,y,t]) * adj_SU2(U[1,x,y,t])
        b = adj_SU2(U[1,xm,yp,t]) * adj_SU2(U[2,xm,y,t]) * U[1,xm,y,t]
        # n  a  a / a  a  n 
        # 3  2  3 / 3  2  3
        # x: const.
        # yp y  y / yp y  y
        # t  tp t / tm tm tm 
        c = U[3,x,yp,t] * adj_SU2(U[2,x,y,tp]) * adj_SU2(U[3,x,y,t])
        d = adj_SU2(U[3,x,yp,tm]) * adj_SU2(U[2,x,y,tm]) * U[3,x,y,tm]
    else #if μ == 3
        # n  a  a / a  a  n 
        # 1  3  1 / 1  3  1 
        # x  xp x / xm xm xm
        # y: const
        # tp t  t / tp t  t 
        a = U[1,x,y,tp] * adj_SU2(U[3,xp,y,t]) * adj_SU2(U[1,x,y,t])
        b = adj_SU2(U[1,xm,y,tp]) * adj_SU2(U[3,xm,y,t]) * U[1,xm,y,t]
        # n a a / a a n 
        # 2 3 2 / 2 3 2
        # x: const 
        # y  yp y / ym ym ym 
        # tp t  t / tp t  t
        c = U[2,x,y,tp] * adj_SU2(U[3,x,yp,t]) * adj_SU2(U[2,x,y,t])
        d = adj_SU2(U[2,x,ym,tp]) * adj_SU2(U[3,x,ym,t]) * U[2,x,ym,t]
    end
    return a + b + c + d
end


# ❗ Not in use! Only for testing purposes
function delta_S_gauge_cube(U, μ, x, y, t, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2, β)
    return β*0.5*(tr((old_coeffs - new_coeffs) * staple_dag_cube(U,μ,x,y,t)))
end

# N_t = N_x = 16
# β = 1.0
# μ = rand(1:3)
# x = rand(1:N_x)
# y = rand(1:N_x)
# t = rand(1:N_t)
# old_field = gaugefield_SU2_cube(N_x, N_t, true)
# new_field = deepcopy(old_field)
# new_field[μ,x,y,t] = ran_SU2(rand())
# delta_S_gauge_cube(old_field, μ, x, y, t, old_field[μ,x,y,t], new_field[μ,x,y,t], β)
# action_cube(old_field) - action_cube(new_field)


# metro! but for D-dim. configs
function metro_cube!(U, μ, x, y, t, step, β, acc)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,x,y,t]
    staple_d = staple_dag_cube(U,μ,x,y,t)
    S_old = β*0.5*tr(U[μ,x,y,t] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,x,y,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end


# 
function chess_metro_cube!(U, step, β, acc)
    NX = size(U,2)
    NT = size(U,4)
    for μ = 1:3
        for trip = 1:2
            for t = 1:NT
                for y = 1:NX
                    for x = 1+mod(t+y+trip,2):2:NX
                        metro_cube!(U, μ, x, y, t, step, β, acc)
                    end
                end
            end
        end
    end
    return nothing
end

# NX = NT = 5
# for trip = 1:2
#     for t = 1:NT
#         for y = 1:NX
#             for x = 1+mod(t+y+trip,2):2:NX
#                 println([x,y,t])
#             end
#         end
#     end
# end

# N_t = 4
# N_x = 4
# list_test = []
# for trip = 1:2
#     for t = 1:NT
#         for y = 1:NX
#             for x = 1+mod(t+y+trip,2):2:NX
#                 push!(list,[x,y,t])
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
