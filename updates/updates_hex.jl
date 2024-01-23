include("updates_head.jl")


#=
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

    # 🐌 More efficient: only use adj_SU2 once 🐌 (but less human-readable, no?)
    if μ == 1
        # n a a a n
        # 23123
        # t, t_pp, t_pp, t_p, t_p
        # x, x_p, x_p, x, x
        a = U[2,t,x] * adj_SU2(U[3,t_pp,x_p]) * adj_SU2(U[1,t_pp,x_p]) * adj_SU2(U[2,t_p,x]) * U[3,t_p,x]
        # n a a a n
        # 32132
        # t, t_mm, t_mm, t_m, t_m
        # x, x_m, x_m, x_m, x_m
        b = U[3,t,x] * adj_SU2(U[2,t_mm,x_m]) * adj_SU2(U[1,t_mm,x_m]) * adj_SU2(U[3,t_m,x_m]) * U[2,t_m,x_m]
    elseif μ == 2
        # n n a a a
        # 13213
        # t_p, t_p, t_m, t_m ,t
        # x_p, x_p, x, x, x
        a = U[1,t_p,x_p] * U[3,t_p,x_p] * adj_SU2(U[2,t_m,x]) * adj_SU2(U[1,t_m,x]) * adj_SU2(U[3,t,x]) 
        # a a a n n
        # 31231
        # t_pp, t_pp, t_p, t_p, t
        # x_p, x_p, x, x, x
        b = adj_SU2(U[3,t_pp,x_p]) * adj_SU2(U[1,t_pp,x_p]) * adj_SU2(U[2,t_p,x]) * U[3,t_p,x] * U[1,t,x]
    elseif μ == 3
        # n n a a a 
        # 12312
        # t_m, t_m, t_p, t_p, t
        # x, x, x_p, x_p, x
        a = U[1,t_m,x] * U[2,t_m,x] * adj_SU2(U[3,t_p,x_p]) * adj_SU2(U[1,t_p,x_p]) * adj_SU2(U[2,t,x])
        # a a a n n
        # 21321
        # t_mm, t_mm, t_m, t_m, t
        # x_m, x_m, x_m, x_m, x
        b = adj_SU2(U[2,t_mm,x_m]) * adj_SU2(U[1,t_mm,x_m]) * adj_SU2(U[3,t_m,x_m]) * U[2,t_m,x_m] * U[1,t,x]
    end

    return a+b
end
=#

# Compute the adjoint (daggered) staple of U[μ,x,t], where U is
# a coeffs_SU2-valued hexagonal config stored in brick-fashion.
function staple_dag_hex(U, μ, x, t)
    # NX = N_x>>1
    NX = size(U,2)
    NT = size(U,3)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    xp = mod1(x+1, NX) # x%NX +1                 
    tp = mod1(t+1, NT) # t%NT +1                 
    xm = mod1(x-1, NX) # (x + NX -2)%NX +1  
    tm = mod1(t-1, NT) # (t + NT -2)%NT +1   
    xpp = mod1(x+2, NX) # (x+1)%NX +1
    tpp = mod1(t+2, NT) # (t+1)%NT +1
    xmm = mod1(x-2, NX) # (x + NX - 3)%NX +1
    tmm = mod1(t-2, NT) # (t + NT - 3)%NT +1

    # coords = half_chess_coords(N_t, N_x)

    # 🐌 More efficient: only use adj_SU2 once 🐌 (but less human-readable, no?)
    if μ == 2
        if iseven(x+t)
            # n  n   a  a  a 
            # 2  1   2  2  1 
            # x  x   xp xp x
            # tp tpp tp t  t
            a = U[2,x,tp] * U[1,x,tpp] * adj_SU2(U[2,xp,tp]) * adj_SU2(U[2,xp,t]) * adj_SU2(U[1,x,t])
            # a  a  a  n  n 
            # 1  2  2  1  2 
            # xm xm xm xm x 
            # tp t  tm tm tm 
            b = adj_SU2(U[1,xm,tp]) * adj_SU2(U[2,xm,t]) * adj_SU2(U[2,xm,tm]) * U[1,xm,tm] * U[2,x,tm]
        else
            # n  a  a  a  n 
            # 1  2  2  1  2 
            # x  xp xp x  x
            # tp t  tm tm tm 
            a = U[1,x,tp] * adj_SU2(U[2,xp,t]) * adj_SU2(U[2,xp,tm]) * adj_SU2(U[1,x,tm]) * U[2,x,tm]
            # n  a   a  a  n 
            # 2  1   2  2  1
            # x  xm  xm xm xm
            # tp tpp tp t  t
            b = U[2,x,tp] * adj_SU2(U[1,xm,tpp]) * adj_SU2(U[2,xm,tp]) * adj_SU2(U[2,xm,t]) * U[1,xm,t]
        end
    else #if μ == 1
        # n  n  a   a  a
        # 2  2  1   2  2
        # xp xp x   x  x
        # t  tp tpp tp t
        a = U[2,xp,t] * U[2,xp,tp] * adj_SU2(U[1,x,tpp]) * adj_SU2(U[2,x,tp]) * adj_SU2(U[2,x,t])
        # a  a   a   n   n 
        # 2  2   1   2   2
        # xp xp  x   x   x 
        # tm tmm tmm tmm tm 
        b = adj_SU2(U[2,xp,tm]) * adj_SU2(U[2,xp,tmm]) * adj_SU2(U[1,x,tmm]) * U[2,x,tmm] * U[2,x,tm]
    end
    # if μ == 1
    #     if iseven(t+x)
    #         # a  a  a  n  n 
    #         # 2  1  1  2  1
    #         # tp t  tm tm tm
    #         # xm xm xm xm x
    #         a = adj_SU2(U[2,tp,xm]) * adj_SU2(U[1,t,xm]) * adj_SU2(U[1,tm,xm]) * U[2,tm,xm] * U[1,tm,x]
    #         # n  n   a  a  a
    #         # 1  2   1  1  2
    #         # tp tpp tp t  t
    #         # x  x   xp xp x
    #         b = U[1,tp,x] * U[2,tpp,x] * adj_SU2(U[1,tp,xp]) * adj_SU2(U[1,t,xp]) * adj_SU2(U[2,t,x])
    #     else
    #         # n  a   a  a  n
    #         # 1  2   1  1  2
    #         # tp tpp tp t  t
    #         # x  xm  xm xm xm 
    #         a = U[1,tp,x] * adj_SU2(U[2,tpp,xm]) * adj_SU2(U[1,tp,xm]) * adj_SU2(U[1,t,xm]) * U[2,t,xm]
    #         # n  a  a  a  n
    #         # 2  1  1  2  1
    #         # tp t  tm tm tm 
    #         # x  xp xp x  x
    #         b = U[2,tp,x] * adj_SU2(U[1,t,xp]) * adj_SU2(U[1,tm,xp]) * adj_SU2(U[2,tm,x]) * U[1,tm,x]
    #     end
    # else #if μ == 2
    #     # n  n  a   a  a
    #     # 1  1  2   1  1
    #     # t  tp tpp tp t
    #     # xp xp x   x  x
    #     a = U[1,t,xp] * U[1,tp,xp] * adj_SU2(U[2,tpp,x]) * adj_SU2(U[1,tp,x]) * adj_SU2(U[1,t,x])
    #     # a  a   a   n   n 
    #     # 1  1   2   1   1 
    #     # tm tmm tmm tmm tm 
    #     # xp xp  x   x   x 
    #     b = adj_SU2(U[1,tm,xp]) * adj_SU2(U[1,tmm,xp]) * adj_SU2(U[2,tmm,x]) * U[1,tmm,x] * U[1,tm,x]
    # end
    return a+b
end

# Just for testing purposes, not to be actually used
function delta_S_gauge_hex(U, μ, x, t, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2, β)
    return β/2*tr((old_coeffs - new_coeffs) * staple_dag_hex(U,μ,x,t))
end

function metro_hex!(U, μ, x, t, step, β, acc)
    new_coeffs = ran_SU2(step) * U[μ,x,t]
    staple_d = staple_dag_hex(U,μ,x,t)
    S_old = β*0.5*tr(U[μ,x,t] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function chess_metro_hex!(U, step, β, acc)
    NX = size(U,2)
    NT = size(U,3)
    μ = 2
    for trip = 1:2
        for start_t = 1:2
            for t = start_t:2:NT
                for x = (1+mod(t+trip,2)):2:NX
                    metro_hex!(U,μ,x,t,step,β,acc)
                    # println(μ, ", ", x, ", ", t)
                end
            end
        end
    end
    μ = 1
    for t = 1:NT
        for x = 2-mod(t,2):2:NX
            metro_hex!(U,μ,x,t,step,β,acc)
            # println(μ, ", ", x, ", ", t)
        end
    end
    return nothing
end

#
function lexico_metro_hex!(U, step, β, acc)
    NX = size(U,2)
    NT = size(U,3)
    for t = 1:NT
        for x = 1:NX
            if x in 1+mod(t+1,2):2:NX
                metro_hex!(U,1,x,t,step,β,acc)
                metro_hex!(U,2,x,t,step,β,acc)
            elseif x in 1+mod(t,2):2:NX
                metro_hex!(U,2,x,t,step,β,acc)
            end
        end
    end
    return nothing
end

# 
function ran_metro_hex!(U, step, β, acc)
    NX = size(U,2)
    NT = size(U,3)
    coords = chess_hex_link_coords(NX,NT)
    r = rand(1:length(coords), length(coords))
    for i = 1:length(coords)
        μ, x, t = coords[r[i]]
        metro_hex!(U,μ,x,t,step,β,acc)
    end    
    return nothing
end

function overrelax_hex!(U, μ, x, t)
    v = proj_SU2(staple_dag_hex(U,μ,x,t))
    U[μ,x,t] = adj_SU2(v *  U[μ,x,t] * v)
    return nothing
end

function chess_overrelax_hex!(U, acc)
    NX = size(U,2)
    NT = size(U,3)
    μ = 2
    for trip = 1:2
        for start_t = 1:2
            for t = start_t:2:NT
                for x = (1+mod(t+trip,2)):2:NX
                    overrelax_hex!(U,μ,x,t)
                    # println(μ, ", ", x, ", ", t)
                end
            end
        end
    end
    μ = 1
    for t = 1:NT
        for x = 2-mod(t,2):2:NX
            overrelax_hex!(U,μ,x,t)
            # println(μ, ", ", x, ", ", t)
        end
    end
    acc[1] += Int(1.5*NX*NT)
    return nothing
end


# NT = 5
# NX = 5
# μ = 1
# for trip = 1:2
#     for trip_x = 1:2
#         for x = trip_x:2:NX
#             for t = (1+mod(x+trip,2)):2:NT
#                 println(μ, ", ", t, ", ", x)
#             end
#         end
#     end
# end
# μ = 2
# trip = 1
# for trip_t = 0:2:2
#     for x = 1:NX
#         for t = trip_t+(1+mod(x+trip,2)):4:NT
#             println(μ, ", ", t, ", ", x)
#         end
#     end
# end
