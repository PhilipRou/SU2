include("updates_head.jl")


#=
function staple_dag_hex(U, μ, t, x)
    # NX = N_x>>1
    NX = size(U,3)
    NX = size(U,2)
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

    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    if μ == 2
        if iseven(x+t)
            # n  n   a  a  a 
            # 2  1   2  2  1 
            # x  x   xp xp x
            # tp tpp tp t  t
            a = U[2,x,tp] * U[1,x,tpp] * adjoint(U[2,xp,tp]) * adjoint(U[2,xp,t]) * adjoint(U[1,x,t])
            # a  a  a  n  n 
            # 1  2  2  1  2 
            # xm xm xm xm x 
            # tp t  tm tm tm 
            b = adjoint(U[1,xm,tp]) * adjoint(U[2,xm,t]) * adjoint(U[2,xm,tm]) * U[1,xm,tm] * U[2,x,tm]
        else
            # n  a  a  a  n 
            # 1  2  2  1  2 
            # x  xp xp x  x
            # tp t  tm tm tm 
            a = U[1,x,tp] * adjoint(U[2,xp,t]) * adjoint(U[2,xp,tm]) * adjoint(U[1,x,tm]) * U[2,x,tm]
            # n  a   a  a  n 
            # 2  1   2  2  1
            # x  xm  xm xm xm
            # tp tpp tp t  t
            b = U[2,x,tp] * adjoint(U[1,xm,tpp]) * adjoint(U[2,xm,tp]) * adjoint(U[2,xm,t]) * U[1,xm,t]
        end
    else #if μ == 1
        # n  n  a   a  a
        # 2  2  1   2  2
        # xp xp x   x  x
        # t  tp tpp tp t
        a = U[2,xp,t] * U[2,xp,tp] * adjoint(U[1,x,tpp]) * adjoint(U[2,x,tp]) * adjoint(U[2,x,t])
        # a  a   a   n   n 
        # 2  2   1   2   2
        # xp xp  x   x   x 
        # tm tmm tmm tmm tm 
        b = adjoint(U[2,xp,tm]) * adjoint(U[2,xp,tmm]) * adjoint(U[1,x,tmm]) * U[2,x,tmm] * U[2,x,tm]
    end
    # if μ == 1
    #     if iseven(t+x)
    #         # a  a  a  n  n 
    #         # 2  1  1  2  1
    #         # tp t  tm tm tm
    #         # xm xm xm xm x
    #         a = adjoint(U[2,tp,xm]) * adjoint(U[1,t,xm]) * adjoint(U[1,tm,xm]) * U[2,tm,xm] * U[1,tm,x]
    #         # n  n   a  a  a
    #         # 1  2   1  1  2
    #         # tp tpp tp t  t
    #         # x  x   xp xp x
    #         b = U[1,tp,x] * U[2,tpp,x] * adjoint(U[1,tp,xp]) * adjoint(U[1,t,xp]) * adjoint(U[2,t,x])
    #     else
    #         # n  a   a  a  n
    #         # 1  2   1  1  2
    #         # tp tpp tp t  t
    #         # x  xm  xm xm xm 
    #         a = U[1,tp,x] * adjoint(U[2,tpp,xm]) * adjoint(U[1,tp,xm]) * adjoint(U[1,t,xm]) * U[2,t,xm]
    #         # n  a  a  a  n
    #         # 2  1  1  2  1
    #         # tp t  tm tm tm 
    #         # x  xp xp x  x
    #         b = U[2,tp,x] * adjoint(U[1,t,xp]) * adjoint(U[1,tm,xp]) * adjoint(U[2,tm,x]) * U[1,tm,x]
    #     end
    # else #if μ == 2
    #     # n  n  a   a  a
    #     # 1  1  2   1  1
    #     # t  tp tpp tp t
    #     # xp xp x   x  x
    #     a = U[1,t,xp] * U[1,tp,xp] * adjoint(U[2,tpp,x]) * adjoint(U[1,tp,x]) * adjoint(U[1,t,x])
    #     # a  a   a   n   n 
    #     # 1  1   2   1   1 
    #     # tm tmm tmm tmm tm 
    #     # xp xp  x   x   x 
    #     b = adjoint(U[1,tm,xp]) * adjoint(U[1,tmm,xp]) * adjoint(U[2,tmm,x]) * U[1,tmm,x] * U[1,tm,x]
    # end
    return a+b
end

# Just for testing purposes, not to be actually used
function delta_S_gauge_hex(U, μ, x, t, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2, β)
    return β/2*tr((old_coeffs - new_coeffs) * staple_dag_hex(U,μ,x,t))
end

function metro_hex!(U, μ, x, t, step, β, acc, group)
    # new_coeffs = ran_SU2(step) * U[μ,x,t]
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    end
    staple_d = staple_dag_hex(U,μ,x,t)
    S_old = β*0.5*real(tr(U[μ,x,t] * staple_d))
    S_new = β*0.5*real(tr(new_coeffs * staple_d))
    if rand() < exp(S_new-S_old)
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function metro_hex_comp!(U, μ, x, t, step, β, acc, group)
    # new_coeffs = ran_SU2(step) * U[μ,x,t]
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    end
    staple_d = staple_dag_hex(U,μ,x,t)
    S_old = β/2*27/4*real(tr(U[μ,x,t] * staple_d))
    S_new = β/2*27/4*real(tr(new_coeffs * staple_d))
    if rand() < exp(S_new-S_old)
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function chess_metro_hex!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    μ = 1
    for trip = 1:2
        for t = trip:2:NT
            # for x = 2-mod(t,2):2:NX
            for x = trip:2:NX
                metro_hex!(U,μ,x,t,step,β,acc, group)
                # println(μ, ", ", x, ", ", t)
            end
        end
    end
    μ = 2
    for trip = 1:2
        for start_t = 1:2
            for t = start_t:2:NT
                # for x = (1+mod(t+trip,2)):2:NX
                for x = trip:2:NX
                    metro_hex!(U,μ,x,t,step,β,acc, group)
                    # println(μ, ", ", x, ", ", t)
                end
            end
        end
    end
    return nothing
end

function chess_metro_hex_comp!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    μ = 1
    for trip = 1:2
        for t = trip:2:NT
            # for x = 2-mod(t,2):2:NX
            for x = trip:2:NX
                metro_hex_comp!(U,μ,x,t,step,β,acc, group)
                # println(μ, ", ", x, ", ", t)
            end
        end
    end
    μ = 2
    for trip = 1:2
        for start_t = 1:2
            for t = start_t:2:NT
                # for x = (1+mod(t+trip,2)):2:NX
                for x = trip:2:NX
                    metro_hex_comp!(U,μ,x,t,step,β,acc, group)
                    # println(μ, ", ", x, ", ", t)
                end
            end
        end
    end
    return nothing
end

#
function lexico_metro_hex!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    for t = 1:NT
        for x = 1:NX
            if x in 1+mod(t+1,2):2:NX
                metro_hex!(U,1,x,t,step,β,acc,group)
                metro_hex!(U,2,x,t,step,β,acc,group)
            elseif x in 1+mod(t,2):2:NX
                metro_hex!(U,2,x,t,step,β,acc, group)
            end
        end
    end
    return nothing
end

# 
function ran_metro_hex!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    coords = hex_links_coords_chess(NX,NT)
    r = rand(1:length(coords), length(coords))
    for i = 1:length(coords)
        μ, x, t = coords[r[i]]
        metro_hex!(U,μ,x,t,step,β,acc, group)
    end    
    return nothing
end

function overrelax_hex!(U, μ, x, t)
    v = proj2man(staple_dag_hex(U,μ,x,t))
    U[μ,x,t] = adjoint(v *  U[μ,x,t] * v)
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

#
function insta_U2_naive_hex(N_x, N_t, Q)
    U = [coeffs_U2(NaN*im,NaN*im,NaN*im,NaN*im) for μ = 1:2, x = 1:N_x, t = 1:N_t]
    for t = 1:N_t
        for x = mod1(t,2):2:N_x
            U[1,x,t] = exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() # Note the factor 2 cf. square instantons
        end
    end
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * coeffs_Id_U2() for x = 1:N_x]
    return U
end

function insta_U2_hex(N_x, N_t, Q)
    U = [coeffs_U2(NaN*im,NaN*im,NaN*im,NaN*im) for μ = 1:2, x = 1:N_x, t = 1:N_t]
    if iseven(Q)
        for t = 1:N_t
            for x = mod1(t,2):2:N_x
                U[1,x,t] = exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() # Note the factor 2 cf. square instantons
            end
        end
        U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * coeffs_Id_U2() for x = 1:N_x]
    else
        for t = 1:N_t
            for x = mod1(t,2):2:N_x
                U[1,x,t] = exp(-im*Q*t*π/N_x/N_t) * (cos(t*π/N_x/N_t)*coeffs_Id_U2() - sin(t*π/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) 
            end
        end
        U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*π/N_x)*coeffs_Id_U2() + sin(x*π/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    end
    return U
end

# function top_charge_U2_hex(U)
#     NX = size(U,2)
#     NT = size(U,3)
#     Q = 0.0
#     for t = 1:NT
#         for x = mod1(t,2):2:NX
#             Q += imag(log(det(hexplaq(U, x, t)))) 
#         end
#     end
#     return Q / 2 / π
# end

# top_charge_U2_hex(insta_U2_hex(32,32,255))

function insta_update_U2_hex!(U,β,acc)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U2_hex(NX,NT,rand([-1,1])) .* U
    if rand() < exp(action_hex(U,β) - action_hex(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += NX*NT*1.5
    end
    return nothing
end

# insta_update() doesn't really work, but cooling the proposal
# only a couple of times (N_insta_cool) seems to do the trick
function insta_cool_update_U2_hex!(U, step, β, N_insta_cool, acc)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U2_naive(NX,NT,rand([-1,1])) .* U
    for i = 1:N_insta_cool
        chess_cool_hex!(U_prop,step,β,[0],"U2")
    end
    if rand() < exp(action_hex(U,β) - action_hex(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end
