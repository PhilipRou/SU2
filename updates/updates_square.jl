include("updates_head.jl")



#
function staple_dag(U, μ, x, t)
    NX = size(U,2)
    NT = size(U,3)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    x_p = mod1(x+1, NX) # x%NX +1                 
    t_p = mod1(t+1, NT) # t%NT +1                 
    x_m = mod1(x-1, NX) # (x + NX -2)%NX +1   
    t_m = mod1(t-1, NT) # (t + NT -2)%NT +1   

    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    if μ == 1
        a = U[2,x_p,t] * adjoint(U[1,x,t_p]) * adjoint(U[2,x,t])
        b = adjoint(U[2,x_p,t_m]) * adjoint(U[1,x,t_m]) * U[2,x,t_m]
    else #if μ == 2
        a = U[1,x,t_p] * adjoint(U[2,x_p,t]) * adjoint(U[1,x,t])
        b = adjoint(U[1,x_m,t_p]) * adjoint(U[2,x_m,t]) * U[1,x_m,t]
    end
    return a + b 
end

# Function to return the non-daggered staple as opposed to staple_dag in SU2_updates.jl
function staple(U, μ, x, t)
    NX = size(U,2)
    NT = size(U,3)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    t_p = mod1(t+1, NT) # t%NT +1                 
    x_p = mod1(x+1, NX) # x%NX +1                 
    t_m = mod1(t-1, NT) # (t + NT -2)%NT +1   
    x_m = mod1(x-1, NX) # (x + NX -2)%NX +1  

    if μ == 1
        # a = U[2,x_p,t] * adjoint(U[1,x,t_p]) * adjoint(U[2,x,t])
        # b = adjoint(U[2,x_p,t_m]) * adjoint(U[1,x,t_m]) * U[2,x,t_m]
        a = U[2,x,t] * U[1,x,t_p] * adjoint(U[2,x_p,t])
        b = adjoint(U[2,x,t_m]) * U[1,x,t_m] * U[2,x_p,t_m]
    else #if μ == 2
        # a = U[1,x,t_p] * adjoint(U[2,x_p,t]) * adjoint(U[1,x,t])
        # b = adjoint(U[1,x_m,t_p]) * adjoint(U[2,x_m,t]) * U[1,x_m,t]
        a = U[1,x,t] * U[2,x_p,t] * adjoint(U[1,x,t_p])
        b = adjoint(U[1,x_m,t]) * U[2,x_m,t] * U[1,x_m,t_p]
    end
    return a + b
end

# ❌ Not in use anymore! Only for testing purposes ❌
function delta_S_gauge(U, μ, x, t, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2, β)
    return β*0.5*real(tr((old_coeffs - new_coeffs) * staple_dag(U,μ,x,t)))
end

function metro!(U, μ, x, t, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    end
    staple_d = staple_dag(U,μ,x,t)
    S_old = -β*0.5*real(tr(U[μ,x,t] * staple_d))
    S_new = -β*0.5*real(tr(new_coeffs * staple_d))
    if rand() < exp(-(S_new-S_old))
        U[μ,x,t] = new_coeffs
        acc[1] += 1/NX/NT/2
    end
    return nothing
end

function lexico_metro!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                metro!(U,μ,x,t,step,β,acc,group)
            end
        end
    end
    return nothing
end

function chess_metro!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    metro!(U,μ,x,t,step, β, acc, group)
                end
            end
        end
    end
    return nothing
end

#
function ran_metro!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    coords = [[rand(1:2), rand(1:NX), rand(1:NT)] for i = 1:2*NX*NT]
    for i = 1:2*NX*NT
        μ, x, t = coords[i]
        metro!(U,μ,x,t,step,β,acc,group)
    end    
    return nothing
end

# the overrelaxation algorithm for one SU(2)-valued link
function overrelax!(U, μ, x, t, acc)
    NX = size(U,2)
    NT = size(U,3)
    v = proj2man(staple_dag(U,μ,x,t))
    # println(typeof(U[μ,x,t]))
    if typeof(U[μ,x,t]) == coeffs_SU2{Float64}
        U[μ,x,t] = adjoint(v * U[μ,x,t] * v)
        acc[1] += 1/2/NX/NT
    elseif typeof(U[μ,x,t]) == coeffs_U2{ComplexF64}
        new_coeffs = adjoint(v * U[μ,x,t] * v)
        staple_d = staple_dag(U,μ,x,t)
        S_old = β*0.5*real(tr(U[μ,x,t] * staple_d))
        S_new = β*0.5*real(tr(new_coeffs * staple_d))
        if rand() < exp(S_old-S_new)
            U[μ,x,t] = new_coeffs
            acc[1] += 1/2/NX/NT
        end
    end
    return nothing
end


# overrelax!() link by link, lexicographically
function lexico_overrelax!(U, acc)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                overrelax!(U,μ,x,t,acc)
            end
        end
    end
    # acc[1] += 1
    return nothing
end

# overrelax!() link by link in a chess pattern
function chess_overrelax!(U, acc)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    overrelax!(U,μ,x,t,acc)
                end
            end
        end
    end
    # acc[1] += 1
    return nothing
end

# Multiply a ±1-instanton config. onto U and accept/reject
# that update with the Metropolis condition
function insta_update_U2!(U,β,acc,Q)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    U_prop = insta_U2(NX,NT,rand([-Q,Q])) .* U
    ΔS = action(U_prop,β) - action(U,β)
    # if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
    if rand() < exp(-ΔS) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 1
    end
    # return nothing
    return action(U_prop, β)
end

# Multiply a ±1-instanton config. onto U and accept/reject
# that update with the Metropolis condition
function insta_update_U2_cheat!(U,β,acc,Q)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    U_prop = insta_U2_naive(NX,NT,rand([-Q,Q])) .* U
    ΔS = action(U_prop,β) - action(U,β)
    # if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
    if rand() < exp(-ΔS) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 1
    end
    # return nothing
    return action(U_prop, β)
end

# Take the logarithm of each link, then link-wise add
# the logarithm of an instanton config. onto it,
# then project back onto the manifold. In other words, an 
# instanton-update in the Lie algebra.
# "bf" stands for brute force!
function insta_update_U2_log_bf!(U,β,acc,Q)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    U_prop = grp2coeffs_U2.(exp.(log.(coeffs2grp.(U)) .+ log.(coeffs2grp.(insta_U2(N_x,N_t,rand([-Q,Q]))))))
    # ΔS = action(U_prop,β) - action(U,β)
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
    # if rand() < exp(-ΔS) # Definitely old minus new!!
        # println("Here we are")
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 1.0
    end
    # return nothing
    return action(U_prop, β)
end

function log_quick(X::coeffs_U2)
    M = coeffs2grp(X)
    E = eigen(M)
    A = E.vectors
    A_inv = adjoint(E.vectors)
    V = diagm(log.(E.values))
    return A * V * A_inv
end

function insta_update_U2_log_bf_quick!(U,β,acc,Q)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    U_prop = grp2coeffs_U2.(exp.(log_quick.(U) .+ log_quick.(insta_U2(N_x,N_t,rand([-Q,Q])))))
    # ΔS = action(U_prop,β) - action(U,β)
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
    # if rand() < exp(-ΔS) # Definitely old minus new!!
        # println("Here we are")
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 1.0
    end
    # return nothing
    return action(U_prop, β)
end

# # Take the logarithm of each link, then link-wise add
# # the logarithm of an instanton config. onto it,
# # then project back onto the manifold. In other words, an 
# # instanton-update in the Lie algebra.
# # No more brute force, but improved versions of log & exp.
# function insta_update_U2_log!(U,β,acc,Q)
#     NX = size(U,2)
#     NT = size(U,3)
#     acc[1] = 0.0
#     U_prop = exp_u2.(log_U2.(U) .+  insta_U2_log(NX,NT,rand([-Q,Q]))) 
#     # ΔS = action(U_prop,β) - action(U,β)
#     if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
#     # if rand() < exp(-ΔS) # Definitely old minus new!!
#         U[:,:,:] = U_prop[:,:,:]
#         acc[1] += 1.0
#     end
#     # return nothing
#     return action(U_prop,β)
# end
function insta_update_U2_log!(U,β,acc,Q)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    U_prop = exp_u2.(log_U2.(U) .+  insta_U2_log(NX,NT,rand([-Q,Q]))) 
    # ΔS = action(U_prop,β) - action(U,β)
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
    # if rand() < exp(-ΔS) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 1.0
    end
    # return nothing
    return action(U_prop,β)
end

function insta_update_U2_log_cheat!(U,β,acc,Q)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    U_prop = exp_u2.(log_U2.(U) .+  insta_U2_log_cheat(NX,NT,rand([-Q,Q]))) 
    # ΔS = action(U_prop,β) - action(U,β)
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
    # if rand() < exp(-ΔS) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 1.0
    end
    # return nothing
    return action(U_prop, β)
end

# 🚧👷 Under construction! 👷🚧
# function insta_flow_update_U2!(U, N_stout_insta, acc, Q)
#     NX = size(U,2)
#     NT = size(U,3)
#     ρ  = 1/N_stout_insta
#     U_prop = stout(U,ρ)
#     for i = 1:N_stout_insta-1
#         U_prop = stout(U_prop,ρ)
#     end
#     U_prop = insta_U2_naive(NX,NT,rand([-Q,Q])) .* U_prop
#     for i = 1:N_stout_insta
#         U_prop = stout(U_prop,-ρ)
#     end
#     if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
#         U[:,:,:] = U_prop[:,:,:]
#         acc[1] += 1
#     end
#     return nothing
# end


# U = gaugefield_U2(N_x, N_t, true);
# U_prop_bf = grp2coeffs_U2.(exp.(log.(coeffs2grp.(U)) .+ log.(coeffs2grp.(insta_U2(N_x,N_t,1)))));
# U_prop_ez = exp_u2.(log_U2.(U) .+  insta_U2_log(N_x,N_t,1)) ;
# bla = isapprox.(U_prop_bf, U_prop_ez)
# # for b in bla
# #     if b == false
# #         bla_count[1] += 1
# #     end
# # end
# # bla_count
# for t = 1:N_t
#     for x = 1:N_x
#         for μ = 1:2
#             if bla[μ,x,t] == false
#                 println(μ, ", ", x, ", ", t )
#                 println(U_prop_bf[μ,x,t])
#                 println(U_prop_ez[μ,x,t])
#                 println(" ")
#             end
#         end
#     end
# end

function metro_mat!(U, μ, x, t, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = coeffs2grp(ran_SU2(step)) * new_coeffs
    elseif group == "U2"
        new_coeffs = coeffs2grp(ran_U2(step)) * new_coeffs
    end
    staple_d = staple_dag(U,μ,x,t)
    S_old = β*0.5*real(tr(U[μ,x,t] * staple_d))
    S_new = β*0.5*real(tr(new_coeffs * staple_d))
    if rand() < exp(S_new-S_old)
        U[μ,x,t] = new_coeffs
        acc[1] += 1/NX/NT/2
    end
    return nothing
end

function chess_metro_mat!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    metro_mat!(U,μ,x,t,step, β, acc, group)
                end
            end
        end
    end
    return nothing
end