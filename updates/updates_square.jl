include("updates_head.jl")



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
    # new_coeffs = ran_SU2(step) * U[μ,x,t]
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    end
    staple_d = staple_dag(U,μ,x,t)
    S_old = β*0.5*real(tr(U[μ,x,t] * staple_d))
    S_new = β*0.5*real(tr(new_coeffs * staple_d))
    if rand() < exp(S_new-S_old)
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function lexico_metro!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
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
    coords = [[rand(1:2), rand(1:NX), rand(1:NT)] for i = 1:2*NX*NT]
    for i = 1:2*NX*NT
        μ, x, t = coords[i]
        metro!(U,μ,x,t,step,β,acc, group)
    end    
    return nothing
end

# the overrelaxation algorithm for one SU(2)-valued link
function overrelax!(U, μ, x, t)
    v = proj_SU2(staple_dag(U,μ,x,t))
    U[μ,x,t] = adjoint(v *  U[μ,x,t] * v)
    return nothing
end


# overrelax!() link by link, lexicographically
function lexico_overrelax!(U, acc)
    NX = size(U,2)
    NT = size(U,3)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                overrelax!(U, μ, x, t)
            end
        end
    end
    acc[1] += 2*NX*NT
    return nothing
end

# overrelax!() link by link in a chess pattern
function chess_overrelax!(U, acc)
    NX = size(U,2)
    NT = size(U,3)
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    overrelax!(U,μ,x,t)
                end
            end
        end
    end
    acc[1] += 2*NX*NT
    return nothing
end

# Create a minimum of the gauge action in the topological
# sector of charge Q. Unfortunately the link values of the 
# odd-Q-instantons are not elements of the center of U(2)
# (i.e. ∉ U(1)), so that 
function insta_U2(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    if iseven(Q)
        U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * coeffs_Id_U2() for x = 1:N_x]
    else
        U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*π/N_x/N_t)*coeffs_Id_U2() - sin(t*π/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*π/N_x)*coeffs_Id_U2() + sin(x*π/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    end
    return U
end

# Create a naive (N_x × N_t)-Q-instanton configuration, in
# the sense that it is not an actual minimum of the gauge
# action for odd Q (see insta_U2() below)
function insta_U2_naive(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * coeffs_Id_U2() for x = 1:N_x]
    return U
end

# Multiply a ±1-instanton config. onto U and accept/reject
# that update with the Metropolis condition
function insta_update_U2!(U,β,acc,Q)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U2(NX,NT,rand([-Q,Q])) .* U
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end

function insta_flow_update_U2!(U, N_stout_insta, acc, Q)
    NX = size(U,2)
    NT = size(U,3)
    ρ  = 1/N_stout_insta
    U_prop = stout(U,ρ)
    for i = 1:N_stout_insta-1
        U_prop = stout(U_prop,ρ)
    end
    U_prop = insta_U2_naive(NX,NT,rand([-Q,Q])) .* U_prop
    for i = 1:N_stout_insta
        U_prop = stout(U_prop,-ρ)
    end
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end