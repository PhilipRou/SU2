# using LinearAlgebra

# include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields.jl")


# ⭕
# ⭕
# ⭕
# ⭕    Implement!!!
# ⭕    Different !!!
# ⭕    Coordinate!!!
# ⭕    System!!!
# ⭕
# ⭕
# ⭕



# # Function to return the non-daggered staple as opposed to staple_dag in SU2_updates.jl
function staple(U, μ, x, t)
    NT = size(U,2)
    NX = size(U,3)
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

# Function to make X traceless and return the exponential of that
function exp_traceless(X::coeffs_SU2)
    X = X - adjoint(X)  # Make X traceless (works 'cause X ∈ SU(2))
    w = sqrt(X.b^2 + X.c^2 + X.d^2)
    s = sin(w)/w
    return coeffs_SU2(cos(w), s*X.b, s*X.c, s*X.d)
end


# ❗  WARNING: this stout smearing uses the fact that for M ∈ SU(2) one can 
#             already construct a traceless matrix by virtue of M - M†.
function stout_single(U, ρ)
    NX = size(U,2)
    NT = size(U,3)
    V = similar(U)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                stap_link = staple(U,μ,x,t) * adjoint(U[μ,x,t])
                stap_link = 0.5*ρ*stap_link
                V[μ,x,t] = exp_traceless(stap_link) * U[μ,x,t]
                # stap_link = coeffs2grp(staple(V,μ,x,t) * adjoint(V[μ,x,t]))
                # V[μ,x,t] = grp2coeffs(exp(α*0.5*(stap_link - adjoint(stap_link)))) * V[μ,x,t]
            end
        end
    end
    return V
end

# ⭕ Definitely reconsider this V = similar(U) bullcrap ⭕
function stout(U, n_stout, ρ)
    V = similar(U)
    if n_stout > 0
        V = stout_single(U,ρ)
        for i = 1:n_stout-1
            V = stout_single(V,ρ)
        end
        return V
    else
        return U
    end
end


# # N_t = 8
# # N_x = 8
# # V = gaugefield_SU2(N_t,N_x,true)
# # μ = rand(1:2)
# # t = rand(1:N_t)
# # x = rand(1:N_x)
# # stap_link = coeffs2grp(staple(V,t,x,μ) * adjoint(V[μ,t,x]))
# # # V[μ,t,x] = grp2coeffs(exp(α*0.5*(stap_link-adjoint(stap_link) - 0.5*tr(stap_link-adjoint(stap_link)))) * coeffs2grp(V[μ,t,x]))
# # V[μ,t,x] = grp2coeffs(exp(0.1*0.5*(stap_link - adjoint(stap_link) - 0.5*tr(stap_link - adjoint(stap_link)) * σ0))) * V[μ,t,x]

# # stap_link - adjoint(stap_link) - 0.5*tr(stap_link - adjoint(stap_link)) * σ0





# ########    Hexagonal Stuff    ########





# function staple_hex(U, μ, t, x)
#     # NX = N_x>>1
#     NT = size(U,2)
#     NX = size(U,3)
#     a = coeffs_SU2(0.0,0.0,0.0,0.0)
#     b = coeffs_SU2(0.0,0.0,0.0,0.0)
#     t_p = mod1(t+1, NT) # t%NT +1                 
#     x_p = mod1(x+1, NX) # x%NX +1                 
#     t_m = mod1(t-1, NT) # (t + NT -2)%NT +1   
#     x_m = mod1(x-1, NX) # (x + NX -2)%NX +1  
#     t_pp = mod1(t+2, NT) # (t+1)%NT +1
#     t_mm = mod1(t-2, NT) # (t + NT - 3)%NT +1

#     # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
#     if μ == 1
#         a = adjoint([U[3,t_p,x]]) * U[2,t_p,x] * U[1,t_pp,x_p] * U[3,t_pp,x_p] * adjoint(U[2,t,x])
#         b = adjoint(U[2,t_m,x_m]) * U[3,t_m,x_m] * U[1,t_mm,x_m] * U[2,t_mm,x_m] * adjoint(U[3,t,x])
#         # a = U[2,t,x] * adjoint(U[3,t_pp,x_p]) * adjoint(U[1,t_pp,x_p]) * adjoint(U[2,t_p,x]) * U[3,t_p,x]
#         # b = U[3,t,x] * adjoint(U[2,t_mm,x_m]) * adjoint(U[1,t_mm,x_m]) * adjoint(U[3,t_m,x_m]) * U[2,t_m,x_m]
#     elseif μ == 2
#         a = U[3,t,x] * U[1,t_m,x] * U[2,t_m,x] * adjoint(U[3,t_p,x_p]) * adjoint(U[1,t_p,x_p]) 
#         b = adjoint(U[1,t,x]) * adjoint(U[3,t_p,x]) * U[2,t_p,x] * U[1,t_pp,x_p] * U[3,t_pp,x_p]
#         # a = U[1,t_p,x_p] * U[3,t_p,x_p] * adjoint(U[2,t_m,x]) * adjoint(U[1,t_m,x]) * adjoint(U[3,t,x]) 
#         # b = adjoint(U[3,t_pp,x_p]) * adjoint(U[1,t_pp,x_p]) * adjoint(U[2,t_p,x]) * U[3,t_p,x] * U[1,t,x]
#     elseif μ == 3
#         a = U[2,t,x] * U[1,t_p,x_p] * U[3,t_p,x_p] * adjoint(U[2,t_m,x]) * adj(U[1,t_m,x])
#         b = adjoint(U[1,t,x]) * adjoint(U[2,t_m,x_m]) * U[3,t_m,x_m] * U[1,t_mm,x_m] * U[2,t_mm,x_m]
#         # a = U[1,t_m,x] * U[2,t_m,x] * adjoint(U[3,t_p,x_p]) * adjoint(U[1,t_p,x_p]) * adjoint(U[2,t,x])
#         # b = adjoint(U[2,t_mm,x_m]) * adjoint(U[1,t_mm,x_m]) * adjoint(U[3,t_m,x_m]) * U[2,t_m,x_m] * U[1,t,x]
#     end

#     return a+b
# end



# # ❗  WARNING: this stout smearing uses the fact that for M ∈ SU(2) one can 
# #             already construct a traceless matrix by virtue of M - M†.
# function stout_single_hex(U, ρ)
#     NT = size(U,2)
#     NX = size(U,3)
#     V = similar(U)
#     for μ = 1:3
#         for t = 1:NT
#             for x = 1:NX
#                 stap_link = staple_hex(U,μ,t,x) * adjoint(U[μ,t,x])
#                 stap_link = 0.5*ρ*stap_link
#                 V[μ,t,x] = exp_traceless(stap_link) * U[μ,t,x]
#                 # stap_link = coeffs2grp(staple(V,μ,t,x) * adjoint(V[μ,t,x]))
#                 # V[μ,t,x] = grp2coeffs(exp(α*0.5*(stap_link - adjoint(stap_link)))) * V[μ,t,x]
#             end
#         end
#     end
#     return V
# end

# function stout_hex(U, n_stout, ρ)
#     V = similar(U)
#     if n_stout > 0
#         V = stout_single(U,ρ)
#         for i = 1:n_stout-1
#             V = stout_single_hex(V,ρ)
#         end
#         return V
#     else
#         return U
#     end
# end





########    3-dimensional Stuff    ########





