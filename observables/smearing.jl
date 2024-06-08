include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_hex.jl")



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

# We cool (🆒) a gauge config by performing a Metropolis
# update, which we only accept if the proposal has lower 
# action.
function cool!(U, μ, x, t, step, β, acc, group)
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    end
    staple_d = staple_dag(U,μ,x,t)
    S_old = β*0.5*real(tr(U[μ,x,t] * staple_d))
    S_new = β*0.5*real(tr(new_coeffs * staple_d))
    if S_old < S_new
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

# Cooling in chequer board pattern
function chess_cool!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    cool!(U,μ,x,t,step,β,acc,group)
                end
            end
        end
    end
    return nothing
end

function ran_cool!(U, cool_degree_in_percent, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    n_cool = round(Int, cool_degree_in_percent*2*NX*NT/100, RoundNearestTiesAway)
    for i = 1:n_cool
        cool!(U,rand(1:2),rand(1:NX),rand(1:NT),step,β,acc,group)
    end
    return nothing
end

# ins = insta_U2(32,32,1);
# ins = gaugefield(32,32,true,"U2","square");
# bla = action(ins,12)
# acc = [0]
# acc_insta = [0]
# for i = 1:100
#     chess_metro!(ins,0.03,12,acc,"U2")
#     insta_update_U2!(ins,β,acc_insta)
#     # chess_cool!(ins,0.03,12,acc,"U2")
# end
# acc_insta[1]/2/32/32/100
# bla - action(ins,12)

# for i = -4:4
#     println(action(insta_U2(32,32,i),12))
# end


################################################################################
# Function to calculate exp(iQ_μ) for the stout smearing procedure (cf. 
# Gattringer/Lang p. 143). The argument is the product of a staple and the 
# corresponding link variable, i.e. a sum of group elements. For SU(2) and U(2)
# the exponential of that quantity can be evaluated without bothering the 
# exponential function, as exp(ianⱼσⱼ) = cos(a)σ₀ + i sin(a) nⱼσⱼ.
function exp_stout(X::coeffs_SU2)
    Y = X - adjoint(X)      # Already traceless for X ∈ SU(2)
    w = sqrt(Y.b^2 + Y.c^2 + Y.d^2)
    s = sin(w)/w
    return coeffs_SU2(cos(w), s*Y.b, s*Y.c, s*Y.d)
end

# # WRONG!!!
# function exp_stout_old(X::coeffs_U2)
#     Q_μ = -1/2*( adjoint(X) - X - tr(adjoint(X) - X)/2 * coeffs_Id_U2() )
#     w = real(sqrt(Q_μ.b^2 + Q_μ.c^2 + Q_μ.d^2))
#     s = sin(w)/w
#     return coeffs_U2(cos(w)+0.0*im, s*Q_μ.b, s*Q_μ.c, s*Q_μ.d)
# end

function exp_stout(X::coeffs_U2)
    Q_μ = -1/2*(adjoint(X) - X) # Not traceless, but intended! It's the Lie algebra of U(2), after all
    return exp_u2(Q_μ)
end

function exp_stout_slow(X::coeffs_U2)
    Q_μ = -1/2*(adjoint(X) - X) 
    return grp2coeffs_U2(exp(coeffs2grp(Q_μ)))
end

# X = ran_U2(rand()) + ran_U2(rand());
# Q_μ = -1/2*( adjoint(X) - X - tr(adjoint(X) - X)/2 * coeffs_Id_U2() );
# ble = grp2coeffs_U2(exp(coeffs2grp(Q_μ))) - exp_traceless_b(X);
# println(abs(ble.a))
# println(abs(ble.b))
# println(abs(ble.c))
# println(abs(ble.d))




# Single stout smearing of a given config U with parameter ρ
function stout(U, ρ)
    NX = size(U,2)
    NT = size(U,3)
    V = similar(U)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                # stap_link = ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])
                V[μ,x,t] = exp_stout(ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])) * U[μ,x,t]
            end
        end
    end
    return V
end

# Single stout smearing of a given config U with parameter ρ
function stout_slow(U, ρ)
    NX = size(U,2)
    NT = size(U,3)
    V = similar(U)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                # stap_link = ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])
                V[μ,x,t] = exp_stout_slow(ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])) * U[μ,x,t]
            end
        end
    end
    return V
end

function stout_midpoint(U, ρ)
    NX = size(U,2)
    NT = size(U,3)
    V = similar(U)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                stap = staple(U,μ,x,t)
                temp = exp_stout(ρ/2 * stap * adjoint(U[μ,x,t])) * U[μ,x,t]
                V[μ,x,t] = exp_stout(ρ * stap * adjoint(temp)) * U[μ,x,t]
            end
        end
    end
    return V
end

# function stout(U, n_stout, ρ)
#     NX = size(U,2)
#     NT = size(U,3)
#     if n_stout > 0
#         V = stout(U,ρ)
#         for i = 1:n_stout-1
#             V = stout(V,ρ)
#         end
#         return V
#     else
#         return U
#     end
# end

        

# # ⭕ Definitely reconsider this V = similar(U) bullcrap ⭕
# function stout(U, n_stout, ρ)
#     V = similar(U)
#     if n_stout > 0
#         V = stout(U,ρ)
#         for i = 1:n_stout-1
#             V = stout(V,ρ)
#         end
#         return V
#     else
#         return U
#     end
# end


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





function cool_hex!(U, μ, x, t, step, β, acc, group)
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    end
    staple_d = staple_dag_hex(U,μ,x,t)
    S_old = β*0.5*real(tr(U[μ,x,t] * staple_d))
    S_new = β*0.5*real(tr(new_coeffs * staple_d))
    if S_old < S_new
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end



function chess_cool_hex!(U, step, β, acc, group)
    NX = size(U,2)
    NT = size(U,3)
    μ = 2
    for trip = 1:2
        for start_t = 1:2
            for t = start_t:2:NT
                for x = (1+mod(t+trip,2)):2:NX
                    cool_hex!(U,μ,x,t,step,β,acc, group)
                end
            end
        end
    end
    μ = 1
    for t = 1:NT
        for x = 2-mod(t,2):2:NX
            cool_hex!(U,μ,x,t,step,β,acc, group)
        end
    end
    return nothing
end

# bla = gaugefield(32,32,true,"U2","hexagon");
# for i = 1:200 chess_metro_hex!(bla,0.2,1,[0],"U2") end
# a1 = action_hex(bla,1)
# cool_hex!(bla,rand(1:2),rand(1:32),rand(1:32),0.2,1.0,[0],"U2")
# a2 = action_hex(bla,1)

# function staple_hex(U, μ, t, x)
#     # NX = N_x>>1
#     NX = size(U,3)
#     NX = size(U,2)
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
# function stout_hex(U, ρ)
#     NX = size(U,3)
#     NX = size(U,2)
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
#         V = stout(U,ρ)
#         for i = 1:n_stout-1
#             V = stout_hex(V,ρ)
#         end
#         return V
#     else
#         return U
#     end
# end





########    3-dimensional Stuff    ########





