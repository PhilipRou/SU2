include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\updates\\updates_hex.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\updates\\updates_cube.jl")



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
    S_old_neg = β*0.5*real(tr(U[μ,x,t] * staple_d))
    S_new_neg = β*0.5*real(tr(new_coeffs * staple_d))
    if S_old_neg < S_new_neg
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
    Y = 1/2*(X - adjoint(X))      # Already traceless for X ∈ SU(2)
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
                # Ω0 = ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])
                # Z0 = -1/2*(adjoint(Ω0) - Ω0)
                # V[μ,x,t] = exp_u2(Z0) * U[μ,x,t]
                V[μ,x,t] = exp_stout(ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])) * U[μ,x,t]
            end
        end
    end
    return V
end

# Single stout smearing of a given config U with parameter ρ
# function stout_slow(U, ρ)
#     NX = size(U,2)
#     NT = size(U,3)
#     V = similar(U)
#     for t = 1:NT
#         for x = 1:NX
#             for μ = 1:2
#                 # stap_link = ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])
#                 V[μ,x,t] = exp_stout_slow(ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])) * U[μ,x,t]
#             end
#         end
#     end
#     return V
# end

function stout(U, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,3)
    if n_stout > 0
        V = stout(U,ρ)
        for i = 1:n_stout-1
            V = stout(V,ρ)
        end
        return V
    else
        return U
    end
end

# function stout_midpoint_old(U, ρ)
#     NX = size(U,2)
#     NT = size(U,3)
#     V = similar(U)
#     for t = 1:NT
#         for x = 1:NX
#             for μ = 1:2
#                 stap = staple(U,μ,x,t)
#                 X1 = exp_stout(ρ/2 * stap * adjoint(U[μ,x,t])) * U[μ,x,t]
#                 V[μ,x,t] = exp_stout(ρ * stap * adjoint(X1)) * U[μ,x,t]
#                 # besser: exp(Z1 - Z0/2) * X1  ̂=  
#             end
#         end
#     end
#     return V
# end

# function stout_midpoint_old(U, n_stout, ρ)
#     NX = size(U,2)
#     NT = size(U,3)
#     if n_stout > 0
#         V = stout_midpoint_old(U,ρ)
#         for i = 1:n_stout-1
#             V = stout_midpoint_old(V,ρ)
#         end
#         return V
#     else
#         return U
#     end
# end

# function stout_midpoint(U, ρ)
#     NX = size(U,2)
#     NT = size(U,3)
#     V = similar(U)
#     for t = 1:NT
#         for x = 1:NX
#             for μ = 1:2
#                 # besser: exp(Z1 - Z0/2) * X1 
#                 stap = staple(U,μ,x,t)
#                 Ω0 = ρ * stap * adjoint(U[μ,x,t])
#                 Z0 = -1/2*(adjoint(Ω0) - Ω0)
#                 X1 = exp_u2(Z0/2) * U[μ,x,t]
#                 Ω1 = ρ * stap * adjoint(X1)
#                 Z1 = -1/2*(adjoint(Ω1) - Ω1)
#                 V[μ,x,t] = exp_u2(Z1 - Z0/2) * X1
#                 # Z0 = ρ * stap * adjoint(U[μ,x,t])
#                 # X1 = exp_stout_slow(1/2 * Z0) * U[μ,x,t]
#                 # Z1 = ρ * stap * adjoint(X1)
#                 # V[μ,x,t] = exp_stout_slow(Z1 - Z0/2) * X1
#             end
#         end
#     end
#     return V
# end

# function stout_midpoint(U, n_stout, ρ)
#     NX = size(U,2)
#     NT = size(U,3)
#     if n_stout > 0
#         V = stout_midpoint(U,ρ)
#         for i = 1:n_stout-1
#             V = stout_midpoint(V,ρ)
#         end
#         return V
#     else
#         return U
#     end
# end

# function stout_midpoint(U, ρ)
#     NX = size(U,2)
#     NT = size(U,3)
#     V = stout(U, 0.5*ρ)     # Inefficient, because staple(U,...) is calculated again below
#     W = similar(U)
#     for t = 1:NT
#         for x = 1:NX
#             for μ = 1:2
#                 # besser: exp(Z1 - Z0/2) * X1 
#                 stap0 = staple(U,μ,x,t)
#                 stap1 = staple(V,μ,x,t)
#                 Ω0 = ρ * stap0 * adjoint(U[μ,x,t])
#                 Ω1 = ρ * stap1 * adjoint(V[μ,x,t])
#                 Z0 = 0.5 * (Ω0 - adjoint(Ω0))
#                 Z1 = 0.5 * (Ω1 - adjoint(Ω1))
#                 W[μ,x,t] = exp_u2(Z1 - Z0/2) * V[μ,x,t]
#             end
#         end
#     end
#     return W
# end

# function stout_midpoint(U, n_stout, ρ)
#     NX = size(U,2)
#     NT = size(U,3)
#     if n_stout > 0
#         V = stout_midpoint(U,ρ)
#         for i = 1:n_stout-1
#             V = stout_midpoint(V,ρ)
#         end
#         return V
#     else
#         return U
#     end
# end

function stout_midpoint(U, ρ)
    NX = size(U,2)
    NT = size(U,3)
    staps = [staple(U,μ,x,t) for μ = 1:2, x = 1:NX, t = 1:NT]
    V = similar(U)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                # besser: exp(Z1 - Z0/2) * X1 
                # stap = staps[μ,x,t]
                V[μ,x,t] = exp_stout(ρ * staps[μ,x,t] * adjoint(U[μ,x,t])) * U[μ,x,t]
            end
        end
    end
    W = similar(U)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                # besser: exp(Z1 - Z0/2) * X1 
                # stap0 = staps[μ,x,t] # staple(U,μ,x,t)
                stap1 = staple(V,μ,x,t)
                # Ω0 = ρ * stap0 * adjoint(U[μ,x,t])
                Ω0 = ρ * staps[μ,x,t] * adjoint(U[μ,x,t])
                Ω1 = ρ * stap1 * adjoint(V[μ,x,t])
                Z0 = 0.5 * (Ω0 - adjoint(Ω0))
                Z1 = 0.5 * (Ω1 - adjoint(Ω1))
                W[μ,x,t] = exp_u2(Z1 - Z0/2) * V[μ,x,t]
            end
        end
    end
    return W
end
        
function stout_midpoint(U, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,3)
    if n_stout > 0
        V = stout_midpoint(U,ρ)
        for i = 1:n_stout-1
            V = stout_midpoint(V,ρ)
        end
        return V
    else
        return U
    end
end



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





function staple_dag_cube_timeslice(U,μ,x,y,t)
    NX = size(U,2)
    NT = size(U,3)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    x_p = mod1(x+1, NX) # x%NX +1                 
    y_p = mod1(y+1, NX) # t%NT +1                 
    x_m = mod1(x-1, NX) # (x + NX -2)%NX +1   
    y_m = mod1(y-1, NX) # (t + NT -2)%NT +1   

    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    if μ == 1
        a = U[2,x_p,y,t] * adjoint(U[1,x,y_p,t]) * adjoint(U[2,x,y,t])
        b = adjoint(U[2,x_p,y_m,t]) * adjoint(U[1,x,y_m,t]) * U[2,x,y_m,t]
    elseif μ == 2
        a = U[1,x,y_p,t] * adjoint(U[2,x_p,y,t]) * adjoint(U[1,x,y,t])
        b = adjoint(U[1,x_m,y_p,t]) * adjoint(U[2,x_m,y,t]) * U[1,x_m,y,t]
    else
        error("μ must be either 1 or 2 for a staple inside of the timeslice")
    end
    return a + b 
end

function stout_cube(U, ρ)
    NX = size(U,2)
    NT = size(U,4)
    V = similar(U)
    for t = 1:NT
        for y = 1:NX
            for x = 1:NX
                for μ = 1:3
                    V[μ,x,y,t] = exp_stout(ρ * adjoint(staple_dag_cube(U,μ,x,y,t)) * adjoint(U[μ,x,y,t])) * U[μ,x,y,t]
                end
            end
        end
    end
    return V
end

function stout_cube(U, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,4)
    if n_stout > 0
        V = stout_cube(U,ρ)
        for i = 1:n_stout-1
            V = stout_cube(V,ρ)
        end
        return V
    else
        return U
    end
end

function stout_cube_timeslice(U,ρ)
    NX = size(U,2)
    NT = size(U,4)
    # V = [coeffs_SU2(0.0, 0.0, 0.0, 0.0) for μ = 1:3, x = 1:NX, y = 1:NX, t = 1:NT]
    V = similar(U)
    for t = 1:NT
        for y = 1:NX
            for x = 1:NX
                for μ = 1:2
                    V[μ,x,y,t] = exp_stout(ρ * adjoint(staple_dag_cube_timeslice(U,μ,x,y,t)) * adjoint(U[μ,x,y,t])) * U[μ,x,y,t]
                end
                V[3,x,y,t] = U[3,x,y,t]
            end
        end
    end
    return V
end

# bla = gaugefield_SU2_cube(8, 8, true);
# bli = stout_cube_timeslice(bla, 0.24);
# # # coeffs_SU2(0.0, 0.0, 0.0, 0.0) in bli
# blub = stout(bla[1:2,:,:,1], 0.24);
# false in isapprox.(bli[1:2,:,:,1], blub)


function stout_cube_timeslice(U, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,4)
    if n_stout > 0
        V = stout_cube_timeslice(U,ρ)
        for i = 1:n_stout-1
            V = stout_cube_timeslice(V,ρ)
        end
        return V
    else
        return U
    end
end

function stout_midpoint_cube(U, ρ)
    NX = size(U,2)
    NT = size(U,4)
    staps = [adjoint(staple_dag_cube(U,μ,x,y,t)) for μ = 1:3, x = 1:NX, y = 1:NX, t = 1:NT]
    V = similar(U)
    for t = 1:NT
        for y = 1:NX
            for x = 1:NX
                for μ = 1:3
                    # besser: exp(Z1 - Z0/2) * X1 
                    # stap = staps[μ,x,t]
                    V[μ,x,y,t] = exp_stout(ρ * staps[μ,x,y,t] * adjoint(U[μ,x,y,t])) * U[μ,x,y,t]
                end
            end
        end
    end
    W = similar(U)
    for t = 1:NT
        for y = 1:NX
            for x = 1:NX
                for μ = 1:3
                    stap1 = adjoint(staple_dag_cube(V,μ,x,y,t))
                    Ω0 = ρ * staps[μ,x,y,t] * adjoint(U[μ,x,y,t])
                    Ω1 = ρ * stap1 * adjoint(V[μ,x,y,t])
                    Z0 = 0.5 * (Ω0 - adjoint(Ω0))
                    Z1 = 0.5 * (Ω1 - adjoint(Ω1))
                    W[μ,x,y,t] = exp_u2(Z1 - Z0/2) * V[μ,x,y,t]
                end
            end
        end
    end
    return W
end

function stout_midpoint_cube(U, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,4)
    if n_stout > 0
        V = stout_midpoint_cube(U,ρ)
        for i = 1:n_stout-1
            V = stout_midpoint_cube(V,ρ)
        end
        return V
    else
        return U
    end
end

function stout_midpoint_cube_timeslice(U,ρ)
    NX = size(U,2)
    NT = size(U,4)
    staps = [adjoint(staple_dag_cube_timeslice(U,μ,x,y,t)) for μ = 1:2, x = 1:NX, y = 1:NX, t = 1:NT]
    V = similar(U)
    for t = 1:NT
        for y = 1:NX
            for x = 1:NX
                for μ = 1:2
                    # besser: exp(Z1 - Z0/2) * X1 
                    # stap = staps[μ,x,t]
                    V[μ,x,y,t] = exp_stout(ρ * staps[μ,x,y,t] * adjoint(U[μ,x,y,t])) * U[μ,x,y,t]
                end
                V[3,x,y,t] = U[3,x,y,t]
            end
        end
    end
    W = similar(U)
    for t = 1:NT
        for y = 1:NX
            for x = 1:NX
                for μ = 1:2
                    stap1 = adjoint(staple_dag_cube_timeslice(V,μ,x,y,t))
                    Ω0 = ρ * staps[μ,x,y,t] * adjoint(U[μ,x,y,t])
                    Ω1 = ρ * stap1 * adjoint(V[μ,x,y,t])
                    Z0 = 0.5 * (Ω0 - adjoint(Ω0))
                    Z1 = 0.5 * (Ω1 - adjoint(Ω1))
                    W[μ,x,y,t] = exp_su2(Z1 - Z0/2) * V[μ,x,y,t]
                end
                W[3,x,y,t] = V[3,x,y,t]
            end
        end
    end
    return W
end

function stout_midpoint_cube_timeslice(U, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,4)
    if n_stout > 0
        V = stout_midpoint_cube_timeslice(U,ρ)
        for i = 1:n_stout-1
            V = stout_midpoint_cube_timeslice(V,ρ)
        end
        return V
    else
        return U
    end
end