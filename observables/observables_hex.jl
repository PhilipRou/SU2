include("observables_head.jl")



#=
# A function that returns the plaquette of a hexagonal lattice. 
# Do we want a volume factor in there❓❓❓
# Mind the coordinates, i.e. the input of x!!!
function hexplaq(U, t, x)
    # NX = N_x>>1     # Mind the coordinates❗
    NX = size(U,3)
    NX = size(U,2)
    t_p = t%NT + 1
    x_p = x%NX + 1
    t_pp = (t+1)%NT + 1
    # n n a a a n
    # 123123
    # t, t, t_pp, t_pp, t_p, t_p
    # x, x, x_p, x_p, x, x 
    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    return U[1,x,t] * U[2,x,t] * adjoint(U[3,t_pp,x_p]) * adjoint(U[1,t_pp,x_p]) * adjoint(U[2,t_p,x]) * U[3,t_p,x]
end
=#

#=
# Action of a hexagonal gauge field. Do we want to sum over every
# plaquette once, or over every space-time point? I think it has
# to be the former, innit?
function action_hex(U, β)
    NX = size(U,3)
    NX = size(U,2)
    S = 2*NT*NX    # later generalization: N_colour * NT * (NX)^d_s
    for t = 1:NT
        for x = 1:NX
            S -= tr(hexplaq(U,x,t))
        end
    end
    return β*S/2    # later generalization: β*S/N_colour
end
=#

# Plaquette on a hexagonal lattice, starting point lowest left corner
function hexplaq(U,x,t)
    NX = size(U,2)
    NT = size(U,3)
    xp = mod1(x+1, NX) # x%NX +1                 
    tp = mod1(t+1, NT) # t%NT +1                 
    # xm = mod1(x-1, NX) # (x + NX -2)%NX +1  
    # tm = mod1(t-1, NT) # (t + NT -2)%NT +1   
    # tmm = mod1(t-2, NT) # (t + NT - 3)%NT +1
    tpp = mod1(t+2, NT) # (t+1)%NT +1
    return U[1,x,t]*U[2,xp,t]*U[2,xp,tp]*adjoint(U[1,x,tpp])*adjoint(U[2,x,tp])*adjoint(U[2,x,t])
end

function action_hex(U,β)
    NX = size(U,2)
    NT = size(U,3)
    S = NX*NT
    for t = 1:NT
        for x = (1+mod(t+1,2)):2:NX
            S -= real(tr(hexplaq(U,x,t)))
        end
    end
    return β*S/2
end


function loop_1x2_hex(U,x,t)
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX)
    # xp2 = mod1(x+2,NX)
    # xp3 = mod1(x+3,NX)
    # xp4 = mod1(x+4,NX)
    tp1 = mod1(t+1,NT)
    tp2 = mod1(t+2,NT)
    tp3 = mod1(t+3,NT)
    tp4 = mod1(t+4,NT)
    # n   n   n   n   n   |   a   a   a   a   a
    # 1   2   2   2   2   |   1   2   2   2   2
    # x   xp1 xp1 xp1 xp1 |   x   x   x   x   x
    # t   t   tp1 tp2 tp3 |   tp4 tp3 tp2 tp1 t
return U[1,x,t] * U[2,xp1,t] * U[2,xp1,tp1] * U[2,xp1,tp2] * U[2,xp1,tp3] * adjoint(U[1,x,tp4]) * adjoint(U[2,x,tp3]) * adjoint(U[2,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
end

# For a spatial Wilson line we'll always start with μ=1 and then μ=2, 
# meaning this way: _/⁻\_/⁻\_/⁻   or   _/⁻\_/⁻\_/⁻\_
function loop_2x1_hex(U,x,t)
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX)
    xp2 = mod1(x+2,NX)
    # xp3 = mod1(x+3,NX)
    # xp4 = mod1(x+4,NX) 
    tp1 = mod1(t+1,NT) 
    tp2 = mod1(t+2,NT) 
    tp3 = mod1(t+3,NT)
    # tp4 = mod1(t+4,NT) 
    # n   n   n   n   n   |   a   a   a   a   a   
    # 1   2   1   2   2   |   1   2   1   2   2 
    # x   xp1 xp1 xp2 xp2 |   xp1 xp1 x   x   x
    # t   t   tp1 tp1 tp2 |   tp3 tp2 tp2 tp1 t
    return U[1,x,t] * U[2,xp1,t] * U[1,xp1,tp1] * U[2,xp2,tp1] * U[2,xp2,tp2] * adjoint(U[1,xp1,tp3]) * adjoint(U[2,xp1,tp2]) * adjoint(U[1,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
end

function loop_2x2_hex(U,x,t)
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX) 
    xp2 = mod1(x+2,NX) 
    # xp3 = mod1(x+3,NX) 
    # xp4 = mod1(x+4,NX) 
    tp1 = mod1(t+1,NT) 
    tp2 = mod1(t+2,NT) 
    tp3 = mod1(t+3,NT) 
    tp4 = mod1(t+4,NT) 
    tp5 = mod1(t+5,NT) 
    # n   n   n   n   n   n   n   |   a   a   a   a   a   a   a 
    # 1   2   1   2   2   2   2   |   1   2   1   2   2   2   2
    # x   xp1 xp1 xp2 xp2 xp2 xp2 |   xp1 xp1 x   x   x   x   x   
    # t   t   tp1 tp1 tp2 tp3 tp4 |   tp5 tp4 tp4 tp3 tp2 tp1 t
    bla =  U[1,x,t] * U[2,xp1,t] * U[1,xp1,tp1] * U[2,xp2,tp1] * U[2,xp2,tp2]  * U[2,xp2,tp3]  * U[2,xp2,tp4] 
    return bla * adjoint(U[1,xp1,tp5]) * adjoint(U[2,xp1,tp4]) * adjoint(U[1,x,tp4]) * adjoint(U[2,x,tp3]) * adjoint(U[2,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
end



# Returns an (Nₜ × Nₓ)-matrix whose entries carry the rectangular (lₜ × lₓ)-loop 
# at the respective (x,t)-points of the lattice. 
function loop_mat(U, l_t, l_x)
    NX = size(U,2)
    NT = size(U,3)
    res = Matrix{coeffs_SU2}(undef, NT, NX)
    for j = 1:NX
        for i = 1:NT
            res[i,j] = coeffs_Id_SU2()
        end
    end
    t_arr = collect(1:NT)
    x_arr = collect(1:NX)
    for i = 1:l_x
        res = res .* U[2,t_arr,x_arr]
        circshift!(x_arr,-(1))    # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    end
    for i = 1:l_t
        res = res .* U[1,t_arr,x_arr]
        circshift!(t_arr,-(1))   
    end
    for i = 1:l_x
        circshift!(x_arr,-(-1))
        res = res .* adjoint.(U[2,t_arr,x_arr])
    end
    for i = 1:l_t
        circshift!(t_arr,-(-1))
        res = res .* adjoint.(U[1,t_arr,x_arr])
    end
    return res
end


################################################################################
# The analogue to loop_mat for hexagonal lattices. Note that we cannot describe
# Wilson loops via their space-time points anymore as here a loop has to 
# encompass the outline of hexagons. Hence in this function the parameters have
# a different meaning:
#       l_x ̂= number of hexagons in x-direction
#       l_t ̂= number of hexagons in t-direction
# As stated above: for a spatial Wilson line, going from left to right, we'll 
# always start with μ=1 and then μ=2, meaning this way:
#        _/⁻\_/⁻\_/⁻   or   _/⁻\_/⁻\_/⁻\_  
# We then go upwards for an even number of links and return in the obvious
# fashion to the starting point.
#=
function loop_mat_hex(U, l_x, l_t)
    NX = size(U,2)
    NT = size(U,3)
    LX = l_x -1
    x_arr = collect(1:NX)
    t_arr = collect(1:NT)

    res = U[1,x_arr,t_arr]
    circshift!(x_arr,-1)
    for i = 1:div(LX,2)
        res = res .* U[2,x_arr,t_arr] .* U[1,x_arr,circshift(t_arr,-1)] .* adjoint.(U[2,circshift(x_arr,-1),t_arr]) .* U[1,circshift(x_arr,-1),t_arr]
        circshift!(x_arr,-2)     # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    end
    for i = 1:Int(LX - 2*div(LX,2))
        res = res .* U[2,x_arr,t_arr] .* U[1,x_arr,circshift(t_arr,-1)]
        circshift!(x_arr,-1)
        circshift!(t_arr,-1)
    end

    for i = 1:l_t
        res = res .* U[2,x_arr,t_arr] .* U[2,x_arr,circshift(t_arr,-1)]
        circshift!(t_arr,-2)
    end

    for i = 1:Int(LX - 2*div(LX,2))
        res = res .* adjoint.(U[1,circshift(x_arr,1),t_arr]) .* adjoint.(U[2,circshift(x_arr,1),circshift(t_arr,1)])
        circshift!(t_arr,1)
        circshift!(x_arr,1)
    end
    for i = 1:div(LX,2)
        res = res .* adjoint.(U[1,circshift(x_arr,1),t_arr]) .* U[2,circshift(x_arr,1),t_arr] .* adjoint.(U[1,circshift(x_arr,2),circshift(t_arr,-1)]) .* adjoint.(U[2,circshift(x_arr,2),t_arr])
        circshift!(x_arr,2)     
    end
    circshift!(x_arr,1)
    res = res .* adjoint.(U[1,x_arr,t_arr])
    
    for i = 1:l_t
        res = res .* adjoint.(U[2, x_arr, circshift(t_arr,1)]) .* adjoint.(U[2, x_arr, circshift(t_arr,2)])
        circshift!(t_arr,2)
    end

    return res
end
=#

# N_x = 32
# N_t = 32
# β   = 1.0
# H = hexfield_SU2(N_x, N_t, true);
# mat1 = loop_mat_hex(H, 3, 5);
# mat2 = loop_mat_hex(temp_gauge_hex(H), 3, 5);
# 0 in isapprox.(mat1, mat2)


# Compute an R×T-Wilson loop at point (x,t), i.e. W(R,T)
function RT_loop_hex(U, R, T, x, t)
    NX = size(U,2)
    NT = size(U,3)
    LX = R-1
    
    res = U[1,x,t]
    x = mod1(x+1,NX)
    for i = 1:div(LX,2)
        xp = mod1(x+1,NX)
        tp = mod1(t+1,NT)
        res = res * U[2,x,t] * U[1,x,tp] * adjoint(U[2,xp,t]) * U[1,xp,t]
        x = mod1(x+2,NX)
        
    end
    for i = 1:Int(LX - 2*div(LX,2))
        tp = mod1(t+1,NT)
        res = res * U[2,x,t] * U[1,x,tp]
        x = mod1(x+1,NX)
        t = tp
    end

    for i = 1:T
        tp = mod1(t+1,NT)
        res = res * U[2,x,t] * U[2,x,tp]
        t = mod1(t+2,NT)
    end

    for i = 1:Int(LX - 2*div(LX,2))
        xm = mod1(x-1,NX)
        tm = mod1(t-1,NT)
        res = res * adjoint(U[1,xm,t]) * adjoint(U[2,xm,tm])
        x = xm
        t = tm
    end
    for i = 1:div(LX,2)
        xm = mod1(x-1,NX)
        xmm = mod1(x-2,NX)
        tp = mod1(t+1,NT)
        res = res * adjoint(U[1,xm,t]) * U[2,xm,t] * adjoint(U[1,xmm,tp]) * adjoint(U[2,xmm,t])
        x = xmm
    end
    x = mod1(x-1,NX)
    res = res * adjoint(U[1,x,t])
    
    for i = 1:T
        tm = mod1(t-1,NT)
        tmm = mod1(t-2,NT)
        res = res * adjoint(U[2,x, tm]) * adjoint(U[2,x, tmm])
        t = tmm
    end
    return res
end    

# ⭕ Revisit when stout implemented! ⭕
#
function measure_loops_hex(U, loops::Array, coords)#, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,3)
    L = length(loops)
    results = [tr(RT_loop_hex(U, loop[1], loop[2], coord[1], coord[2])) for loop in loops, coord in coords]
    mean_vals = [sum(results[i,:]) for i = 1:L] ./(NX*0.5*NT)        
    return mean_vals
end

#
function rhomb_half_loop(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX)
    xp2 = mod1(x+2,NX)
    xp3 = mod1(x+3,NX)
    # xp4 = mod1(x+4,NX) 
    tp1 = mod1(t+1,NT) 
    tp2 = mod1(t+2,NT) 
    tp3 = mod1(t+3,NT)
    tp4 = mod1(t+4,NT) 
    # n   n   n   n   n   | a   n   | a   a   a   a   a
    # 1   2   1   2   2   | 1   2   | 1   2   2   2   2  
    # x   xp1 xp1 xp2 xp2 | xp1 xp1 | x   x   x   x   x
    # t   t   tp1 tp1 tp2 | tp3 tp3 | tp4 tp3 tp2 tp1 t
    bla = U[1,x,t] * U[2,xp1,t] * U[1,xp1,tp1] * U[2,xp2,tp1] * U[2,xp2,tp2] 
    bla = bla * adjoint(U[1,xp1,tp3]) * U[2,xp1,tp3]
    return bla * adjoint(U[1,x,tp4]) * adjoint(U[2,x,tp3]) * adjoint(U[2,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
end

#
function edge_loop_hex(U, x, t)    
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX)
    xp2 = mod1(x+2,NX)
    # xp3 = mod1(x+3,NX)
    # xp4 = mod1(x+4,NX) 
    tp1 = mod1(t+1,NT) 
    tp2 = mod1(t+2,NT) 
    tp3 = mod1(t+3,NT)
    tp4 = mod1(t+4,NT) 
    tm  = mod1(t-1,NT)
    # n   a   n   n   n   | a   n   n   n   | a   a   a   a   a
    # 1   2   1   2   2   | 1   2   2   2   | 1   2   2   2   2 
    # x   xp1 xp1 xp2 xp2 | xp1 xp1 xp1 xp1 | x   x   x   x   x 
    # t   tm  tm  tm  t   | tp1 tp1 tp2 tp3 | tp4 tp3 tp2 tp1 t
    bla = U[1,x,t] * adjoint(U[2,xp1,tm]) * U[1,xp1,tm] * U[2,xp2,tm] * U[2,xp2,t]
    bla = bla * adjoint(U[1,xp1,tp1]) * U[2,xp1,tp1] * U[2,xp1,tp2] * U[2,xp1,tp3]
    return bla * adjoint(U[1,x,tp4]) * adjoint(U[2,x,tp3]) * adjoint(U[2,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
end

#
function rhomb_loop(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX)
    xp2 = mod1(x+2,NX)
    xp3 = mod1(x+3,NX)
    # xp4 = mod1(x+4,NX) 
    tp1 = mod1(t+1,NT) 
    tp2 = mod1(t+2,NT) 
    tp3 = mod1(t+3,NT)
    # tp4 = mod1(t+4,NT) 
    tm  = mod1(t-1,NT)
    # n   a   n   n   n   n   n   | a   n   | a   a   a   a   a 
    # 1   2   1   2   1   2   2   | 1   2   | 1   2   1   2   2
    # x   xp1 xp1 xp2 xp2 xp3 xp3 | xp2 xp2 | xp1 xp1 x   x   x 
    # t   tm  tm  tm  t   t   tp1 | tp2 tp2 | tp3 tp2 tp2 tp1 t
    bla = U[1,x,t] * adjoint(U[2,xp1,tm]) * U[1,xp1,tm] * U[2,xp2,tm] * U[1,xp2,t] * U[2,xp3,t] * U[2,xp3,tp1]
    bla = bla * adjoint(U[1,xp2,tp2]) * U[2,xp2,tp2]
    return bla * adjoint(U[1,xp1,tp3]) * adjoint(U[2,xp1,tp2]) * adjoint(U[1,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
end

#
function L_loop_hex(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX)
    xp2 = mod1(x+2,NX)
    # xp3 = mod1(x+3,NX)
    # xp4 = mod1(x+4,NX) 
    tp1 = mod1(t+1,NT) 
    tp2 = mod1(t+2,NT) 
    tp3 = mod1(t+3,NT)
    tp4 = mod1(t+4,NT) 
    tp5 = mod1(t+5,NT) 
    tp6 = mod1(t+6,NT) 
    tm  = mod1(t-1,NT)
    # n   a   n   n   n   | a   n   n   n   n   n   | a   a   a   a   a   a   a
    # 1   2   1   2   2   | 1   2   2   2   2   2   | 1   2   2   2   2   2   2
    # x   xp1 xp1 xp2 xp2 | xp1 xp1 xp1 xp1 xp1 xp1 | x   x   x   x   x   x   x
    # t   tm  tm  tm  t   | tp1 tp1 tp2 tp3 tp4 tp5 | tp6 tp5 tp4 tp3 tp2 tp1 t
    bla = U[1,x,t] * adjoint(U[2,xp1,tm]) * U[1,xp1,tm] * U[2,xp2,tm] * U[2,xp2,t]
    bla = bla * adjoint(U[1,xp1,tp1]) * U[2,xp1,tp1] * U[2,xp1,tp2] * U[2,xp1,tp3] * U[2,xp1,tp4] * U[2,xp1,tp5]
    return bla * adjoint(U[1,x,tp6]) * adjoint(U[2,x,tp5]) * adjoint(U[2,x,tp4]) * adjoint(U[2,x,tp3]) * adjoint(U[2,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
end

#= 
# For debugging purposes
function spatial_ind_loop_mat_hex(l_x)
    L = l_x-1
    big_loop = div(L,2)
    smol_lop = L - 2*div(L,2)
    return big_loop, smol_lop
end
=#

# loop_mat_hex_mike

function top_charge_U2_hex(U)
    NX = size(U,2)
    NT = size(U,3)
    Q = 0.0
    for t = 1:NT
        for x = mod1(t,2):2:NX
            Q += imag(log(det(hexplaq(U, x, t)))) 
        end
    end
    return Q / 2 / π
end

