include("observables_head.jl")

# function mywrite(path, obs)
#     bla = open(path, "a")
#     write(bla, "$obs\n")
#     close(bla)
#     return nothing
# end

# function mywrite(path, obs::Vector)
#     bla = open(path, "a")
#     for i = 1:length(obs)
#         blu = obs[i]
#         write(bla, "$blu\t")
#     end
#     write(bla, "\n")
#     close(bla)
#     return nothing
# end

# function mywrite(path, obs)
#     bla = open(path,"a")
#     writedlm(bla, obs)
#     # writedlm(bla, "\n")
#     close(bla)
#     return nothing
# end

# function mywrite(path, obs::Array)
#     bla = open(path,"a")
#     writedlm(bla, transpose(obs))
#     # writedlm(bla, "\n")
#     close(bla)
#     return nothing
# end

function mywrite(path, obs::Array{coeffs_SU2}, N_x, N_t)
    mat = [get_array(obs[μ,x,t]) for μ = 1:2, x = 1:N_x, t = 1:N_t ]
    bla = open(path,"w")
    writedlm(bla, mat)
    close(bla)
    return nothing
end

function read_last_conf(path, N_x, N_t)
    V1 = readdlm(path)
    V1 = V1[1:size(V1,1)-1, 1:size(V1,2)-1]
    return reshape([coeffs_SU2(V1[i,1], V1[i,2], V1[i,3], V1[i,4]) for i = 1:size(V1,1)], (2,N_x,N_t))
end


function plaq(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    x_p = mod1(x+1, NX) # x%NX + 1
    t_p = mod1(t+1, NT) # t%NT + 1
    # return mult_SU2(U.U[2,t,x], mult_SU2(U.U[1,t,x_p], mult_SU2(adjoint(U.U[2,t_p,x]), adjoint(U.U[1,t,x]))))
    return U[1,x,t] * U[2,x_p,t] * adjoint(U[1,x,t_p]) * adjoint(U[2,x,t])
end

function action(U, β)
    NX = size(U,2)
    NT = size(U,3)
    S = 2*NX*NT   # later generalization: N_colour * NT * (NX)^d_s
    for t = 1:NT
        for x = 1:NX
            S -= real(tr(plaq(U,x,t)))
        end
    end
    return β*S/2    # later generalization: β*S/N_colour
end

# A (2×3)-Wilson loop written by hand for debugging purposes
function loop_2x3_square(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX)
    xp2 = mod1(x+2,NX)
    tp1 = mod1(t+1,NT)
    tp2 = mod1(t+2,NT)
    tp3 = mod1(t+3,NT)
    
    # n   n   n   n   n   |  a   a   a   a   a 
    # 1   1   2   2   2   |  1   1   2   2   2 
    # x   xp1 xp2 xp2 xp2 |  xp1 x   x   x   x 
    # t   t   t   tp1 tp2 |  tp3 tp3 tp2 tp1 t 

    res = U[1,x,t] * U[1,xp1,t] * U[2,xp2,t] * U[2,xp2,tp1] * U[2,xp2,tp2]
    res = res * adjoint(U[1,xp1,tp3]) * adjoint(U[1,x,tp3]) * adjoint(U[2,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
    return res
end

# Compute an (R×T)- Wilson loop at the coordinate (x,t)
function RT_loop(U, R, T, x, t)
    NX = size(U,2)
    NT = size(U,3)
    # loop = coeffs_Id_SU2()
    loop = U[1,x,t]
    x = mod1(x+1,NX)
    for i = 2:R
        loop *= U[1,x,t]
        x = mod1(x+1, NX)
    end
    for i = 1:T
        loop *= U[2,x,t]
        t = mod1(t+1,NT)
    end
    for i = 1:R
        x = mod1(x-1, NX)
        loop *= adjoint(U[1,x,t])
    end
    for i = 1:T
        t = mod1(t-1,NT)
        loop *= adjoint(U[2,x,t])
    end
    return loop
end


#=
# Returns an (Nₓ × Nₜ)-matrix whose entries carry the rectangular (R×T)-loop 
# at the respective (x,t)-points of the lattice. 
function loop_mat(U, R, T)
    NX = size(U,2)
    NT = size(U,3)

    res = [coeffs_Id_SU2() for x = 1:NX, t = 1:NT] # coeffs of the identity in (NT×NX)-matrix
    x_arr = collect(1:NX)
    t_arr = collect(1:NT)
    for i = 1:R
        res = res .* U[1,x_arr,t_arr]
        circshift!(x_arr,-1)    # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    end
    for i = 1:T
        res = res .* U[2,x_arr,t_arr]
        circshift!(t_arr,-1)   
    end
    for i = 1:R
        circshift!(x_arr,1)
        res = res .* adjoint.(U[1,x_arr,t_arr])
    end
    for i = 1:T
        circshift!(t_arr,1)
        res = res .* adjoint.(U[2,x_arr,t_arr])
    end
    return res
end
=#

# Returns an (Nₓ × Nₜ)-matrix whose entries carry the rectangular (R×T)-loop 
# at the respective (x,t)-points of the lattice. 
function loop_mat(U, R, T)
    NX = size(U,2)
    NT = size(U,3)
    return [RT_loop(U,R,T,x,t) for x = 1:NX, t = 1:NT]
end

#
function RT_loop_mike(U, R, T, x, t, avg_U)
    NX = size(U,2)
    NT = size(U,3)

    res = U[1,x,t]
    x = mod1(x+1,NX)
    for i = 2:R-1
        res *= avg_U[1,x,t]
        x = mod1(x+1, NX)
    end
    res *= U[1,x,t]
    x = mod1(x+1,NX)
    
    res *= U[2,x,t]
    t = mod1(t+1,NT)
    for i = 2:T-1
        res *= avg_U[2,x,t]
        t = mod1(t+1,NT)
    end
    res *= U[2,x,t]
    t = mod1(t+1,NT)

    x = mod1(x-1,NX)
    res *= adjoint(U[1,x,t])
    for i = 2:R-1
        x = mod1(x-1, NX)
        res *= adjoint(avg_U[1,x,t])
    end
    x = mod1(x-1,NX)
    res *= adjoint(U[1,x,t])

    t = mod1(t-1,NT)
    res *= adjoint(U[2,x,t])
    for i = 2:T-1
        t = mod1(t-1,NT)
        res *= adjoint(avg_U[2,x,t])
    end
    t = mod1(t-1,NT)
    res *= adjoint(U[2,x,t])

    return res
end

# Returns an (Nₜ × Nₓ)-matrix just like loop_mat(), but uses 
# C.Michael, NPB 259, 58, eq.(6)
# for noise reduction. This can only be done for those links
# of the loop which are not part of the corners. Hence we will
# have to get an if-statement for each Wilson-line (i.e. the 
# four "if T/x > 2"-blocks). If the Wilson line is long enough,
# we can drag the first and last link of the Wilson line out of
# the inner for-loop and replace the links inside with avg_U[...].
# Regarding tests: only possible by comparing simulation results.
function loop_mat_mike(U, R, T, β)
    NX = size(U,2)
    NT = size(U,3)

    if R<=2 && T<=2
        return loop_mat(U,R,T)
    end
    
    stap_field = [staple(U,μ,x,t) for μ = 1:2, x = 1:NX, t = 1:NT]
    # d_field contains the variable d of the paper, evaluated at each coordinate 
    d_field = sqrt.(det.(stap_field))
    avg_U = [(besseli(2,β*d_field[μ,x,t]) / (besseli(1,β*d_field[μ,x,t]) * d_field[μ,x,t])) * stap_field[μ,x,t] for μ = 1:2, x = 1:NX, t = 1:NT]

    return [RT_loop_mike(U,R,T,x,t,avg_U) for x = 1:NX, t = 1:NT]
end


# A function which measures everything one can measure (yet) using Wilson loops.
# The loops are specified in the array 'loops' in which tuples of [n_t, n_x] 
# are specified, i.e. loops = [[1,1], [1,2], [1,4], [3,4], ...]
function measure_loops_corrs(U, loops::Array, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,3)
    L = length(loops)
    t_arr = collect(1:NT)
    # Transpose loop_mat because Julia is columnn-major, see below
    results = Vector{Matrix{Float64}}(undef,L)
    for i = 1:L
        # "results" conatins matrices corresponding to different resp. loops;
        # these matrices contain the trace (2× first entry of quaternion array) of
        # the loop at the resp. space-time point
        mat = loop_mat(stout(U,n_stout,ρ),loops[i][1],loops[i][2])
        results[i] = tr.(mat)
    end
    # Now for each loop we want to obtain a column in "summed":
    # the position in that column is equal to the t-index of the time slice 
    # over which we sum our loop-observable
    summed = Matrix{Float64}(undef,NT,L)
    for i = 1:L
        summed[:,i] = [sum(results[i][:,t]) for t = 1:NT] 
    end

    # A vector containing correlation matrices
    corrs = Array{Float64}(undef,L,L,NT)
    for τ = 1:NT
        corrs[:,:,τ] = [sum(summed[:,i] .* summed[circshift(t_arr,-τ),j]) for i = 1:L, j = 1:L] # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    end
    corrs ./= NX^2*NT

    mean_vals = Vector{Float64}(undef,L)
    for i = 1:L
        mean_vals[i] = sum(summed[:,i])
    end
    mean_vals ./= (NX*NT)

    return corrs, mean_vals
end
# So the result of measure_loops is an array containing corrs and mean_vals:
#   corrs[:,:,t] is the correlation matrix at phys. time t
#   mean_vals[i] is the mean value of the i-th loop (average of the lattice)


# A function just to get the mean values 'cause it's faster
function measure_loops(U, loops::Array, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,3)
    L = length(loops)
    results = [tr.(loop_mat(stout(U,n_stout,ρ), loop[1], loop[2])) for loop in loops]
    mean_vals = [sum(results[i]) for i = 1:L] ./(NX*NT)
    return mean_vals
end

# Same as measure_loops(), but using loop_mat_mike() instead of loop_mat()
function measure_loops_mike(U, loops::Array, β)
    NT = size(U,3)
    NX = size(U,2)
    L = length(loops)
    results = [tr.(loop_mat_mike(U, loops[i][1], loops[i][2], β)) for i = 1:L]
    mean_vals = [sum(results[i]) for i = 1:L] ./(NX*NT)
    return mean_vals
end

# By edge we mean Wilson loops of the shape below (uncomment the shape for more clearness)
#   _
#  | |_
#  |_._|
function edge_loop_1(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    # n   n   n   a   n   a   a   a  
    # 1   1   2   1   2   1   2   2
    # x   xp  xpp xp  xp  x   x   x
    # t   t   t   tp  tp  tpp tp  t
    xp  = mod1(x+1, NX) # x%NX + 1
    tp  = mod1(t+1, NT) # t%NT + 1
    xpp = mod1(x+2, NX) 
    tpp = mod1(t+2, NT) 
    return U[1,x,t] * U[1,xp,t] * U[2,xpp,t] * adjoint(U[1,xp,tp]) * U[2,xp,tp] * adjoint(U[1,x,tpp]) * adjoint(U[2,x,tp]) * adjoint(U[2,x,t])
end

# test_field = gaugefield_SU2(32,32,true)
# X = rand(1:32)
# T = rand(1:32)
# edge_loop_1(test_field,X,T)

# By "L" we mean Wilson loops of the shape below (uncomment the shape for more clearness)
#   _
#  | |
#  | |_
#  |_._|
function L_loop_1(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    # n   n   n   a   n   n   |  a   a   a   a
    # 1   1   2   1   2   2   |  1   2   2   2
    # x   xp  xpp xp  xp  xp  |  x   x   x   x 
    # t   t   t   tp  tp  tpp |  tp3 tpp tp  t
    xp  = mod1(x+1, NX) # x%NX + 1
    xpp = mod1(x+2, NX) 
    tp  = mod1(t+1, NT) # t%NT + 1
    tpp = mod1(t+2, NT) 
    tp3 = mod1(t+3, NT) 
    return U[1,x,t] * U[1,xp,t] * U[2,xpp,t] * adjoint(U[1,xp,tp]) * U[2,xp,tp] * U[2,xp,tpp] * adjoint(U[1,x,tp3]) * adjoint(U[2,x,tpp]) * adjoint(U[2,x,tp]) * adjoint(U[2,x,t])
end

#
function rhomb_half_loop_square(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp  = mod1(x+1, NX) # x%NX + 1
    tp  = mod1(t+1, NT) # t%NT + 1
    # xpp = mod1(x+2, NX) 
    tpp = mod1(t+2, NT) 
    a = U[1,x,t]*U[2,xp,t] + U[2,x,t]*U[1,x,tp]
    b = adjoint(U[1,x,tp])*U[2,x,tp] + U[2,xp,tp]*adjoint(U[1,x,tpp])
    # c = a*b/det(a*b)
    # return c*adjoint(U[2,x,tp])*adjoint(U[2,x,t])
    return proj_SU2(a*b*adjoint(U[2,x,tp])*adjoint(U[2,x,t]))
end

#
function rhomb_loop_square(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp  = mod1(x+1, NX) # x%NX + 1
    tp  = mod1(t+1, NT) # t%NT + 1
    xpp = mod1(x+2, NX) 
    tpp = mod1(t+2, NT) 
    a = adjoint(U[2,x,t])*U[1,x,t] + U[1,x,tp]*adjoint(U[2,xp,t])
    b = U[1,xp,t]*U[2,xpp,t] + U[2,xp,t]*U[1,xp,tp]
    c = adjoint(U[1,xp,tp])*U[2,xp,tp] + U[2,xpp,tp]*adjoint(U[1,xp,tpp])
    d = adjoint(U[2,xp,tp])*adjoint(U[1,x,tp]) + adjoint(U[1,x,tpp])*adjoint(U[2,x,tp])
    # return a*b*c*d/det(a*b*c*d)
    return proj_SU2(a*b*c*d)
end



# clover action
# topological charge (even though non-existent in 2D SU(2))
# Polyakov loop???





########    U(2) Shenanigans    ########





function top_charge_U2(U)
    NX = size(U,2)
    NT = size(U,3)
    return imag(sum([log(det(plaq(U, x, t))) for x = 1:NX, t = 1:NT])) / 2 / π
end
