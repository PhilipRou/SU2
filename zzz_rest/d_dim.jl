
####    D-dimensional shenanigans   #### 🚧👷 under construction 👷🚧 





# Given a coordinate coord = [t,x_1, x_2, ...] we want to find out the one-dimensional
# index of the lattice:
#       ind = 1 + (t-1) + (x_1 - 1)⋅ N_t + (x_2 - 1)⋅ N_t⋅ N_x + ... (x_i - 1)⋅ N_t⋅ N_x^(i-1)
function ind_D_2_one(coord::Vector, N_t, N_x, D)
    ind = 1
    ind += coord[1] - 1
    ind += (coord[2] - 1) * N_t     #❗❗❗ Here we see that we MUST HAVE D ≥ 2 ❗❗❗
    for i = 3:D
        ind += (coord[i]-1) * N_t * N_x^(i-2)
    end
    return ind
end

# Given a one-dimensional index we want to find the d-dimensional coordinate again
function ind_one_2_D(ind::Int, N_t, N_x, D)
    num = deepcopy(ind)
    coord = Vector{Int}(undef, D)
    for i = 0:D-2               #❗❗❗ Here we also see that we MUST HAVE D ≥ 2 ❗❗❗
        pos = D-i
        coord[pos] = div(num-1, N_t * N_x^(pos-2) ) +1
        num = (num-1) % (N_t * N_x^(pos-2)) +1
    end
    coord[1] = num
    return coord
end


# Construct a D-dimensional gaugefield with SU(2)-valued links
function gaugefield_SU2_D(N_t::Int, N_x::Int, D::Int, hot::Bool)
    N = N_x .* ones(Int, D)
    N[1] = N_t
    U = Array{coeffs_SU2}(undef, D, N...)
    indices = [ind_one_2_D(i,N_t,N_x,D) for i = 1:N_t*N_x^(D-1)]
    if hot
        for ind in indices
            for μ = 1:D
                U[μ,ind...] = ran_SU2(rand())
            end
        end
    else
        for ind in indices
            for μ = 1:D
                U[μ,ind...] = coeffs_SU2(1.0,0.0,0.0,0.0)
            end
        end
    end
    return U
end

# function to get a list of all D-dim. coordinates that a lattice
# of size N_t × N_x^(D-1) has to offer. The first index is incremented
# the fastest.
function D_dim_coords(N_t, N_x, D)
    return [ind_one_2_D(i, N_t, N_x, D) for i = 1:N_t*N_x^(D-1)]
end





####    Updates    ####






# staple_dag but for D-dimensional configs, D arbitrary
function staple_dag_D(U, μ, coord::Vector)
    D  = size(U,1)  # number of space-time dimensions
    NT = size(U,2)
    NX = size(U,3)
    # N: vector containing number of sites in resp. directeion, i.e. N = [NT, NX, NX, NX, ...]
    N    = NX .* ones(Int,D)
    N[1] = NT
    # news: all the directional indices ν (nu/new) to be summed over later
    news = collect(1:D)
    popat!(news, μ)
    # vecs: we want vecs[:,ρ] to be the unit vector in ρ-direction
    vecs = Matrix{Int}(I,D,D)
    # wings: will store the (D-1) summands of the staple, i.e. the P_i in [Gattr./Lang Eq.(4.20)]
    wings = Array{coeffs_SU2}(undef, D-1)

    for i = 1:D-1
        arg_1 = mod1.(coord+vecs[:,μ], N[μ])
        arg_2 = mod1.(coord+vecs[:,news[i]], N[news[i]])
        arg_3 = coord
        arg_4 = mod1.(mod1.(coord+vecs[:,μ], N[μ]) - vecs[:,news[i]], N[news[i]])
        arg_5 = mod1.(coord-vecs[:,news[i]], N[news[i]])
        arg_6 = mod1.(coord-vecs[:,news[i]], N[news[i]])

        wings[i] = U[news[i], arg_1...] * adjoint(U[μ, arg_2...]) * adjoint(U[news[i], arg_3...])
        wings[i] += adjoint(U[news[i], arg_4...]) * adjoint(U[μ, arg_5...]) * U[news[i], arg_6...]
    end
    
    return sum(wings)
end

# metro! but for D-dim. configs
function metro_D!(U, μ, coord::Vector, step)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,coord...]
    staple_d = staple_dag_D(U,μ,coord)
    S_old = β*0.5*tr(U[μ,coord...] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,coord...] = new_coeffs
        acc[1] += 1
    end
    return nothing
end


# Returns an array containing all 1-dim. indices of an (N_T × N_X^{D-1}) - lattice,
# but in a chess-pattern. Not intended for regular use, just once per simulation
# due to inefficiency (see e.g. comment on chess_metro_D! )
function indices_chess(N_t, N_x, D)
    indices = []
    for trip = 1:2
        # i runs over the number of (consecutive) blocks of length N_t, of which there are N_x^(D-1)-many
        for i = 1:N_x^(D-1) 
            # "pos" determines whether vec (below) starts on the first or second index of the i-th block
            # and has to respect all dimensions (hence the j-part)
            pos = i-1 + sum([div(i-1, N_x^j) for j = 1:D])
            # vec counts every other index in the i-th block of length N_t
            vec = Vector(1 + mod(trip+pos,2) + (i-1)*N_t : 2 : i*N_t)
            push!(indices, vec...)
        end
    end
    return indices
end


# test_list = indices_chess(4,4,3)

# for i = 1:64
#     @assert i in test_list
# end

function D_dim_coords_chess(N_t, N_x, D)
    indices = indices_chess(N_t, N_x, D)
    return ind_one_2_D.(indices)
end

# chess_metro!, but in D dimensions
# ❗ Requires an array "coords" in advance, obtained via D_dim_coords_chess ❗
# (Because initiating this list everytime would be too costly, just do it 
# once per simulation is doable, albeit a little ugly)
function chess_metro_D!(U, step, coords)
    D = size(U,1)
    for μ = 1:D
        for coord in coords
            metro_D!(U,μ,coord,step)
        end
    end
    return nothing
end


# # overrelax!, but in D dimensions (D arbitrary)
# function overrelax_D!(U, μ, coord::Vector)
#     v = proj2man!(staple_dag_D(U,μ,coord))
#     U[μ,coord...] = adjoint(v *  U[μ,coord...] * v)
#     return nothing
# end

# overrelax!, but in D dimensions (D arbitrary)
function overrelax_D!(U, μ, coord::Vector)
    v = proj2man!(staple_dag_D(U,μ,coord))
    U[μ,coord...] = adjoint(v *  U[μ,coord...] * v)
    return nothing
end


# 
function chess_overrelax_3d!(U)
    NT = size(U,2)
    NX = size(U,3)

    # μ = 1   # t-direction
    # for trip = 1:2
    #     for y = 1:N_x
    #         for x = (1+mod(y+trip,2)):2:N_x
    #             for t = 1:N_t
    #                 metro_D!(U, μ, [t,x,y], step)
    #             end
    #         end
    #     end
    # end

    # for μ = 2:3   # x-direction
    #     for trip = 1:2
    #         for y = 1:N_x
    #             for x = 1:N_x
    #                 for t = (1+mod(y+trip,2)):2:N_t
    #                     metro_D!(U, μ, [t,x,y], step)
    #                 end
    #             end
    #         end
    #     end
    # end

    for μ = 1:3
        for trip = 1:2
            for y = 1:NX
                for x = 1:NX
                    for t = (1+mod(x+y+trip,2):2:NT)
                        overrelax_D!(U, μ, [t,x,y])
                    end
                end
            end
        end
    end
    return nothing
end


# chess_overrelax!, but in D dimensions
# ❗ Requires an array "coords" in advance, obtained via D_dim_coords_chess ❗
# (Because initiating this list everytime would be too costly, just do it 
# once per simulation is doable, albeit a little ugly)
function chess_overrelax_D!(U, coords)
    D  = size(U,1)
    for μ = 1:D
        for coord in coords
            overrelax_D!(U,μ,coord)
        end
    end
    acc[1] += D*NT*NX^(D-1)
    return nothing
end



#=
using BenchmarkTools

N_t = N_x = 16
μ = rand(1:3)
t = rand(1:N_t)
x = rand(1:N_x)
y = rand(1:N_x)
test_field = gaugefield_SU2_3d(N_t, N_x, true)

@benchmark staple_dag_D(test_field, μ, [t,x,y])     # (9 ± 21) μs
@benchmark staple_dag_3d(test_field, μ, t, x, y)    # (125 ± 132) ns
=#
