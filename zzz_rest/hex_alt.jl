include("C:\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")

# We want to get only one half of the chequer board, that being 
# the vertices on which two non-trivial links are defined
# (i.e. the "starting points" of the plaquettes)
function half_chess_coords(N_t, N_x)
    coords = []
    trip = 1  # Only want half the indices
    for x = 1:N_x
        for t = (1+mod(x+trip,2)):2:N_t
            push!(coords, [t,x])
        end
    end
    return coords
end

# function other_half_chess_coords(N_t, N_x)
#     indices = []
#     trip = 2  # Now we want the other half
#     for x = 1:N_x
#         for t = (1+mod(x+trip,2)):2:N_t
#             push!(indices, [t,x])
#         end
#     end
#     return indices
# end

function hex_links_coords_chess(N_t, N_x)
    coords = []
    trip = 1  # Only want half the indices
    for x = 1:N_x
        for t = (1+mod(x+trip,2)):2:N_t
            push!(coords, [1,t,x], [2,t,x])
            push!(coords, [1,mod1(t+1, N_t),x])
        end
    end
    return coords
end

# hex_links_coords_chess(10,10)


# test_inds = half_chess_indices(10,10)
# for bla in test_inds
#     println(bla[1], bla[2])
# end

function hexfield_SU2_alt(N_t::Int64, N_x::Int64, hot::Bool)
    # @assert iseven(N_t) "Input N_t must be even (recall PBC for hex. lattices)"
    # @assert iseven(N_x) "Input N_x must be even (recall PBC for hex. lattices)"
    # U = Array{coeffs_SU2}(undef, 2, N_t, N_x)
    # U = [NaN for μ = 1:2, t = 1:N_t, x = 1:N_x]
    U = [coeffs_SU2(NaN, NaN, NaN, NaN) for μ = 1:2, t = 1:N_t, x = 1:N_x]
    # coords = half_chess_coords(N_t, N_x)
    # other_coords = other_half_chess_coords(N_t, N_x)
    coords = hex_links_coords_chess(N_t, N_x)
    if hot
        for μ = 1:2
            for coord in coords
                U[μ,coord[1], coord[2]] = ran_SU2(rand())
            end
        end
    else 
        for μ = 1:2
            for coord in coords
                U[μ,coord[1], coord[2]] = coeffs_SU2(1.0,0.0,0.0,0.0) # coeffs of 1∈ SU(2)
            end
        end
    end
    return U
end


function staple_dag_hex_alt(U, μ, t, x)
    # NX = N_x>>1
    NX = size(U,2)
    NT = size(U,3)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    tp = mod1(t+1, NT) # t%NT +1                 
    xp = mod1(x+1, NX) # x%NX +1                 
    tm = mod1(t-1, NT) # (t + NT -2)%NT +1   
    xm = mod1(x-1, NX) # (x + NX -2)%NX +1  
    tpp = mod1(t+2, NT) # (t+1)%NT +1
    tmm = mod1(t-2, NT) # (t + NT - 3)%NT +1

    # coords = half_chess_coords(N_t, N_x)

    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    if μ == 1
        if iseven(t+x)
            # a  a  a  n  n 
            # 2  1  1  2  1
            # tp t  tm tm tm
            # xm xm xm xm x
            a = adjoint(U[2,tp,xm]) * adjoint(U[1,t,xm])* adjoint(U[1,tm,xm]) * U[2,tm,xm] * U[1,tm,x]
            # n  n   a  a  a
            # 1  2   1  1  2
            # tp tpp tp t  t
            # x  x   xp xp x
            b = U[1,tp,x] * U[2,tpp,x] * adjoint(U[1,tp,xp]) * adjoint(U[1,t,xp]) * adjoint(U[2,t,x])
        else
            # n  a   a  a  n
            # 1  2   1  1  2
            # tp tpp tp t  t
            # x  xm  xm xm xm 
            a = U[1,tp,x] * adjoint(U[2,tpp,xm]) * adjoint(U[1,tp,xm]) * adjoint(U[1,t,xm]) * U[2,t,xm]
            # n  a  a  a  n
            # 2  1  1  2  1
            # tp t  tm tm tm 
            # x  xp xp x  x
            b = U[2,tp,x] * adjoint(U[1,t,xp]) * adjoint(U[1,tm,xp]) * adjoint(U[2,tm,x]) * U[1,tm,x]
        end
    else #if μ == 2
        # n  n  a   a  a
        # 1  1  2   1  1
        # t  tp tpp tp t
        # xp xp x   x  x
        a = U[1,t,xp] * U[1,tp,xp] * adjoint(U[2,tpp,x]) * adjoint(U[1,tp,x]) * adjoint(U[1,t,x])
        # a  a   a   n   n 
        # 1  1   2   1   1 
        # tm tmm tmm tmm tm 
        # xp xp  x   x   x 
        b = adjoint(U[1,tm,xp]) * adjoint(U[1,tmm,xp])* adjoint(U[2,tmm,x]) * U[1,tmm,x] * U[1,tm,x]
    end
    return a+b
end

# ❗ Just for testing purposes, not to be actually used
function delta_S_gauge_hex_alt(U, μ, t, x, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2)
    return β/2*tr((old_coeffs - new_coeffs) * staple_dag_hex_alt(U,μ,t,x))
end

function action_hex_alt(U,β)
    NX = size(U,2)
    NT = size(U,3)
    S = 2*NT*NX
    for coord in half_chess_coords(NT, NX)
        t = coord[1]
        x = coord[2]
        tp = mod1(t+1, NT) # t%NT +1                 
        xp = mod1(x+1, NX) # x%NX +1                 
        # tm = mod1(t-1, NT) # (t + NT -2)%NT +1   
        # xm = mod1(x-1, NX) # (x + NX -2)%NX +1  
        tpp = mod1(t+2, NT) # (t+1)%NT +1
        # tmm = mod1(t-2, NT) # (t + NT - 3)%NT +1
        S -= tr(U[2,t,x]*U[1,t,xp]*U[1,tp,xp]*adjoint(U[2,tpp,x])*adjoint(U[1,tp,x])*adjoint(U[1,t,x]))
    end
    S *= 0.5*β
end


# β = 1.0
# N_t = N_x = 10
# test_hex_alt = hexfield_SU2(N_t, N_x, true)
# link_coords = hex_links_coords_chess(N_t, N_x)
# coord = rand(link_coords)
# # old_coord = deepcopy(coord)
# # coord = [1,2,3]
# test_hex_neu = deepcopy(test_hex_alt)
# test_hex_neu[coord...] = ran_SU2(rand())
# delta_S_gauge_hex_alt(test_hex_alt,coord...,test_hex_alt[coord...],test_hex_neu[coord...])
# action_hex_alt(test_hex_alt) - action_hex_alt(test_hex_neu)
# # test_hex_alt[coord...]
# # test_hex_neu[coord...]

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

function metro_hex!(U, μ, t, x, step)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,t,x]
    staple_d = staple_dag_hex(U,μ,t,x)
    S_old = β*0.5*tr(U[μ,t,x] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function metro_hex_alt!(U, μ, t, x, step)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,t,x]
    staple_d = staple_dag_hex_alt(U,μ,t,x)
    S_old = β*0.5*tr(U[μ,t,x] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x] = new_coeffs
        acc[1] += 1
    end
    return nothing
end


# metro_hex_alt!(test_hex_alt,coord...,0.1)

function lexico_metro_hex_alt!(U,step)
    NX = size(U,2)
    NT = size(U,3)
    for μ = 1:2
        for t = 1:NT
            for x = 1:NX
                metro_hex_alt!(U, μ, t, x, step)
            end
        end
    end
    return nothing
end

function lexico_metro_hex!(U,step)
    NT = size(U,2)
    NX = size(U,2)
    for μ = 1:3
        for t = 1:NT
            for x = 1:NX
                metro_hex!(U, μ, t, x, step)
            end
        end
    end
    return nothing
end



#=
using BenchmarkTools

N_t = 32
N_x = 32
β   = 1.0
ϵ   = 0.2
test_field_hex = hexfield_SU2(N_t>>1, N_x, true)
test_field_hex_alt = hexfield_SU2_alt(N_t, N_x, true)
acc = [0]

# coords = hex_links_coords_chess(N_t,N_x)

@benchmark lexico_metro_hex!(test_field_hex, ϵ) # (700 ± 300) μs
@benchmark lexico_metro_hex_alt!(test_field_hex_alt, ϵ) # (950 ± 350) μs

coord = rand(hex_links_coords_chess(N_t>>1,N_x))
@benchmark metro_hex!(test_field_hex, coord[1],coord[2],coord[3], ϵ)  # (0.5±0.5)μs
@benchmark metro_hex_alt!(test_field_hex_alt, coord[1],coord[2],coord[3], ϵ) # (0.5±0.5)μs

# size(test_field_hex)
# size(test_field_hex_alt)



# Bottom line: the old way works just fine!
# Edit: For N_t = N_x = 64 the alternative method actually has less allocs, approx. same runtime
=#


