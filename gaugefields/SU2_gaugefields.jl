################################################################################
using LinearAlgebra


# Use of symbols:
#   ❗ Attention, keep this comment in mind w.r.t. the object defined below
#   ⭕ Note to self: improve at next convenience
#   Any other emoji: desperate attempt to be funny


# Pauli matrices:
σ0 = [1 0; 0 1]     # not really a Pauli matrix 🙄
σ1 = [0 1; 1 0]
σ2 = [0 -im; im 0]
σ3 = [1 0; 0 -1]
Σ  = [σ0, σ1, σ2, σ3]
Σ_im = [σ0, im*σ1, im*σ2, im*σ3]

# Construct a Lie group element out of a linear combination of Lie algebra elements
function alg2grp_SU2(coeffs::Array)
    @assert length(coeffs) == 3 "lie2mat needs an array of three coeffs for σ₁ to σ₃"
    ϕ = zeros(ComplexF64, 2, 2)
    for i = 1:3
        ϕ += coeffs[i] * Σ[i+1]
    end
    return exp(im*ϕ)
end

# One can represent any SU(2)-matrix M with the help of four coeffs a,b,c,d
# and the Pauli matrices, provided that a²+b²+c²+d² = 1. Then M is given as
#               M = aσ₀ + i(bσ₁ + cσ₂ + dσ₃).
mutable struct coeffs_SU2{T <: Real}
    a::T
    b::T
    c::T
    d::T
    function coeffs_SU2(
        a::T,
        b::T,
        c::T,
        d::T
        ) where {T <: Real}
        return new{T}(a,b,c,d)
    end
end

# SU(2) coeffs have to square up to one! But sometimes we need that, e.g. staple:
function Base.:+(X::coeffs_SU2, Y::coeffs_SU2)
    return coeffs_SU2(X.a+Y.a, X.b+Y.b, X.c+Y.c, X.d+Y.d)
end

function Base.:-(X::coeffs_SU2, Y::coeffs_SU2)
    return coeffs_SU2(X.a-Y.a, X.b-Y.b, X.c-Y.c, X.d-Y.d)
end

# Multiply two SU(2)-matrices whose coefficients are given in X and Y
function Base.:*(X::coeffs_SU2, Y::coeffs_SU2)
    a = X.a*Y.a - X.b*Y.b - X.c*Y.c - X.d*Y.d 
    b = X.a*Y.b + X.b*Y.a - X.c*Y.d + X.d*Y.c
    c = X.a*Y.c + X.b*Y.d + X.c*Y.a - X.d*Y.b
    d = X.a*Y.d - X.b*Y.c + X.c*Y.b + X.d*Y.a
    return coeffs_SU2(a,b,c,d)
end

function Base.:*(α, X::coeffs_SU2)
    return coeffs_SU2(α*X.a, α*X.b, α*X.c, α*X.d)
end

function Base.:/(X::coeffs_SU2, α)
    return coeffs_SU2(X.a/α, X.b/α, X.c/α, X.d/α)
end

function Base.:(==)(X::coeffs_SU2, Y::coeffs_SU2)
    return [X.a,X.b,X.c,X.d] == [Y.a,Y.b,Y.c,Y.d]
end

function Base.isapprox(X::coeffs_SU2, Y::coeffs_SU2)
    # if isapprox(X.a, Y.a)
    #     if isapprox(X.b, Y.b)
    #         if isapprox(X.c, Y.c)
    #             if isapprox(X.d, Y.d)
    #                 return true
    #             end
    #         end
    #     end
    # else 
    #     return false
    # end
    return isapprox([X.a,X.b,X.c,X.d], [Y.a,Y.b,Y.c,Y.d])
end

function LinearAlgebra.tr(X::coeffs_SU2)
    return 2*X.a
end

function LinearAlgebra.det(X::coeffs_SU2)
    return X.a^2 + X.b^2 + X.c^2 + X.d^2 
end

# Take the adjoint of coeffs_SU2, i.e. the adjoint of the matrix 
# corresponding to said coeffs_SU2
function LinearAlgebra.adjoint(X::coeffs_SU2)
    return coeffs_SU2(X.a, -X.b, -X.c, -X.d)
end

# Take the matrix logarithm of coeffs_SU2
# (note: staying in matrix-form or even returning a matrix
# is significantly less efficient)
function log_SU2(X::coeffs_SU2)
    return acos(X.a)/sqrt(1-X.a^2) * coeffs_SU2(0.0, X.b, X.c, X.d)
end

# Map an su(2)-element (from the Lie-algebra) onto the
# manifold SU(2) via the exponential function
function exp_su2(X::coeffs_SU2)
    absX = sqrt(X.b^2+X.c^2+X.d^2)
    A = cos(absX)
    B = sin(absX) / absX
    return coeffs_SU2(A, B*X.b, B*X.c, B*X.d)
end

# ❗ Very inefficient, only for debugging purposes ❗
function get_array(X::coeffs_SU2)
    return [X.a, X.b, X.c, X.d]
end

# Generate coeffs_SU2(x₀,x₁,x₂,x₃) such that 
# x₀σ₀ + i ∑ₖ xₖσₖ ∈ SU(2) and close to the identity:
# the smaller ϵ, the closer, and ϵ = 0 gives the identity
function ran_SU2(ϵ)
    # r1 = 2 * (rand()-0.5)
    # r2 = 2 * (rand()-0.5)
    # r3 = 2 * (rand()-0.5)
    r1, r2, r3 = 2 .* (rand(3) .- 0.5)
    absr = sqrt(r1^2 + r2^2 + r3^2)
    return coeffs_SU2(sqrt(1-ϵ^2), ϵ*r1/absr, ϵ*r2/absr, ϵ*r3/absr)
end

# Create a random element of su(2), i.e. the Lie algebra (mind the
# lower case letters!)
function ran_su2(ϵ)
    r = ϵ .* (2 .* rand(3) .- 1)
    return coeffs_SU2(0.0, r[1], r[2], r[3])
end

# Quickly get the coefficients of the identity element
function coeffs_Id_SU2()
    return coeffs_SU2(1.0, 0.0, 0.0, 0.0)
end

# Project coeffs_SU2 onto SU2 (since addition is allowed it may happen that
# some coeffs_SU2 do not describe an SU2 element anymore)
function proj2man(X::coeffs_SU2)
    return X/sqrt(det(X))
end

# Given coeffs_SU2 create the corresponding SU2-matrix
function coeffs2grp(X::coeffs_SU2)
    # @assert isapprox(sum(coeffs.^2), 1.0) "coeffs2grp needs coeffs which square up to 1.0"
    # return sum(coeffs .* Σ_im)
    return X.a*σ0 + im*X.b*σ1 + im*X.c*σ2 + im*X.d*σ3
end

# Given an SU2-matrix return the corresponding coeffs_SU2
function grp2coeffs(mat)
    return coeffs_SU2(real(mat[1,1]), imag(mat[1,2]), real(mat[1,2]), imag(mat[1,1]))
end

# # Take the adjoint of an SU(2) matrix, where the input is an array of
# # four coefficients needed in the quaternionic representation
# function adj_SU2(X::coeffs_SU2)
#     return coeffs_SU2(X.a, -X.b, -X.c, -X.d)
# end

# mutable struct gaugefield_SU2
#     U::Array{Vector{Float64}, 3}
#     N_t::Int64
#     N_x::Int64
#     V::Int64
#     hot::Bool
#     acc_count::Int64

#     function gaugefield_SU2(
#         N_t::Int64,
#         N_x::Int64,
#         hot::Bool
#         )
#         acc_count = 0
#         V = N_t*N_x
#         # U = zeros(ComplexF64, 2, 2, 2, N_t, N_x)
#         U = Array{Vector{Float64}, 3}(undef, 2, N_t, N_x)
#         if hot
#             for μ = 1:2
#                 for t = 1:N_t
#                     for x = 1:N_x
#                         U[μ,t,x] = ran_SU2(rand())
#                     end
#                 end
#             end
#         else 
#             for μ = 1:2
#                 for t = 1:N_t
#                     for x = 1:N_x
#                         U[μ,t,x] = [1,0,0,0]
#                     end
#                 end
#             end       
#         end
#         return new(U, N_t, N_x, V, hot, acc_count)
#     end
# end                                                                                                   


# Construct a square gauge field: an (2 × N_x × N_t)-Array with 
# entries which are SU(2)-valued, or rather in our case, coeffs_SU2-valued 
# (see struct "coeffs_SU2" above). This means that in order to access the
# link in μ-direction at space-time point n = (x,t), we need U[μ,x,t], where
# μ = 1 corresponds to the x-direction ("sideways") and
# μ = 2 corresponds to the t-direction ("upwards").
function gaugefield_SU2(N_x::Int64, N_t::Int64, hot::Bool)
    U = Array{coeffs_SU2}(undef, 2, N_x, N_t)
    if hot
        for t = 1:N_t
            for x = 1:N_x
                for μ = 1:2
                    U[μ,x,t] = ran_SU2(rand())
                end
            end
        end
    else 
        for t = 1:N_t
            for x = 1:N_x
                for μ = 1:2
                    U[μ,x,t] = coeffs_Id_SU2() # coeffs of identity ∈ SU(2)
                end
            end
        end       
    end
    return U
end

#=
# Apply comb gauge onto a square config, except it doesn't work for some reason...
function comb_gauge(U)
    NX = size(U,2)
    NT = size(U,3)
    V = gaugefield_SU2(NX, NT, false)
    V[1,NX,1] = reduce(*,U[1,:,1])
    Ω_slice = U[2,:,1]
    Ω_slice[2:NX] = [reduce(*, U[1,1:x,1]) for x = 1:NX-1] .* Ω_slice[2:NX]
    for t = 2:NT-1
        V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice,-1))
        Ω_slice = Ω_slice .* U[2,:,t]
    end
    V[2,:,NT] = Ω_slice
    for t = NT:NT
        V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice,-1))
        Ω_slice = Ω_slice .* U[2,:,t]
    end
    return V
end
=#



####    Hexagonal shenanigans   #### 





#=
# Construct a hexagonal gauge field with SU(2)-valued links.
# ❗❗❗ Here N_t and N_x denote the number of hexagonal cells in resp. direction ❗❗❗
function hexfield_SU2(N_t::Int64, N_x::Int64, hot::Bool)
    # The order of the μ index is:
    #   μ   vec dir
    #__________________
    #   1   ̂μ   ➡
    #   2   ̂ν   ↗
    #   3   ̃λ   ↘
    @assert iseven(N_t) "Input N_t must be even (recall PBC for hex. lattices)"
    @assert iseven(N_x) "Input N_x must be even (recall PBC for hex. lattices)"
    U = Array{coeffs_SU2}(undef, 3, N_t, N_x)
    if hot
        for μ = 1:3
            for t = 1:N_t
                for x = 1:N_x
                    U[μ,t,x] = ran_SU2(rand())
                end
            end
        end
    else 
        for μ = 1:3
            for t = 1:N_t
                for x = 1:N_x
                    U[μ,t,x] = coeffs_SU2(1.0,0.0,0.0,0.0) # coeffs of 1∈ SU(2)
                end
            end
        end       
    end
    return U
end
=#

# Produce an array containing the indices to be used on a hexagonal config
# (see hexfield_SU2() below). Each element is of the form [μ,x,t] with
# μ ∈ {1,2},  x ∈ {1,...,N_x}  and  t ∈ {1,...,N_t}. 
# The way to store a hexagonal config is the following: we leave out every 
# other link in 2-direction, creating a "shifted brick-structure". The 
# remaining links to be used are the ones created by the function below.
function hex_links_coords_chess(N_x, N_t)
    coords = []
    μ = 1
    for t = 1:N_t
        for x = 2-mod(t,2):2:N_x
            push!(coords, [μ,x,t])
            # println(μ, ", ", x, ", ", t)
        end
    end
    μ = 2
    for trip = 1:2
        for start_t = 1:2
            for t = start_t:2:N_t
                for x = (1+mod(t+trip,2)):2:N_x
                    push!(coords, [μ,x,t])
                    # println(μ, ", ", x, ", ", t)
                end
            end
        end
    end
    return coords
end

# Produce an array of all legal [μ,x,t] indices in lexicographical order
function hex_link_coords_lex(N_x, N_t)
    coords = []
    for t = 1:N_t
        for x = 1:N_x
            if x in 1+mod(t+1,2):2:N_x
                push!(coords, [1,x,t], [2,x,t])
            elseif x in 1+mod(t,2):2:N_x
                push!(coords, [1,x,t])
            end
        end
    end
    return coords
end

# @benchmark hex_links_coords_chess(10,10)     # (7±12) μs
# @benchmark chess_hex_link_coords_lex(10,10) # (9±17) µs

# Produce an array of coords [x,t] on which there are two non-trivial 
# gauge links, i.e. where also U[1,x,t] is defined
function half_chess_coords(N_x, N_t)
    coords = []
    # trip = 1  # Only want half the indices
    for t = 1:N_t
        for x = mod1(t,2):2:N_x
            push!(coords, [x,t])
        end
    end
    return coords
end

function hex_link_coords(N_x, N_t)
    coords = []
    μ = 1
    for t = 1:N_t
        for x = mod1(t,2):2:N_x
            push!(coords, [μ,x,t])
        end
    end
    μ = 2
    for t = 1:N_t
        for x = 1:N_x
            push!(coords, [μ,x,t])
        end
    end
    return coords
end

# Construct a hexagonal gauge field. See comment on hex_links_coords_chess() 
# above. Every taboo link is constructed as (coeffs_SU2 of) a NaN matrix.
function hexfield_SU2(N_x::Int64, N_t::Int64, hot::Bool)
    @assert iseven(N_x)
    @assert iseven(N_t)
    U = [coeffs_SU2(NaN,NaN,NaN,NaN) for μ = 1:2, x = 1:N_x, t = 1:N_t]
    # coords = hex_links_coords_chess(N_x, N_t)
    coords = hex_link_coords(N_x, N_t)
    if hot
        for coord in coords
            U[coord[1], coord[2], coord[3]] = ran_SU2(rand())
        end
    else 
        for coord in coords
            U[coord[1], coord[2], coord[3]] = coeffs_Id_SU2()
        end
    end
    return U
end

# Apply temporal gauge onto a square config
# ❗❗❗ Works but is significantly less efficient than temp_gauge(),
# which also works on hexagonal configs
function temp_gauge_hex(U)
    NX = size(U,2)
    NT = size(U,3)
    V = hexfield_SU2(NX, NT, false)
    Ω_slice = [coeffs_Id_SU2() for x = 1:NX] 
    for t = 1:NT
        # V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice,-1))
        V[1,1+mod(t+1,2):2:NX,t] = Ω_slice[1+mod(t+1,2):2:NX] .* U[1,1+mod(t+1,2):2:NX,t] .* adjoint.(circshift(Ω_slice,-1)[1+mod(t+1,2):2:NX])
        Ω_slice = Ω_slice .* U[2,:,t]
    end
    V[2,:,NT] = Ω_slice
    return V
end





####    3-dimensional shenanigans    ####





# Construct a D-dimensional gaugefield with SU(2)-valued links
function gaugefield_SU2_cube(N_x::Int, N_t::Int, hot::Bool)
    U = Array{coeffs_SU2}(undef, 3, N_x, N_x, N_t)
    if hot
        for t = 1:N_t
            for y = 1:N_x
                for x = 1:N_x
                    for μ = 1:3
                        U[μ,x,y,t] = ran_SU2(rand())
                    end
                end
            end 
        end
    else
        for t = 1:N_t
            for y = 1:N_x
                for x = 1:N_x
                    for μ = 1:3
                        U[μ,x,y,t] = coeffs_Id_SU2()
                    end
                end
            end 
        end
    end
    return U
end


