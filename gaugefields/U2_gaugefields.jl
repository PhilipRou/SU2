################################################################################
# include("SU2_gaugefields.jl")
using LinearAlgebra



# Essentially the same as coeffs_SU2, just for the group U(2). For M ∈ U(2)
# there is one additional parameter: a phase φ, such that det(M) = exp(iφ).
# The most efficient way is to just copy coeffs_SU2 and allow Complex types as
# as values for a,...,d. 
mutable struct coeffs_U2{T <: Number}
    a::T
    b::T
    c::T
    d::T
    function coeffs_U2(
        a::T,
        b::T,
        c::T,
        d::T
        ) where {T <: Number}
        return new{T}(a,b,c,d)
    end
end

# U(2) coeffs have to square up to exp(iφ)! But sometimes we need to add 
# matrices and leave the group, e.g. staple.
function Base.:+(X::coeffs_U2, Y::coeffs_U2)
    return coeffs_U2(X.a+Y.a, X.b+Y.b, X.c+Y.c, X.d+Y.d)
end

function Base.:-(X::coeffs_U2, Y::coeffs_U2)
    return coeffs_U2(X.a-Y.a, X.b-Y.b, X.c-Y.c, X.d-Y.d)
end

# Multiply two U(2)-matrices whose coefficients are given in X and Y
function Base.:*(X::coeffs_U2, Y::coeffs_U2)
    a = X.a*Y.a - X.b*Y.b - X.c*Y.c - X.d*Y.d 
    b = X.a*Y.b + X.b*Y.a - X.c*Y.d + X.d*Y.c
    c = X.a*Y.c + X.b*Y.d + X.c*Y.a - X.d*Y.b
    d = X.a*Y.d - X.b*Y.c + X.c*Y.b + X.d*Y.a
    return coeffs_U2(a,b,c,d)
end

function Base.:*(α, X::coeffs_U2)
    return coeffs_U2(α*X.a, α*X.b, α*X.c, α*X.d)
end

function Base.:/(X::coeffs_U2, α)
    return coeffs_U2(X.a/α, X.b/α, X.c/α, X.d/α)
end

function Base.:(==)(X::coeffs_U2, Y::coeffs_U2)
    return [X.a,X.b,X.c,X.d] == [Y.a,Y.b,Y.c,Y.d]
end

function Base.isapprox(X::coeffs_U2, Y::coeffs_U2)
    return isapprox([X.a,X.b,X.c,X.d], [Y.a,Y.b,Y.c,Y.d])
end

function LinearAlgebra.tr(X::coeffs_U2)
    return 2*X.a
end

function LinearAlgebra.det(X::coeffs_U2)
    return X.a^2 + X.b^2 + X.c^2 + X.d^2 
end

# Take the adjoint of coeffs_U2, i.e. the adjoint of the matrix 
# corresponding to said coeffs_U2
function LinearAlgebra.adjoint(X::coeffs_U2)
    return coeffs_U2(adjoint(X.a), -adjoint(X.b), -adjoint(X.c), -adjoint(X.d))
end

# ❗ Very inefficient, only for debugging purposes ❗
function get_array(X::coeffs_U2)
    return [X.a, X.b, X.c, X.d]
end

# Essentially the same as ran_SU2, but with the additional phase factor 
# due to det(M) = exp(iφ), φ ∈ [0,2π). To that end generate a random 
# SU(2)-matrix and multiply it with exp(iφ/2).
function ran_U2(ϵ)
    # r1 = 2 * (rand()-0.5)
    # r2 = 2 * (rand()-0.5)
    # r3 = 2 * (rand()-0.5)
    r1, r2, r3 = 2 .* (rand(3) .- 0.5)
    absr = sqrt(r1^2 + r2^2 + r3^2)
    ph_fac = exp(im*ϵ*π*rand()) # Decidedly not exp(2*im*ϵ*π*rand()), see above
    return coeffs_U2(ph_fac*sqrt(1-ϵ^2), ph_fac*ϵ*r1/absr, ph_fac*ϵ*r2/absr, ph_fac*ϵ*r3/absr)
end

# Quickly get the coefficients of the identity element
function coeffs_Id_U2()
    return coeffs_U2(1.0+0.0*im, 0.0*im, 0.0*im, 0.0*im)
end

# Project coeffs_U2 onto U2 (since addition is allowed it may happen that
# some coeffs_U2 do not describe a U2 element anymore)
function proj_U2(X::coeffs_U2)
    return X/sqrt(abs(det(X)))
end

# Given coeffs_U2 create the corresponding U2-matrix
function coeffs2grp(X::coeffs_U2)
    return X.a*σ0 + im*X.b*σ1 + im*X.c*σ2 + im*X.d*σ3
end

# Given a U2-matrix return the corresponding coeffs_U2
function grp2coeffs_U2(mat)
    d = sqrt(det(mat))
    M = mat/d
    return coeffs_U2(d*real(M[1,1]), d*imag(M[1,2]), d*real(M[1,2]), d*imag(M[1,1]))
end

# Construct a square gauge field: an (2 × N_t × N_x)-Array with 
# entries which are SU(2)-valued, or rather in our case, coeffs_SU2-valued 
# (see struct "coeffs_SU2" above). This means that in order to access the
# link in μ-direction at space-time point n = (t,x), we need U[μ,t,x], where
# μ = 1 corresponds to the t-direction ("upwards") and 
# μ = 2 corresponds to the x-direction ("sideways").
function gaugefield_U2(N_x::Int64, N_t::Int64, hot::Bool)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    if hot
        for t = 1:N_t
            for x = 1:N_x
                for μ = 1:2
                    U[μ,x,t] = ran_U2(rand())
                end
            end
        end
    else 
        for t = 1:N_t
            for x = 1:N_x
                for μ = 1:2
                    U[μ,x,t] = coeffs_Id_U2() # coeffs of identity ∈ SU(2)
                end
            end
        end       
    end
    return U
end

# Apply temporal gauge onto a square config
function temp_gauge_U2(U)
    NX = size(U,2)
    NT = size(U,3)
    V = gaugefield_U2(NX, NT, false)
    Ω_slice = [coeffs_Id_U2() for x = 1:NX] 
    for t = 1:NT
        V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice,-1))
        Ω_slice = Ω_slice .* U[2,:,t]
    end
    V[2,:,NT] = Ω_slice
    return V
end





####    Hexagonal shenanigans   #### 





# Construct a hexagonal gauge field. See comment on chess_hex_link_coords() in
# gaugefeilds.jl. Every taboo link is constructed as (coeffs_U2 of) a NaN matrix.
function hexfield_U2(N_x::Int64, N_t::Int64, hot::Bool)
    @assert iseven(N_x)
    @assert iseven(N_t)
    U = [coeffs_U2(NaN*(1+im),NaN*(1+im),NaN*(1+im),NaN*(1+im)) for μ = 1:2, x = 1:N_x, t = 1:N_t]
    coords = chess_hex_link_coords(N_x, N_t)
    if hot
        for coord in coords
            U[coord[1], coord[2], coord[3]] = ran_U2(rand())
        end
    else 
        for coord in coords
            U[coord[1], coord[2], coord[3]] = coeffs_Id_U2()
        end
    end
    return U
end

# Apply temporal gauge onto a square config
# ❗❗❗ Works but is significantly less efficient than temp_gauge(),
# which also works on hexagonal configs
function temp_gauge_hex_U2(U)
    NX = size(U,2)
    NT = size(U,3)
    V = hexfield_U2(NX, NT, false)
    Ω_slice = [coeffs_Id_U2() for x = 1:NX] 
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
function gaugefield_U2_cube(N_x::Int, N_t::Int, hot::Bool)
    U = Array{coeffs_U2}(undef, 3, N_x, N_x, N_t)
    if hot
        for t = 1:N_t
            for y = 1:N_x
                for x = 1:N_x
                    for μ = 1:3
                        U[μ,x,y,t] = ran_U2(rand())
                    end
                end
            end 
        end
    else
        for t = 1:N_t
            for y = 1:N_x
                for x = 1:N_x
                    for μ = 1:3
                        U[μ,x,y,t] = coeffs_Id_U2()
                    end
                end
            end 
        end
    end
    return U
end


