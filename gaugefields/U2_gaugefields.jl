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

function Base.:*(α::Number, X::coeffs_U2)
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

# # Take the matrix logarithm of coeffs_U2
# # (note: staying in matrix-form or even returning a matrix
# # is significantly less efficient)
# function log_U2(Y::coeffs_U2)
#     z = sqrt(det(Y))
#     X = Y/z
#     ϕ = imag(log(z))
#     res = acos(X.a)/sqrt(1-X.a^2) * coeffs_U2(0.0*im, X.b, X.c, X.d)
#     return res + coeffs_U2(im*ϕ, 0.0*im, 0.0*im, 0.0*im)
# end

function grp2coeffs_u2(M::Matrix)
    x0 = im*imag(M[1,1]+M[2,2])/2
    x1 = Complex(imag(M[1,2]))
    x2 = Complex(real(M[1,2]))
    x3 = Complex(imag(M[1,1]-M[2,2])/2)
    return coeffs_U2(x0, x1, x2, x3)
end

function log_U2(X::coeffs_U2)
    M = coeffs2grp(X)
    E = eigen(M)
    A = E.vectors
    A_inv = adjoint(E.vectors)
    V = diagm(log.(E.values))
    return grp2coeffs_u2(A * V * A_inv)
end 

const ϵ_for_exp = 10^(-15) # for exp_u2() below
# Map a u(2)-element (from the physicists' Lie-algebra) onto the
# manifold U(2) via the exponential function
function exp_u2(Y::coeffs_U2)
    # X = coeffs_U2(Y.a, im*Y.b, im*Y.c, im*Y.d)
    ϵ = 10^(-15)
    A = exp(Y.a) * cos(sqrt(Y.b^2+Y.c^2+Y.d^2))
    B = exp(Y.a) * sin(sqrt(Y.b^2+Y.c^2+Y.d^2 + ϵ)) / sqrt(Y.b^2+Y.c^2+Y.d^2 + ϵ)
    return coeffs_U2(A, B*Y.b, B*Y.c, B*Y.d)
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
    r1, r2, r3, r4 = 2 .* (rand(4) .- 0.5)
    absr = sqrt(r1^2 + r2^2 + r3^2)
    ph_fac = exp(im*ϵ*π*r4) # Decidedly not exp(2*im*ϵ*π*r4), see above
    # ⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕
    return ph_fac*coeffs_U2(sqrt(1-ϵ^2), ϵ*r1/absr, ϵ*r2/absr, ϵ*r3/absr)
    # return sum([ph_fac*sqrt(1-ϵ^2), ph_fac*ϵ*r1/absr, ph_fac*ϵ*r2/absr, ph_fac*ϵ*r3/absr] .* Σ_im)
end

# Create a random element of u(2), i.e. the Lie algebra (mind the
# lower case letter!)
function ran_u2(ϵ)
    r = ϵ .* (2 .* rand(4) .- 1) .+ 0.0im
    r[1] = im*r[1]
    return coeffs_U2(r[1], r[2], r[3], r[4])
end

# Quickly get the coefficients of the identity element
function coeffs_Id_U2()
    return coeffs_U2(1.0+0.0*im, 0.0*im, 0.0*im, 0.0*im)
end

# Given coeffs_U2 create the corresponding U2-matrix
function coeffs2grp(X::coeffs_U2)::Matrix{ComplexF64}
    return X.a*σ0 + im*X.b*σ1 + im*X.c*σ2 + im*X.d*σ3
end

# function proj_U2(A::Matrix)
#     SVD = svd(A)
#     return SVD.U * SVD.Vt / sqrt(det(A))
# end

function proj2man(U::coeffs_U2)
    A = coeffs2grp(U)
    SVD = svd(A)
    A = SVD.U * SVD.Vt 
    return grp2coeffs_U2(A)
end

function proj2man_mat_U2(A)
    SVD = svd(A)
    return SVD.U * SVD.Vt 
end

# Given a U2-matrix return the corresponding coeffs_U2
function grp2coeffs_U2(mat)
    d = sqrt(det(mat))
    M = mat/d
    return coeffs_U2(d*real(M[1,1]), d*imag(M[1,2]), d*real(M[1,2]), d*imag(M[1,1]))
end

# Construct a square gauge fielC: an (2 × N_t × N_x)-Array with 
# entries which are SU(2)-valued, or rather in our case, coeffs_SU2-valued 
# (see struct "coeffs_SU2" above). This means that in order to access the
# link in μ-direction at space-time point n = (t,x), we need U[μ,t,x], where
# μ = 1 corresponds to the t-direction ("upwards") and 
# μ = 2 corresponds to the x-direction ("sideways").
function gaugefield_U2(N_x::Int64, N_t::Int64, hot::Bool)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t) # ⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕⭕
    # U = Array{Matrix}(undef, 2, N_x, N_t)
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

# gaugefield_U2(8, 8, true);

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





# Construct a hexagonal gauge field. See comment on hex_links_coords_chess() in
# gaugefeilds.jl. Every taboo link is constructed as (coeffs_U2 of) a NaN matrix.
function hexfield_U2(N_x::Int64, N_t::Int64, hot::Bool)
    @assert iseven(N_x)
    @assert iseven(N_t)
    U = [coeffs_U2(NaN*im,NaN*im,NaN*im,NaN*im) for μ = 1:2, x = 1:N_x, t = 1:N_t]
    coords = hex_links_coords_chess(N_x, N_t)
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


