include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")


# Create a minimum of the gauge action in the topological
# sector of charge Q. Unfortunately the link values of the 
# odd-Q-instantons are not elements of the center of U(2)
# (i.e. ∉ U(1)), so that 
function insta_U2(N_x, N_t, Q)
    # U = gaugefield_U2(N_x,N_t,false)
    # if iseven(Q)
    #     U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x, t = 1:N_t]
    #     U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    #     U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * coeffs_Id_U2() for x = 1:N_x]
    # else
    #     # U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*π/N_x/N_t)*coeffs_Id_U2() - sin(t*π/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
    #     # U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    #     # U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*π/N_x)*coeffs_Id_U2() + sin(x*π/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    #     ### For more on the below, see "insta_half_int" somewhere else in this pile of code
    #     c1 = convert.(ComplexF64, π/2*[1,0,0])
    #     c2 = convert.(ComplexF64, π/2*[0,1,0])
    #     U[1,1:N_x-1,:] = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x-1, t = 1:N_t]
    #     U[1,N_x,:] = [exp(-im*Q*t*π/N_x/N_t) * exp_u2(coeffs_U2(0.0im, c1[1], c1[2], c1[3])) for t = 1:N_t]
    #     U[2,:,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1:N_x]
    # end
    # return U
    U = gaugefield_U2(N_x,N_t,false)
    x_fac = exp_u2(coeffs_U2(0.0im, complex(π/(1+mod(Q,2))), 0.0im, 0.0im))
    t_fac = exp_u2(coeffs_U2(0.0im, 0.0im, complex(π/(1+mod(Q,2))), 0.0im))
    U[1,1:N_x-1,:] = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x-1, t = 1:N_t]
    U[1,N_x,:] = [exp(-im*Q*t*π/N_x/N_t) * x_fac for t = 1:N_t]
    U[2,:,N_t] = [exp(im*Q*x*π/N_x) * t_fac for x = 1:N_x]
    return U
end

# Create a naive (N_x × N_t)-Q-instanton configuration, in
# the sense that it is not an actual minimum of the gauge
# action for odd Q (see insta_U2() below)
function insta_U2_naive(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * coeffs_Id_U2() for x = 1:N_x]
    return U
end

function def_not_insta_U2(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*π/N_x/N_t)*coeffs_Id_U2() - sin(t*π/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*π/N_x)*coeffs_Id_U2() + sin(x*π/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    return U
end

function insta_U2_split(N_x, N_t, Q, v)
    # M = grp2coeffs_U2(exp(sum(im .* v .* [σ1,σ2,σ3])))
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * grp2coeffs_U2(exp(sum(-(im*t*π/N_x/N_t) .* v .* [σ1,σ2,σ3]))) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * grp2coeffs_U2(exp(sum((im*x*π/N_x) .* v .* [σ1,σ2,σ3]))) for x = 1:N_x]
    return U
end

# The same as insta_U2(), only that instead of the Lie group
# elements now Lie algebra elements are produced (i.e. the
# logarithm of each link of insta_U2())
function insta_U2_log(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    δ = Complex(1.0e-15)    # To ensure that exp_u2() does not produce NaN's
    if iseven(Q)
        U[1,:,:]       = [coeffs_U2(-im*Q*t*π/N_x/N_t, δ, δ, δ) for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [coeffs_U2(δ, δ, δ, δ) for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [coeffs_U2(im*Q*x*π/N_x, δ, δ, δ) for x = 1:N_x]
    else
        U[1,:,:]       = [coeffs_U2(-im*Q*t*π/N_x/N_t, δ,δ, -t*π/N_x/N_t + 0.0im)  for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [coeffs_U2(δ,δ,δ,δ) for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [coeffs_U2(im*Q*x*π/N_x, δ,δ, x*π/N_x + 0.0im) for x = 1:N_x]
    end
    return U
end

function insta_U2_log_cheat(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    δ = Complex(1.0e-15)    # To ensure that exp_u2() does not produce NaN's
    # if iseven(Q)
        U[1,:,:]       = [coeffs_U2(-im*Q*t*π/N_x/N_t, δ, δ, δ) for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [coeffs_U2(δ, δ, δ, δ) for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [coeffs_U2(im*Q*x*π/N_x, δ, δ, δ) for x = 1:N_x]
    # else
    #     U[1,:,:]       = [coeffs_U2(-im*Q*t*π/N_x/N_t, δ,δ, -t*π/N_x/N_t + 0.0im)  for x = 1:N_x, t = 1:N_t]
    #     U[2,:,1:N_t-1] = [coeffs_U2(δ,δ,δ,δ) for x = 1:N_x, t = 1:N_t-1]
    #     U[2,:,N_t]     = [coeffs_U2(im*Q*x*π/N_x, δ,δ, x*π/N_x + 0.0im) for x = 1:N_x]
    # end
    return U
end