include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")

using Plots
using LsqFit
using LaTeXStrings
using LinearAlgebra
using Roots
using DelimitedFiles
using Statistics
using Optim

function insta_action(β, N_c, N_x, N_t, Q, z)
    ReTrP = (N_c-1)*cos(2*π*z/N_x/N_t) + cos(2*π*(Q - (N_c-1)*z)/N_x/N_t ) 
    return β*N_x*N_t*(1 - ReTrP/N_c)
end

# function insta_U2_w(N_x, N_t, Q, z)
#     w = π*(2*z-Q)
#     U = Array{coeffs_U2}(undef, 2, N_x, N_t)
#     U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*w/N_x/N_t)*coeffs_Id_U2() - sin(t*w/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
#     U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
#     U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*w/N_x)*coeffs_Id_U2() + sin(x*w/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
#     return U
# end

function insta_U2_z(N_x, N_t, Q, z)
    # w = π*(2*z-Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * exp_u2(-t*2*π/N_x/N_t * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * exp_u2(x*2*π/N_x * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x]
    return U
end

# function insta_U2_slow(N_x, N_t, Q, z)
#     # w = π*(2*z-Q)
#     U = Array{coeffs_U2}(undef, 2, N_x, N_t)
#     U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * grp2coeffs_U2(exp(coeffs2grp(-t*π/N_x/N_t * (z-Q) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))))) for x = 1:N_x, t = 1:N_t]
#     U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
#     U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * grp2coeffs_U2(exp(coeffs2grp(x*π/N_x * (z-Q) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))))) for x = 1:N_x]
#     return U
# end

function insta_U2_comb(N_x, N_t, Q)
    insta = gaugefield_U2(N_x,N_t,false)
    x_vec = complex.([π/(1+mod(Q,2)), 0.0, 0.0])
    t_vec = complex.([0.0, π/(1+mod(Q,2)), 0.0])
    x_fac = exp_u2(coeffs_U2(0.0im, x_vec[1], x_vec[2], x_vec[3]))
    t_fac = exp_u2(coeffs_U2(0.0im, t_vec[1], t_vec[2], t_vec[3]))
    insta[1,1:N_x-1,:] = [exp(-im*Q*π*(t-1)/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x-1, t = 1:N_t]
    insta[1,N_x,:] = [exp(-im*Q*π*(t-1+N_x)/N_x/N_t) * x_fac for t = 1:N_t]
    insta[2,:,N_t] = [exp(im*Q*x*π/N_x) * t_fac for x = 1:N_x]
    return insta
end

# action(insta_U2_comb(N_x,N_t,1),1) - action(insta_U2(N_x,N_t,1),1)

function insta_U2_comb_odd(N_x, N_t, Q, M_rot)
    insta = gaugefield_U2(N_x,N_t,false)
    x_vec = M_rot * complex.([π/(1+mod(Q,2)), 0.0, 0.0])
    t_vec = M_rot * complex.([0.0, π/(1+mod(Q,2)), 0.0])
    x_fac = exp_u2(coeffs_U2(0.0im, x_vec[1], x_vec[2], x_vec[3]))
    t_fac = exp_u2(coeffs_U2(0.0im, t_vec[1], t_vec[2], t_vec[3]))
    insta[1,1:N_x-1,:] = [exp(-im*Q*π*(t-1)/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x-1, t = 1:N_t]
    insta[1,N_x,:] = [exp(-im*Q*π*(t-1+N_x)/N_x/N_t) * x_fac for t = 1:N_t]
    insta[2,:,N_t] = [exp(im*Q*x*π/N_x) * t_fac for x = 1:N_x]
    return insta
end

function insta_U2_comb_even(N_x, N_t, Q, coeffs_α_and_vec)
    ### coeffs_α_and_vec: first entry is alpha, second, third and fourth the vector. Needs to be
    ### one vector for optim() to work
    insta = gaugefield_U2(N_x,N_t,false)
    x_vec = complex.(coeffs_α_and_vec[2:4])
    t_vec = coeffs_α_and_vec[1] .* x_vec
    x_fac = exp_u2(coeffs_U2(0.0im, x_vec[1], x_vec[2], x_vec[3]))
    t_fac = exp_u2(coeffs_U2(0.0im, t_vec[1], t_vec[2], t_vec[3]))
    insta[1,1:N_x-1,:] = [exp(-im*Q*π*(t-1)/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x-1, t = 1:N_t]
    insta[1,N_x,:] = [exp(-im*Q*π*(t-1+N_x)/N_x/N_t) * x_fac for t = 1:N_t]
    insta[2,:,N_t] = [exp(im*Q*x*π/N_x) * t_fac for x = 1:N_x]
    return insta
end

# function insta_U2_comb_even(N_x, N_t, Q, α, x_vector)
#     ### coeffs_α_and_vec: first entry is alpha, second, third and fourth the vector. Needs to be
#     ### one vector for optim() to work
#     insta = gaugefield_U2(N_x,N_t,false)
#     x_vec = complex.(x_vector)
#     t_vec = α .* x_vec
#     x_fac = exp_u2(coeffs_U2(0.0im, x_vec[1], x_vec[2], x_vec[3]))
#     t_fac = exp_u2(coeffs_U2(0.0im, t_vec[1], t_vec[2], t_vec[3]))
#     insta[1,1:N_x-1,:] = [exp(-im*Q*π*(t-1)/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x-1, t = 1:N_t]
#     insta[1,N_x,:] = [exp(-im*Q*π*(t-1+N_x)/N_x/N_t) * x_fac for t = 1:N_t]
#     insta[2,:,N_t] = [exp(im*Q*x*π/N_x) * t_fac for x = 1:N_x]
#     return insta
# end

function insta_U2_z_comb(N_x, N_t, Q, z)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,1:N_x-1,:] = [exp(-im*Q*π*(t-1)/N_x/N_t) * exp_u2(-(t-1)*2*π/N_x/N_t * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x-1, t = 1:N_t]
    U[1,N_x,:]     = [exp(-im*Q*π*(t-1+N_x)/N_x/N_t) * exp_u2(-(t-1+N_x)*2*π/N_x/N_t * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * exp_u2(x*2*π/N_x * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x]
    return U
end

# q,z = rand(1:10,2)
# @assert isapprox(action(insta_U2_z(32,32,q,z),1), action(insta_U2_z_comb(32,32,q,z),1))

function insta_U2_z_comb(N_x, N_t, Q, z, coeffs)
    ### coeffs: 2 real numbers telling how to stretch the outer slices with
    ### the factor exp(i*coeffs[1]*σ_3) resp. exp(i*coeffs[2]*σ_3)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,1:N_x-1,:] = [exp(-im*Q*π*(t-1)/N_x/N_t) * exp_u2(-(t-1)*2*π/N_x/N_t * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x-1, t = 1:N_t]
    U[1,N_x,:]     = [exp(-im*Q*π*(t-1+N_x)/N_x/N_t) * exp_u2(-(t-1+N_x)*2*π/N_x/N_t * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * exp_u2(x*2*π/N_x * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x]
    for t = 1:N_t
        U[1,N_x,t] = exp_u2(coeffs[1]*coeffs_U2(0.0im,0.0im,0.0im,complex(1.0))) * U[1,N_x,t]
    end
    for x = 1:N_x
        U[2,x,N_t] = exp_u2(coeffs[2]*coeffs_U2(0.0im,0.0im,0.0im,complex(1.0))) * U[2,x,N_t]
    end
    return U
end

function two_metric(M::Matrix, N::Matrix)
    O = M-N
    return real(sqrt(sum(O.*conj(O))))
end

function two_metric(X::coeffs_U2, Y::coeffs_U2)
    Z = X-Y
    return real(sqrt(2)*sqrt(conj(Z.a)*Z.a + conj(Z.b)*Z.b + conj(Z.c)*Z.c + conj(Z.d)*Z.d))
end

function two_metric_squared(X::coeffs_U2, Y::coeffs_U2)
    Z = X-Y
    return real(2* (conj(Z.a)*Z.a + conj(Z.b)*Z.b + conj(Z.c)*Z.c + conj(Z.d)*Z.d))
end

function two_metric_field(U, V)
    NX = size(U,2)
    NT = size(U,3)
    return sqrt(sum([two_metric_squared(U[μ,x,t], V[μ,x,t]) for μ = 1:2, x = 1:NX, t = 1:NT ])) / 2*NX*NT
end

# Needs M to be ∈ SU(2)!
function Pauli_coeffs(M::Matrix)
    N  = -im*log(M)
    v1 = real(N[2,1])
    v2 = imag(N[2,1])
    v3 = real(N[1,1])
    return [v1,v2,v3]
    # abs_v = acos(real(M[1,1]))
    # fac = abs_v/sin(abs_v)
    # v1 = fac*imag(M[1,2])
    # v2 = fac*real(M[1,2])
    # v3 = fac*imag(M[1,1])
    # return [v1,v2,v3]
end

# # Needs M to be ∈ SU(2)!
# function Pauli_coeffs(X::coeffs_U2)
#     M = coeffs2grp(X)
#     Pauli_coeffs(M)
# end

# Needs M to be ∈ SU(2)!
function Pauli_coeffs(X::coeffs_U2)
    Y = log_U2(X)
    return [Y.b, Y.c, Y.d]
end

# M1 = coeffs2grp(ran_U2(rand()))
# M1 = M1/sqrt(det(M1))
# isapprox(exp(im*sum(Pauli_coeffs(M1) .* [σ1,σ2,σ3])), M1)


function find_rot_mat(a,b)
    v = cross(a,b)
    s = sqrt(v'*v)
    c = a'*b
    M_v = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]
    return I(3) + M_v + (1-c)/s^2 * M_v^2
end

# a1 = rand(3)
# a2 = rand(3)
# a1 = a1/sqrt(a1'*a1)
# a2 = a2/sqrt(a2'*a2)
# M1 = find_rot_mat(a1,a2)
# isapprox(M1 * a1, a2)


function rot_mat2quat(M::Matrix)
    x0x1 = (M[2,3]-M[3,2])/4
    x0x2 = (M[3,1]-M[1,3])/4
    x0x3 = (M[1,2]-M[2,1])/4
    x1x3 = (M[1,3]+M[3,1])/4
    x0_over_x1 = x0x3/x1x3
    x0 = sqrt(x0x1 * x0_over_x1)
    x1 = x0x1/x0
    x2 = x0x2/x0
    x3 = x0x3/x0
    x0, x1, x2, x3 = complex.([x0, x1, x2, x3])
    return coeffs_U2(x0, x1, x2, x3)
end

function vec2Pauli(v)
    return im*sum(v.*[σ1,σ2,σ3])
end

function Pauli2vec(M::Matrix)
    v1 = imag(M[1,2])
    v2 = real(M[1,2])
    v3 = imag(M[1,1])
    return [v1,v2,v3]
end

function Pauli2vec(X::coeffs_U2)
    return [X.b, X.c, X.d]
end

# heatmap([real(tr(plaq(insta_U2_comb_odd(N_x, N_t, 1, λ0),x,t))) for x = 1:N_x, t = 1:N_t])
# action(insta_U2_comb_odd(N_x, N_t, 1, λ0),1) - action(insta_max,1)
# action(insta_U2_comb_odd(N_x, N_t, 1, λ0),1) - action(insta_U2(N_x,N_t,1),1)
# action(insta_U2_comb_odd(N_x, N_t, 2, λ0),1) - action(insta_U2(N_x,N_t,2),1)

function two_metric_field_insta_rot(U, M_rot_Lie_coeffs)
    NX = size(U,2)
    NT = size(U,3)
    q = round(Int,top_charge_U2(U))
    L1 = [0 0 0; 0 0 -1; 0 1 0]
    L2 = [0 0 1; 0 0 0; -1 0 0]
    L3 = [0 -1 0; 1 0 0; 0 0 0]
    M_rot = exp(sum(M_rot_Lie_coeffs.*[L1,L2,L3]))
    insta = insta_U2_comb_odd(NX,NT,q,M_rot)
    return sqrt(sum([two_metric_squared(U[μ,x,t], insta[μ,x,t]) for μ = 1:2, x = 1:NX, t = 1:NT ])) / 2*NX*NT
end

# two_metric_field_insta_rot(blabla, [0.0,0.0,0.0])

# a = Pauli2vec(bla_max[1,N_x,1]/sqrt(det(bla_max[1,N_x,1])))
# b = Pauli_coeffs(bla_max[1,N_x,1]/sqrt(det(bla_max[1,N_x,1])))

function minimum_insta_metric(U)
    NX = size(U,2)
    NT = size(U,3)
    V = max_gauge(U,"U2")
    q = round(Int,top_charge_U2(U))
    if iseven(q)
        U1_fac_ratio_outer = exp(-im*q*π*(1+NX)/NX/NT) / sqrt(det(V[1,NX,2])) 
        for t = 1:NT
            V[1,NX,t] = U1_fac_ratio_outer * V[1,NX,t]
        end
        U1_fac_ratio_outer = exp(im*q*π/NX) / sqrt(det(V[2,1,NT]))
        for x = 1:NX
            V[2,x,NT] = U1_fac_ratio_outer * V[2,x,NT]
        end
        v_x = Pauli_coeffs(V[1,NX,1]/sqrt(det(V[1,NX,1])))
        optim_metric(α) = two_metric_field(V,insta_U2_comb_even(NX,NT,q,[α,v_x[1],v_x[2],v_x[3]]))
        return minimum([optimize(optim_metric,[α],NelderMead()).minimum for α in -1.0:1:1.0])
        # V[1,NX,:] = [-exp(-im*q*π*(t-1+NX)/NX/NT) * coeffs_Id_U2() for t = 1:NT]
        # V[2,:,NT] = [-exp(im*q*x*π/NX) * coeffs_Id_U2() for x = 1:NX]
    else # if isodd(q)
        v_x = Pauli2vec(V[1,NX,1]/sqrt(det(V[1,NX,1])))
        v_t = Pauli2vec(V[2,1,NT]/sqrt(det(V[2,1,NT])))
        v_x_insta = [1.0, 0.0, 0.0]
        v_t_insta = [0.0, 1.0, 0.0]
        M1 = find_rot_mat(v_x, v_x_insta)
        M2 = find_rot_mat(M1*v_t, v_t_insta)
        rot_quat = rot_mat2quat(M2*M1)
        U1_fac_ratio_outer = exp(-im*q*π*(1+NX)/NX/NT) / sqrt(det(V[1,NX,2])) 
        for t = 1:NT
            V[1,NX,t] = U1_fac_ratio_outer * rot_quat * V[1,NX,t] * adjoint(rot_quat)
        end
        U1_fac_ratio_outer = exp(im*q*π/NX) / sqrt(det(V[2,1,NT]))
        for x = 1:NX
            V[2,x,NT] = U1_fac_ratio_outer * rot_quat * V[2,x,NT] * adjoint(rot_quat)
        end
        return two_metric_field(V,insta_U2_comb(NX,NT,q))
    end
end

# bla = gaugefield_U2(32,32,true);
# top_charge_U2(bla)
# for i = 1:200 chess_metro!(bla,0.1,6.0,[0.0],"U2") end
# bla = stout_midpoint(bla, 5e3, 0.1)
# minimum_insta_metric(bla)
optimize_insta_metric(bla)
# action(bla,1)
# insta_action(1,2,32,32,-4,-2)

function optimize_insta_metric(U)
    NX = size(U,2)
    NT = size(U,3)
    q = round(Int,top_charge_U2(U))
    
    V = max_gauge(U,"U2")
    U1_fac_ratio_outer_x = exp(-im*q*π*(1+NX)/NX/NT)/sqrt(det(V[1,NX,2])) 
    for t = 1:NT
        V[1,NX,t] = U1_fac_ratio_outer_x * V[1,NX,t]
    end
    U1_fac_ratio_outer_t = exp(im*q*π/NX) / sqrt(det(V[2,1,NT]))
    for x = 1:NX
        V[2,x,NT] = U1_fac_ratio_outer_t * V[2,x,NT]
    end
    if iseven(q)
        optim_metric_even(coeffs) = two_metric_field(V, insta_U2_comb_even(NX,NT,q,coeffs))
        return minimum([optimize(optim_metric_even,start_coeffs,NelderMead()).minimum for start_coeffs in [zeros(4), ones(4), 2 .* rand(4) .-1]])
    else
        optim_metric_odd(coeffs) = two_metric_field_insta_rot(V, coeffs)
        return minimum([optimize(optim_metric_odd,start_coeffs,NelderMead()).minimum for start_coeffs in [zeros(3), ones(3), 2 .* rand(3) .-1]])
    end
end

# bla = gaugefield_U2(32,32,true);
# for i = 1:5000 chess_metro!(bla, 0.1, 5.0, [0.0], "U2") end
# q_temp = top_charge_U2(bla)
# for i = 1:5000 bla = stout_midpoint(bla,0.1) end
# q_temp = top_charge_U2(bla)
# two_metric_field(max_gauge(bla, "U2"), insta_U2_comb(32,32,q_temp))
# minimum_insta_metric(bla)
# minimum([optimize_insta_metric(bla,zeros(4)),optimize_insta_metric(bla,ones(4)),optimize_insta_metric(bla,rand(4))]) 
# # coeffs2grp(bla[1,32,1]/sqrt(det(bla[1,32,1])))

function optimize_special_metric(U, start_coeffs, z)
    @assert length(start_coeffs) == 2 "start_coeffs must number 2, for each outer slice to be optimized by exp(i*start_coeffs[1/2]*σ_3)"
    NX = size(U,2)
    NT = size(U,3)
    q = round(Int,top_charge_U2(U))
    
    V = max_gauge(U,"U2")
    U1_fac_ratio_outer_x = exp(-im*q*π*(1+NX)/NX/NT)/sqrt(det(V[1,NX,2])) 
    for t = 1:NT
        V[1,NX,t] = U1_fac_ratio_outer_x * V[1,NX,t]
    end
    U1_fac_ratio_outer_t = exp(im*q*π/NX) / sqrt(det(V[2,1,NT]))
    for x = 1:NX
        V[2,x,NT] = U1_fac_ratio_outer_t * V[2,x,NT]
    end
    
    optim_metric(coeffs) = two_metric_field(V, insta_U2_z_comb(NX,NT,q,z,coeffs))
    return optimize(optim_metric,start_coeffs,NelderMead()).minimum
end





################################################################################





# @assert 1==0 "Do we really want to start a smearing run???"
β   = 6.0
L   = 32
N_x = L
N_t = L
hot = true
ρ   = 0.1
N_therm   = 500
N_meas    = 20
N_sepa    = 10
N_metro   = 1
N_over    = 3
N_smear   = Int(5e4)
acc_wish  = 0.8
ϵ         = 0.1
m_smear_inds = collect(2000:2000:N_smear)

base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data_new"    
actions_path = string(base_path,"\\non_smeared_actions.txt")
params_path = string(base_path, "\\params.txt")


params = "Square Simulation with:
β            = $β
N_t          = $N_t
N_x          = $N_x 
hot          = $hot

N_therm      = $N_therm
N_metro      = $N_metro
N_over       = $N_over
N_sepa       = $N_sepa
N_meas       = $N_meas
acc_wish     = $acc_wish

N_smear      = $N_smear
ρ            = $ρ"

bla = open(params_path, "a")
write(bla, params)
close(bla)


acc_metro = [0.0]
acc_over  = [0.0]
# actions_therm = []
U = gaugefield(N_x, N_t, hot, "U2", "square")
for therm = 1:N_therm
    chess_metro!(U,ϵ,β,acc_metro,"U2")
    ϵ *= sqrt(acc_metro[1] / acc_wish) # only update ϵ acc. to Metropolis
    # push!(actions_therm, action(U,β))
end
# plot(actions_therm)
actions = []
for meas = 1:N_meas
    for sepa = 1:N_sepa
        for met = 1:N_metro
            chess_metro!(U,ϵ,β,acc_metro,"U2")
        end
        for over = 1:N_over
            chess_overrelax!(U,acc_over)
        end
        push!(actions, action(U,β))
    end
    smeared_actions = [action(U,β)]
    smeared_charges = [top_charge_U2(U)]
    smeared_m_opt   = [optimize_insta_metric(U)]
    smeared_m_anal  = [minimum_insta_metric(U)]
    V = stout_midpoint(U,ρ)
    push!(smeared_actions,action(V,β))
    push!(smeared_charges,top_charge_U2(V))
    count = 0
    for smear = 1:N_smear
        V = stout_midpoint(V,ρ)
        push!(smeared_actions,action(V,β))
        push!(smeared_charges,top_charge_U2(V))
        if smear in m_smear_inds
            push!(smeared_m_opt, optimize_insta_metric(V))
            push!(smeared_m_anal, minimum_insta_metric(V))
        end
        if smear%Int(N_smear/20) == 0
            count += 5
            println("Measurement Nr.: $meas, Smearing Progress: $count%")
        end
    end # smear
    # push!(smeared_m_opt, optimize_insta_metric(V))
    # push!(smeared_m_anal, minimum_insta_metric(V))
    S_path = string(base_path,"\\sms_measnr_$meas.txt")
    Q_path = string(base_path,"\\smq_measnr_$meas.txt")
    smeared_m_opt_path = string(base_path, "\\smeared_m_opt_measnr_$meas.txt")
    smeared_m_anal_path = string(base_path, "\\smeared_m_anal_measnr_$meas.txt")
    writedlm(S_path, smeared_actions)
    writedlm(Q_path, smeared_charges)
    writedlm(smeared_m_opt_path, smeared_m_opt)
    writedlm(smeared_m_anal_path, smeared_m_anal)
end # meas
writedlm(actions_path, actions)
println("We're done here!")






base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data_new\\prel_smears"
# actions_path = string(base_path, "non_smeard_actions.txt")    
# actions = readdlm(actions_path)
# plot(actions)

for meas = 1:1
    show_s      = true
    show_z      = true
    show_q      = false
    show_m_opt  = false
    show_m_anal = false

    z = -1
    ### in prel_smears interesting configs:
    # meas 1    q -3   z -1/-2
    # meas 4    q  9   z  4/ 5
    # meas 10   q -1   z  0/-1

    S_path = string(base_path,"\\sms_measnr_$meas.txt")
    Q_path = string(base_path,"\\smq_measnr_$meas.txt")
    smeared_m_opt_path = string(base_path, "\\smeared_m_opt_measnr_$meas.txt")
    smeared_m_anal_path = string(base_path, "\\smeared_m_anal_measnr_$meas.txt")
    smeared_actions = readdlm(S_path)
    smeared_charges = readdlm(Q_path)
    smeared_m_opt   = readdlm(smeared_m_opt_path)
    smeared_m_anal  = readdlm(smeared_m_anal_path) 
    flow_times   = ρ.*collect(0:N_smear+1)
    flow_times_m = ρ.*vcat([0],m_smear_inds)
    start_τ_s = 10
    start_τ_q = 0
    start_τ_m = 600
    start_ind_s = findall(x->x==start_τ_s, flow_times)[1]
    start_ind_q = findall(x->x==start_τ_q, flow_times)[1]
    start_ind_m = findall(x->x==start_τ_m, flow_times_m)[1]
    skip_s = 5

    q = round(Int,last(smeared_charges))
    
    if show_s
        image_s = plot(
            flow_times[start_ind_s:skip_s:end], 
            smeared_actions[start_ind_s:skip_s:end] ./ β,
            title = latexstring("\$S/\\beta\$ during smearing\n\$\\beta = $β, L = $L, \\rho = $ρ, q = $q,\\, \\mathrm{Nr. meas.} = $meas\$"),
            xlabel = latexstring("flow time \$\\tau\$"),
            label = :false,
            yaxis = :log
        )
        if show_z
            hline!([insta_action(β,2,L,L,q,z)/β], label = latexstring("\$S/\\beta\$ of spec. conf. for \$(q,z) = ($q,$z)\$"))
        end
        display(image_s)
    end

    if show_q
        image_q = plot(
            flow_times[start_ind_q:end], 
            smeared_charges[start_ind_q:end],
            title = latexstring("\$q\$ during smearing\n\$\\beta = $β, L = $L, \\rho = $ρ, \\mathrm{Nr. meas.} = $meas\$"),
            xlabel = latexstring("flow time \$\\tau\$"),
            label = :false
        )
        display(image_q)
    end

    if show_m_opt
        image_m_opt = scatter(
            flow_times_m[start_ind_m:end], 
            smeared_m_opt[start_ind_m:end-1],
            title = latexstring("optim. \$||U-\\mathrm{inst.}||\$ during smearing\n\$\\beta = $β, L = $L, \\rho = $ρ, q = $q,\\, \\mathrm{Nr. meas.} = $meas\$"),
            xlabel = latexstring("flow time \$\\tau\$"),
            label = :false,
            yaxis = :log
        )
        display(image_m_opt)
    end

    if show_m_anal
        image_m_anal = scatter(
            flow_times_m[start_ind_m:end], 
            smeared_m_anal[start_ind_m:end-1],
            title = latexstring("anal. \$||U-\\mathrm{inst.}||\$ during smearing\n\$\\beta = $β, L = $L, \\rho = $ρ, q = $q,\\, \\mathrm{Nr. meas.} = $meas\$"),
            xlabel = latexstring("flow time \$\\tau\$"),
            label = :false,
            yaxis = :log
        )
        display(image_m_anal)
    end
end