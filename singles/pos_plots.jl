include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\observables\\smearing.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\analyze\\SU2_analyze_head.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\analyze\\SU2_jackknives.jl")


using Plots
using LsqFit
using LaTeXStrings
using LinearAlgebra
# using Roots
using DelimitedFiles
using Statistics
using Optim
using QuadGK
using Measures

data_path ="C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Lattice_projects\\Lattice2024\\PoS" 
fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Lattice_projects\\Lattice2024\\PoS"

cb_hexes  = [
    "#377eb8", # blue
    "#ff7f00", # orange
    "#4daf4a", # green
    "#984ea3", # purple
    "#e41a1c", # red
    "#f781bf", # pink
    "#999999", # grey
    "#a65628", # brown
]
cb_colors = parse.(Colorant, cb_hexes)
cb_blue, cb_orange, cb_green, cb_purple, cb_red, cb_pink, cb_grey, cb_brown  = cb_colors

function insta_action_U2_min(β, N_x, N_t, q)
    return β*N_x*N_t*(1-cos(q*π/N_x/N_t))
end

function insta_action_U2(β, N_x, N_t, q, z)
    N_c = 2
    ReTrP = (N_c-1)*cos(2*π*z/N_x/N_t) + cos(2*π*(q - (N_c-1)*z)/N_x/N_t ) 
    return β*N_x*N_t*(1 - ReTrP/N_c)
end

function insta_action_Nc(β, N_c, N_x, N_t, Q, z)
    ReTrP = (N_c-1)*cos(2*π*z/N_x/N_t) + cos(2*π*(Q - (N_c-1)*z)/N_x/N_t ) 
    return β*N_x*N_t*(1 - ReTrP/N_c)
end

function analytic_susc_U1(β)
    kern(ϕ) = ϕ^2 * exp(β*cos(ϕ))
    return quadgk(kern, -π, π)[1] / besseli(0,β) / (2*π)^3
end

function analytic_plaq_U2(β)
    # b = big(β)
    b = β
    numer(α) = besseli(0,b*cos(α)) + besseli(2,b*cos(α))
    denom(α) = 2*besseli(1,b*cos(α))/cos(α)
    return quadgk(numer,0,π/2)[1]/quadgk(denom,0,π/2)[1] - 1/b
end

function analytic_plaq_U2_approx(β)
    b = big(β)
    numer(α) = exp(b*cos(α)) / sqrt(cos(α)) * 2
    denom(α) = 2 * exp(b*cos(α)) / sqrt(cos(α))^3
    return quadgk(numer,0,π/2)[1]/quadgk(denom,0,π/2)[1] - 1/b
end

function analytic_plaq_U2_approx_2(β)
    b = big(β)
    numer(α) = exp(b*cos(α)) / sqrt(cos(α)) * (2 - 14/(8*b*cos(α)))
    denom(α) = 2 * exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)))
    return quadgk(numer,0,π/2)[1]/quadgk(denom,0,π/2)[1] - 1/b
end

function analytic_plaq_U2_approx_3(β)
    b = big(β)
    numer(α) = exp(b*cos(α)) / sqrt(cos(α)) * (2 - 14/(8*b*cos(α)) + 21/(32*(b*cos(α))^2))
    denom(α) = 2 * exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)) -15/(128*(b*cos(α))^2))
    return quadgk(numer,0,π/2)[1]/quadgk(denom,0,π/2)[1] - 1/b
end

function analytic_susc_U2(β)
    b = β
    denom(α)   = besseli(1,b*cos(α))/cos(α)
    numer(α) = α^2 * besseli(1,b*cos(α))/cos(α)
    return quadgk(numer,0,π/2)[1] / quadgk(denom,0,π/2)[1] / π^2
end

function analytic_susc_U2_approx(β)
    b = big(β)
    denom(α)   = exp(b*cos(α)) / sqrt(cos(α))^3
    numer(α) = α^2 * exp(b*cos(α)) / sqrt(cos(α))^3
    return quadgk(numer,0,π/2)[1] / quadgk(denom,0,π/2)[1] / π^2
end

function analytic_susc_U2_approx_2(β)
    b = big(β)
    denom(α)   = exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)))
    numer(α) = α^2 * exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)))
    return quadgk(numer,0,π/2)[1] / quadgk(denom,0,π/2)[1] / π^2
end

function analytic_susc_U2_approx_3(β)
    b = big(β)
    denom(α)   = exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)) - 15/(128*(b*cos(α))^2))
    numer(α) = α^2 * exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)) - 15/(128*(b*cos(α))^2))
    return quadgk(numer,0,π/2)[1] / quadgk(denom,0,π/2)[1] / π^2
end

function insta_U2_z(N_x, N_t, q, z)
    # w = π*(2*z-q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*q*t*π/N_x/N_t) * exp_u2(-t*2*π/N_x/N_t * (z-q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*q*x*π/N_x) * exp_u2(x*2*π/N_x * (z-q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x]
    return U
end

function small_U2_single_lie_direction(lie_dir::Int, ϵ)
    ### "ϵ" is short for "hitsize"
    if lie_dir == 0
        return coeffs_U2(exp(im*ϵ/2), complex(0.0), complex(0.0), complex(0.0))
    elseif lie_dir == 1
        return coeffs_U2(complex(cos(ϵ/2)), complex(sin(ϵ/2)), complex(0.0),      complex(0.0))
    elseif lie_dir == 2
        return coeffs_U2(complex(cos(ϵ/2)), complex(0.0),      complex(sin(ϵ/2)), complex(0.0))
    elseif lie_dir == 3
        return coeffs_U2(complex(cos(ϵ/2)), complex(0.0),      complex(0.0),      complex(sin(ϵ/2)))
    else
        error("Group directions 'lie_dir' are labeled from 0 to 3")
    end
end

function small_U2_lie_direction_vec(hitsize, v::Vector)
    w = hitsize * sqrt(1/(v'*v)) .* v
    return grp2coeffs_U2(exp(im*sum(w.*Σ)/2))
end

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

function insta_U2_z_comb(N_x, N_t, Q, z)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,1:N_x-1,:] = [exp(-im*Q*π*(t-1)/N_x/N_t) * exp_u2(-(t-1)*2*π/N_x/N_t * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x-1, t = 1:N_t]
    U[1,N_x,:]     = [exp(-im*Q*π*(t-1+N_x)/N_x/N_t) * exp_u2(-(t-1+N_x)*2*π/N_x/N_t * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * exp_u2(x*2*π/N_x * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x]
    return U
end

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
    return sqrt(sum([two_metric_squared(U[μ,x,t], V[μ,x,t]) for μ = 1:2, x = 1:NX, t = 1:NT ])) / (2*NX*NT)
end

function find_rot_mat(a,b)
    v = cross(a,b)
    s = sqrt(v'*v)
    c = a'*b
    M_v = [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]
    return I(3) + M_v + (1-c)/s^2 * M_v^2
end

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

function Pauli2vec(M::Matrix)
    v1 = imag(M[1,2])
    v2 = real(M[1,2])
    v3 = imag(M[1,1])
    return [v1,v2,v3]
end

function Pauli2vec(X::coeffs_U2)
    return [X.b, X.c, X.d]
end

function two_metric_field_insta_rot(U, M_rot_Lie_coeffs)
    NX = size(U,2)
    NT = size(U,3)
    q = round(Int,top_charge_U2(U))
    L1 = [0 0 0; 0 0 -1; 0 1 0]
    L2 = [0 0 1; 0 0 0; -1 0 0]
    L3 = [0 -1 0; 1 0 0; 0 0 0]
    M_rot = exp(sum(M_rot_Lie_coeffs.*[L1,L2,L3]))
    insta = insta_U2_comb_odd(NX,NT,q,M_rot)
    return sqrt(sum([two_metric_squared(U[μ,x,t], insta[μ,x,t]) for μ = 1:2, x = 1:NX, t = 1:NT ])) / (2*NX*NT)
end

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
    else
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
    end
    return two_metric_field(V,insta_U2_comb(NX,NT,q))
end

function optimize_insta_metric(U, start_coeffs)
    @assert length(start_coeffs) == 4 "start_coeffs must number 4, the last 3 of which are the vector"
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
        return optimize(optim_metric_even, start_coeffs,NelderMead()).minimum
    else
        optim_metric_odd(coeffs) = two_metric_field_insta_rot(V, coeffs)
        return optimize(optim_metric_odd, start_coeffs[2:4],NelderMead()).minimum
    end
end

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





########  Let's disturb some more locals  ########





L        = 32
q        = 1
z        = 5
N_smear  = Int(1e4)
ρ        = 0.1
pertsize = 1e-1

dist_actions = Array{Float64}(undef,N_smear+1,4)
dist_charges = Array{Float64}(undef,N_smear+1,4)
U = insta_U2_z(L, L, q, z)
μ = rand(1:2)
x = rand(1:L)
t = rand(1:L)

#=
for lie_dir = 0:3
    count = 0
    println("    Perturbation Lie direction: $lie_dir")
    V = deepcopy(U)
    V[μ,x,t] = small_U2_single_lie_direction(lie_dir, pertsize) * V[μ,x,t]
    dist_actions[1,lie_dir+1] = action(V,1)
    dist_charges[1,lie_dir+1] = top_charge_U2(V)
    for smear = 1:N_smear
        V = stout_midpoint(V,ρ)
        dist_actions[smear+1,lie_dir+1] = action(V,1)
        dist_charges[smear+1,lie_dir+1] = top_charge_U2(V)
        if smear%(N_smear/20) == 0
            count += 5
            println("        Smearing progress: $count%")
        end
    end
end
=#

dist_actions_path = string(data_path,"\\dist_actions\\dist_actions_q_$(q)_z_$(z)_pertsize_$pertsize.txt")
# # writedlm(dist_actions_path, dist_actions)
dist_actions = readdlm(dist_actions_path)

let
    # l_styles = [:solid, :dash, :dashdot, :dashdotdot]
    l_styles = [:solid, :solid, :dash, :dash]
    l_widths = [2, 2, 1, 1]
    alphas   = [0.7, 0.7, 1, 1]
    τ0 = 2
    start_ind = Int(τ0/ρ)
    image_dist = plot(
        # title  = latexstring("\$L = $L, q = $q, z = $z, \\mathrm{pert.\\ size} = $pertsize\$"),
        xlabel = latexstring("flow time \$\\tau/a^2\$"),
        ylabel = latexstring("\$S/\\beta\$"),
        tickfontsize = 11,
        labelfontsize = 16,
        legendfontsize = 11,
        legend = :top,
        # foreground_color_legend = :false,
        background_color_legend = RGBA(1.0,1.0,1.0,0.8),
        # left_margin = 2mm,
    )
    # greys = (:grey84, :grey77, :grey70, :grey63, :grey56, :grey49, :grey42)
    greys = (:grey77, :grey70, :grey63, :grey56, :grey49, :grey42, :grey30)
    for lie_dir = 0:3
        image_dist = plot!(
            ρ.*Array(0:N_smear)[start_ind:end],
            dist_actions[start_ind:end,lie_dir+1],
            label = latexstring("\$j = $lie_dir\$"),
            # color = greys[2*lie_dir+1],
            color = cb_colors[lie_dir+1],
            # markershape = m_shapes[lie_dir+1],
            # markersize = 2,
            linewidth = l_widths[lie_dir+1], # 1.2,
            alpha = alphas[lie_dir+1],
            linestyle = l_styles[lie_dir+1]
        )
    end
    for z_prime = 1:z-1
        j = z-z_prime
        image_dist = hline!(
            [insta_action_U2(1,L,L,q,j)],
            label = latexstring("\$S\$ for \$q=$q, z=$j\$"),
            color = greys[j+length(greys)-z+1],
            linestyle = :dot,
            legend = :right,
            linewidth = 1.5,
            alpha = 0.8,
        )
    end
    image_dist = hline!(
        [insta_action_U2_min(1,L,L,q)],
        label = latexstring("lower bound for \$q=1\$"),
        color = :black, # cb_red,
        linestyle = :dot,
        legend = :right,
        linewidth = 2.5,
        alpha = 0.8,
    )
    # image_dist = plot!(ylim = (0.0, 0.02))
    image_dist = plot!(yaxis = :log)
    display(image_dist)
end # let
# end # z
# end # q

image_dist_path = string(fig_path,"\\dist_actions\\dist_actions_q_$(q)_z_$(z)_pertsize_$pertsize.pdf")
# savefig(image_dist_path)



####################################################################################################################################



L        = 32
q        = 5
z        = 3
N_smear  = Int(1e4)
ρ        = 0.1
pertsize = 1e-1

dist_actions = Array{Float64}(undef,N_smear+1)
dist_charges = Array{Float64}(undef,N_smear+1)
U = insta_U2_z(L, L, q, z)
μ = rand(1:2)
x = rand(1:L)
t = rand(1:L)
# v = 2 .* rand(4) .- 0.5
v = [1,0,0,-1]

let
    count = 0
    println("Perturbation Lie direction: $(round.(v, digits = 2))")
    V = deepcopy(U)
    V[μ,x,t] = small_U2_lie_direction_vec(pertsize, v) * V[μ,x,t]
    dist_actions[1] = action(V,1)
    dist_charges[1] = top_charge_U2(V)
    for smear = 1:N_smear
        V = stout_midpoint(V,ρ)
        dist_actions[smear+1] = action(V,1)
        dist_charges[smear+1] = top_charge_U2(V)
        if smear%(N_smear/20) == 0
            count += 5
            println("    Smearing progress: $count%")
        end
    end
end



# dist_actions_path = string(data_path,"\\dist_actions_VEC_q_$(q)_z_$(z)_pertsize_$(pertsize).txt")
# # writedlm(dist_actions_path, dist_actions)
# dist_actions = readdlm(dist_actions_path)



let
    τ0 = 2
    start_ind = Int(τ0/ρ)
    image_dist = plot(
        title  = latexstring("\$L = $L, q = $q, z = $z, \\mathrm{pert.\\, size} = $pertsize\$\n Lie dir.: \$$(round.(v,digits = 2))\$"),
        xlabel = latexstring("flow time \$\\tau/a^2\$"),
        ylabel = latexstring("\$S/\\beta\$"),
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 12,
        legend = :top,
        # foreground_color_legend = :false,
        background_color_legend = RGBA(1.0,1.0,1.0,0.8),
        # left_margin = 2mm,
    )
    # greys = (:grey84, :grey77, :grey70, :grey63, :grey56, :grey49, :grey42)
    greys = (:grey77, :grey70, :grey63, :grey56, :grey49, :grey42, :grey30)
    image_dist = plot!(
        ρ.*Array(0:N_smear)[start_ind:end],
        dist_actions[start_ind:end],
        label = "Measurements",
        color = cb_colors[1],
        linewidth = 2,
        alpha = 1,
    )
    image_dist = hline!(
        [insta_action_U2_min(1,L,L,q)],
        label = latexstring("lower bound for \$q=1\$"),
        color = cb_red,
        linestyle = :dot,
        legend = :right,
        linewidth = 2.5,
        alpha = 0.8,
    )
    # for z_prime = 1:z
    #     j = z-z_prime+1
    #     image_dist = hline!(
    #         [insta_action_U2(1,L,L,q,j)],
    #         label = latexstring("\$S\$ for \$q=1, z=$j\$"),
    #         color = greys[j],
    #         linestyle = :dot,
    #         legend = :right,
    #         linewidth = 2.5,
    #         alpha = 0.8,
    #     )
    # end
    # image_dist = plot!(ylim = (0.0, 0.02))
    display(image_dist)
end


image_dist_path = string(fig_path,"\\disturbed_locals_SINGLE_LINK.pdf")
# savefig(image_dist_path)





#### special plots ####

L        = 32
q        = 1
z        = 5
N_smear  = Int(1e4)
ρ        = 0.1
pertsize = 1e-1
dist_actions_path = string(data_path,"\\dist_actions\\dist_actions_q_$(q)_z_$(z)_pertsize_$pertsize.txt")
dist_actions = readdlm(dist_actions_path)

special_colors = [cb_pink, cb_orange, cb_purple, cb_blue]
let
    # l_styles = [:solid, :dash, :dashdot, :dashdotdot]
    l_styles = [:solid, :solid, :dash, :dash]
    l_widths = [4.5, 4.5, 1.7, 1.5]
    color_inds = [1,2,3,4]
    alphas   = [0.6, 0.45, 1, 1]
    τ0 = 2
    start_ind = Int(τ0/ρ)
    image_dist = plot(
        # title  = latexstring("\$L = $L, q = $q, z = $z, \\mathrm{pert.\\ size} = $pertsize\$"),
        xlabel = latexstring("flow time \$\\tau/a^2\$"),
        ylabel = latexstring("\$S/\\beta\$"),
        tickfontsize = 11,
        labelfontsize = 16,
        legendfontsize = 11,
        legend = :right,
        # foreground_color_legend = :false,
        background_color_legend = RGBA(1.0,1.0,1.0,0.8),
        # left_margin = 2mm,
    )
    greys = (:grey84, :grey77, :grey70, :grey63, :grey56, :grey49, :grey42, :grey35)
    # greys = (:grey77, :grey70, :grey63, :grey56, :grey49, :grey42, :grey30)
    for z_prime = 1:z-1
        j = z-z_prime
        image_dist = hline!(
            [insta_action_U2(1,L,L,q,j)],
            label = :false,# latexstring("\$S\$ for \$q=$q, z=$j\$"),
            color = greys[2*j],
            # color = greys[j+length(greys)-z+1],
            linestyle = :dot,
            linewidth = 2.5,
            alpha = 0.8,
        )
    end
    for lie_dir = 0:3
        image_dist = plot!(
            ρ.*Array(0:N_smear)[start_ind:end],
            dist_actions[start_ind:end,lie_dir+1],
            label = latexstring("\$j = $lie_dir\$"),
            # color = greys[2*lie_dir+1],
            color = special_colors[color_inds[lie_dir+1]],
            # markershape = m_shapes[lie_dir+1],
            # markersize = 2,
            linewidth = l_widths[lie_dir+1], # 1.2,
            alpha = alphas[lie_dir+1],
            linestyle = l_styles[lie_dir+1]
        )
    end

    annotate!(300-55, 0.28, text(latexstring("\$(q=1, z=4)\$"), :black, 11))
    annotate!(355-55, 0.145, text(latexstring("\$(q=1, z=3)\$"), :black, 11))
    annotate!(370-55, 0.0555, text(latexstring("\$(q=1, z=2)\$"), :black, 11))
    annotate!(415-55, 0.0112, text(latexstring("\$(q=1, z=1)\$"), :black, 11))
    image_dist = hline!(
        [insta_action_U2_min(1,L,L,q)],
        label = latexstring("lower bound for \$q=1\$"),
        color = :black, # cb_red,
        linestyle = :dot,
        linewidth = 2.5,
        alpha = 0.8,
    )
    # image_dist = plot!(ylim = (0.0, 0.02))
    image_dist = plot!(yaxis = :log, legend = :right)
    display(image_dist)
end # let

image_dist_path = string(fig_path,"\\dist_actions\\dist_actions_q_$(q)_z_$(z)_pertsize_$pertsize.pdf")
# savefig(image_dist_path)





########  *How* special are our special configs?  ########





#=
function small_U2_single_lie_direction(lie_dir, ϵ)
    v = [0.0, 0.0, 0.0]
    if lie_dir == 0
        return coeffs_U2(exp(im*ϵ/2), complex(0.0), complex(0.0), complex(0.0))
    elseif lie_dir == 1
        return coeffs_U2(complex(cos(ϵ/2)), complex(sin(ϵ/2)), complex(0.0),    complex(0.0))
    elseif lie_dir == 2
        return coeffs_U2(complex(cos(ϵ/2)), complex(0.0),    complex(sin(ϵ/2)), complex(0.0))
    elseif lie_dir == 3
        return coeffs_U2(complex(cos(ϵ/2)), complex(0.0),    complex(0.0),    complex(sin(ϵ/2)))
    else
        error("Group directions 'lie_dir' are labeled from 0 to 3")
    end
end

# abs(det(small_U2_single_lie_direction(rand(0:3),rand())))
# bla = small_U2_single_lie_direction(rand(0:3),rand())
# coeffs2grp(bla * adjoint(bla))
# for r = 0:3
    # reps = rand()
    # @assert isapprox(coeffs2grp(small_U2_single_lie_direction(r,reps)), exp(im*reps*Σ[r+1]))
# end

function vary_config_U2(U, ϵ)
    NX = size(U,2)
    NT = size(U,3)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                stap_dag = staple_dag(U, μ, x, t)
                old_link = U[μ,x,t]
                s_old = -real(tr(old_link*stap_dag))
                for lie_dir in [0,1,2,3]
                    new_link_p = small_U2_single_lie_direction(lie_dir,+ϵ) * old_link
                    new_link_m = small_U2_single_lie_direction(lie_dir,-ϵ) * old_link
                    delta_s_p = -real(tr(new_link_p*stap_dag)) - s_old
                    delta_s_m = -real(tr(new_link_m*stap_dag)) - s_old
                    if delta_s_p < 0
                        error("Action decrease in μ = $μ, x = $x, t = $t, ϵ = $ϵ in pos. group direction $lie_dir")
                    elseif delta_s_m < 0
                        error("Action decrease in μ = $μ, x = $x, t = $t, ϵ = $ϵ in neg. group direction $lie_dir")
                    end
                end
            end
        end
    end
    return nothing
end

function vary_config_ran_U2(U, ϵ, N_prop)
    NX = size(U,2)
    NT = size(U,3)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                stap_dag = staple_dag(U, μ, x, t)
                old_link = U[μ,x,t]
                s_old = -real(tr(old_link*stap_dag))
                for N = 1:N_prop
                    # r = ran_U2(ϵ*(2*rand()-1))
                    r = ran_U2(ϵ)
                    new_link = r * old_link
                    s_new = -real(tr(new_link*stap_dag))
                    if s_new-s_old < 0 
                        error("Action decrease in μ = $μ, x = $x, t = $t, ϵ = $ϵ for random element $(coeffs2grp(r))")
                    end
                end
            end
        end
    end
end


N_x = N_t = 32
for q = -10:10
    for z = -10:10
        insta = insta_U2_z(N_x, N_t, q, z)
        for pow = 8:8
            vary_config_U2(insta,10.0^(-pow))
            # vary_config_ran_U2(insta,10.0^(-pow),10)
        end
    end
end

let
    N_x = N_t = 32
    q = 2
    z = 1
    ϵ = 1e-7
    insta = insta_U2_z(N_x, N_t, q, z)
    vary_config_U2(insta,ϵ)
end

# bla = small_U2_single_lie_direction(rand(0:3),1e-3); coeffs2grp(bla)[1,1]

function vary_special_config_U2(tol, L, q_vals, z_vals, ϵ_powers, file_name)
    mac_prec = eps()
    for q in q_vals
    for z in z_vals; insta = insta_U2_z(L, L, q, z)
    for pow in ϵ_powers; ϵ = 10.0^(-pow)
        for t = 1:L
        for x = 1:L
        for μ = 1:2
            stap_dag = staple_dag(insta, μ, x, t)
            old_link = insta[μ,x,t]
            s_old    = -real(tr(old_link*stap_dag))
            for lie_dir in [0,1,2,3]
                new_link_p = small_U2_single_lie_direction(lie_dir,+ϵ) * old_link
                new_link_m = small_U2_single_lie_direction(lie_dir,-ϵ) * old_link
                delta_s_p = -real(tr(new_link_p*stap_dag)) - s_old
                delta_s_m = -real(tr(new_link_m*stap_dag)) - s_old
                if delta_s_p < tol * mac_prec
                    open(file_name, "a") do io
                        write(io, "q = $q, z = $z, μ = $μ, x = $x, t = $t, pow = -$pow, lie_dir = +$lie_dir, ΔS/eps = $(delta_s_p/mac_prec)\n")
                    end
                end
                if delta_s_p < tol * mac_prec
                    open(file_name, "a") do io
                        write(io, "q = $q, z = $z, μ = $μ, x = $x, t = $t, pow = -$pow, lie_dir = -$lie_dir, ΔS/eps = $(delta_s_m/mac_prec) \n")
                    end
                end
            end
        end # t
        end # x
        end # μ
    end # q
    end # z
    end # pow
    return nothing
end

for hitsize_pow = 1:7
    # hitsize = 11
    vary_file_name = string(data_path, "\\vary_hitsize_$hitsize_pow.txt")
    vary_special_config_U2(10, 32, 0:5, 0:5, hitsize_pow:hitsize_pow, vary_file_name)
end
=#





### Take a "special config" of size L², q, z. The action difference under
### variation of the link at [μ,x,t] along the Lie-direction lie_dir is measured for
### various hitsizes contained in hitsize_range.
function sdiff_min_spec_conf(L, q, z, μ, x, t, hitsize_range, lie_dir)
    insta    = insta_U2_z(L, L, q, z)
    old_link = insta[μ,x,t]
    stap_dag = staple_dag(insta,μ,x,t)
    s_old    = -real(tr(old_link*stap_dag))
    s_diff  = []
    for hitsize in hitsize_range
        new_link = small_U2_single_lie_direction(lie_dir, hitsize) * old_link
        s_new    = -real(tr(new_link*stap_dag))
        push!(s_diff, s_new-s_old)
    end
    return s_diff
end

last_fit_converged = true
let
    L = 32
    resol = 1e-7
    magnif = 1/resol
    model_sq(x, p) = p[1] .* (magnif .* (x .- p[2])) .^2 .+ p[3]
    p_sdiff = [1.0, 0.0, 0.0]
    hitsize_range  = Vector(-3.5*resol:resol:3.5*resol)
    fit_plot_xvals = Vector(-4.0*resol:resol/10:4.0*resol)
    for q = 1:1
    for z = 5:5
    μ_vals = rand(1:2, 5)
    x_vals = [rand(1:L-1), rand(1:L-1), L,           L, rand(1:L)]
    t_vals = [rand(1:L-1), L,           rand(1:L-1), L, rand(1:L)]
    for link_ind = 5:5
        μ = μ_vals[link_ind]
        x = x_vals[link_ind]
        t = t_vals[link_ind]
        for lie_dir = 0:3
            sdiff     = sdiff_min_spec_conf(L, q, z, μ, x, t, hitsize_range, lie_dir) ./ eps()
            fit_sdiff = curve_fit(model_sq, hitsize_range, sdiff, p_sdiff)
            a = round.(fit_sdiff.param[1]*magnif^2, sigdigits = 3)
            b = round.(fit_sdiff.param[2], sigdigits = 3)
            c = round.(fit_sdiff.param[3], sigdigits = 3)
            last_fit_converged = fit_sdiff.converged
            image_sdiff = plot(
                # title = latexstring("\$ \\Delta S \$ of the special config at \$(q,z) = ($q,$z)\$ \n \$(\\mu, x, t) = ($μ,$x,$t)\$, Lie-dir.: \$ $lie_dir\$"),
                xlabel = latexstring("step size \$\\epsilon_1\$"),
                ylabel = latexstring("\$\\Delta S / \\left(\\beta \\cdot \\epsilon_\\textrm{machine} \\right)\$"),
                rightmargin = 5mm,
                labelfontsize = 15,
                tickfontsize = 10,
                legendfontsize = 10,
                legend = :top
            )
            image_sdiff = plot!(
                fit_plot_xvals,
                model_sq(fit_plot_xvals, fit_sdiff.param),
                label = latexstring("\$\\mathrm{fit}(x) = a\\cdot (x-b)^2 + c\$ \n \$a = $a\$ \n \$b = $b\$ \n \$c = $c\$"),
                color = cb_orange
            )
            image_sdiff = scatter!(
                hitsize_range,
                sdiff,
                label = "Measurements",
                color = cb_blue
            )
            display(image_sdiff)
        end # lie_dir
    end # link_ind
    end # z
    end # q
end # let

image_sdiff_liedir_path = string(fig_path, "\\sdiff_q_1_z_5_lie_dir_3.pdf")
# savefig(image_sdiff_liedir_path)





### Take a "special config" of size L², q, z. The action difference under
### variation of the link at [μ,x,t] along the Lie-direction lie_dir is measured for
### various hitsizes contained in hitsize_range.
function sdiff_min_spec_conf_vec(L, q, z, μ, x, t, hitsize_range, v)
    insta = insta_U2_z(L, L, q, z)
    old_link = insta[μ,x,t]
    stap_dag = staple_dag(insta,μ,x,t)
    s_old    = -real(tr(old_link*stap_dag))
    s_diff  = []
    for hitsize in hitsize_range
        new_link = small_U2_lie_direction_vec(hitsize, v) * old_link
        s_new    = -real(tr(new_link*stap_dag))
        push!(s_diff, s_new-s_old)
    end
    return s_diff
end

last_fit_converged = true
let
    L = 32
    resol = 1e-7
    magnif = 1/resol
    model_sq(x, p) = p[1] .* (magnif .* (x .- p[2])) .^2 .+ p[3]
    p_sdiff = [1.0, 0.0, 0.0]
    hitsize_range  = Vector(-3.5*resol:resol:3.5*resol)
    fit_plot_xvals = Vector(-4.0*resol:resol/10:4.0*resol)
    for q = 1:1
    for z = 5:5
        μ_vals = rand(1:2, 5)
        x_vals = [rand(1:L-1), rand(1:L-1), L,           L, rand(1:L)]
        t_vals = [rand(1:L-1), L,           rand(1:L-1), L, rand(1:L)]
        for link_ind = 5:5
            μ = μ_vals[link_ind]
            x = x_vals[link_ind]
            t = t_vals[link_ind]
            for num_pert = 1:1
                v = 2 .* rand(4) .- 1
                println("v = $(round.(v,digits=3))")
                # v = [rand().-0.5, 0.0, 0.0, rand().-0.5]
                # v = [0.0, rand().-0.5, rand().-0.5, 0.0]
                sdiff     = sdiff_min_spec_conf_vec(L, q, z, μ, x, t, hitsize_range, v) ./ eps()
                fit_sdiff = curve_fit(model_sq, hitsize_range, sdiff, p_sdiff)
                a = round.(fit_sdiff.param[1]*magnif^2, sigdigits = 3)
                b = round.(fit_sdiff.param[2], sigdigits = 3)
                c = round.(fit_sdiff.param[3], sigdigits = 3)
                last_fit_converged = fit_sdiff.converged
                image_sdiff = plot(
                    # title = latexstring("\$ \\Delta S \$ of the special config at \$(q,z) = ($q,$z)\$ \n \$(\\mu, x, t) = ($μ,$x,$t)\$ \n Lie vec.: \$ $(round.(v,digits=2))\$"),
                    xlabel = latexstring("step size \$\\epsilon\$"),
                    ylabel = latexstring("\$\\Delta S / \\left(\\beta \\cdot \\epsilon_\\textrm{machine} \\right)\$"),
                    rightmargin = 5mm,
                    labelfontsize = 15,
                    tickfontsize = 10,
                    legendfontsize = 10,
                    # legend = :inside,
                    legend = :top
                )
                image_sdiff = plot!(
                    fit_plot_xvals,
                    model_sq(fit_plot_xvals, fit_sdiff.param),
                    label = latexstring("\$\\mathrm{fit}(x) = a\\cdot (x-b)^2 + c \$ \n \$a = $a\$ \n \$b = $b\$ \n \$c = $c\$"),
                    color = cb_orange
                )
                image_sdiff = scatter!(
                    hitsize_range,
                    sdiff,
                    label = "Measurements",
                    color = cb_blue
                )
                display(image_sdiff)
            end # num_pert
        end # link_ind
    end # z
    end # q
end # let


image_sdiff_lievec_path = string(fig_path, "\\sdiff_q_1_z_5_lie_vec.pdf")
# savefig(image_sdiff_lievec_path)




### Take a "special config" of size L², q, z and perturb a link at [μ,x,t] with hitsize
### pert_size in the Lie-direction pert_dir. The action difference under
### variation of the link at [μ,x,t] along the Lie-direction lie_dir is measured for
### various hitsizes contained in hitsize_range.
function sdiff_min_pert_spec_conf(L, q, z, μ, x, t, hitsize_range, lie_dir, pert_dir, pert_size)
    insta = insta_U2_z(L, L, q, z)
    insta[μ,x,t] = small_U2_single_lie_direction(pert_dir, pert_size) * insta[μ,x,t]
    old_link = insta[μ,x,t]
    stap_dag = staple_dag(insta,μ,x,t)
    s_old    = -real(tr(old_link*stap_dag))
    s_diff  = []
    for hitsize in hitsize_range
        new_link = small_U2_single_lie_direction(lie_dir, hitsize) * old_link
        s_new    = -real(tr(new_link*stap_dag))
        push!(s_diff, s_new-s_old)
    end
    return s_diff
end

last_fit_converged = true
let
    L = 32
    pert_size = 1e-3
    resol = 1e-7
    magnif = 1/resol
    model_sq_pert(x, p) = p[1] .* (magnif .* (x .- p[2])) .^2 .+ p[3]
    p_sdiff_pert = [1.0, 1e-8, -1e-2]
    hitsize_range  = Vector(-3.5*resol:resol:3.5*resol)
    fit_plot_xvals = Vector(-4.0*resol:resol/10:4.0*resol)
    for q = 1:1
    for z = 5:5
        μ_vals = rand(1:2, 5)
        x_vals = [rand(1:L-1), rand(1:L-1), L,           L, rand(1:L)]
        t_vals = [rand(1:L-1), L,           rand(1:L-1), L, rand(1:L)]
        for link_ind = 5:5
            μ = μ_vals[link_ind]
            x = x_vals[link_ind]
            t = t_vals[link_ind]
            for pert_dir = 3:3
                lie_directions = [0,1,2,3]
                popat!(lie_directions, pert_dir+1)
                for lie_dir = 0:0
                    sdiff     = sdiff_min_pert_spec_conf(L, q, z, μ, x, t, hitsize_range, lie_dir, pert_dir, pert_size) ./ eps()
                    fit_sdiff = curve_fit(model_sq_pert, hitsize_range, sdiff, p_sdiff_pert)
                    a = round.(fit_sdiff.param[1]*magnif^2, sigdigits = 3)
                    b = round.(fit_sdiff.param[2], sigdigits = 3)
                    c = round.(fit_sdiff.param[3], sigdigits = 3)
                    last_fit_converged = fit_sdiff.converged
                    image_sdiff = plot(
                        # title = latexstring("\$ \\Delta S \$ of the special config at \$(q,z) = ($q,$z)\$ \n \$(\\mu, x, t) = ($μ,$x,$t)\$, pert. size: \$ $pert_size\$ \n pert. dir.: \$ $pert_dir\$, Lie-dir.: \$ $lie_dir\$"),
                        xlabel = latexstring("step size \$\\epsilon_2\$"),
                        ylabel = latexstring("\$\\Delta S / \\left(\\beta \\cdot \\epsilon_\\textrm{machine} \\right)\$"),
                        rightmargin = 5mm,
                        labelfontsize = 15,
                        tickfontsize = 10,
                        legendfontsize = 10,
                        legend = :top
                    )
                    image_sdiff = plot!(
                        fit_plot_xvals,
                        model_sq_pert(fit_plot_xvals, fit_sdiff.param),
                        label = latexstring("\$\\mathrm{fit}(x) = a\\cdot (x-b)^2 + c\$ \n \$a = $a\$ \n \$b = $b\$ \n \$c = $c\$"),
                        color = cb_orange
                    )
                    image_sdiff = scatter!(
                        hitsize_range,
                        sdiff,
                        label = "Measurements",
                        color = cb_blue
                    )
                    display(image_sdiff)
                end # lie_dir
            end # pert_dir
        end # link_ind
    end # z
    end # q
end # let

image_sdiff_lievec_path = string(fig_path, "\\sdiff_q_1_z_5_second_pert.pdf")
# savefig(image_sdiff_lievec_path)

### Observations:
### pert. size 1e-3, z = 0:3 
### indepently of (μ,x,t)
### q = 0
###     parabolas with b, c not sign. different from zero
### q = 1
###     same and:
###     for pert_dir, lie_dir ∈ {0,3}
###     z = 0: b ~ -8.9e-9
###     z = 1: b ~ +9.3e-9
###     z = 2: b ~ +2.8e-8
###     z = 3: b ~ +4.7e-8
### q = 2
###     for pert_dir, lie_dir ∈ {0,3}
###     z = 0: b ~ -3.8e-8
###     z = 1: b ~ +/-0.0   (O(1e-23), embedding of instanton config)
###     z = 2: b ~ +3.8e-8
###     z = 3: b ~ +7.5e-8
### q = 3
###     for pert_dir, lie_dir ∈ {0,3}
###     z = 0: b ~ -8.5e-8
###     z = 1: b ~ -2.8e-8
###     z = 2: b ~ +2.8e-8
###     z = 3: b ~ +8.5e-8
### q = 4
###     for pert_dir, lie_dir ∈ {0,3}
###     z = 0: b ~ -1.5e-7
###     z = 1: b ~ -7.5e-8
###     z = 2: b ~ +/-0.0   (O(1e-10), embedding of instanton config)
###     z = 3: b ~ +7.5e-8
### q = 5
###     for pert_dir, lie_dir ∈ {0,3}
###     z = 0: b ~ -2.4e-7
###     z = 1: b ~ -1.4e-7
###     z = 2: b ~ -4.6e-8
###     z = 3: b ~ +4.6e-8





### Take a "special config" of size L², q, z and perturb a link at [μ,x,t] with hitsize
### pert_size in the Lie-direction given by the VECTOR v. The action difference under
### variation of the link at [μ,x,t] along the Lie-direction lie_dir is measured for
### various hitsizes contained in hitsize_range.
function sdiff_min_pert_spec_conf_vec(L, q, z, μ, x, t, hitsize_range, lie_dir, v, pert_size)
    insta = insta_U2_z(L, L, q, z)
    insta[μ,x,t] = small_U2_lie_direction_vec(pert_size, v) * insta[μ,x,t]
    old_link = insta[μ,x,t]
    stap_dag = staple_dag(insta,μ,x,t)
    s_old    = -real(tr(old_link*stap_dag))
    s_diff  = []
    for hitsize in hitsize_range
        new_link = small_U2_single_lie_direction(lie_dir, hitsize) * old_link
        s_new    = -real(tr(new_link*stap_dag))
        push!(s_diff, s_new-s_old)
    end
    return s_diff
end

last_fit_converged = true
let
    L = 32
    pert_size = 1e-3
    resol = 1e-7
    magnif = 1/resol
    model_sq_pert(x, p) = p[1] .* (magnif .* (x .- p[2])) .^2 .+ p[3]
    p_sdiff_pert = [1.0, 1e-8, -1e-2]
    hitsize_range  = Vector(-3.5*resol:resol:3.5*resol)
    fit_plot_xvals = Vector(-4.0*resol:resol/10:4.0*resol)
    for q = 1:1
    for z = 5:5
        μ_vals = rand(1:2, 5)
        x_vals = [rand(1:L-1), rand(1:L-1), L,           L, rand(1:L)]
        t_vals = [rand(1:L-1), L,           rand(1:L-1), L, rand(1:L)]
        for link_ind = 1:4
            μ = μ_vals[link_ind]
            x = x_vals[link_ind]
            t = t_vals[link_ind]
            for num_pert = 1:1
                w = 2 .* rand(4) .- 1
                for lie_dir = 0:3
                    v = copy(w)
                    v[lie_dir+1] = 0.0
                    sdiff     = sdiff_min_pert_spec_conf_vec(L, q, z, μ, x, t, hitsize_range, lie_dir, v, pert_size) ./ eps()
                    fit_sdiff = curve_fit(model_sq_pert, hitsize_range, sdiff, p_sdiff_pert)
                    a = round.(fit_sdiff.param[1]*magnif^2, sigdigits = 3)
                    b = round.(fit_sdiff.param[2], sigdigits = 3)
                    c = round.(fit_sdiff.param[3], sigdigits = 3)
                    last_fit_converged = fit_sdiff.converged
                    image_sdiff = plot(
                        title = latexstring("\$ \\Delta S \$ of the special config at \$(q,z) = ($q,$z)\$ \n \$(\\mu, x, t) = ($μ,$x,$t)\$, pert. size: \$ $pert_size\$ \n pert. vec.: \$ $(round.(v,digits=2))\$ Lie-dir.: \$ $lie_dir\$"),
                        xlabel = "hitsize",
                        ylabel = latexstring("\$\\Delta S / \\left(\\beta \\cdot \\epsilon_\\textrm{machine} \\right)\$"),
                        rightmargin = 5mm,
                        labelfontsize = 15,
                        tickfontsize = 10,
                        legend = :inside
                    )
                    image_sdiff = plot!(
                        fit_plot_xvals,
                        model_sq_pert(fit_plot_xvals, fit_sdiff.param),
                        label = latexstring("\$\\mathrm{fit}(x) = a\\cdot (x-b)^2 + c\$ \n \$a = $a\$ \n \$b = $b\$ \n \$c = $c\$"),
                        color = cb_orange
                    )
                    image_sdiff = scatter!(
                        hitsize_range,
                        sdiff,
                        label = "Measurements",
                        color = cb_blue
                    )
                    display(image_sdiff)
                end # lie_dir
            end # num_pert
        end # link_ind
    end # z
    end # q
end # let

### Observations:
### pert. size 1e-3, z = 0:3 
### indepently of (μ,x,t)
### Parameters b,c significantly different from 0 only for Lie direction ∈ {0,3} 
###     But in those cases b,c significantly larger than before, when pert. was
###     not carried out in linear combination of Lie directions





### Take a "special config" of size L², q, z and perturb a link which is PART OF THE PLAQUETTE
### with hitsize pert_size in the Lie-direction pert_dir. The action difference under
### variation of the link at [μ,x,t] along the Lie-direction lie_dir is measured for
### various hitsizes contained in hitsize_range.
function sdiff_min_pert_spec_conf_vec_neigh(L, q, z, μ, x, t, hitsize_range, lie_dir, pert_dir, pert_size)
    insta = insta_U2_z(L, L, q, z)
    insta[μ,x,t] = small_U2_single_lie_direction(pert_dir, pert_size) * insta[μ,x,t]
    ### Now the old link is the one on the other side of one of the plaquettes
    old_link = coeffs_Id_U2()
    stap_dag = coeffs_Id_U2()
    if μ == 1
        old_link = insta[μ,x,mod1(t+1,L)]
        stap_dag = staple_dag(insta,μ,x,mod1(t+1,L))
    elseif μ == 2
        old_link = insta[μ,mod1(x+1,L),t]
        stap_dag = staple_dag(insta,μ,mod1(x+1,L),t)
    end
    s_old    = -real(tr(old_link*stap_dag))
    s_diff  = []
    for hitsize in hitsize_range
        new_link = small_U2_single_lie_direction(lie_dir, hitsize) * old_link
        s_new    = -real(tr(new_link*stap_dag))
        push!(s_diff, s_new-s_old)
    end
    return s_diff
end

last_fit_converged = true
let
    L = 32
    pert_size = 1e-3
    resol = 1e-4
    magnif = 1/resol
    model_sq_pert(x, p) = p[1] .* (magnif .* (x .- p[2])) .^2 .+ p[3]
    p_sdiff_pert = [1.0, 1e-8, -1e-2]
    hitsize_range  = Vector(-3.5*resol:resol:3.5*resol)
    fit_plot_xvals = Vector(-4.0*resol:resol/10:4.0*resol)
    for q = 1:1
    for z = 5:5
        μ_vals = rand(1:2, 5)
        x_vals = [rand(1:L-1), rand(1:L-1), L,           L, rand(1:L)]
        t_vals = [rand(1:L-1), L,           rand(1:L-1), L, rand(1:L)]
        for link_ind = 5:5
            μ = μ_vals[link_ind]
            x = x_vals[link_ind]
            t = t_vals[link_ind]
            for pert_dir = 0:3
                lie_directions = [0,1,2,3]
                popat!(lie_directions, pert_dir+1)
                for lie_dir in lie_directions
                    sdiff     = sdiff_min_pert_spec_conf_vec_neigh(L, q, z, μ, x, t, hitsize_range, lie_dir, pert_dir, pert_size) ./ eps()
                    fit_sdiff = curve_fit(model_sq_pert, hitsize_range, sdiff, p_sdiff_pert)
                    a = round.(fit_sdiff.param[1]*magnif^2, sigdigits = 3)
                    b = round.(fit_sdiff.param[2], sigdigits = 3)
                    c = round.(fit_sdiff.param[3], sigdigits = 3)
                    last_fit_converged = fit_sdiff.converged
                    image_sdiff = plot(
                        title = latexstring("\$ \\Delta S \$ of the special config at \$(q,z) = ($q,$z)\$ \n NEIGH \$(\\mu, x, t) = ($μ,$x,$t)\$, pert. size: \$ $pert_size\$ \n pert. dir.: \$ $pert_dir\$, Lie-dir.: \$ $lie_dir\$"),
                        xlabel = "hitsize",
                        ylabel = latexstring("\$\\Delta S / \\left(\\beta \\cdot \\epsilon_\\textrm{machine} \\right)\$"),
                        rightmargin = 5mm,
                        labelfontsize = 15,
                        tickfontsize = 10,
                        legend = :inside
                    )
                    image_sdiff = plot!(
                        fit_plot_xvals,
                        model_sq_pert(fit_plot_xvals, fit_sdiff.param),
                        label = latexstring("\$\\mathrm{fit}(x) = a\\cdot (x-b)^2 + c\$ \n \$a = $a\$ \n \$b = $b\$ \n \$c = $c\$"),
                        color = cb_orange
                    )
                    image_sdiff = scatter!(
                        hitsize_range,
                        sdiff,
                        label = "Measurements",
                        color = cb_blue
                    )
                    display(image_sdiff)
                end # lie_dir
            end # pert_dir
        end # link_ind
    end # z
    end # q
end # let


### Observations:
### pert. size 1e-3, q = 5, z = 3 
### Parameters b,c significantly different from 0 only for Lie direction ∈ {1,2} 
###     Especially for [μ,x,t] = [1,24,32]





######## Global perturbation and measure metric distance ########





L = 32
q = 1
z = 5
N_smear = Int(1e4)
ρ = 0.1
epsilons = round.(vcat([0.0], [10.0^(-i) for i = 2:7]), digits = 8)

dist_actions = Array{Float64}(undef,N_smear+1,length(epsilons))
metric_inds  = Int.(N_smear .* collect(0.99:0.002:1.0))
dist_metrics_opt = []
dist_metrics_min = []

#=
for i in eachindex(epsilons)
    println("Started smearing at i = $i")
    push!(dist_metrics_opt, [])
    # push!(dist_metrics_min, [])
    count = 0
    ϵ = epsilons[i]
    U = [ran_U2(rand()*ϵ) for μ = 1:2, x = 1:L, t = 1:L] .* insta_U2_z(L, L, q, z)
    dist_actions[1,i]   = action(U,1)
    V = stout_midpoint(U,ρ)
    for smear = 1:N_smear
        dist_actions[smear+1,i] = action(V,1)
        V = stout_midpoint(V,ρ)
        if smear in metric_inds
            push!(dist_metrics_opt[i], minimum([optimize_insta_metric(V, coeffs) for coeffs in [zeros(4), ones(4), 2 .* rand(4) .-1]]))
            # push!(dist_metrics_min[i], minimum_insta_metric(V))
        end
        if smear%(N_smear/20) == 0
            count += 5
            println("    Smearing progress: $count%")
        end
    end
end

dist_actions_path = string(data_path,"\\glob_dist_actions\\glob_dist_actions.txt")
# writedlm(dist_actions_path, dist_actions)
for i = 1:length(epsilons)
    dist_metric_opt_path = string(data_path,"\\glob_dist_actions\\glob_dist_metric_opt_epsnr_$i.txt")
    dist_metric_min_path = string(data_path,"\\glob_dist_actions\\glob_dist_metric_min_epsnr_$i.txt")
    # writedlm(dist_metric_opt_path, dist_metrics_opt[i])
    # writedlm(dist_metric_min_path, dist_metrics_min[i])
end
=#

dist_actions_path = string(data_path,"\\glob_dist_actions\\glob_dist_actions.txt")
dist_actions = readdlm(dist_actions_path)
dist_metrics_opt = []
dist_metrics_min = []
for i = 1:length(epsilons)
    dist_metric_opt_path = string(data_path,"\\glob_dist_actions\\glob_dist_metric_opt_epsnr_$i.txt")
    dist_metric_min_path = string(data_path,"\\glob_dist_actions\\glob_dist_metric_min_epsnr_$i.txt")
    push!(dist_metrics_opt, readdlm(dist_metric_opt_path))
    push!(dist_metrics_min, readdlm(dist_metric_min_path))
end



let
    τ0 = 2
    τf = 750
    start_ind = Int(τ0/ρ)
    final_ind = Int(τf/ρ)
    image_dist = plot(
        xlabel = latexstring("flow time \$\\tau/a^2\$"),
        ylabel = L"S/\beta",
        xticks = 0:250:750,
        tickfontsize = 10,
        labelfontsize = 15,
        legendfontsize = 11,
        legend = :top,
        background_color_legend = RGBA(1.0,1.0,1.0,0.8),
    )
    greys = (:grey77, :grey70, :grey63, :grey56, :grey49, :grey42, :grey30)
    for i = 1:length(epsilons)-1 
        j = length(epsilons)-i+1
        image_dist = plot!(
            ρ.*Array(0:N_smear)[start_ind:final_ind],
            dist_actions[start_ind:final_ind,j],
            label = latexstring("\$\\epsilon = $(epsilons[j])\$"),
            color = greys[i],
            linewidth = 2.5,
            alpha = 1
        )
    end
    image_dist = plot!(
        ρ.*Array(0:N_smear)[start_ind:final_ind],
        dist_actions[start_ind:final_ind,1],
        label = latexstring("\$\\epsilon = 0.0\$"),
        color = :black,
        linewidth = 2,
        alpha = 0.8
    )
    image_dist = hline!(
        [insta_action_U2_min(1,L,L,q)],
        label = latexstring("lower bound for \$q=1\$"),
        color = cb_red,
        linestyle = :dot,
        legend = :right,
        linewidth = 2.5,
        alpha = 0.8,
    )
    display(image_dist)
end

let
    image_insta_metrics = plot(
        title = latexstring("Optim. norm of glob. dist. smeared configs\n\$(q,z)=($q,$z) compared to instanton"),
        xlabel = latexstring("flow time \$\\tau\$"),
    )
    for i = 2:length(epsilons)
        image_insta_metrics = scatter!(
            ρ.*metric_inds,
            dist_metrics_opt[i],
            label = latexstring("\$\\epsilon = $(epsilons[i])\$"),
            yaxis = :log,
            xticks = ρ.*metric_inds,
            legend = :outerright,
        )
    end
    display(image_insta_metrics)
end

let
    image_insta_metrics = plot(
        title = latexstring("Anal. norm of glob. dist. smeared configs\n\$(q,z)=($q,$z) compared to instanton"),
        xlabel = latexstring("flow time \$\\tau\$"),
    )
    for i = 2:length(epsilons)
        image_insta_metrics = scatter!(
            ρ.*metric_inds,
            dist_metrics_min[i],
            label = latexstring("\$\\epsilon = $(epsilons[i])\$"),
            yaxis = :log,
            xticks = ρ.*metric_inds,
            legend = :outerright,
        )
    end
    display(image_insta_metrics)
end





######## Local perturbation and measure metric distance ########





L        = 32
q        = 1
z        = 5
N_smear  = Int(1e4)
ρ        = 0.1
pertsize = 1e-1

lie_dirs = [1,2]
dist_actions = Array{Float64}(undef,N_smear+1,2)
dist_insta_metric_z1 = []
dist_insta_metric_z4 = []
dist_insta_metric    = []
U = insta_U2_z(L, L, q, z)
μ = rand(1:2)
x = rand(1:L)
t = rand(1:L)

metric_inds_z4 = collect(1000:100:1600)
metric_inds_z1 = collect(2300:100:3500)
metric_inds    = collect(9500:100:10000)

for lie_dir_ind = 1:length(lie_dirs)
    lie_dir = lie_dirs[lie_dir_ind]
    count = 0
    println("Perturbation Lie direction: $lie_dir")
    V = deepcopy(U)
    V[μ,x,t] = small_U2_single_lie_direction(lie_dir, pertsize) * V[μ,x,t]
    dist_actions[1,lie_dir_ind] = action(V,1)

    for smear = 1:N_smear
        V = stout_midpoint(V,ρ)
        dist_actions[smear+1,lie_dir_ind] = action(V,1)

        if smear in metric_inds_z4
            push!(dist_insta_metric_z4, minimum([optimize_special_metric(V, coeffs, 4) for coeffs in [zeros(2), ones(2), 2 .* rand(2) .-1]]))
        end
        if smear in metric_inds_z1
            push!(dist_insta_metric_z1, minimum([optimize_special_metric(V, coeffs, 1) for coeffs in [zeros(2), ones(2), 2 .* rand(2) .-1]]))
        end
        if smear in metric_inds
            push!(dist_insta_metric, minimum([optimize_insta_metric(V, coeffs) for coeffs in [zeros(4), ones(4), 2 .* rand(4) .-1]]))
        end

        if smear%(N_smear/20) == 0
            count += 5
            println("    Smearing progress: $count%")
        end
    end
end

dist_actions_path = string(data_path,"\\loc_dist_actions\\loc_dist_actions_q_$(q)_z_$(z)_pertsize_$(pertsize).txt")
dist_insta_metric_z4_path = string(data_path,"\\loc_dist_actions\\loc_dist_metric_z4_q_$(q)_z_$(z)_pertsize_$(pertsize).txt")
dist_insta_metric_z1_path = string(data_path,"\\loc_dist_actions\\loc_dist_metric_z1_q_$(q)_z_$(z)_pertsize_$(pertsize).txt")
dist_insta_metric_path = string(data_path,"\\loc_dist_actions\\loc_dist_metric_insta_q_$(q)_z_$(z)_pertsize_$(pertsize).txt")
# # writedlm(dist_actions_path, dist_actions)
# # writedlm(dist_insta_metric_z4_path, dist_insta_metric_z4)
# # writedlm(dist_insta_metric_z1_path, dist_insta_metric_z1)
# # writedlm(dist_insta_metric_path, dist_insta_metric)
# dist_actions = readdlm(dist_actions_path)
# dist_insta_metric_z4 = readdlm(dist_insta_metric_z4_path)
# dist_insta_metric_z1 = readdlm(dist_insta_metric_z1_path)
# dist_insta_metric = readdlm(dist_insta_metric_path)


let
    scatter(ρ.*metric_inds_z4, dist_insta_metric_z4[1:length(metric_inds_z4)], label = latexstring("\$z=4, j=1\$"))
    scatter!(ρ.*metric_inds_z4, dist_insta_metric_z4[length(metric_inds_z4)+1:end], label = latexstring("\$z=4, j=2\$"))
    plot!(title = latexstring("Half-optim. norm of loc. dist. smeared configs \n\$(q,z)=(1,5)\$ compared to \$(1,4)\$"),
        xlabel = latexstring("flow time \$\\tau\$"),
        yaxis = :log
    )
end

let
    scatter(ρ.*metric_inds_z1, dist_insta_metric_z1[1:length(metric_inds_z1)], label = latexstring("\$z=1, j=1\$"))
    scatter!(ρ.*metric_inds_z1, dist_insta_metric_z1[length(metric_inds_z1)+1:end], label = latexstring("\$z=1, j=2\$"))
    plot!(title = latexstring("Half-optim. norm of loc. dist. smeared configs \n\$(q,z)=(1,5)\$ compared to \$(1,1)\$"),
        xlabel = latexstring("flow time \$\\tau\$"),
        yaxis = :log
    )
end

let
    scatter(ρ.*metric_inds, dist_insta_metric[1:length(metric_inds)], label = latexstring("\$j=1\$"))
    scatter!(ρ.*metric_inds, dist_insta_metric[length(metric_inds)+1:end], label = latexstring("\$j=2\$"))
    plot!(title = latexstring("Optim. norm of loc. dist. smeared configs \n\$(q,z)=(1,5)\$ compared to instanton"),
        xlabel = latexstring("flow time \$\\tau\$"),
        yaxis = :log
    )
end
