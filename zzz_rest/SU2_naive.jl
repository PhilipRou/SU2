N_t = 128
N_x = 128
β   = 5.0
ϵ   = 1.0
N_therm = 2000
acc_count = [0]

using LinearAlgebra
# using BenchmarkTools
################################################################################
# Pauli matrices:
σ0 = [1 0; 0 1]
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

function plaq(U, t, x)
    t_plus = t%N_t + 1
    x_plus = x%N_x + 1
    return U[2,t,x] * U[1,t,x_plus] * adjoint(U[2,t_plus,x]) * adjoint(U[1,t,x])
end

function action(U)
    S = 2*N_t*N_x  # later generalization: N_dim * U.V
    for t = 1:N_t
        for x = 1:N_x
            S -= real(tr(plaq(U,t,x)))
        end
    end
    return β*S/2
end

function staple_dag(U, t, x, dir)
    a = zeros(ComplexF64, 2, 2)
    b = zeros(ComplexF64, 2, 2)
    t_plu = t%N_t +1                 
    x_plu = x%N_x +1                 
    t_min = (t + N_t -2)%N_t +1   
    x_min = (x + N_x -2)%N_x +1   

    if dir == 1
        a = U[2,t_plu,x] * adjoint(U[1,t,x_plu]) * adjoint(U[2,t,x])
        b = adjoint(U[2,t_plu,x_min]) * adjoint(U[1,t,x_min]) * U[2,t,x_min]
    elseif dir == 2
        a = U[1,t,x_plu] * adjoint(U[2,t_plu,x]) * adjoint(U[1,t,x])
        b = adjoint(U[1,t_min,x_plu]) * adjoint(U[2,t_min,x]) * U[1,t_min,x]
    end
    return a + b
end

function delta_S_gauge(U, t::Int64, x::Int64, dir::Int64, old_link, new_link)
    return β*0.5*real(tr((old_link - new_link) * staple_dag(U,t,x,dir)))
end

function local_metro!(U, t, x, dir)
    old_link = deepcopy(U[dir,t,x])
    new_link = alg2grp_SU2(ϵ.* (rand(3).-0.5)) * old_link
    if rand() < exp(-delta_S_gauge(U, t, x, dir, old_link, new_link))
        U[dir,t,x] = new_link
        acc_count[1] += 1
    end
    return nothing
end
    
function lexiko_metro!(U)
    for dir = 1:2
        for t = 1:N_t
            for x = 1:N_x
                local_metro!(U,t,x,dir)
            end
        end
    end
    return nothing
end





actions = []
U = Array{Array,3}(undef,2,N_t,N_x)
for μ = 1:2
    for t = 1:N_t
        for x = 1:N_x
            ran_coeffs = ϵ .*2 .* (rand(3) .- 0.5)
            U[μ,t,x] = alg2grp_SU2(ran_coeffs)
        end
    end
end

count = 0

for i = 1:N_therm
    lexiko_metro!(U)
    push!(actions, action(U))
    if i%(N_therm/10) == 0
        count += 1
        println(10*count, "% done")
    end
end

auto_corr_time(actions[1:2000])
acc_count[1] / (2 * 128^2 * N_therm)
mean(actions[100:2000])/β/128^2
bootstrap(actions[100:2000]./β./128^2, 10, 50)

plot(actions[1:200]./β./128^2)

