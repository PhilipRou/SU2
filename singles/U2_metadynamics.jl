include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\analyze\\SU2_analyze_head.jl")


using Plots
using LaTeXStrings
using LinearAlgebra
using DelimitedFiles
using Statistics

# function top_charge_U2(U)
#     NX = size(U,2)
#     NT = size(U,3)
#     return sum([imag(log(det(plaq(U, x, t)))) for x = 1:NX, t = 1:NT]) / 2 / π
# end

function meta_charge(U)
    NX = size(U,2)
    NT = size(U,3)
    return sum([imag(det(plaq(U, x, t))) for x = 1:NX, t = 1:NT]) / 2 / π
end

function ind2metaq(q_max, δq, i)
    # return Array(-q_max:δq:q_max)[i]
    return -q_max + (i-1)*δq
end

#=
function lower_ind(q_max, δq, q)
    return round(Int, (q+q_max)/δq, RoundDown) +1
end

# ind2metaq(q_max, δq, 3)
# lower_ind(q_max, δq, -15)
# [lower_ind(q_max, δq, i) for i = -14.9999:δq:15.0001] == Array(1:3001)
=#

function naive_ind(q_max, δq, q)
    return (q+q_max)/δq +1
end

# isapprox([round(Int, naive_ind(q_max, δq, i), RoundDown) for i = -14.9999:δq:15.0001], Array(1:3001))

function update_bias!(bias_pot, q_max, δq, w, q)
    naive_ind = naive_ind(q_max, δq, q)
    lower_ind = round(Int, naive_ind, RoundDown)
    bla = (q - ind2metaq(q_max, δq, lower_ind)) / δq
    bias_pot[lower_ind]     += w * (1 - bla)
    bias_pot[lower_ind + 1] += w * bla
    return nothing
end

function read_bias(bias_pot, q_max, δq, q)
    naive_ind = naive_ind(q_max, δq, q)
    lower_ind = round(Int, naive_ind, RoundDown)
    bla = (q - ind2metaq(q_max, δq, lower_ind)) / δq
    return w*bias_pot[lower_ind] + (1-w)*bias_pot[lower_ind+1]
end

function metro_bias!(U, μ, x, t, step, β, acc, group, bias_pot, q_max, δq)
    NX = size(U,2)
    NT = size(U,3)
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    end
    U_new = deepcopy(U)             # ⛔ we can do better...
    U_new[μ,x,t] = new_coeffs       # ⛔ 
    new_metaq = meta_charge(U_new)  # ⛔ 
    old_metaq = meta_charge(U)
    # new_metaq = old_metaq
    # if μ == 1
    #     t_m = mod1(t-1,NT)
    #     new_plaq = 
    #     new_metaq -= (imag(det(plaq(U,x,t)) + det(plaq(U,x,t_m)))) / (2*π)
    # elseif μ == 2

    # end
    staple_d = staple_dag(U,μ,x,t)
    S_old = -β*0.5*real(tr(U[μ,x,t] * staple_d)) + read_bias(bias_pot, q_max, δq, metaq_old)
    S_new = -β*0.5*real(tr(new_coeffs * staple_d)) + read_bias(bias_pot, q_max, δq, metaq_new)
    if rand() < exp(-(S_new-S_old))
        U[μ,x,t] = new_coeffs
        acc[1] += 1/NX/NT/2
    end
    return nothing
end

function chess_metro_bias!(U, step, β, acc, group, bias_pot, q_max, δq, w)
    NX = size(U,2)
    NT = size(U,3)
    acc[1] = 0.0
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    metro_bias!(U,μ,x,t,step, β, acc, group, bias_pot, q_max, δq)
                    update_bias!(bias_pot, q_max, δq, w, meta_charge(U))
                end
            end
        end
    end
    return nothing
end



β = 12.0
L = 32
N_x = L
N_t = L
N_therm = 250
N_meas  = 1000
N_skip  = 10

ϵ = 0.1
acc_wish = 0.8
acc_therm = [0.0]
acc_metro = [0.0]

U = gaugefield_U2(N_x, N_t, true)
actions_therm = []
for therm = 1:N_therm
    chess_metro!(U, ϵ, β, acc_therm, "U2")
    ϵ *= sqrt(acc_therm[1] / acc_wish)
    # push!(actions_therm, action(U,β))
end
# plot(actions_therm)

charges = []
meta_charges = []
for meas = 1:N_meas
    for skip = 1:N_skip
        chess_metro!(U, ϵ, β, acc_metro, "U2")
    end
    push!(charges, top_charge_U2(U))
    push!(meta_charges, meta_charge(U))
end
# plot(charges)
# plot!(meta_charges)

histogram2d(round.(Int,charges), meta_charges, bins = 100)







β = 12.0
L = 32
N_x = L
N_t = L
N_therm = 250
N_meas  = 1000
N_skip  = 10

δq = 1e-2
w  = 1e-3
q_max = 15
# length(Array(-q_max:δq:q_max))
bias_pot = Vector{Float64}(undef, Int(2*q_max/δq+1))

ϵ = 0.1
acc_wish = 0.8
acc_therm = [0.0]
acc_metro = [0.0]

U = gaugefield_U2(N_x, N_t, true)
actions_therm = []
for therm = 1:N_therm
    chess_metro!(U, ϵ, β, acc_therm, "U2")
    ϵ *= sqrt(acc_therm[1] / acc_wish)
    # push!(actions_therm, action(U,β))
end
# plot(actions_therm)

charges = []
meta_charges = []
for meas = 1:N_meas
    for skip = 1:N_skip
        chess_metro_bias!(U, ϵ, β, acc_therm, "U2", bias_pot, q_max, δq, w)
    end
    push!(charges, top_charge_U2(U))
    push!(meta_charges, meta_charge(U))
end
# plot(charges)
# plot!(meta_charges)

histogram2d(round.(Int,charges), meta_charges, bins = 100)