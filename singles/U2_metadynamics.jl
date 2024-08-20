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

function meta_charge(U, N_smear_meta, ρ_meta)
    NX = size(U,2)
    NT = size(U,3)
    V = stout_midpoint_fast(U, N_smear_meta, ρ_meta)
    return sum([imag(det(plaq(V, x, t))) for x = 1:NX, t = 1:NT]) / 2 / π
end

function wilson_charges(U, N_smear_meta, N_smear, ρ)
    V = stout_midpoint_fast(U, N_smear_meta, ρ)
    for smear = 0:N_smear-1
        V = stout_midpoint_fast(V, N_smear_meta, ρ)
        # println("1st loop")
    end
    q_meta = top_charge_U2_wil(V)
    for smear = 0:N_smear-N_smear_meta-1
        V = stout_midpoint_fast(V, N_smear_meta, ρ)
        # println("2nd loop")
    end
    q = top_charge_U2_wil(V)
    return q, q_meta
end

bla = gaugefield_U2(8, 8, true);
for i = 1:100 chess_metro!(bla, 0.1, 5.0, [0.0], "U2") end
wilson_charges(bla, 1, 2, 0.1)
top_charge_U2_wil(stout_midpoint_fast(bla, 2, 0.1))

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
    ind = naive_ind(q_max, δq, q)
    lower_ind = round(Int, ind, RoundDown)
    bla = (q - ind2metaq(q_max, δq, lower_ind)) / δq
    bias_pot[lower_ind]     += w * (1 - bla)
    bias_pot[lower_ind + 1] += w * bla
    return nothing
end

function read_bias(bias_pot, q_max, δq, q)
    ind = naive_ind(q_max, δq, q)
    lower_ind = round(Int, ind, RoundDown)
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
    metaq_new = meta_charge(U_new, N_smear_meta, ρ_meta)  # ⛔ 
    metaq_old = meta_charge(U, N_smear_meta, ρ_meta)
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
        metaq_old = metaq_new
    end
    update_bias!(bias_pot, q_max, δq, w, metaq_old)
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
                    # @show x,t,trip,μ
                    metro_bias!(U,μ,x,t,step, β, acc, group, bias_pot, q_max, δq)
                end
            end
        end
    end
    return nothing
end



for N_smear = 0:10
    β = 8.0
    L = 32
    N_x = L
    N_t = L
    N_therm = 250
    N_meas  = 1000
    N_skip  = 10
    N_smear_meta = N_smear
    ρ_meta       = 0.1

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

    count = 0
    charges = []
    meta_charges = []
    for meas = 1:N_meas
        for skip = 1:N_skip
            chess_metro!(U, ϵ, β, acc_metro, "U2")
        end
        push!(charges, top_charge_U2(U))
        push!(meta_charges, meta_charge(U, N_smear_meta, ρ_meta))
        if meas%Int(N_meas/20) == 0
            count += 5
            println("Progress: $count%")
        end
    end
    # plot(charges)
    # plot!(meta_charges)

    hist_q = histogram2d(round.(Int,charges), meta_charges, bins = 100, title = "β = $β, L = $L, N_smear = $N_smear_meta")
    display(hist_q)
end




β = 12.0
L = 32
N_x = L
N_t = L
N_therm = 250
N_build = 1e4
N_smear_meta = 3
ρ_meta       = 0.1
# N_meas  = 1000
# N_skip  = 1

δq = 1e-2
w  = 1e-3
q_max = 15
# length(Array(-q_max:δq:q_max))
bias_pot = zeros(Int(2*q_max/δq+1))

ϵ = 0.1
acc_wish = 0.8
acc_therm = [0.0]
acc_metro = [0.0]

U = gaugefield_U2(N_x, N_t, false)
actions_therm = []
for therm = 1:N_therm
    chess_metro!(U, ϵ, β, acc_therm, "U2")
    ϵ *= sqrt(acc_therm[1] / acc_wish)
    # push!(actions_therm, action(U,β))
end
# plot(actions_therm)

count = 0
acc_bias_count = 0
for build = 1:N_build
    # chess_metro_bias!(U, ϵ, β, acc_therm, "U2", bias_pot, q_max, δq, w)
    metaq_old = meta_charge(U, N_smear_meta, ρ_meta)
    U_copy = deepcopy(U)
    chess_metro!(U_copy, ϵ, β, [0.0], "U2")
    metaq_new = meta_charge(U_copy, N_smear_meta, ρ_meta)
    S_bias_old = read_bias(bias_pot, q_max, δq, metaq_old)
    S_bias_new = read_bias(bias_pot, q_max, δq, metaq_new)
    if rand() < exp(-(S_bias_new - S_bias_old))
        U = U_copy
        metaq_old = metaq_new
        acc_bias_count += 1
    end
    update_bias!(bias_pot, q_max, δq, w, metaq_old)
    bias_pot = (bias_pot + reverse(bias_pot))/2
    if build%Int(N_build/20) == 0
        count += 5
        println("Progress building: $count%")
    end
end

# histogram2d(round.(Int,charges), meta_charges, bins = 100, title = "N_smear = $N_smear_meta")

cut = 900
# bias_pot_sym = (bias_pot + reverse(bias_pot))/2
display(plot(-q_max:δq:q_max, bias_pot, xticks = -14:2:14))
display(plot(Array(-q_max:δq:q_max)[cut+1:end-cut], bias_pot[cut+1:end-cut], xticks = -14:2:14))




