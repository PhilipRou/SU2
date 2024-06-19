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

function insta_action(β, N_c, N_x, N_t, Q, z)
    ReTrP = (N_c-1)*cos(2*π*z/N_x/N_t) + cos(2*π*(Q - (N_c-1)*z)/N_x/N_t ) 
    return β*N_x*N_t*(1 - ReTrP/N_c)
end

function insta_U2_w(N_x, N_t, Q, z)
    w = π*(2*z-Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*w/N_x/N_t)*coeffs_Id_U2() - sin(t*w/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*w/N_x)*coeffs_Id_U2() + sin(x*w/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    return U
end

function insta_U2_z(N_x, N_t, Q, z)
    # w = π*(2*z-Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * exp_u2(-t*2*π/N_x/N_t * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * exp_u2(x*2*π/N_x * (z-Q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x]
    return U
end

function insta_U2_slow(N_x, N_t, Q, z)
    # w = π*(2*z-Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * grp2coeffs_U2(exp(coeffs2grp(-t*π/N_x/N_t * (z-Q) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))))) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * grp2coeffs_U2(exp(coeffs2grp(x*π/N_x * (z-Q) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))))) for x = 1:N_x]
    return U
end


for i = 1:100
    N_t = 32
    N_x = 32
    β = 8.0
    global U = gaugefield_U2(N_t, N_x, true)
    actions = []
    for i = 1:300 chess_metro!(U, 0.05, β, [0.0], "U2"); push!(actions, action(U,β)) end
    Q = round(Int, top_charge_U2(U))
    U = U.*insta_U2(N_x, N_t, 1-Q)
    for i = 1:20 chess_metro!(U, 0.05, β, [0.0], "U2"); push!(actions, action(U,β)) end
    # display(plot(actions))
    V = stout(U, 0.01)
    for i = 1:10
        V = stout(V, 0.01)
    end
    Q = round(Int, top_charge_U2(V))
    println("Q = $Q")
    if abs(round(Int,Q)) == 1
        break
    end
end



ρ = 0.01
smear_pot = 5
N_smear = 10^smear_pot
smeared_actions = [action(U,1.0)/N_x/N_t]
count = 0
Queues = []
V = stout(U,ρ);
for i = 1:N_smear
    push!(smeared_actions, action(V,1.0)/N_x/N_t)
    V = stout(V,ρ)
    if i%Int(N_smear/100) == 0
        count += 1
        println("Smearing Progress: $count%")
        push!(Queues, top_charge_U2(V))
    end
end
println("Done!")

plot(smeared_actions_10[10^4:100:10^5], label = :false, legend = :right)
plot!(smeared_actions_8[10^4:100:10^5], label = :false)
plot!(smeared_actions_9[10^4:100:10^5], label = :false)
hline!([insta_action(1.0, 2, N_x, N_t, 0, 0.0)/N_x/N_t], label = "z = 0")
hline!([insta_action(1.0, 2, N_x, N_t, 0, 0.5)/N_x/N_t], label = "z = 0.5")
hline!([insta_action(1.0, 2, N_x, N_t, 0, 1.0)/N_x/N_t], label = "z = 1")
# plot!(xaxis = :log)


# # # @assert 1==0 "Bro, just check it twice"; writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\sms\\sms_12.txt", smeared_actions)
# using DelimitedFiles
smeared_actions_8 = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\old_sms\\sms_8.txt")
smeared_actions_9 = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\old_sms\\sms_9.txt")
smeared_actions_10 = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\old_sms\\sms_11.txt")
# Q = 1




first_not_ind = 0
for i = 1:length(smeared_actions_2)
    if isnan(smeared_actions_2[i])
        # first_not_ind = i
        println(i)
        break
    end
end
last_ind = first_not_ind - 1


function find_z(actions, ind)
    f(z) = insta_action(1.0, 2, N_x, N_t, 0, z) / N_x / N_t - actions[ind]
    # f(2000)
    # z_smear = find_zero(f, (0,10^5))
    return find_zero(f, (0,ind))
end
z_not_int = find_z(smeared_actions_9,10^5)
# z_low = round(Int, z_smear, RoundDown)
# z_up = round(Int, z_smear, RoundUp)
z_int = round(Int, 2*find_z(smeared_actions_9,10^5))/2
# action_insta_up = insta_action(1.0, 2, N_x, N_t, Q, z_up) / N_x / N_t
# action_insta_low = insta_action(1.0, 2, N_x, N_t, Q, z_low) / N_x / N_t
action_insta_int = insta_action(1.0, 2, N_x, N_t, 0.0, z_int) /N_x /N_t
# action_diff = round(smeared_actions[last_ind] - action_insta_low, sigdigits = 3)
# smeared_actions[10^5] - smeared_actions[last_ind]

smeared_actions_sub = smeared_actions .- minimum(smeared_actions)
upper_end = 10^5
lower_end = 10^2
step_size = 10
rho_vals = Array(lower_end:step_size:upper_end)
image_smeared_actions = plot(
    ρ .* rho_vals,
    smeared_actions_sub[rho_vals],
    title = latexstring("\$ S/\\beta V \$ for a Stout smeared 2D U(2) field \n initial \$ q = $Q, \\rho = $ρ \$" ),
    xlabel = latexstring("Smearing time \$\\tau\$"),
    label = "Smeared actions minus minimum",
    yaxis = :log,
    # xticks = [10.0^i for i = 0:0.5:3]
)




# z_max = N_x*N_t/2 + Int(round(Q/2, RoundDown))
# action_insta_max = insta_action(1.0, 2, N_x, N_t, Q, z_max) / N_x / N_t

# N_blocks = 100
# L = round(Int, length(smeared_actions)/N_blocks)
L = 100
image_test = plot(
    ρ .* Array(1:L),
    smeared_actions[1:L],
    # color = :blue,
    label = :false
)
# for i = 2:40
#     image_test = plot!(
#         ρ .* Array(1+(i-1)*L:i*L),
#         smeared_actions[1+(i-1)*L:i*L],
#         color = :blue,
#         label = :false
#     )
# end
image_test = plot!(
    ρ .* Array(1+L:10000:last_ind),
    smeared_actions[1+L:10000:last_ind],
    # color = :blue,
    label = :false
)
image_test = plot!(xaxis = :log)
display(image_test)

image_smeared_actions = plot(
    ρ .* Array(1:10^4),
    smeared_actions[1:10^4],
    title = latexstring("\$ S/\\beta V \$ for a smeared 2D U(2) field of top. charge \$ q = $Q\$"),
    xlabel = latexstring("Smearing time \$\\tau\$"),
    label = latexstring("Smeared action, \$\\rho = $ρ\$"),
    xaxis = :log,
    xticks = [10.0^i for i = -2:2]
)
# image_smeared_actions = hline!([action(insta_U2(N_x,N_t,Q),1.0)/N_x/N_t], label = "Min. instanton action")
# image_smeared_actions = hline!([action_insta_max], label = "Max. instanton action")
image_smeared_actions = hline!([action_insta_low], label = latexstring("Instanton action with \$ z = $z_low\$"))
image_smeared_actions = hline!([action_insta_up], label = latexstring("Instanton action with \$ z = $z_up\$"))

# image_smeared_action = plot!(xaxis = :log)#, xticks = [10.0^i for i = -2:6] )
# action(insta_U2(N_x,N_t,Q),1.0)/N_x/N_t

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\smeared_actions\\smeared_actions_q_$Q.3.pdf")









let
    # @assert 1==0 "Do we really want to start a smearing run???"
    β   = 3.0
    L   = 16
    N_x = L
    N_t = L
    ρ   = 0.01
    N_therm   = 500
    N_meas    = 100
    N_sepa    = 100
    N_smear   = 10^5
    acc_wish  = 0.8
    ϵ         = 0.1
    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_data_6"
    acc_metro = [0.0]
    actions_therm = []
    U = gaugefield(N_x, N_t, true, "U2", "square")
    for therm = 1:N_therm
        chess_metro!(U,ϵ,β,acc_metro,"U2")
        ϵ *= sqrt(acc_metro[1] / acc_wish) # only update ϵ acc. to Metropolis
        # push!(actions_therm, action(U,β))
    end
    # plot(actions_therm)
    actions = []
    for meas = 1:N_meas
        for sepa = 1:N_sepa
            chess_metro!(U,ϵ,β,acc_metro,"U2")
            push!(actions, action(U,β))
        end
        smeared_actions = [action(U,β)]
        smeared_charges = [top_charge_U2(U)]
        V = stout(U,ρ)
        count = 0
        for smear = 1:N_smear
            push!(smeared_actions,action(V,β))
            push!(smeared_charges,top_charge_U2(V))
            V = stout(V,ρ)
            if smear%Int(N_smear/100) == 0
                count += 1
                println("Measurement Nr.: $meas, Smearing Progress: $count%")
            end
        end
        S_path = string(base_path,"\\sms_$meas.txt")
        Q_path = string(base_path,"\\smq_$meas.txt")
        writedlm(S_path,smeared_actions)
        writedlm(Q_path,smeared_charges)
    end
    actions_path = string(base_path,"\\non_smeared_actions.txt")
    params_path = string(base_path, "\\params.txt")
    writedlm(actions_path, actions)
    writedlm(params_path, "L = $L\n β = $β")
    println("We're done here!")
end





β       = 16.0
L       = 32
ρ       = 0.01
N_sepa  = 100
N_smear = 10^5

base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_6\\sms_data_6"
base_fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_6"
nsms_path = string(base_path, "\\non_smeared_actions.txt")
nsms_fig_path = string(base_fig_path,"\\non_smeared_actions.pdf")

nsms = readdlm(nsms_path) ./ (β*L^2)

image_actions = plot(
    nsms,
    legend = :false,
    title  = latexstring("\$S/\\beta V \$ Time Series for Unsmeared Configs. \n 2D U(2), \$\\beta = $β, L = $L, N_{sepa} = $N_sepa\$ "),
    xlabel = "Monte Carlo Time"
)
display(image_actions)
# savefig(nsms_fig_path)





β       = 3.0
L       = 16
ρ       = 0.01
N_sepa  = 100
N_smear = 10^5
N_meas  = 100

base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_6\\sms_data_6"
base_fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_6"

queues = []

for i = 1:N_meas
# i = 15
    smq_path  = string(base_path, "\\smq_$i.txt")
    sms_path  = string(base_path, "\\sms_$i.txt")
    smq_fig_path  = string(base_fig_path, "\\smq_$i.pdf")
    sms_fig_path  = string(base_fig_path, "\\sms_$i.pdf")

    smq  = readdlm(smq_path)
    sms  = readdlm(sms_path) ./ (β*L^2)

    Q = round(Int,last(smq))
    push!(queues, Q)
    sms = sms .- minimum(sms)

    #=
    cut = 100
    sep = 100
    taus = ρ .* Array(1:N_smear+1)
    smq_plot = vcat(smq[1:cut], smq[cut+1:sep:end])
    sms_plot = vcat(sms[1:cut], sms[cut+1:sep:end])
    tau_plot = vcat(taus[1:cut], taus[cut+1:sep:end])

    q_window = 1:300
    image_smq = plot(
        taus[q_window],
        smq[q_window],
        legend = :false,
        title  = latexstring("Top. Charge \$q\$ of Smeared Config Nr. $i\n 2D U(2), \$\\beta = $β, L = $L, \\rho = $ρ\$"),
        xlabel = latexstring("Smearing Time \$\\tau\$")
    )
    display(image_smq)
    # savefig(smq_fig_path) # 🟥🟥🟥

    s_window = cut:length(sms_plot)-50
    image_sms = plot(
        tau_plot[s_window],
        sms_plot[s_window],
        legend = :false,
        title  = latexstring("Subtracted \$S/\\beta V \$ of Smeared Config Nr. $i\n 2D U(2), \$\\beta = $β, L = $L, \\rho = $ρ\$, final \$q=$Q\$"),
        xlabel = latexstring("Smearing Time \$\\tau\$"),
        yaxis = :log,
        # xticks = 0:100:900
    )
    display(image_sms)
    # savefig(sms_fig_path) # 🟥🟥🟥
    =#
end

q_max = maximum(abs.(queues))
configs_by_q_list = []
for q = 0:q_max
    q_inds = []
    for i = 1:N_meas
        if abs(queues[i]) == q
            push!(q_inds, i)
        end
    end
    push!(configs_by_q_list, q_inds)
end
configs_by_q_list[1]
[length(configs_by_q_list[q+1]) for q = 0:q_max]



# the measurements with final charge q=1:
# top_one_runs = [2,5,6,8,12] # 16, 3.0, 15
# top_one_runs = [8,11,15] # 32, 8.0, 15
# top_one_runs = Array(1:15) # 32, 16.0, 15
let
    q = 3
    confs_with_that_q = configs_by_q_list[q+1]

    β       = 3.0
    L       = 16
    ρ       = 0.01
    N_sepa  = 100
    N_smear = 10^5
    # N_smear = 5*10^4

    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_6\\sms_data_6"
    base_fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms"

    cut = 100
    sep = 100
    plot_length = cut + length(cut+1:sep:N_smear) + 1
    plot_sms = Array{Float64}(undef, plot_length, length(confs_with_that_q))
    for i = 1:length(confs_with_that_q)
    # i = 1
        conf = confs_with_that_q[i]
        sms_path = string(base_path, "\\sms_$conf.txt")
        sms = readdlm(sms_path) / (β*L^2)
        plot_sms[:,i] =  vcat(sms[1:cut], sms[cut+1:sep:end])
    end

    taus = ρ .* Array(1:N_smear+1)
    tau_plot = vcat(taus[1:cut], taus[cut+1:sep:end])
    s_window = 200:plot_length

    global image_top_ones = plot(
        tau_plot[s_window],
        plot_sms[s_window,1],
        # legend = :false,
        label = "Conf. Nr. $(confs_with_that_q[1])",
        title  = latexstring("\$S/\\beta V \$ of Various Smeared Configs.\n 2D U(2), \$\\beta = $β, L = $L, \\rho = $ρ\$, final \$q=\\pm $q \$"),
        xlabel = latexstring("Smearing Time \$\\tau\$"),
        # yaxis = :log,
        xticks = 0:100:900
    )
    for i = 2:length(confs_with_that_q)
        image_top_ones = plot!(
            tau_plot[s_window],
            plot_sms[s_window, i],
            label = "Conf. Nr. $(confs_with_that_q[i])",
        )
    end
    display(image_top_ones)
end

let
    q = 3
    z = q/2
    image_top_ones = hline!(
        [insta_action(β,2,L,L,q,z) / (β*L^2)],
        linestyle = :dash,
        color = :red,
        label = latexstring("Insta: \$q = $q, z = $z\$")
    )
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\smeared_actions_q_11.pdf")





q = 12
let
    confs_with_that_q = configs_by_q_list[q+1]

    β       = 3.0
    L       = 16
    ρ       = 0.01
    N_sepa  = 100
    N_smear = 10^5
    # N_smear = 5*10^4

    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_6\\sms_data_6"
    base_fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms"

    cut = 100
    sep = 100
    plot_length = cut + length(cut+1:sep:N_smear) + 1
    plot_sms = Array{Float64}(undef, plot_length, length(confs_with_that_q))
    long_sms = Array{Float64}(undef, N_smear+1, length(confs_with_that_q))
    last_sms = Array{Float64}(undef, length(confs_with_that_q))
    for i = 1:length(confs_with_that_q)
    # i = 1
        conf = confs_with_that_q[i]
        sms_path = string(base_path, "\\sms_$conf.txt")
        sms = readdlm(sms_path) / (β*L^2)
        long_sms[:,i] = sms
        plot_sms[:,i] =  vcat(sms[1:cut], sms[cut+1:sep:end])
        last_sms[i] = last(sms)
    end

    z = q/2
    last_sms .-= insta_action(β, 2, L, L, q, z)

    # z = q/2
    plot_sms = plot_sms .- [insta_action(β,2,L,L,q,z) / (β*L^2)]

    taus = ρ .* Array(1:N_smear+1)
    tau_plot = vcat(taus[1:cut], taus[cut+1:sep:end])
    s_window = 200:plot_length

    global image_top_ones = plot(
        tau_plot[s_window],
        plot_sms[s_window,1],
        # legend = :false,
        label = "Conf. Nr. $(confs_with_that_q[1])",
        title  = latexstring("\$S/\\beta V \$ minus Insta action at \$z=$z\$  \n 2D U(2), \$\\beta = $β, L = $L, \\rho = $ρ\$, final \$q=\\pm $q \$"),
        xlabel = latexstring("Smearing Time \$\\tau\$"),
        yaxis = :log,
        xticks = 0:100:900,
        yticks = [10.0^(-i) for i = 4:15],
        # legendfontsize = 5
    )
    for i = 2:length(confs_with_that_q)
        image_top_ones = plot!(
            tau_plot[s_window],
            plot_sms[s_window, i],
            label = "Conf. Nr. $(confs_with_that_q[i])",
        )
    end
    display(image_top_ones)



    # display(histogram(last_sms, title = "Last smeared action of configs with q = $q", bins = 100))#
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_6\\smeared_actions_q_$q.pdf")







#=

# the measurements with final q=2:
# top_two_runs = [3,4,7,10,13,14] # 16, 3.0, 15
top_two_runs = [2,4,12] # 32, 8.0, 15

cut = 100
sep = 100
plot_length = cut + length(cut+1:sep:N_smear) + 1
top_two_sms = Array{Float64}(undef, plot_length, length(top_two_runs))
for i = 1:length(top_two_runs)
    run = top_two_runs[i]
    sms_path = string(base_path, "\\sms_$run.txt")
    sms = readdlm(sms_path) / (β*L^2)
    top_two_sms[:,i] =  vcat(sms[1:cut], sms[cut+1:sep:end])
end

taus = ρ .* Array(1:N_smear+1)
tau_plot = vcat(taus[1:cut], taus[cut+1:sep:end])
s_window = 400:plot_length 

let
    global image_top_twos = plot(
        tau_plot[s_window],
        top_two_sms[s_window,1],
        # legend = :false,
        label = "Conf. Nr. $(top_two_runs[1])",
        title  = latexstring("\$S/\\beta V \$ of Various Smeared Configs.\n 2D U(2), \$\\beta = $β, L = $L, \\rho = $ρ\$, final \$q=\\pm 2 \$"),
        xlabel = latexstring("Smearing Time \$\\tau\$"),
        # yaxis = :log,
        xticks = 0:100:900
    )
    for i = 2:length(top_two_runs)
        image_top_twos = plot!(
            tau_plot[s_window],
            top_two_sms[s_window, i],
            label = "Conf. Nr. $(top_two_runs[i])",
        )
    end
    display(image_top_twos)
end

image_top_twos = hline!(
    [insta_action(β,2,L,L,2,1) / (β*L^2)],
    linestyle = :dash,
    color = :red,
    label = latexstring("Insta: \$q = 2, z = 1\$")
)

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\smeared_actions_q_2.pdf")







# the measurements with different final q's:
# top_dif_runs = [2,3,15,1,9]
# q_vals = [1,2,3,4,6]
top_dif_runs = [14,8,2,10,6,5,7]
q_vals = [0,1,2,3,4,7,8]

cut = 100
sep = 100
plot_length = cut + length(cut+1:sep:N_smear) + 1
top_two_sms = Array{Float64}(undef, plot_length, length(top_dif_runs))
for i = 1:length(top_dif_runs)
    run = top_dif_runs[i]
    sms_path = string(base_path, "\\sms_$run.txt")
    sms = readdlm(sms_path) / (β*L^2)
    top_two_sms[:,i] =  vcat(sms[1:cut], sms[cut+1:sep:end])
end

taus = ρ .* Array(1:N_smear+1)
tau_plot = vcat(taus[1:cut], taus[cut+1:sep:end])
s_window = 110:plot_length -50

let
    global image_top_difs = plot(
        tau_plot[s_window],
        top_two_sms[s_window,1],
        # legend = :false,
        label = latexstring("final \$q=$(q_vals[1])\$"),
        title  = latexstring("\$S/\\beta V \$ of Various Smeared Configs. of Different \$q\$ \n 2D U(2), \$\\beta = $β, L = $L, \\rho = $ρ\$"),
        xlabel = latexstring("Smearing Time \$\\tau\$"),
        # yaxis = :log,
        xticks = 0:100:900
    )
    for i = 2:length(top_dif_runs)
        image_top_difs = plot!(
            tau_plot[s_window],
            top_two_sms[s_window, i],
            label = latexstring("final \$q=$(q_vals[i])\$"),
        )
    end
    display(image_top_difs)
end

# let
#     image_top_difs = hline!(
#         [insta_action(β,2,L,L,1,0.5) / (β*L^2)],
#         linestyle = :dash,
#         color = palette(:default)[1],
#         label = latexstring("Insta: \$q = 1, z = 0.5\$")
#     )
#     image_top_difs = hline!(
#         [insta_action(β,2,L,L,2,1) / (β*L^2)],
#         linestyle = :dash,
#         color = palette(:default)[2],
#         label = latexstring("Insta: \$q = 2, z = 1\$")
#     )
#     image_top_difs = hline!(
#         [insta_action(β,2,L,L,3,1.5) / (β*L^2)],
#         linestyle = :dash,
#         color = palette(:default)[3],
#         label = latexstring("Insta: \$q = 3, z = 1.5\$")
#     )
#     image_top_difs = hline!(
#         [insta_action(β,2,L,L,4,2) / (β*L^2)],
#         linestyle = :dash,
#         color = palette(:default)[4],
#         label = latexstring("Insta: \$q = 4, z = 2\$")
#     )
#     image_top_difs = hline!(
#         [insta_action(β,2,L,L,6,3) / (β*L^2)],
#         linestyle = :dash,
#         color = palette(:default)[5],
#         label = latexstring("Insta: \$q = 6, z = 3\$")
#     )
#     display(image_top_difs)
# end

let
    image_top_difs = hline!(
        [insta_action(β,2,L,L,0,0) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[1],
        label = latexstring("Insta: \$q = 0, z = 0\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,1,0.5) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[2],
        label = latexstring("Insta: \$q = 1, z = 0.5\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,2,1) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[3],
        label = latexstring("Insta: \$q = 2, z = 1\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,3,1.5) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[4],
        label = latexstring("Insta: \$q = 3, z = 1.5\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,4,2) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[5],
        label = latexstring("Insta: \$q = 4, z = 2\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,7,3.5) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[6],
        label = latexstring("Insta: \$q = 7, z = 8\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,8,4) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[7],
        label = latexstring("Insta: \$q = 8, z = 4\$")
    )
    display(image_top_difs)
end


# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\smeared_actions_q_dif.pdf")

=#
q = 1
for z = 0.5:0.5:2.5
# z = 0.5
# z = 1.5
    L = 16
    N_x = L
    N_t = L
    β = 3.0
    ρ = 0.1
    N_smear = 10^4

    q = 1
    # z = 5
    U = insta_U2_w(N_x, N_t, q, z);
    # W = [ran_U2(0.0001) for μ = 1:2, x = 1:N_x, t = 1:N_t] .* U;
    # action(U,β)
    # action(W,β)
    # insta_action(β,2,N_x,N_t,1,5) - action(W,β)
    # insta_action(β,2,N_x,N_t,1,5.5) - insta_action(β,2,N_x,N_t,1,5.0)

    num_eps = 10
    epsilons = [10.0^(-i) for i = 4:4+num_eps]
    push!(epsilons, 0.0)
    num_eps = length(epsilons)
    actions = Array{Float64}(undef,N_smear+1,num_eps);
    charges = Array{Float64}(undef,N_smear+1,num_eps);

    for i = 1:num_eps
        ϵ = epsilons[i]
        W = [ran_U2(ϵ) for μ = 1:2, x = 1:N_x, t = 1:N_t] .* U
        V = stout_midpoint(W,ρ);
        insta_smeared_actions = [action(V,β)]
        insta_smeared_charges = [top_charge_U2(V)]
        count = 0
        for smear = 1:N_smear
            V = stout_midpoint(V,ρ)
            push!(insta_smeared_actions,action(V,β))
            push!(insta_smeared_charges,top_charge_U2(V))
            if smear%Int(N_smear/100) == 0
                count += 1
                println("q = $q, z= $z, Conf. Nr. $i, Smearing Progress: $count%")
            end
        end
        actions[:,i] = insta_smeared_actions
        charges[:,i] = insta_smeared_charges
    end

    cut = 100
    sep = 100
    taus = ρ .* Array(1:N_smear+1)
    # smq_plot = vcat(smq[1:cut], smq[cut+1:sep:end])
    # sms_plot = vcat(sms[1:cut], sms[cut+1:sep:end])
    # tau_plot = vcat(taus[1:cut], taus[cut+1:sep:end])

    image_queue = plot(
        taus,
        round.(Int,charges[:,1])
    )
    # image_actions = plot(
    #     insta_smeared_actions[1:2000]
    # )
    # image_actions = hline!([insta_action(β,2,N_x,N_t,1,5)])

    # image_actions = plot(
    #     insta_smeared_actions,
    #     label = :false,
    #     title = latexstring("Smeared Actions of an Instanton with \$q = $q, z = $z\$ \n 2D U(2), \$ L = $L, \\beta = $β, \\rho = $ρ \$ ")
    # )
    # image_actions = hline!(
    #     [insta_action(β,2,N_x,N_t,1,5)],
    #     label = latexstring("Insta: \$ q = 1, z=5\$")
    # )
    # image_actions = hline!(
    #     [insta_action(β,2,N_x,N_t,1,0.5)],
    #     label = latexstring("Insta: \$ q = 1, z=0.5\$")
    # )
    # for z = 1:0.5:4.5
    #     image_actions = hline!([insta_action(β,2,N_x,N_t,1,z)], color = :grey)
    # end
    # display(image_actions)


    image_actions = plot(
        taus,
        actions[:,1],
        label = latexstring("\$\\epsilon = $(round(epsilons[1], sigdigits = 1))\$"),
        title = latexstring("Smeared Actions of pert. Insta. with q = $q, z = $z\n 2D U(2), \$ L = $L, \\beta = $β, \\rho = $ρ \$ "),
        legend = :right,
        xlabel = latexstring("Smearing Time \$\\tau\$")
    )
    for i = 2:num_eps
        image_actions = plot!(
            taus,
            actions[:,i],
            label = latexstring("\$\\epsilon = $(round(epsilons[i], sigdigits = 1))\$"),
            )
    end
    image_actions = hline!(
        [insta_action(β,2,N_x,N_t,q,z)],
        label = latexstring("Insta. Action: \$ q = $q, z=$z\$"),
        linestyle = :dash,
        color = :red
    )
    image_actions = hline!(
        [insta_action(β,2,N_x,N_t,q,q/2)],
        label = latexstring("Insta. Action: \$ q = $q, z=$(q/2)\$"),
        linestyle = :dash,
        color = :red
    )
    if mod(z,1) != 0
        image_actions = hline!(
            [insta_action(β,2,N_x,N_t,q,z-0.5)],
            label = latexstring("Insta. Action: \$ q = $q, z=$(z-0.5)\$"),
            linestyle = :dash,
            color = :green
        )
    end
    display(image_actions)


    # image_actions = hline!(
    #     # [action(insta_U2_z(N_x,N_t,1,0.5), β)],
    #     [insta_action(β,2,N_x,N_t,1,0.5)],
    #     label = :false, # latexstring("Insta. Action: \$ q = 1, z=$(q/2)\$"),
    #     linestyle = :dash,
    #     color = :green
    # )

    # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\smeared_insta_q_$q.z_$z.pdf")
end



Q = 1
Z = 1.0;
# action(insta_U2_slow(N_x,N_t,Q,Z), 1)
action(insta_U2_z(N_x,N_t,Q,Z/2), 1)
action(insta_U2_w(N_x,N_t,Q,1), 1)
insta_action(1,2,N_x,N_t,Q,1/2)*2



L = 16
N_x = L
N_t = L
β = 12.0
ρ = 0.1
N_smear = 10^3

U = gaugefield_U2(N_x, N_t, true);

let
    actions_therm = []
    for therm = 1:200
        chess_metro!(U,0.1,β,[0.0],"U2") 
        push!(actions_therm, action(U,β))
    end
    display(plot(actions_therm))
end

V_normal = stout(U,ρ)
V_mid = stout_midpoint(U,ρ)
actions_smear_normal = [action(U,β), action(V_normal,β)]
actions_smear_mid = [action(U,β), action(V_mid,β)]
count = 0
for smear = 1:N_smear
    V_normal = stout(V_normal,ρ)
    V_mid = stout_midpoint(V_mid,ρ)
    push!(actions_smear_normal,action(V_normal,β))
    push!(actions_smear_mid,action(V_mid,β))
    if smear%Int(N_smear/100) == 0
        count += 1
        println("Smearing Progress: $count%")
    end
end

plot(actions_smear_normal[1:100]./(β*L^2))
plot!(actions_smear_mid[1:100]./(β*L^2))

plot!(Array(1:10:100),actions_smear_normal[1:10]./(β*L^2))
plot!(Array(1:10:100),actions_smear_mid[1:10]./(β*L^2))
plot((actions_smear_normal .- actions_smear_mid)./(β*L^2))

# savefig("C:\\Users\\proue\\Downloads\\test.pdf")


L = 16
N_x = L
N_t = L
ρ = 0.1

u = insta_U2_z(N_x,N_t,1,1);
v = [ran_U2(0.01) for μ = 1:2, x = 1:N_x, t = 1:N_t] .* u
# v = stout_midpoint(u,0.1);
s = [];
for i = 1:10000
    push!(s,action(v,1))
    v = stout_midpoint(v,0.1)
end

# plot(vcat(s[1:100], s[101:100:end]))
# plot(s[5000:end])
plot(0.1*Array(100:length(s)), s[100:end])
hline!([insta_action(1,2,N_x,N_t,1,1)])
hline!([insta_action(1,2,N_x,N_t,1,0.5)])
# hline!([insta_action(1,2,N_x,N_t,1,1.5)])

heatmap([real(tr(plaq(v,x,t))) for x = 1:N_x, t = 1:N_t])#,xticks = 1:16, yticks = 1:16)
heatmap([real(tr(plaq(temp_gauge(v,"U2"),x,t))) for x = 1:N_x, t = 1:N_t])#,xticks = 1:16, yticks = 1:16)
heatmap([real(tr(plaq(max_gauge(v,"U2"),x,t))) for x = 1:N_x, t = 1:N_t])#,xticks = 1:16, yticks = 1:16)
# plot!(s[7000:end] .- insta_action(1,2,N_x,N_t,1,0.5))
last(s) - insta_action(1,2,N_x,N_t,1,0.5)

heatmap([real(tr(plaq(insta_U2_z(N_x,N_t,1,0.5),x,t))) for x = 1:N_x, t = 1:N_t])#,xticks = 1:16, yticks = 1:16)
15/16*L^2-sum([real(tr(plaq(insta_U2_z(N_x,N_t,1,0.5),x,t))) for x = 1:N_x, t = 1:N_t-1])/2
15/16*insta_action(1,2,L,L,1,0.5)

uu = temp_gauge(v, "U2");
uuu = max_gauge(v, "U2");
# action(v,1) - action(uuu,1)
log_U2(uu[1,1,16])
log_U2(uu[1,1,5]/sqrt(det(u[1,1,5])))
log_U2(uu[2,9,16]) / (π/L)

# log_U2(uu[1,1,2])
log_U2(uuu[1,16,3])
log_n_squared(uuu[1,1,16])
coeffs2grp(uuu[2,9,16])^2
uuu[1,1,16].a * [1 0; 0 1]
uuu[1,1,15].a * [1 0; 0 1]
imag(log_U2(uu[1,1,3]).a)/π*L^2

isapprox(coeffs2grp(uuu[1,16,16]), uuu[1,16,16].a * [1 0; 0 1])
coeffs2grp(uuu[2,8,16]) * adjoint(coeffs2grp(uuu[2,9,16]))

action(uuu,1)
for x = N_x>>1:N_x
    X = uuu[2,x,L]
    uuu[2,x,L] = coeffs_U2(X.a, -X.b, -X.c, -X.d)
end

log_U2(uuu[1,L,L])
log_U2(uuu[2,7,L])
log_U2(uuu[1,1,1])
log_U2(uuu[2,1,16])

sum([log_U2(uuu[1,L,L]).b , log_U2(uuu[1,L,L]).c , log_U2(uuu[1,L,L]).d] .* [log_U2(uuu[2,1,L]).b, log_U2(uuu[2,1,L]).c, log_U2(uuu[2,1,L]).d])

sqrt(log_n_squared(uuu[2,7,16])) / π
sqrt(log_n_squared(uuu[1,16,1])) / π

sqrt(log_n_squared(uuu[1,L,L])) / π
sqrt(log_n_squared(uuu[1,L,L]*uuu[2,1,L])) / π
sqrt(log_n_squared(uuu[2,1,L]*adjoint(uuu[1,L,1]))) / π
sqrt(log_n_squared(uuu[1,L,L]*uuu[2,1,L]*adjoint(uuu[1,L,1]))) / π
sqrt(log_n_squared(uuu[2,1,L]*adjoint(uuu[1,L,1])*adjoint(uuu[2,L,L]))) / π
sqrt(log_n_squared(uuu[1,L,L]*uuu[2,1,L]*adjoint(uuu[1,L,1])*adjoint(uuu[2,L,L]))) / π
sqrt(log_n_squared(plaq(uuu,L,L))) / π


coeffs2grp(uuu[1,L,L])
coeffs2grp(uuu[1,L,L]*uuu[2,1,L])
coeffs2grp(uuu[1,L,L]*uuu[2,1,L]*adjoint(uuu[1,L,1]))
coeffs2grp(uuu[1,L,L]*uuu[2,1,L]*adjoint(uuu[1,L,1])*adjoint(uuu[2,L,L]))

uuu[1,L,L]*uuu[2,1,L] - uuu[2,1,L]*uuu[1,L,L]
log_U2(uuu[1,L,L])
log_U2(uuu[2,1,L])
log_U2(uuu[1,L,L]*uuu[2,1,L])

exp_u2(log_U2(uuu[1,L,L]) + log_U2(uuu[2,1,L]))



isapprox(plaq(uuu,8,16), uuu[1,8,16] * uuu[2,9,16] * adjoint(uuu[2,8,16]) )
uuu[2,9,16] * adjoint(uuu[2,8,16])


isapprox(exp_u2(2*log_U2(uuu[2,9,16])), uuu[2,9,16]*uuu[2,9,16])
log_U2(exp_u2(2*log_U2(uuu[2,9,16])))


[real(log_n_squared(uuu[1,x,1])) for x = 1:L]
log_U2(uuu[2,7,16] * adjoint(uuu[2,8,16]))
coeffs2grp(uuu[2,7,16]) * adjoint(coeffs2grp(uuu[2,8,16]))

σ0 = [1 0; 0 1]     # not really a Pauli matrix 🙄
σ1 = [0 1; 1 0]
σ2 = [0 -im; im 0]
σ3 = [1 0; 0 -1]

for xxx = 4:10
    b1, c1, d1 = uuu[2,xxx,16].b, uuu[2,xxx,16].c, uuu[2,xxx,16].d
    b2, c2, d2 = uuu[2,xxx+1,16].b, uuu[2,xxx+1,16].c, uuu[2,xxx+1,16].d
    A1 = sum(im.*[b1,c1,d1] .* [σ1, σ2, σ3])
    A2 = sum(im.*[b2,c2,d2] .* [σ1, σ2, σ3])
    println(A1*adjoint(A2))
end


function log_n_squared(X::coeffs_U2)
    bla = log_U2(X)
    return bla.b^2 + bla.c^2 + bla.d^2
end

bla = plot([real(log_n_squared(uuu[1,x,1])) for x = 1:L])
for t = 2:N_t
    bla = plot!([real(log_n_squared(uuu[1,x,t])) for x = 1:L])
end
display(bla)

# bli = plot([real(log_n_squared(uu[1,1,t])) for t = 1:16])
# for x = 2:N_x
#     bli = plot!([real(log_n_squared(uu[1,x,t])) for t = 1:16])
# end
# display(bli)

ble = plot([real(log_n_squared(uuu[2,x,1])) for x = 1:16])
for t = 2:N_t
    ble = plot!([real(log_n_squared(uuu[2,x,t])) for x = 1:16])
end
display(ble)
coeffs2grp(uu[2,5,16])
log_n_squared(uu[2,2,16])
log_U2(uu[2,3,16])
log_U2(uu[1,16,16])

log_n_squared(plaq(uu,1,1))
coeffs2grp(plaq(uu,16,16))
[imag(log(plaq(uu,x,t).a)) for x = 1:N_x, t = 1:N_t]
imag(log(plaq(uu,1,1).a))*16^2/π
plaq(uu,1,1).a

# top_charge_U2(uuu)

isapprox([coeffs2grp(plaq(uuu,x,t)) for x = 1:N_x, t = 1:N_t], [plaq(uuu,1,1).a*[1 0; 0 1] for x = 1:N_x, t = 1:N_t])
isapprox([coeffs2grp(uuu[1,x,t]) for x = 1:N_x-1, t = 1:N_t], [uuu[1,x,t].a*[1 0; 0 1] for x = 1:N_x-1, t = 1:N_t])


coeffs2grp(uu[1,1,16] * uu[2,2,16])# * adjoint(uu[1,1,1]) * adjoint(uu[2,1,16]))
# log_U2(uu[1,1,16] * uu[2,2,16] * adjoint(uu[1,1,1]) * adjoint(uu[2,1,16]))
log_U2(uu[1,3,16] * uu[2,4,16])
log_U2(uu[2,2,16])


# blu = plot([real(log_n_squared(uu[2,1,t])) for t = 1:16])
# for x = 2:N_x
#     blu = plot!([real(log_n_squared(uu[2,x,t])) for t = 1:16])
# end
# display(blu)


blab = scatter([imag(log_U2(uuu[2,x,16]).a)/π*N_x for x = 1:N_x], yticks = -L/2:L/2)
# blab = scatter!([x/L for x = 1:7])
# blab = scatter!(8:16, [x/L - 1 for x = 8:16])


blub = scatter([imag(log_U2(uuu[1,x,1]).a)/π for x = 1:N_x] .* L^2);
for t = 2:N_t
    blub = scatter!([imag(log_U2(uuu[1,x,t]).a)/π for x = 1:N_x] .* L^2)
end
blub = plot!(yticks = -32:2)
display(blub)

X = 1
blub_ar = [imag(log_U2(uuu[1,X,mod1(t+1,N_t)]).a)/π - imag(log_U2(uuu[1,X,t]).a)/π for t = 1:N_t]
sum(blub_ar[1:end-1])
blub_ar[1] * L^2



# UU = gaugefield_U2(16, 16, true);
# ρ = 0.2
# s_stout = []
# s_stout_mid_old = []
# s_stout_mid = []
# v = stout(UU,ρ)
# v_mid_old = stout_midpoint_old(UU,ρ)
# v_mid = stout_midpoint(UU,ρ)
# for i = 1:100
#     v = stout(v,ρ)
#     v_mid_old = stout_midpoint_old(v_mid_old,ρ)
#     v_mid = stout_midpoint(v_mid,ρ)
#     push!(s_stout, action(v,1))
#     push!(s_stout_mid_old, action(v_mid_old,1))
#     push!(s_stout_mid, action(v_mid,1))
# end

# plot(s_stout)
# plot!(s_stout_mid_old)
# plot!(s_stout_mid)


function insta_half_int(Q, N_x, N_t)
    U = gaugefield_U2(N_x, N_t, false)
    r1 = rand(3) .- 0.5
    r2 = rand(3) .- 0.5
    r3 = cross(r1,r2)
    c1 = complex(π/2 * r1/sqrt(sum(r1.^2)))
    c2 = complex(π/2 * r3/sqrt(sum(r3.^2)))
    # c2 = c1
    U[1,1:N_x-1,:] = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x-1, t = 1:N_t]
    U[1,N_x,:] = [exp(-im*Q*t*π/N_x/N_t) * exp_u2(coeffs_U2(0.0im, c1[1], c1[2], c1[3])) for t = 1:N_t]
    U[2,:,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1:N_x]
    # U[2,1:N_x>>1,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1:N_x>>1]
    # U[2,N_x>>1+1:N_x,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1+N_x>>1:N_x]
    return U
end

function insta_full_int(Q, N_x, N_t)
    U = gaugefield_U2(N_x, N_t, false)
    r1 = rand(3) .- 0.5
    r2 = rand(3) .- 0.5
    r3 = cross(r1,r2)
    c1 = complex(π * r1/sqrt(sum(r1.^2)))
    c2 = complex(π * r3/sqrt(sum(r3.^2)))
    # c2 = c1
    U[1,1:N_x-1,:] = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x-1, t = 1:N_t]
    U[1,N_x,:] = [exp(-im*Q*t*π/N_x/N_t) * exp_u2(coeffs_U2(0.0im, c1[1], c1[2], c1[3])) for t = 1:N_t]
    U[2,:,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1:N_x]
    # U[2,1:N_x>>1,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1:N_x>>1]
    # U[2,N_x>>1+1:N_x,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1+N_x>>1:N_x]
    return U
end

# iii = insta_half_int(1,L,L);
# isapprox([coeffs2grp(plaq(iii,x,t)) for x = 1:N_x, t = 1:N_t], [plaq(iii,1,1).a*[1 0; 0 1] for x = 1:N_x, t = 1:N_t])

action(insta_half_int(-1,L,L),1)
action(insta_full_int(4,L,L),1)
insta_action(1,2,L,L,4,2)

# log_U2(uuu[2,8,16])

r1 = rand(3) .- 0.5
r2 = rand(3) .- 0.5
r3 = cross(r1,r2)
r1 = r1/sqrt(sum(r1.^2))
r3 = r3/sqrt(sum(r3.^2))

r1 = π/2 .* r1
r3 = π/2 .* r3

r4 = cross(r1,r3)
r4./(π/2)^2

sum(r1.*[σ1,σ2,σ3]) * sum(r3.*[σ1,σ2,σ3])
im*sum(r4.*[σ1,σ2,σ3])

exp_u2(coeffs_U2(0.0im, complex(r1)...)) * exp_u2(coeffs_U2(0.0im, complex(r3)...))
coeffs2grp(exp_u2(coeffs_U2(0.0im, complex(r3)...)) * exp_u2(coeffs_U2(0.0im, complex(r1)...)) * exp_u2(coeffs_U2(0.0im, complex(-2 .*r4)...)))



bla = gaugefield_U2(L, L, true);
for i = 1:200 chess_metro!(bla, 0.1, 6.0, [0.0], "U2") end


top_charge_U2(stout(bla,200,0.1))
action_before = action(bla,1)
s = [action(bla,1)]
for i = 1:10000
    push!(s,action(insta_half_int(1,L,L).*bla,1))
end
histogram(s)

action_after = action(exp_u2.(log_U2.(bla) .+ log_U2.(insta_half_int(-1,L,L))),1)








λ1 = Complex.([0 1 0; 1 0 0; 0 0 0])
λ2 = Complex.([0 -im 0; im 0 0; 0 0 0])
λ3 = Complex.([1 0 0; 0 -1 0; 0 0 0])
λ4 = Complex.([0 0 1; 0 0 0; 1 0 0])
λ5 = Complex.([0 0 -im; 0 0 0; im 0 0])
λ6 = Complex.([0 0 0; 0 0 1; 0 1 0])
λ7 = Complex.([0 0 0; 0 0 -im; 0 im 0])
λ8 = Complex.([1 0 0; 0 1 0; 0 0 -2]) / sqrt(3)

Λ = [λ1, λ2, λ3, λ4, λ5, λ6, λ7, λ8]

λ0 = Complex.([1 0 0; 0 1 0; 0 0 1])
λ_zero = Complex.([0 0 0; 0 0 0; 0 0 0])

function ran_U3(ϵ)
    return exp(ϵ * im * sum((rand(8) .- 0.5).*Λ))
end

function insta_U3(N_x, N_t, Q)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    w = 0.0
    if mod(Q%3,3) == 1
        w = -2*π/sqrt(3)
    elseif mod(Q%3,3) == 2
        w = 2*π/sqrt(3)
    end
    U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((-im*t*w)/(N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp((im*x*w)/(N_x) * λ8) for x = 1:N_x]
    return U
end

function insta_U3_w(N_x, N_t, Q, z)
    w = sqrt(3) * 2 * pi * (z-Q/3)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((-im*t*w)/(N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp((im*x*w)/(N_x) * λ8) for x = 1:N_x]
    return U
end

function stout_U3(U, ρ)
    NX = size(U,2)
    NT = size(U,3)
    V = similar(U)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                Ω0 = ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])
                Z0 = -1/2*(adjoint(Ω0) - Ω0)
                V[μ,x,t] = exp(Z0) * U[μ,x,t]
                # V[μ,x,t] = exp_stout(ρ * staple(U,μ,x,t) * adjoint(U[μ,x,t])) * U[μ,x,t]
            end
        end
    end
    return V
end

function action_U3(U, β)
    NX = size(U,2)
    NT = size(U,3)
    S = 3*NX*NT   # later generalization: N_colour * NT * (NX)^d_s
    for t = 1:NT
        for x = 1:NX
            S -= real(tr(plaq(U,x,t)))
        end
    end
    return β*S/3    # later generalization: β*S/N_colour
end


function write_conf_U3(U, path)
    NX = size(U,2)
    NT = size(U,3)
    bla = []
    for μ = 1:2
        for x = 1:NX
            for t = 1:NT
                push!(bla, U[μ,x,t])
            end
        end
    end
    writedlm(path, bla, ',')
    return nothing
end

function read_config_U3(path)
    blub = readdlm(path, ',', ComplexF64)
    L = Int(sqrt(size(blub,1)/2))
    N = Int(sqrt(size(blub,2)))
    return [reshape(blub[t+(x-1)*L+(μ-1)*L^2,:], N, N) for μ = 1:2, x = 1:L, t = 1:L] 
end

function temp_gauge_U3(U)
    NX = size(U,2)
    NT = size(U,3)
    V = [convert.(ComplexF64, λ0) for μ = 1:2, x = 1:NX, t = 1:NT]
    Ω_slice = [convert.(ComplexF64, λ0) for x = 1:NX] 
    for t = 1:NT
        V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice,-1))
        Ω_slice = Ω_slice .* U[2,:,t]
    end
    V[2,:,NT] = Ω_slice
    return V
end

function max_gauge_U3(U)
    NX = size(U,2)
    NT = size(U,3)
    V = [convert.(ComplexF64, λ0) for μ = 1:2, x = 1:NX, t = 1:NT]
    Ω_slice = [convert.(ComplexF64, λ0)]
    for x = 1:NX-1
        next_el = last(Ω_slice) * U[1,x,1]
        push!(Ω_slice, next_el)
    end
    Ω_slice_copy = Ω_slice
    V[1,NX,1] = last(Ω_slice) * U[1,NX,1]
    for t = 2:NT
        Ω_slice = Ω_slice .* U[2,:,t-1]
        V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice, -1))
    end
    V[2,:,NT] = Ω_slice .* U[2,:,NT] .* adjoint.(Ω_slice_copy)
    return V
end

function Gell_coeffs(M::Matrix)
    x8 = -sqrt(3)/2*M[3,3]
    x3 = M[1,1] - x8/sqrt(3)
    x1 = real(M[2,1])
    x2 = imag(M[2,1])
    x4 = real(M[3,1])
    x5 = imag(M[3,1])
    x6 = real(M[3,2])
    x7 = imag(M[3,2])
    return [x1, x2, x3, x4, x5, x6, x7, x8]
end

function Gell_coeffs_log(M::Matrix{ComplexF64})
    phase_fac = det(M)^(1/3)
    A = M/phase_fac
    return Gell_coeffs(real.(-im*log(A)))
end

# conf_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\smeared_U3_conf_q_1.txt"
# write_conf_U3(V, conf_path)
# v = read_config_U3(conf_path);






L = 16
ρ = 0.1
U = [ran_U3(rand()) for μ = 1:2, x = 1:L, t = 1:L] .* insta_U3_w(L, L, 1, 2);
N_smear = 5*10^4
count = 0

# v = stout_midpoint(u,0.1);
V = stout_U3(U, ρ)
s = [action_U3(V,1)];
for i = 1:N_smear
    push!(s,action_U3(V,1))
    V = stout_U3(V,ρ)
    if i%Int(N_smear/100) == 0
        count += 1
        println("Smearing Progress U3: $count%")
    end
end

plot(s)
insta_action(1,3,L,L,1,1/3)
last(s)


# 🚧👷 Still under construction 👷🚧
function insta_third_int_U3(Q, N_x, N_t)
    U = [complex(convert.(Float64, λ0)) for μ = 1:2, x = 1:N_x, t = 1:N_t]
    r1 = rand(8) .- 0.5
    r2 = rand(8) .- 0.5
    # r3 = cross(r1,r2)
    c1 = complex(π/3 * r1/sqrt(sum(r1.^2)))
    c2 = complex(π/3 * r2/sqrt(sum(r2.^2)))
    # c2 = c1
    U[1,1:N_x-1,:] = [exp(-im*2*Q*t*π/3/N_x/N_t) * complex(λ0) for x = 1:N_x-1, t = 1:N_t]
    U[1,N_x,:] =     [exp(-im*2*Q*t*π/3/N_x/N_t) * exp(im*sum(r1.*Λ)) for t = 1:N_t]
    U[2,:,N_t] =     [exp(im*2*Q*x*π/3/N_x) * exp(im*sum(r2.*Λ)) for x = 1:N_x]
    # U[2,1:N_x>>1,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1:N_x>>1]
    # U[2,N_x>>1+1:N_x,N_t] = [exp(im*Q*x*π/N_x) * exp_u2(coeffs_U2(0.0im, c2[1], c2[2], c2[3])) for x = 1+N_x>>1:N_x]
    return U
end

insta_third_int_U3(1,L,L);

action_U3(insta_U3_w(L,L,1,3/3),1)
action_U3(insta_third_int_U3(1,L,L),1)
insta_action(1,3,L,L,1,1/3)

insta_action(1,2,16,16,5,2.5)/256





L = 16
N_x = L
N_t = L
heatmap([real(tr(plaq(v,x,t))) for x = 1:N_x, t = 1:N_t])#,xticks = 1:16, yticks = 1:16)
heatmap([real(tr(plaq(temp_gauge_U3(v),x,t))) for x = 1:N_x, t = 1:N_t])#,xticks = 1:16, yticks = 1:16)
heatmap([real(tr(plaq(max_gauge_U3(v),x,t))) for x = 1:N_x, t = 1:N_t])#,xticks = 1:16, yticks = 1:16)


# r = rand(8)
# M = sum(r.*Λ);
# isapprox(Gell_coeffs(M), r)


uu = temp_gauge_U3(v);
uuu = max_gauge_U3(v);


imag(log(plaq(uuu,1,1)[1,1]))/π*256
isapprox([plaq(uuu,x,t) for x = 1:N_x, t = 1:N_t], [exp(im*2*π/3/N_x/N_t) * λ0 for x = 1:N_x, t = 1:N_t])


let
    println("Absolute value squared of Gell-Mann-coefficients IN X DIRECTION:")
    for x = 1:N_x
        for t = 1:N_t
            bla = sum(Gell_coeffs_log(uuu[1,x,t]).^2)
            if bla > 10.0^(-10)
                println("x = $x, t = $t, |r|²: $bla")
            end
        end
    end
end

let
    println("Absolute value squared of Gell-Mann-coefficients IN T DIRECTION:")
    for x = 1:N_x
        for t = 1:N_t
            bla = sum(Gell_coeffs_log(uuu[2,x,t]).^2)
            if bla > 10.0^(-10)
                println("x = $x, t = $t, |r|²: $bla")
            end
        end
    end
end


sqrt(sum(Gell_coeffs_log(uuu[1,N_x,1]).^2))
sqrt(sum(Gell_coeffs_log(uuu[2,9,N_t]).^2))


Gell_coeffs_log(uuu[2,1,N_t])

let
    μ = 2
    x = 15
    t = N_t
    # uuu[μ,x,t] * uuu[μ,x,t] * uuu[μ,x,t]
    println(imag(log((uuu[μ,x,t] * uuu[μ,x,t] * uuu[μ,x,t])[1,1]))/π)
    # U_t(x,N_t) -> x/8*π for x = 1:N_x/2,   above it is -π + x/8*π
end

let
    μ = 1
    x = 16
    t = 1
    # uuu[μ,x,t] * uuu[μ,x,t] * uuu[μ,x,t]
    println(imag(log((uuu[μ,x,t] * uuu[μ,x,t] * uuu[μ,x,t])[1,1]))/π*L^2)
end


uuu[1,N_x,1] * uuu[2,1,1]    #* adjoint(uuu[1,N_x,2]) * adjoint(uuu[2,N_x,1])
uuu[1,N_x,1] * uuu[2,1,1] * adjoint(uuu[1,N_x,2]) * adjoint(uuu[2,1,1])


for x = 1:N_x
    μ = 2
    t = N_t
    println(transpose(Gell_coeffs_log(uuu[μ,x,t])))
end

for t = 1:N_t
    μ = 1
    x = N_x
    println(transpose(Gell_coeffs_log(uuu[μ,x,t])))
end


