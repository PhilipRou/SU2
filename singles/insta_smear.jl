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
    # @assert 1==0 "We don't want to erase that which has already been obtained!"
    β   = 8.0
    L   = 32
    N_x = L
    N_t = L
    ρ   = 0.01
    N_therm   = 500
    N_meas    = 50
    N_sepa    = 100
    N_smear   = 10^5
    acc_wish  = 0.8
    ϵ         = 0.1
    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_data_3"
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





β       = 3.0
L       = 16
ρ       = 0.01
N_sepa  = 100
N_smear = 10^5

base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_data"
base_fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms"
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
N_meas  = 15

base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_data"
base_fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms"

for i = 15:15 # 1:N_meas
# i = 15
    smq_path  = string(base_path, "\\smq_$i.txt")
    sms_path  = string(base_path, "\\sms_$i.txt")
    smq_fig_path  = string(base_fig_path, "\\smq_$i.pdf")
    sms_fig_path  = string(base_fig_path, "\\sms_$i.pdf")

    smq  = readdlm(smq_path)
    sms  = readdlm(sms_path) ./ (β*L^2)

    Q = round(Int,last(smq))
    sms = sms .- minimum(sms)

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
        xticks = 0:100:900
    )
    display(image_sms)
    # savefig(sms_fig_path) # 🟥🟥🟥
end





# the measurements with final charge q=1:
top_one_runs = [2,5,6,8,12]

β       = 3.0
L       = 16
ρ       = 0.01
N_sepa  = 100
N_smear = 10^5

base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_1\\sms_data"
base_fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms"

cut = 100
sep = 100
plot_length = cut + length(cut+1:sep:N_smear) + 1
top_one_sms = Array{Float64}(undef, plot_length, length(top_one_runs))
for i = 1:length(top_one_runs)
    run = top_one_runs[i]
    sms_path = string(base_path, "\\sms_$run.txt")
    sms = readdlm(sms_path) / (β*L^2)
    top_one_sms[:,i] =  vcat(sms[1:cut], sms[cut+1:sep:end])
end

taus = ρ .* Array(1:N_smear+1)
tau_plot = vcat(taus[1:cut], taus[cut+1:sep:end])
s_window = 400:plot_length-300

let
    global image_top_ones = plot(
        tau_plot[s_window],
        top_one_sms[s_window,1],
        # legend = :false,
        label = "Run Nr. $(top_one_runs[1])",
        title  = latexstring("\$S/\\beta V \$ of Various Smeared Configs.\n 2D U(2), \$\\beta = $β, L = $L, \\rho = $ρ\$, final \$q=\\pm 1 \$"),
        xlabel = latexstring("Smearing Time \$\\tau\$"),
        # yaxis = :log,
        xticks = 0:100:900
    )
    for i = 2:length(top_one_runs)
        image_top_ones = plot!(
            tau_plot[s_window],
            top_one_sms[s_window, i],
            label = "Run Nr. $(top_one_runs[i])",
        )
    end
    display(image_top_ones)
end

image_top_ones = hline!(
    [insta_action(β,2,L,L,0,0.5) / (β*L^2)]
)





# the measurements with final q=2:
top_two_runs = [3,4,7,10,13,14]

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
s_window = 400:plot_length -50

let
    global image_top_twos = plot(
        tau_plot[s_window],
        top_two_sms[s_window,1],
        # legend = :false,
        label = "Run Nr. $(top_two_runs[1])",
        title  = latexstring("\$S/\\beta V \$ of Various Smeared Configs.\n 2D U(2), \$\\beta = $β, L = $L, \\rho = $ρ\$, final \$q=\\pm 2 \$"),
        xlabel = latexstring("Smearing Time \$\\tau\$"),
        # yaxis = :log,
        xticks = 0:100:900
    )
    for i = 2:length(top_two_runs)
        image_top_twos = plot!(
            tau_plot[s_window],
            top_two_sms[s_window, i],
            label = "Run Nr. $(top_two_runs[i])",
        )
    end
    display(image_top_twos)
end

image_top_ones = hline!(
    [insta_action(β,2,L,L,2,0.0) / (β*L^2)]
)








# the measurements with different final q's:
top_dif_runs = [2,3,15,1,9]
q_vals = [1,2,3,4,6]

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
s_window = 105:plot_length -50

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

let
    image_top_difs = hline!(
        [insta_action(β,2,L,L,1,0.5) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[1],
        label = latexstring("Insta: \$q = 1, z = 0.5\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,2,1) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[2],
        label = latexstring("Insta: \$q = 2, z = 1\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,3,1.5) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[3],
        label = latexstring("Insta: \$q = 3, z = 1.5\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,4,2) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[4],
        label = latexstring("Insta: \$q = 4, z = 2\$")
    )
    image_top_difs = hline!(
        [insta_action(β,2,L,L,6,3) / (β*L^2)],
        linestyle = :dash,
        color = palette(:default)[5],
        label = latexstring("Insta: \$q = 6, z = 3\$")
    )
    display(image_top_difs)
end

image_top_difs = hline!(
    [insta_action(β,2,L,L,1,0) / (β*L^2)],
    linestyle = :dash,
    color = :red,
    label = latexstring("Insta: \$q = 1, z = 0/1\$")
)
