using LinearAlgebra
using StatsBase 
using DelimitedFiles

include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")



################################################################################



####    Observable params    ####
for beta in [6.0]
    for L in [16]
        # for q_insta in [1]
        for N_stout in [10^5]  #[0, 10^2, 10^3, 10^4]

            # @assert 1==0 "Do you really want L = $L?"

            comment = "Cross check with Stephan using SUPER stout"
            
            ####    Update params   ####
            global group    = "U2"  # For "group" choose from: "SU2" or "U2"
            global β        = beta
            global N_t      = L
            global N_x      = L
            global ϵ        = 0.2   # Starting value for step size in generation of random elements
            global hot      = true  # true: hot start, false: cold start
            global read_last_config = false
            global N_metro  = 1                         # N_metro-many Metropolois sweeps followed by...
            global N_over   = 3                         # ...N_over-many overrelaxation sweeps, followed by...
            global N_insta  = 0                         # ...N_insta instanton updates
            global N_sepa   = 10
            global N_therm  = 500                       # For each therm.-sweep (N_metro+N_over+N_insta) sweeps will be performed
            # global N_meas   = Int( 100* (round(300 * (128/N_t)^2 / (N_metro+N_over+N_insta) /100, RoundNearestTiesAway))) 
            global N_meas   = 800
            # N_meas is similar to N_therm, but now we measure after (N_metro+N_over+N_insta) sweeps
            global N_stout_insta = 100
            global Q_start   = 0
            global Q_insta   = 0 # q_insta

            # if group != "SU2" && N_over != 0
            #     error("Overrelaxation is only implemented for SU(2) yet, so please set N_over to 0")
            # end
            if group != "U2" && N_insta != 0
                error("Instanton updates are only implemented for U(2) yet, so please set N_insta to 0")
            end

            ####    Measurement params    ####
            global n_stout   = N_stout      # 🟥🟥🟥
            global ρ         = 0.1
            loops     = [[1,1], [1,2], [2,1], [2,2]]#, [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
            # loops   = [[2^i,R] for i = 0:Int(log2(N_t)), R = 1:4:N_x ]


            ####    Accpetance and Progress    ####
            global counter        = 0
            global acc_metro      = [0.0]
            global acc_over       = [0.0]
            global acc_insta      = [0.0]
            global acc_wish       = 0.85
            # global acc_wish_insta = 0.45
            # global cool_deg       = 3.5


            ####    Handling directories    ####
            base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\"
            last_base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\"
            base_path = string(base_path, group)
            last_base_path = string(last_base_path, group)
            base_path = string(base_path,"_data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ")
            last_base_path = string(last_base_path,"_data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ")

            count_path = string(base_path, "\\sim_count.txt")


            # Have we already simulated with these parameters (excluding "loops")? If so,
            # make a new directory and keep track of how often via "sim_count.txt"
            if !isdir(base_path)    
                mkpath(base_path)
                writedlm(count_path, 1)
                base_path = string(base_path, "\\sim_count_1")
            else 
                sim_count = Int(readdlm(count_path)[1])
                sim_count += 1
                writedlm(count_path, sim_count)
                last_sim_count = sim_count -1
                last_base_path = string(base_path, "\\sim_count_$last_sim_count")
                base_path = string(base_path, "\\sim_count_$sim_count")
            end
            mkdir(base_path)
            
            # actions_path     = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\para_data\\actions_eps_$ϵ._beta_$β._L_$N_t._Nr_$run.txt"
            insta_delta_s_plus_path = string(base_path, "\\insta_delta_s_plus.txt") #⭕⭕⭕⭕⭕

            actions_path = string(base_path,"\\actions.txt")
            actions_clover_path = string(base_path,"\\actions_clover.txt")
            actions_unsm_path = string(base_path,"\\actions_unsm.txt")
            actions_clover_unsm_path = string(base_path,"\\actions_clover_unsm.txt")
            params_path = string(base_path,"\\params.txt") #
            acceptances_path = string(base_path,"\\acceptances.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
            last_conf_path   = string(base_path,"\\last_config.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\last_config_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
            top_charge_path = string(base_path,"\\top_charge.txt")
            top_charge_wil_path = string(base_path,"\\top_charge_wil.txt")
            corr_mat_paths = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
            mean_vals_path = string(base_path,"\\mean_vals.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
            mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")
            edge_loop_means_path = string(base_path,"\\edge_loop_means.txt")
            L_loop_means_path = string(base_path,"\\L_loop_means.txt")
            last_conf_path = string(last_base_path,"\\last_config.txt")
            rhomb_means_path = string(base_path,"\\rhomb_means.txt")
            half_rhomb_means_path = string(base_path,"\\half_rhomb_means.txt")


            # Write down the parameters of the simulation
            params = "Square Simulation with:
            N_t          = $N_t
            N_x          = $N_x 
            β            = $β
            hot          = $hot

            N_metro      = $N_metro
            N_over       = $N_over
            N_insta      = $N_insta
            N_stout_insta = $N_stout_insta
            N_therm      = $N_therm
            N_meas       = $N_meas
            N_sepa       = $N_sepa
            acc_wish     = $acc_wish
            Q_start      = $Q_start
            Q_insta      = $Q_insta

            n_stout      = $n_stout
            ρ            = $ρ
            loops        = $loops
            
            comment      = $comment"

            bla = open(params_path, "a")
            write(bla, params)
            close(bla)




            ################################################################################




            # U = gaugefield_SU2(N_t, N_x, hot)
            U = gaugefield(N_x, N_t, hot, group, "square")      # "lattice" manually set to "square"
            if read_last_config 
                U = read_last_conf(last_conf_path, N_t, N_x)    # Watch out ❗❗❗❗❗
            end

            for i = 1:N_therm
                chess_metro!(U,ϵ,β,acc_metro,group)
                ϵ *= sqrt(acc_metro[1] / acc_wish) # only update ϵ acc. to Metropolis
                global epsilon = deepcopy(ϵ)
                println("Acceptance: $(round(acc_metro[1],digits = 3)), ϵ: $epsilon ")
            end
            bla = open(params_path, "a")
            write(bla, "\n  ϵ = $epsilon")
            close(bla)

            U = insta_U2(N_x, N_t, -round(Int, top_charge_U2(U)) + Q_start) .* U
            for i = 1:50 chess_metro!(U,ϵ,β,acc_metro,group) end

            acc_metro    = [0.0]
            acc_over     = [0.0]
            acc_insta    = [0.0]

            insta_delta_s_plus = [0.0] #⭕⭕⭕⭕⭕

            for meas = 1:N_meas 
                if meas%(Int(N_meas/100)) == 1
                    println(" ")
                    println("We're already ", counter, "% deep in the simulation with N_t = $N_t, N_x = $N_x, β = $β and ϵ = $ϵ !")
                    # mywrite(last_conf_path, U, N_x, N_t)
                    # bla = open(last_conf_path, "a")
                    # write(bla, "Progress in simulation: $counter %")
                    # close(bla)
                    counter += 1
                end
                for sepa = 1:N_sepa
                    for metro = 1:N_metro
                        chess_metro!(U,ϵ,β,acc_metro,group)
                    end
                    for over = 1:N_over
                        chess_overrelax!(U,acc_over)
                    end
                    for j = 1:N_insta
                        # # insta_cool_update_U2!(U,ϵ,β,N_insta_cool,acc_insta)
                        # # insta_ran_cool_update_U2!(U,ϵ,β,N_insta_cool,cool_deg,acc_insta)
                        # insta_delta_s_plus[1] = insta_update_U2_log!(U,β,acc_insta,Q_insta)
                        insta_delta_s_plus[1] = insta_update_U2_log_bf_quick!(U,β,acc_insta,Q_insta) # ⭕⭕⭕⭕⭕
                        # # insta_flow_update_U2!(U, N_stout_insta, acc_insta, Q_insta)
                    end
                end

                U = proj2man.(U)

                
                # mywrite(insta_delta_s_plus_path, insta_delta_s_plus) # ⭕⭕⭕⭕⭕
                
                mywrite(acceptances_path, [acc_metro[1], acc_over[1], acc_insta[1]])
                mywrite(actions_unsm_path, action(U,β))
                mywrite(actions_clover_unsm_path, action_clover(U,β))

                V = stout(U, n_stout, ρ)
                
                mywrite(actions_path, action(V,β))
                mywrite(actions_clover_path, action_clover(V,β))
                mywrite(top_charge_path, top_charge_U2(V))
                mywrite(top_charge_wil_path, top_charge_U2_wil(V))
                # mywrite(mean_vals_path, results[2])
                
                # results = measure_RT_loops_corrs(U,loops,n_stout,ρ)
                # for t = 1:N_t
                #     mywrite(corr_mat_paths[t], results[1][:,:,t])
                # end

            end
            println(" ")
            println("We're done!")
        end
    end
end


# acceptances = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\N_t_32.N_x_32._beta_6.0\\n_stout_0._rho_0.12\\sim_count_1\\acceptances.txt")
# last(acceptances)/2/N_t/N_x/(N_metro+N_over)/(N_meas+N_therm)
# mod_acc = [(acceptances[i+1]-acceptances[i])/2/N_t/N_x/(N_metro+N_over) for i = 1:length(acceptances)-1]
# plot(mod_acc)

# loop_means = readdlm(mean_vals_mike_path)
# scatter([mean(loop_means[:,i]) for i = 1:length(loops)], yerror = [std(loop_means[:,i]) for i = 1:length(loops)])
# edge_means = readdlm(edge_loop_means_path)
# L_means = readdlm(L_loop_means_path)
# scatter!([17],[mean(edge_means)], yerror = [std(edge_means)])
# scatter!([18],[mean(L_means)], yerror = [std(L_means)])



#=
for L in [16,32,64,96]
# L = 16
    N_x = L
    N_t = L
    acc_chess = [0]
    acc_lex = [0]
    acc_ran = [0]
    β = 6.0
    ϵ_chess = 0.2
    ϵ_lex = 0.2
    ϵ_ran = 0.2
    acc_wish = 0.9
    n_meas = 10000
    counter = 0
    # for i = 1:100
    #     chess_metro_hex!(test_hex, 0.2)
    # end
    test_lex = gaugefield_SU2(N_x,N_t,false)

    test_chess = deepcopy(test_lex)
    test_ran = deepcopy(test_lex)
    actions_lex = []
    actions_chess = []
    actions_ran = []
    acceptances_lex = []
    acceptances_chess = []
    acceptances_ran = []

    for i = 1:200
        acc_copy_chess = acc_chess[1]
        acc_copy_lex = acc_lex[1]
        acc_copy_ran = acc_ran[1]
        # for j = 1:N_metro
        lexico_metro!(test_lex, ϵ_lex, β, acc_lex)
        chess_metro!(test_chess, ϵ_chess, β, acc_chess)
        ran_metro!(test_ran, ϵ_ran, β, acc_ran)
        # end
        # mywrite(acceptances_path, acc[1])
        ϵ_chess *= sqrt((acc_chess[1]-acc_copy_chess) /2/N_t/N_x / acc_wish)
        ϵ_lex *= sqrt((acc_lex[1]-acc_copy_lex) /2/N_t/N_x / acc_wish)
        ϵ_ran *= sqrt((acc_ran[1]-acc_copy_ran) /2/N_t/N_x / acc_wish)



        global ϵ_chess_write = deepcopy(ϵ_chess)
        global ϵ_lex_write = deepcopy(ϵ_lex)
        global ϵ_ran_write = deepcopy(ϵ_ran)
    end


    for i = 1:n_meas
        if i%(Int(n_meas/100)) == 1
            println(" ")
            println("We're already ", counter, "% deep in the SQUARE simulation with N_t = $N_t, N_x = $N_x and β = $β !")
        #     # mywrite(last_conf_path, U)
        #     # bla = open(last_conf_path, "a")
        #     # write(bla, "Progress in simulation: $counter %")
        #     # close(bla)
            counter += 1
        end
        lexico_metro!(test_lex, ϵ_lex, β, acc_lex)
        chess_metro!(test_chess, ϵ_chess, β, acc_chess)
        ran_metro!(test_ran, ϵ_ran, β, acc_ran)
        push!(actions_lex, action(test_lex, β))
        push!(actions_chess, action(test_chess, β))
        push!(actions_ran, action(test_ran, β))
        push!(acceptances_lex, acc_lex[1])
        push!(acceptances_chess, acc_chess[1])
        push!(acceptances_ran, acc_ran[1])
    end
    # image = plot(actions_chess[1:200], label = "chess")
    # image = plot!(actions_lex[1:200], label = "lexicog.")
    # image = plot!(actions_ran[1:200], label = "random")
    # image = plot!(title = "Action Time Series, Different Metrop. Updates, Hex.
    # N_t = N_x = $N_t,  β = $β",
    # xlabel = "MC time")
    # display(image)
    # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\actions_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.pdf")
    
    # auto_chess = []
    # auto_lex = []
    # auto_ran = []
    # counter = 0
    # for start = 50:10:n_meas-50
    #     if start%(Int(n_meas/50)) == 1
    #         # println(" ")
    #         println("Progress for autos: $counter%")
    #         # mywrite(last_conf_path, U)
    #         # bla = open(last_conf_path, "a")
    #         # write(bla, "Progress in simulation: $counter %")
    #         # close(bla)
    #         counter += 2
    #     end
    #     push!(auto_chess, auto_corr_time(actions_chess[start:n_meas]))
    #     push!(auto_lex, auto_corr_time(actions_lex[start:n_meas]))
    #     push!(auto_ran, auto_corr_time(actions_ran[start:n_meas]))
    # end
    # image = plot(50:10:n_meas-50, auto_chess, label = "chess")
    # image = plot!(50:10:n_meas-50, auto_lex, label = "lexicog.")
    # image = plot!(50:10:n_meas-50, auto_ran, label = "ran.")
    # image = plot!(title = "Int. Autoc. Time, Different Metrop. Updates, Hex.
    # N_t = N_x = $N_t,  β = $β,  ϵ = $ϵ",
    # xlabel = "MC time at which summation began")
    # display(image)

    # auto_chess = []
    # auto_lex = []
    # auto_ran = []
    # for i = 50:50:n_meas-50
    #     push!(auto_chess, auto_corr_norm(actions_chess,i))
    #     push!(auto_lex, auto_corr_norm(actions_lex,i))
    #     push!(auto_ran, auto_corr_norm(actions_ran,i))
    # end

    # image_auto = plot(50:50:n_meas-50, auto_chess, label = "chess")
    # image_auto = plot!(50:50:n_meas-50, auto_lex, label = "lexicog.")
    # image_auto = plot!(50:50:n_meas-50, auto_ran, label = "ran.")
    # image_auto = plot!(title = "Norm. Autocorr. Γ(t), Different Metrop. Updates, Hex.
    # N_t = N_x = $N_t,  β = $β",
    # xlabel = "MC time t")
    # display(image_auto)
    # mkpath("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β")
    # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\auto_norm_N_t_$N_t._N_x_$N_x.pdf")
    # println("chess: ", mean(auto_chess), " ± ", std(auto_chess)/sqrt(length(auto_chess)))
    # println("lexicog.: ", mean(auto_lex),  " ± ", std(auto_lex)/sqrt(length(auto_lex)))
    # println("random: ", mean(auto_ran),  " ± ", std(auto_ran)/sqrt(length(auto_ran)))
    # println(" ")

    # mkpath("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β")

    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\eps_chess_$N_t._N_x_$N_x.txt", ϵ_chess_write)
    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\eps_lex_$N_t._N_x_$N_x.txt", ϵ_lex_write)
    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\eps_ran_$N_t._N_x_$N_x.txt", ϵ_ran_write)

    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\actions_lex_N_t_$N_t._N_x_$N_x.txt", actions_lex)
    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\actions_chess_N_t_$N_t._N_x_$N_x.txt", actions_chess)
    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\actions_ran_N_t_$N_t._N_x_$N_x.txt", actions_ran)

    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\acceptances_lex_N_t_$N_t._N_x_$N_x.txt", acceptances_lex)
    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\acceptances_chess_N_t_$N_t._N_x_$N_x.txt", acceptances_chess)
    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\auto_corr\\beta_$β\\acceptances_ran_N_t_$N_t._N_x_$N_x.txt", acceptances_ran)

end



# N_x = 128
# N_t = 128
# β   = 8.0
# ϵ   = 0.2
# acc = [0]
# acc_wish = 0.9
# N_metro = 3
# N_over  = 1

# U   = gaugefield_SU2(N_x, N_t, true)

# N_therm = 200
# N_meas  = 1000

# for i = 1:N_therm
#     acc_copy = acc[1]
#     for j = 1:N_metro
#         chess_metro!(U,ϵ,β,acc)
#     end
#     for j = 1:N_over
#         chess_overrelax!(U,acc)
#     end
#     # mywrite(acceptances_path, acc[1])
#     ϵ *= sqrt((acc[1]-acc_copy)  /2/N_t/N_x/(N_metro+N_over) / acc_wish)
#     epsilon = deepcopy(ϵ)
#     println(epsilon)
# end

# loops = [[2,3]]
# loops_norm = []
# loops_mike = []
# loops_hand = []
# for i = 1:N_meas
#     for j = 1:N_metro
#         chess_metro!(U,ϵ,β,acc)
#     end
#     for j = 1:N_over
#         chess_overrelax!(U,acc)
#     end
#     push!(loops_norm, measure_RT_loops(U,loops,0,0.0))
#     push!(loops_mike, measure_RT_loops_mike(U,loops))
#     push!(loops_hand, mean(tr.([loop_2x3_square(U,x,t) for x = 1:N_x, t = 1:N_t])))
# end

# loops_norm

# b_size = Int(round(2*auto_corr_time([loops_norm[i][1] for i = 1:N_meas]) + 1, RoundUp))    
# jackknife([loops_norm[i][1] for i = 1:N_meas], b_size)
# b_size = Int(round(2*auto_corr_time([loops_mike[i][1] for i = 1:N_meas]) + 1, RoundUp))    
# jackknife([loops_mike[i][1] for i = 1:N_meas], b_size)
# b_size = Int(round(2*auto_corr_time(loops_hand) + 1, RoundUp))    
# jackknife(loops_hand, b_size)

# # Results for normal:   0.6057(4)
# # Results for mike:     0.6054(3)
# # Results for händisch: 0.6507(4)


# # println("Normal : ",  mean(loops_norm)[1])
# # println("Mike   : ",  mean(loops_mike)[1])
# # println("Hand   : ",  mean(loops_hand)[1])
=#
