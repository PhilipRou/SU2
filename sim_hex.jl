using LinearAlgebra
using StatsBase
using DelimitedFiles

include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\observables\\observables_hex.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\updates\\updates_hex.jl")



################################################################################



####    Observable params    ####
for beta in [12] # [2.0,4,0,6.0,8.0]  
    for L = 32:32:32
        #beta = 3.0
        # L = 32
        ####    Update params   ####
        global group    = "U2"   # For "group" choose from: "SU2" or "U2"
        global β        = beta
        global N_t      = 2*L
        global N_x      = L
        global ϵ        = 0.2   # Step size for generation of random elements
        global hot      = true  # true: hot start, false: cold start
        global read_last_config = false
        global N_metro  = 3                         # N_metro-many Metropolois sweeps followed by...
        global N_over   = 0                         # ...N_over-many overrelaxation sweeps, followed by...
        global N_insta  = 1                         # ...N_insta insta_cool updates, using...
        global N_insta_cool = 1                     # ...N_insta_cool cooling steps
        global N_therm  = Int(500/(N_metro+N_over+N_insta))   # For each therm.-sweep (N_metro+N_over+N_insta) sweeps will be performed
        global N_meas   = Int( 100* (round(1600 * (2*128^2/N_t/N_x) / (N_metro+N_over+N_insta) /100, RoundNearestTiesAway))) 
        # N_meas is similar to N_therm, but now we measure after (N_metro+N_over+N_insta) sweeps

        if group != "SU2" && N_over != 0
            error("Overrelaxation is only implemented for SU(2) yet, so please set N_over to 0")
        end
        if group != "U2" && N_insta != 0
            error("Instanton updates are only implemented for U(2) yet, so please set N_insta to 0")
        end

        ####    Measurement params    ####
        global n_stout   = 0
        global ρ         = 0.12
        loops     = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
        # loops = [[1,1]]
        # loops   = [[2^i,R] for i = 0:Int(log2(N_t)), R = 1:4:N_x ]
        # loops   = [[N_t, 0]]    # Polyakov


        ####    Accpetance and Progress    ####
        global counter      = 0
        global acc_metro    = [0]
        global acc_over     = [0]
        global acc_insta    = [0]
        global acc_wish     = 0.9


        ####    Handling directories    ####
        base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\"
        last_base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\"
        base_path = string(base_path, group)
        last_base_path = string(last_base_path, group)
        base_path = string(base_path,"_data\\hex_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ")
        last_base_path = string(last_base_path,"_data\\hex_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ")

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
        
        # actions_path     = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\para_data\\actions_eps_$ϵ._beta_$β._L_$N_t._Nr_$run.txt"
        actions_path = string(base_path,"\\actions.txt")        
        params_path = string(base_path,"\\params.txt") #
        acceptances_path = string(base_path,"\\acceptances.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
        last_conf_path   = string(base_path,"\\last_config.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\last_config_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
        top_charge_path = string(base_path,"\\top_charge.txt")
        corr_mat_paths = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
        mean_vals_path = string(base_path,"\\mean_vals.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
        # mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")
        rhomb_half_loop_means_path = string(base_path,"\\rhomb_half_loop_means.txt")
        edge_loop_means_path = string(base_path,"\\edge_loop_means.txt")
        L_loop_means_path = string(base_path,"\\L_loop_means.txt")
        rhomb_loop_means_path = string(base_path,"\\rhomb_loop_means.txt")
        # last_conf_path = string(last_base_path,"\\last_config.txt")
        
        
        # Write down the parameters of the simulation
        params = "Hexagonal Simulation with:
        N_t          = $N_t
        N_x          = $N_x 
        β            = $β
        hot          = $hot

        N_metro      = $N_metro
        N_over       = $N_over
        N_insta      = $N_insta
        N_insta_cool = $N_insta_cool
        N_therm      = $N_therm
        N_meas       = $N_meas
        acc_wish     = $acc_wish

        n_stout      = $n_stout
        ρ            = $ρ
        loops        = $loops"

        bla = open(params_path, "a")
        write(bla, params)
        close(bla)




        ################################################################################


        

        # U = gaugefield_SU2(N_t, N_x, true)
        H = gaugefield(N_x, N_t, hot, group, "hexagon")      # "lattice" manually set to "hexagon"
        legal_coords = half_chess_coords(N_x,N_t)

        for i = 1:N_therm
            acc_copy = acc_metro[1]
            for j = 1:N_metro
                chess_metro_hex!(H,ϵ,β,acc_metro,group)
            end
            for j = 1:N_over
                chess_overrelax_hex!(H,acc_over)
            end
            for j = 1:N_insta
                insta_cool_update_U2_hex!(H,ϵ,β,N_insta_cool,acc_insta)
            end
            # mywrite(acceptances_path, acc[1])
            ϵ *= sqrt((acc_metro[1]-acc_copy)  /1.5/N_t/N_x/N_metro / acc_wish)
            global epsilon = deepcopy(ϵ)
            println(epsilon)
        end
        bla = open(params_path, "a")
        write(bla, "\n  ϵ = $epsilon")
        close(bla)

        # actions     = Vector{Float64}(undef, N_meas)
        # actions_hex = Vector{Float64}(undef, N_meas)
        # means_1x1   = Vector{Float64}(undef, N_meas)
        # means_1x2   = Vector{Float64}(undef, N_meas)
        # means_2x1   = Vector{Float64}(undef, N_meas)
        # means_2x2   = Vector{Float64}(undef, N_meas)

        acc_metro    = [0]
        acc_over     = [0]
        acc_insta    = [0]

        for i = 1:N_meas
            if i%(Int(N_meas/100)) == 1
                println(" ")
                println("We're already ", counter, "% deep in the HEXAGONAL simulation with N_t = $N_t, N_x = $N_x, β = $β and ϵ = $ϵ !")
                # mywrite(last_conf_path, U)
                # bla = open(last_conf_path, "a")
                # write(bla, "Progress in simulation: $counter %")
                # close(bla)
                counter += 1
            end
            for j = 1:N_metro
                chess_metro_hex_comp!(H,ϵ,β,acc,group)  #❗❗❗❗❗
            end
            for j = 1:N_over
                chess_overrelax_hex!(H,acc)
            end
            for j = 1:N_insta
                insta_cool_update_U2_hex!(H,ϵ,β,N_insta_cool,acc_insta)
            end

            # loops_summed_1x1 = 0.0
            # loops_summed_1x2 = 0.0
            # loops_summed_2x1 = 0.0
            # loops_summed_2x2 = 0.0
            # for x = 1:N_x
            #     for t = 1:N_t
            #         loops_summed_1x1 += tr(hexplaq(H,t,x))
            #         loops_summed_1x2 += tr(loop_1x2_hex(H,t,x))
            #         loops_summed_2x1 += tr(loop_2x1_hex(H,t,x))
            #         loops_summed_2x2 += tr(loop_2x2_hex(H,t,x))
            #     end
            # end

            # means_1x1[i] = loops_summed_1x1/N_t/N_x
            # means_1x2[i] = loops_summed_1x2/N_t/N_x
            # means_2x1[i] = loops_summed_2x1/N_t/N_x
            # means_2x2[i] = loops_summed_2x2/N_t/N_x
            # # actions[i] = action(U)
            # actions_hex[i] = action_hex(H)


            mywrite(acceptances_path, [acc_metro[1], acc[1]])
            mywrite(actions_path, action_hex_comp(H,β))         #❗❗❗❗❗
            mywrite(top_charge_path, top_charge_U2_hex_comp(H)) #❗❗❗❗❗
            loop_means = measure_RT_loops_hex(H,loops,legal_coords)#,n_stout,ρ)
            mywrite(mean_vals_path, loop_means)
            edge_loop_mean = sum([tr(edge_loop_hex(H,coord[1],coord[2])) for coord in half_chess_coords(N_x,N_t)])/N_x/(0.5*N_t)
            # rhomb_half_loop_mean = sum([tr(rhomb_half_loop(H,coord[1],coord[2])) for coord in half_chess_coords(N_x,N_t)])/N_x/(0.5*N_t)
            L_loop_mean = sum([tr(L_loop_hex(H,coord[1],coord[2])) for coord in half_chess_coords(N_x,N_t)])/N_x/(0.5*N_t)
            # rhomb_loop_mean = sum([tr(rhomb_loop(H,coord[1],coord[2])) for coord in half_chess_coords(N_x,N_t)])/N_x/(0.5*N_t)
            mywrite(edge_loop_means_path, edge_loop_mean)
            # mywrite(rhomb_half_loop_means_path, rhomb_half_loop_mean)
            mywrite(L_loop_means_path, L_loop_mean)
            # mywrite(rhomb_loop_means_path, rhomb_loop_mean)
        end
        # push!(actions_hex, action_hex(H))


        # image = plot(1:N_meas, actions, label = "square")
        # image = plot!(1:N_meas, actions_hex, label = "hexagonal")
        # image = plot!(mean(actions).*ones(N_meas), color = :cornflowerblue, ribbon = std(actions), label = :false)
        # image = plot!(mean(actions_hex).*ones(N_meas), color = :orange, ribbon = std(actions), label = :false)
        # image = plot!(legend = :left)

        # actions_hex_double = 2 .*actions_hex
        # # image = plot!(actions_hex_double, color = :red, label = "hexagonal × 2")
        # # image = plot!(mean(actions_hex_double).*ones(N_meas), ribbon = std(actions_hex_double), color = :red, label = :false)

        # display(image)


        # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\hex_data\\actions_hex_$N_t.txt", actions_hex)
        # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\hex_data\\means_1x1_$N_t.txt", means_1x1)
        # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\hex_data\\means_1x2_$N_t.txt", means_1x2)
        # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\hex_data\\means_2x1_$N_t.txt", means_2x1)
        # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\hex_data\\means_2x2_$N_t.txt", means_2x2)
        
        println(" ")
        println("We're done!")
    end
end






# acceptances = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\beta_6.0\\N_t_32.N_x_32\\n_stout_0._rho_0.12\\sim_count_1\\acceptances.txt")
# last(acceptances)/1.5/N_t/N_x/(N_metro+N_over)/(N_meas+N_therm)
# mod_acc = [(acceptances[i+1]-acceptances[i])/0.75/2/N_t/N_x/(N_metro+N_over) for i = 1:length(acceptances)-1]
# plot(mod_acc)






#=
for L in [16,32,64,96]
# L = 16
    N_x = L
    N_t = 2*L
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
    test_lex = hexfield_SU2(N_x,N_t,false)

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
        lexico_metro_hex!(test_lex, ϵ_lex, β, acc_lex)
        chess_metro_hex!(test_chess, ϵ_chess, β, acc_chess)
        ran_metro_hex!(test_ran, ϵ_ran, β, acc_ran)
        # end
        # mywrite(acceptances_path, acc[1])
        ϵ_chess *= sqrt((acc_chess[1]-acc_copy_chess) /2/0.75/N_t/N_x / acc_wish)
        ϵ_lex *= sqrt((acc_lex[1]-acc_copy_lex) /2/0.75/N_t/N_x / acc_wish)
        ϵ_ran *= sqrt((acc_ran[1]-acc_copy_ran) /2/0.75/N_t/N_x / acc_wish)


        global ϵ_chess_write = deepcopy(ϵ_chess)
        global ϵ_lex_write = deepcopy(ϵ_lex)
        global ϵ_ran_write = deepcopy(ϵ_ran)
    end


    for i = 1:n_meas
        if i%(Int(n_meas/100)) == 1
            println(" ")
            println("We're already ", counter, "% deep in the HEXAGONAL simulation with N_t = $N_t, N_x = $N_x and β = $β !")
        #     # mywrite(last_conf_path, U)
        #     # bla = open(last_conf_path, "a")
        #     # write(bla, "Progress in simulation: $counter %")
        #     # close(bla)
            counter += 1
        end
        lexico_metro_hex!(test_lex, ϵ_lex, β, acc_lex)
        chess_metro_hex!(test_chess, ϵ_chess, β, acc_chess)
        ran_metro_hex!(test_ran, ϵ_ran, β, acc_ran)
        push!(actions_lex, action_hex(test_lex, β))
        push!(actions_chess, action_hex(test_chess, β))
        push!(actions_ran, action_hex(test_ran, β))
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
    # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\actions_beta_$β._N_t_$N_t._N_x_$N_x._epsilon_$ϵ.pdf")
    
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
    # mkpath("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β")
    # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\auto_norm_N_t_$N_t._N_x_$N_x.pdf")
    # println("chess: ", mean(auto_chess), " ± ", std(auto_chess)/sqrt(length(auto_chess)))
    # println("lexicog.: ", mean(auto_lex),  " ± ", std(auto_lex)/sqrt(length(auto_lex)))
    # println("random: ", mean(auto_ran),  " ± ", std(auto_ran)/sqrt(length(auto_ran)))
    # println(" ")

    mkpath("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β")

    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\eps_chess_$N_t._N_x_$N_x.txt", ϵ_chess_write)
    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\eps_lex_$N_t._N_x_$N_x.txt", ϵ_lex_write)
    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\eps_ran_$N_t._N_x_$N_x.txt", ϵ_ran_write)

    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\actions_lex_N_t_$N_t._N_x_$N_x.txt", actions_lex)
    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\actions_chess_N_t_$N_t._N_x_$N_x.txt", actions_chess)
    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\actions_ran_N_t_$N_t._N_x_$N_x.txt", actions_ran)

    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\acceptances_lex_N_t_$N_t._N_x_$N_x.txt", acceptances_lex)
    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\acceptances_chess_N_t_$N_t._N_x_$N_x.txt", acceptances_chess)
    writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\acceptances_ran_N_t_$N_t._N_x_$N_x.txt", acceptances_ran)

    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\auto_norm_lex_N_t_$N_t._N_x_$N_x.txt", auto_lex)
    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\auto_norm_chess_N_t_$N_t._N_x_$N_x.txt", auto_chess)
    # writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2_data\\hex_data\\auto_corr\\beta_$β\\auto_norm_ran_N_t_$N_t._N_x_$N_x.txt", auto_ran)
end
=#



