using LinearAlgebra
using StatsBase 
using DelimitedFiles

include("D:\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("D:\\Physik Uni\\julia_projects\\SU2\\observables\\observables_cube.jl")
include("D:\\Physik Uni\\julia_projects\\SU2\\updates\\updates_cube.jl")



################################################################################



####    Observable params   ####    
for β in [8.0]
    # for N_stout in [7]
    comment = "Serious run for plaq- and clover self-correlation (\"loops\" unneccessary)"
    ####    Update params   ####
    global N_t = N_x = 8
    global hot       = true

    global N_metro   = 1        # N_metro-many Metropolois sweeps followed by...
    global N_over    = 3        # ...N_over-many overrelaxation sweeps will be performed...
    global N_therm   = 100      # ...for N_therm times,
    global N_meas    = 100
    global N_sepa    = 10
    global ϵ         = 0.2
    global acc_wish  = 0.85
    # N_meas    = Int( 100* (round(3200 * (128/N_t)^2 / (N_metro+N_over) /100, RoundNearestTiesAway))) # ...plus N_meas times,
    #                                                       # i.e. (N_therm + N_meas) ⋅ (N_metro + N_over) sweeps in total
    # global n_stout   = 7
    global ρ         = 0.24
    global smear_nums= [0,1,3,7,15]
    # global loops     = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    # global loops     = [[1,1], [2,2], [3,3], [4,4], [5,5], [6,6]]
    global loops     = [[1,1]]
    # loops   = [[2^i,R] for i = 0:Int(log2(N_t)), R = 1:4:N_x ]
    # loops   = [[N_t, 0]]    # Polyakov
    
    global counter = 0
    global acc = [0.0]


    ####    Handling directories    ####
    base_path = "D:\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\smear_nums_$smear_nums._rho_$ρ"
    # last_base_path = "D:\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\smear_nums_$smear_nums._rho_$ρ"
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
    
    # actions_path     = "D:\\Physik Uni\\julia_projects\\SU2_data\\para_data\\actions_eps_$ϵ._beta_$β._L_$N_t._Nr_$run.txt"
    params_path = string(base_path,"\\params.txt") #
    acceptances_path = string(base_path,"\\acceptances.txt") #"D:\\Physik Uni\\julia_projects\\SU2_data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
    last_conf_path   = string(base_path,"\\last_config.txt") #"D:\\Physik Uni\\julia_projects\\SU2_data\\last_config_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
    corr_mat_paths = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["D:\\Physik Uni\\julia_projects\\SU2_data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
    mean_vals_path = string(base_path,"\\mean_vals.txt") #"D:\\Physik Uni\\julia_projects\\SU2_data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
    corr_mat_clover_paths = [string(base_path,"\\corrs_clover_t_$t.txt") for t = 1:N_t] #["D:\\Physik Uni\\julia_projects\\SU2_data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
    mean_vals_clover_path = string(base_path,"\\mean_vals_clover.txt") #"D:\\Physik Uni\\julia_projects\\SU2_data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
    # mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")
    # last_conf_path = string(last_base_path,"\\last_config.txt")


    params = "3d Simulation with:
    N_t     = $N_t
    N_x     = $N_x 
    β       = $β
    hot     = $hot
    ϵ       = $ϵ

    N_metro = $N_metro
    N_over  = $N_over
    N_therm = $N_therm
    N_meas  = $N_meas
    N_sepa  = $N_sepa

    smear_nums = $smear_nums
    ρ          = $ρ
    loops      = $loops
    
    comment = $comment"

    # writedlm(params_path,params)
    # mkpath(params_path)
    bla = open(params_path, "a")
    write(bla, params)
    close(bla)





    ################################################################################





    U = gaugefield_SU2_cube(N_x, N_t, hot)

    for i = 1:N_therm
        chess_metro_cube!(U,ϵ,β,acc)
        ϵ *= sqrt(acc[1] / acc_wish)
        println("Acceptance: $(round(acc[1],digits = 3)), ϵ: $ϵ ")
    end

    for meas = 1:N_meas
        if meas%(Int(N_meas/100)) == 1
            println(" ")
            println("We're already ", counter, "% deep in the 3-dim. simulation with N_x = $N_x, N_t = $N_t, β = $β and ϵ = $ϵ !")
            counter += 1
        end
        for sepa = 1:N_sepa
            for metro = 1:N_metro
                chess_metro_cube!(U,ϵ,β,acc)
            end
            for over = 1:N_over
                chess_overrelax_cube!(U)
            end
        end
        
        U = proj2man.(U)
        parity_hit_cube!(U)

        # results = measure_RT_loops_corrs_cube(U, loops, n_stout, ρ)
        # results = measure_RT_loops_corrs_cube_selfonly(U, loops, n_stout, ρ)
        # results_plaq = measure_plaq_12(U, n_stout, ρ)
        results = measure_s_wil(U, smear_nums, ρ)
        # results_clover = measure_clover(U, n_stout, ρ)
        corr_mats = results[1]
        mean_vals = results[2]
        # corr_mats_plaq = results_plaq[1]
        # mean_vals_plaq = results_plaq[2]
        # corr_mats_clover = results_clover[1]
        # mean_vals_clover = results_clover[2]
        mywrite(acceptances_path, acc[1])
        mywrite(mean_vals_path, mean_vals)
        # mywrite(mean_vals_clover_path, mean_vals_clover)
        for t = 1:N_t
            mywrite(corr_mat_paths[t], corr_mats[:,:,t])
            # mywrite(corr_mat_paths[t], corr_mats[:,t])
            # mywrite(corr_mat_paths[t], corr_mats_plaq[t])
            # mywrite(corr_mat_clover_paths[t], corr_mats_clover[t])
        end
        # mywrite(mean_vals_path, loop_means)
        # mywrite(mean_vals_mike_path, loop_means_mike)
    end
    println(" ")
    println("We're done!")
end
# end

