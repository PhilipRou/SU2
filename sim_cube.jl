using LinearAlgebra
using StatsBase 
# using Plots
# using BenchmarkTools
using DelimitedFiles

include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_cube.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_cube.jl")



################################################################################



####    Basic params    ####
D           = 3
N_t         = 8
N_x         = 8
β           = 6.0
hot         = true
ϵ           = 0.2
acc_wish    = 0.85


####    Update params   ####
N_metro = 3              # N_metro-many Metropolois sweeps followed by...
N_over  = 1              # ...N_over-many overrelaxation sweeps will be performed...
N_therm = Int(200/(N_metro+N_over))  # ...for N_therm times,
N_meas  = Int(800/(N_metro+N_over))  # ...plus N_meas times,
#                               # i.e. (N_therm + N_meas) ⋅ (N_metro + N_over) sweeps in total

####    We've gotta track acceptance somehow...    ####
acc = [0]

####    Observable params   ####
n_stout = 0
for β in [4.0]
# for L = 32:32:128
# L = 32
    ####    Update params   ####
    N_t = N_x = 32
    N_metro   = 3                           # N_metro-many Metropolois sweeps followed by...
    N_over    = 1                           # ...N_over-many overrelaxation sweeps will be performed...
    N_therm   = Int(200/(N_metro+N_over))   # ...for N_therm times,
    N_meas    = Int(2400/(N_metro+N_over))
    # N_meas    = Int( 100* (round(3200 * (128/N_t)^2 / (N_metro+N_over) /100, RoundNearestTiesAway))) # ...plus N_meas times,
    #                                                       # i.e. (N_therm + N_meas) ⋅ (N_metro + N_over) sweeps in total
    ρ         = 0.12
    loops     = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    # loops   = [[2^i,R] for i = 0:Int(log2(N_t)), R = 1:4:N_x ]
    # loops   = [[N_t, 0]]    # Polyakov


    ####    Dirty tricks    ####
    counter = 0
    acc = [0]


    ####    Handling directories    ####
    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ"
    last_base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ"
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
    params_path = string(base_path,"\\params.txt") #
    acceptances_path = string(base_path,"\\acceptances.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
    last_conf_path   = string(base_path,"\\last_config.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\last_config_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
    corr_mat_paths = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
    mean_vals_path = string(base_path,"\\mean_vals.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")
    last_conf_path = string(last_base_path,"\\last_config.txt")


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

    n_stout = $n_stout
    ρ       = $ρ
    loops   = $loops"

    # writedlm(params_path,params)
    # mkpath(params_path)
    bla = open(params_path, "a")
    write(bla, params)
    close(bla)





    ################################################################################





    U = gaugefield_SU2_D(N_t, N_x, D, true)
    # U = gaugefield_SU2(N_t, N_x, hot)
    # U = read_last_conf(last_conf_path, N_t, N_x)   # Watch out ❗❗❗❗❗

    for i = 1:N_therm
        acc_copy = acc[1]
        for j = 1:N_metro
            chess_metro_3d!(U,ϵ,β,acc)
        end
        for j = 1:N_over
            chess_overrelax_3d!(U,acc)
        end
        ϵ *= sqrt((acc[1]-acc_copy)  /3/N_t/N_x/(N_metro+N_over) / acc_wish)
    end

    for i = 1:N_meas 
        if i%(Int(N_meas/100)) == 1
            println(" ")
            println("We're already ", counter, "% deep in the $D-dim. simulation with N_t = $N_t, N_x = $N_x, β = $β and ϵ = $ϵ !")
            # mywrite(last_conf_path, U)
            # bla = open(last_conf_path, "a")
            # write(bla, "Progress in simulation: $counter %")
            # close(bla)
            counter += 1
        end
        for j = 1:N_metro
            chess_metro_3d!(U,ϵ,β,acc)
        end
        for j = 1:N_over
            chess_overrelax_3d!(U,acc)
        end

        results = measure_loops_corrs_3d(U, loops)
        corr_mats = (results[1] .+ results[2] .+ results[3]) ./ 3 # Only possible if N_t = N_x
        mean_vals = (results[4] .+ results[5] .+ results[6]) ./ 3 # Only possible if N_t = N_x
        mywrite(acceptances_path, acc[1])
        mywrite(mean_vals_path, mean_vals)
        for t = 1:N_t
            mywrite(corr_mat_paths[t], corr_mats[:,:,t])
        end


        # mywrite(mean_vals_path, loop_means)
        # mywrite(mean_vals_mike_path, loop_means_mike)
    end
    println(" ")
    println("We're done!")
end
# end
