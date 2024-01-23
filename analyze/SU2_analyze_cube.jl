include("SU2_analyze_head.jl")




β           = 4.0
D           = 3
N_t = N_x   = 32
ϵ           = 0.2
n_stout     = 0
ρ           = 0.12
sim_count   = 3
loops     = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
N_metro   = 3
N_over    = 1

base_path           = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
params_path         = string(base_path,"\\params.txt") #
acceptances_path    = string(base_path,"\\acceptances.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
last_conf_path      = string(base_path,"\\last_config.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\last_config_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
corr_mat_paths      = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
mean_vals_path      = string(base_path,"\\mean_vals.txt") #"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
# mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")
# last_conf_path      = string(last_base_path,"\\last_config.txt")

mean_vals = readdlm(mean_vals_path)
n_meas = size(mean_vals,1)
L = length(loops)

acc = readdlm(acceptances_path)
# last(acc) /( (N_metro+N_over) * D * N_t * N_x^(D-1) * n_meas )
for i = 1:n_meas
    acc[i] /= (N_metro+N_over) * D * N_t * N_x^(D-1) * i
end

last(acc)
# plot(acc)



mean_vals = readdlm(mean_vals_path)

b_sizes = [Int(round(2*auto_corr_time(mean_vals[:,i])+1, RoundNearestTiesAway)) for i = 1:L]
jacks = [jackknife(mean_vals[:,i], b_sizes[i]) for i = 1:L]
loop_means = [jacks[i][1] for i = 1:L]
loop_errs = [jacks[i][2] for i = 1:L]

scatter(string.(loops), loop_means, yerror = loop_errs)



# For each t we want corr_mat_ar[t] to contain the time series of the correlation
# matrix at Euclidean time t
corr_mat_ar = Array{Array{Matrix{Float64}}}(undef, N_t)
for t = 1:N_t
    corr_mats = readdlm(corr_mat_paths[t])
    corr_mat_ar[t] = [corr_mats[Int((i-1)*L+1) : Int(i*L), 1:L] for i = 1:n_meas]
end

e_val_means = Matrix{Float64}(undef, L, N_t)
e_val_errs = Matrix{Float64}(undef, L, N_t)
for t = 1:N_t
    e_val_mat = Matrix{Float64}(undef, n_meas, L)
    for i = 1:n_meas
        e_val_mat[i,:] = abs.(eigen(corr_mat_ar[t][i]).values)
    end
    b_sizes = [Int(round(2*auto_corr_time(e_val_mat[:,i])+1, RoundNearestTiesAway)) for i = 1:L]
    jacks = [jackknife(e_val_mat[:,i], b_sizes[i]) for i = 1:L]
    e_val_jack_means = [jacks[i][1] for i = 1:L]
    e_val_jack_errs = [jacks[i][2] for i = 1:L]
    e_val_means[:,t] = e_val_jack_means
    e_val_errs[:,t] = e_val_jack_errs
end

image_ev = scatter(1:N_t, e_val_means[1,:], yerror = e_val_errs[1,:])
for l = 2:L
    image_ev = scatter!(1:N_t, e_val_means[l,:], yerror = e_val_errs[l,:])
end
display(image_ev)



scatter(1:10, e_val_means[16,1:10], yerror = e_val_errs[16,1:10])
