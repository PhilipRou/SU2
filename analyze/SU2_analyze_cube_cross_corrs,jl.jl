include("SU2_analyze_head.jl")



β           = 8.0
N_t = N_x   = 8
smear_nums  = [0,1,3,7,15]
ρ           = 0.24
sim_count   = 1


base_path           = "D:\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\smear_nums_$smear_nums._rho_$ρ\\sim_count_$sim_count"
acceptances_path    = string(base_path,"\\acceptances.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
corr_mat_paths      = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["D:\\Physik Uni\\julia_projects\\SU2\\data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
mean_vals_path      = string(base_path,"\\mean_vals.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"

mean_vals = readdlm(mean_vals_path)
n_meas  = size(mean_vals,1)
n_smear = size(mean_vals,2)

evs      = Array{Float64}(undef, N_t, n_smear)
evs_errs = Array{Float64}(undef, N_t, n_smear)
for τ = 1:N_t
    # τ = 1
    raw_corrs = readdlm(corr_mat_paths[τ])
    corr_mats_for_b_size = Array{Float64}(undef,n_smear,n_smear,n_meas)
    for i = 1:N_meas
        corr_mats_for_b_size[:,:,i] = raw_corrs[(i-1)*n_smear+1:i*n_smear , :]
    end
    b_size    = maximum([round(Int, 2*auto_corr_time(corr_mats_for_b_size[s1,s2,:]), RoundNearestTiesAway) for s1 in eachindex(smear_nums), s2 in eachindex(smear_nums) ])
    corr_mats = [(raw_corrs[(i-1)*n_smear+1:i*n_smear , :] + transpose(raw_corrs[(i-1)*n_smear+1:i*n_smear , :])) / 2 for i = 1:n_meas]
    jack = jack_corr_mat_ev(corr_mats,b_size)
    evs[τ,:]      = reverse(jack[1])
    evs_errs[τ,:] = reverse(jack[2])
end

start_ind = 3
image_evs = plot_corrs(N_t, evs[:,start_ind], evs_errs[:,start_ind], "EV nr. $start_ind", palette(:default)[start_ind])
for smear = start_ind+1:n_smear
    image_evs = plot_corrs!(N_t, evs[:,smear], evs_errs[:,smear], "EV nr. $smear", image_evs, palette(:default)[smear])
end
image_evs = plot!(
    title = "EV's of corr. mats. of s_wil\n β = $β, L = $N_t, n_smear = $smear_nums, ρ = $ρ",
    xlabel = latexstring("\$ t \$")
)
display(image_evs)


