include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")



β           = 8.0
N_t = N_x   = 32
smear_nums  = [0,1,3,7,15]
ρ           = 0.24
sim_count   = 2

base_path           = "D:\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\smear_nums_$smear_nums._rho_$ρ\\sim_count_$sim_count"
acceptances_path    = string(base_path,"\\acceptances.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
corr_mat_paths      = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["D:\\Physik Uni\\julia_projects\\SU2\\data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
mean_vals_path      = string(base_path,"\\mean_vals.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"

mean_vals = readdlm(mean_vals_path)
n_meas  = size(mean_vals,1)
n_smear = size(mean_vals,2)

corr_mats_list = []
corr_mats_small_list = []
for t = 1:N_t
    raw_corrs = readdlm(corr_mat_paths[t])
    corr_mats       = [(raw_corrs[(i-1)*n_smear+1:i*n_smear , :] + transpose(raw_corrs[(i-1)*n_smear+1:i*n_smear , :])) / 2 for i = 1:n_meas]
    corr_mats_small = [(raw_corrs[(i-1)*n_smear+2:i*n_smear , 2:n_smear] + transpose(raw_corrs[(i-1)*n_smear+2:i*n_smear , 2:n_smear])) / 2 for i = 1:n_meas]
    push!(corr_mats_list, corr_mats)
    push!(corr_mats_small_list, corr_mats_small)
end




b_sizes = [Int(round(2*auto_corr_time(mean_vals[:,i])+1, RoundNearestTiesAway)) for i = 1:n_smear]
jacks = [jackknife(mean_vals[:,i], b_sizes[i]) for i = 1:n_smear]
loop_means = [jacks[i][1] for i = 1:n_smear]
loop_errs = [jacks[i][2] for i = 1:n_smear]

scatter(smear_nums, loop_means, yerror = loop_errs, label = "Smeared s_wil loop exp. values", xticks = smear_nums)






small    = true
evs      = Array{Float64}(undef, N_t, n_smear)
evs_errs = Array{Float64}(undef, N_t, n_smear)
if small
    evs      = Array{Float64}(undef, N_t, n_smear-1)
    evs_errs = Array{Float64}(undef, N_t, n_smear-1)
end
for τ = 1:N_t
    # τ = 1
    corr_mats = corr_mats_list[τ]
    if small
        corr_mats = corr_mats_small_list[τ]
    end
    L = size(corr_mats[1],2)
    corr_mats_for_b_size = Array{Float64}(undef,L,L,n_meas)
    for i = 1:n_meas
        corr_mats_for_b_size[:,:,i] = corr_mats[i]
    end
    b_size = maximum([round(Int, 2*auto_corr_time(corr_mats_for_b_size[s1,s2,:]), RoundNearestTiesAway) for s1 = 1:L, s2 = 1:L ])
    jack = jack_corr_mat_ev(corr_mats,b_size)
    evs[τ,:]      = reverse(jack[1])
    evs_errs[τ,:] = reverse(jack[2])
end

start_ind = 2
end_ind   = n_smear
image_evs = plot_corrs(N_t, evs[:,start_ind], evs_errs[:,start_ind], "EV nr. $start_ind", palette(:default)[start_ind])
for smear = start_ind+1:end_ind
    image_evs = plot_corrs!(N_t, evs[:,smear], evs_errs[:,smear], "EV nr. $smear", image_evs, palette(:default)[smear])
end
image_evs = plot!(
    title = "EV's of FULL corr. mats. of s_wil\n β = $β, L = $N_t, n_smear = $smear_nums, ρ = $ρ",
    xlabel = latexstring("\$ t \$"),
    right_margin = 10mm,
    yaxis = :log,
    ylims = (1e-11, 1e-3)
)
display(image_evs)





small         = true
conn_evs      = Array{Float64}(undef, N_t, n_smear)
conn_evs_errs = Array{Float64}(undef, N_t, n_smear)
if small
    conn_evs      = Array{Float64}(undef, N_t, n_smear-1)
    conn_evs_errs = Array{Float64}(undef, N_t, n_smear-1)
end
for τ = 1:N_t
    # τ = 1
    corr_mats = corr_mats_list[τ]
    if small
        corr_mats = corr_mats_small_list[τ]
    end
    L = size(corr_mats[1],2)
    corr_mats_for_b_size = Array{Float64}(undef,L,L,n_meas)
    for i = 1:n_meas
        corr_mats_for_b_size[:,:,i] = corr_mats[i]
    end
    b_size_1  = maximum([round(Int, 2*auto_corr_time(corr_mats_for_b_size[s1,s2,:]), RoundNearestTiesAway) for s1 = 1:L, s2 = 1:L])
    b_size_2  = maximum([round(Int, 2*auto_corr_time(mean_vals[:,s]), RoundNearestTiesAway) for s = 1:L])    
    b_size    = maximum([b_size_1, b_size_2])
    jack = jack_conn_corr_mat_ev(corr_mats, mean_vals, b_size)
    conn_evs[τ,:]      = reverse(jack[1])
    conn_evs_errs[τ,:] = reverse(jack[2])
end

start_ind = 1
end_ind   = 1
image_conn_evs = plot_corrs(N_t, conn_evs[:,start_ind], conn_evs_errs[:,start_ind], "EV nr. $start_ind", palette(:default)[start_ind])
for smear = start_ind+1:end_ind
    image_conn_evs = plot_corrs!(N_t, conn_evs[:,smear], conn_evs_errs[:,smear], "EV nr. $smear", image_evs, palette(:default)[smear])
end
image_conn_evs = plot!(
    title = "EV's of connected corr. mats. of s_wil\n β = $β, L = $N_t, n_smear = $smear_nums, ρ = $ρ",
    xlabel = latexstring("\$ t \$"),
    right_margin = 10mm,
    yaxis = :log,
    # ylim = (1e-12, 1e-3)
)
display(image_conn_evs)





ev_nr = 2
small = true
mass_conn_ev_2pt_mean = []
mass_conn_ev_2pt_err  = []
for τ = 1:N_t>>1
    # τ = 1
    corr_mats_t1 = corr_mats_list[τ]
    corr_mats_t2 = corr_mats_list[mod1(τ+2,N_t)]  
    if small
        corr_mats_t1 = corr_mats_small_list[τ]
        corr_mats_t2 = corr_mats_small_list[mod1(τ+2,N_t)]  
    end
    L = size(corr_mats_t1[1], 2)
    corr_mats_for_b_size_t1 = Array{Float64}(undef,L,L,n_meas)
    corr_mats_for_b_size_t2 = Array{Float64}(undef,L,L,n_meas)
    for i = 1:n_meas
        corr_mats_for_b_size_t1[:,:,i] = corr_mats_t1[i]
        corr_mats_for_b_size_t2[:,:,i] = corr_mats_t2[i]
    end
    b_size_1  = maximum([round(Int, 2*auto_corr_time(mean_vals[:,s]), RoundNearestTiesAway) for s = 1:L])    
    b_size_2  = maximum([round(Int, 2*auto_corr_time(corr_mats_for_b_size_t1[s1,s2,:]), RoundNearestTiesAway) for s1 = 1:L, s2 = 1:L ])
    b_size_3  = maximum([round(Int, 2*auto_corr_time(corr_mats_for_b_size_t2[s1,s2,:]), RoundNearestTiesAway) for s1 = 1:L, s2 = 1:L ])
    b_size    = maximum([b_size_1, b_size_2, b_size_3])
    jack = jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1,corr_mats_t2,mean_vals,b_size,ev_nr) ./2
    push!(mass_conn_ev_2pt_mean, jack[1])
    push!(mass_conn_ev_2pt_err, jack[2])
end

t_vals_2pt = collect(1:length(mass_conn_ev_2pt_mean)) .- 0.5
image_mass_conn_ev_2pt = scatter(
    t_vals_2pt, circshift(mass_conn_ev_2pt_mean,1), yerror = circshift(mass_conn_ev_2pt_err,1), label = "EV nr. $ev_nr"
)





ev_nr = 1
small = true
mass_conn_ev_2pt_2_mean = []
mass_conn_ev_2pt_2_err  = []
for τ = 1:N_t>>1
    # τ = 1
    corr_mats_t1 = corr_mats_list[τ]
    corr_mats_t2 = corr_mats_list[mod1(τ+2,N_t)]  
    if small
        corr_mats_t1 = corr_mats_small_list[τ]
        corr_mats_t2 = corr_mats_small_list[mod1(τ+2,N_t)]  
    end
    L = size(corr_mats_t1[1], 2)
    corr_mats_for_b_size_t1 = Array{Float64}(undef,L,L,n_meas)
    corr_mats_for_b_size_t2 = Array{Float64}(undef,L,L,n_meas)
    for i = 1:n_meas
        corr_mats_for_b_size_t1[:,:,i] = corr_mats_t1[i]
        corr_mats_for_b_size_t2[:,:,i] = corr_mats_t2[i]
    end
    b_size_1  = maximum([round(Int, 2*auto_corr_time(mean_vals[:,s]), RoundNearestTiesAway) for s = 1:L])    
    b_size_2  = maximum([round(Int, 2*auto_corr_time(corr_mats_for_b_size_t1[s1,s2,:]), RoundNearestTiesAway) for s1 = 1:L, s2 = 1:L ])
    b_size_3  = maximum([round(Int, 2*auto_corr_time(corr_mats_for_b_size_t2[s1,s2,:]), RoundNearestTiesAway) for s1 = 1:L, s2 = 1:L ])
    b_size    = maximum([b_size_1, b_size_2, b_size_3])
    jack = jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1,corr_mats_t2,mean_vals,b_size,ev_nr) ./2
    push!(mass_conn_ev_2pt_2_mean, jack[1])
    push!(mass_conn_ev_2pt_2_err, jack[2])
end

t_vals_2pt_2 = collect(1:length(mass_conn_ev_2pt_2_mean))
image_mass_conn_ev_2pt_2 = scatter(
    t_vals_2pt_2, circshift(mass_conn_ev_2pt_2_mean,1), yerror = circshift(mass_conn_ev_2pt_2_err,1), label = "EV nr. $ev_nr"
)
