include("SU2_analyze_head.jl")



β           = 8.0
N_t = N_x   = 32
n_stout     = 7
ρ           = 0.24
sim_count   = 1
# loops       = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
loops       = [[1,1], [2,2], [3,3], [4,4], [5,5], [6,6]]
# loops       = [[1,1]]
# N_metro     = 1
# N_over      = 3

base_path           = "D:\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
acceptances_path    = string(base_path,"\\acceptances.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
last_conf_path      = string(base_path,"\\last_config.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\last_config_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
corr_mat_paths      = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["D:\\Physik Uni\\julia_projects\\SU2\\data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
mean_vals_path      = string(base_path,"\\mean_vals.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
# corr_mat_paths      = [string(base_path,"\\corrs_clover_t_$t.txt") for t = 1:N_t] #["D:\\Physik Uni\\julia_projects\\SU2\\data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
# mean_vals_path      = string(base_path,"\\mean_vals_clover.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"

mean_vals = readdlm(mean_vals_path)
# ses = s_wil.(mean_vals)
# b_size = round(Int, 2*auto_corr_time(mean_vals)+1)
# jackknife(mean_vals, b_size)
n_meas = size(mean_vals,1)
L = size(mean_vals,2)

for l = 1:L
    mean_vals[:,l] ./= l^2
end





actions = s_from_plaq.(1, mean_vals[:,1]);
plot(actions)

b_sizes = [Int(round(2*auto_corr_time(mean_vals[:,i])+1, RoundNearestTiesAway)) for i = 1:L]
jacks = [jackknife(mean_vals[:,i], b_sizes[i]) for i = 1:L]
loop_means = [jacks[i][1] for i = 1:L]
loop_errs = [jacks[i][2] for i = 1:L]

scatter(string.(loops), loop_means, yerror = loop_errs, label = "Normed Wilson loop exp. values")

plaq_corr_mean = Vector{Float64}(undef,N_t);
plaq_corr_err = Vector{Float64}(undef,N_t);
for t = 1:N_t
    plaq_corrs = readdlm(corr_mat_paths[t])[:,1]
    b_size = round(Int, 2*auto_corr_time(plaq_corrs)+1)
    plaq_corr_mean[t], plaq_corr_err[t] = jackknife(plaq_corrs, b_size)
end
# image_full = plot_corrs(N_t, plaq_corr_mean, plaq_corr_err, :false)
image_full = plot(0:N_t-1, plaq_corr_mean[circshift(1:N_t,1)], yerror = plaq_corr_err[circshift(1:N_t,1)], label = :false)
image_full = plot!(
    title = "Full Plaquette Corr. \n β = $β, L = $N_t, n_smear = $n_stout, ρ = $ρ",
    yaxis = :log
)
display(image_full)

plaq_corr_con_mean = Vector{Float64}(undef,N_t);
plaq_corr_con_err = Vector{Float64}(undef,N_t);
b_sizes = Vector{Float64}(undef,N_t);
for t = 1:N_t
    # t = 1
    plaq_corrs = readdlm(corr_mat_paths[t])[:,1]
    plaq_means = mean_vals[:,1]
    b_size = maximum([round(Int, 2*auto_corr_time(plaq_corrs)+1), round(Int, 2*auto_corr_time(plaq_means)+1)])
    b_sizes[t] = b_size
    plaq_corr_con_mean[t], plaq_corr_con_err[t] = jack_conn_corr_self(plaq_corrs, plaq_means, b_size)
end
# image_con = plot_corrs(N_t, plaq_corr_con_mean, plaq_corr_con_err, :false)
image_con = plot(0:N_t-1, plaq_corr_con_mean[circshift(1:N_t,1)], yerror = plaq_corr_con_err[circshift(1:N_t,1)], label = :false)
image_con = plot!(
    title = "Connected Plaquette Corr. \n β = $β, L = $N_t, n_smear = $n_stout, ρ = $ρ",
    ylims = (1e-10, 1e-7),
    yaxis = :log,
)
display(image_con)
#=
loop_corr_con_mean = Array{Float64}(undef,N_t,L);
loop_corr_con_err = Array{Float64}(undef,N_t,L);
# b_sizes = Vector{Float64}(undef,N_t);
for t = 1:N_t
    # t = 1
    loop_corrs = readdlm(corr_mat_paths[t])
    loop_means = mean_vals
    b_size = [maximum([round(Int, 2*auto_corr_time(loop_corrs[:,l])+1), round(Int, 2*auto_corr_time(loop_means[:,l])+1)]) for l = 1:L]
    # b_sizes[t] = b_size
    for l = 1:L
        loop_corr_con_mean[t,l], loop_corr_con_err[t,l] = jack_conn_corr_self(loop_corrs[:,l], loop_means[:,l], b_size[l])
    end
end
image_con = plot_corrs(N_t, loop_corr_con_mean[:,1], loop_corr_con_err[:,1], latexstring("\$(1\\times 1)\$"))
for l = 2:L
    plot_corrs!(N_t, loop_corr_con_mean[:,l], loop_corr_con_err[:,l], latexstring("\$($l\\times $l)\$"), image_con)
end
image_con = plot!(
    title = "Connected Corr. of normed Wilson loops \n β = $β, n_meas = $n_meas, n_smear = $n_stout, ρ = $ρ"
)

# savefig("D:\\Physik Uni\\Master_Thesis\\glueballs\\plaq_corrs_n_smear_$n_stout.n_meas_$n_meas.pdf")
=#

mass_conn_2pt_mean = []
mass_conn_2pt_err = []
for t = 1:N_t
    loop_corrs_t1 = readdlm(corr_mat_paths[t])[:,1]
    loop_corrs_t2 = readdlm(corr_mat_paths[mod1(t+1,N_t)])[:,1]
    loop_means = mean_vals[:,1]
    τ_corrs_t1 = auto_corr_time(loop_corrs_t1)
    τ_corrs_t2 = auto_corr_time(loop_corrs_t2)
    τ_means = auto_corr_time(loop_means)
    τ_int = maximum([τ_corrs_t1, τ_corrs_t2, τ_means])
    b_size = round(Int, 2*τ_int+1)
    jack = jack_mass_conn_corr_self_2pt(loop_corrs_t1, loop_corrs_t2, loop_means, b_size)
    push!(mass_conn_2pt_mean, jack[1])
    push!(mass_conn_2pt_err, jack[2])
end
t_vals_2pt = collect(1:length(mass_conn_2pt_mean)) .- 0.5
image_mass_conn_2pt = scatter(
    t_vals_2pt, 
    mass_conn_2pt_mean, 
    yerror = mass_conn_2pt_err, 
    xticks = t_vals_2pt,
    title = "2-pt. Masses of conn. plaq. corr.\n β = $β, L = $N_t, n_smear = $n_stout, ρ = $ρ",
    label = latexstring("\$m_{\\textrm{eff}}\\left(t+\\frac{a}{2}\\right) = \\log\\left(\\frac{C(t)}{C(t+a)} \\right) \$"),
    xlabel = latexstring("\$ t\$"),
    ylim = (0.5, 1.5),
    legend = :topleft
    # legendfontsize = 10
)

mass_conn_2pt_2_mean = []
mass_conn_2pt_2_err = []
for t = 1:N_t
    loop_corrs_t1 = readdlm(corr_mat_paths[t])[:,1]
    loop_corrs_t2 = readdlm(corr_mat_paths[mod1(t+2,N_t)])[:,1]
    loop_means = mean_vals[:,1]
    τ_corrs_t1 = auto_corr_time(loop_corrs_t1)
    τ_corrs_t2 = auto_corr_time(loop_corrs_t2)
    τ_means = auto_corr_time(loop_means)
    τ_int = maximum([τ_corrs_t1, τ_corrs_t2, τ_means])
    b_size = round(Int, 2*τ_int+1)
    jack = jack_mass_conn_corr_self_2pt(loop_corrs_t1, loop_corrs_t2, loop_means, b_size) ./2
    push!(mass_conn_2pt_2_mean, jack[1])
    push!(mass_conn_2pt_2_err, jack[2])
end
t_vals_2pt_2 = collect(1:length(mass_conn_2pt_2_mean))
image_mass_conn_2pt_2 = scatter(
    t_vals_2pt_2,
    mass_conn_2pt_2_mean, 
    yerror = mass_conn_2pt_2_err, 
    xticks = t_vals_2pt_2,
    title = "2-pt. Masses of conn. plaq. corr.\n β = $β, L = $N_t, n_smear = $n_stout, ρ = $ρ",
    label = latexstring("\$m_{\\textrm{eff}}(t) = \\frac{1}{2} \\log\\left(\\frac{C(t-a)}{C(t+a)} \\right) \$"),
    xlabel = latexstring("\$ t\$"),
    ylim = (0.5, 1.5),
    legend = :topleft
)

mass_conn_3pt_mean = [];
mass_conn_3pt_err = [];
for t = 1:N_t
    loop_corrs_t1 = readdlm(corr_mat_paths[t])[:,1]
    loop_corrs_t2 = readdlm(corr_mat_paths[mod1(t+1,N_t)])[:,1]
    loop_corrs_t3 = readdlm(corr_mat_paths[mod1(t+2,N_t)])[:,1]
    loop_means = mean_vals[:,1]
    τ_corrs_t1 = auto_corr_time(loop_corrs_t1)
    τ_corrs_t2 = auto_corr_time(loop_corrs_t2)
    τ_corrs_t3 = auto_corr_time(loop_corrs_t3)
    τ_means = auto_corr_time(loop_means)
    τ_int = maximum([τ_corrs_t1, τ_corrs_t2, τ_means])
    b_size = round(Int, 2*τ_int+1)
    jack = jack_mass_conn_corr_self_3pt(loop_corrs_t1, loop_corrs_t2, loop_corrs_t3, loop_means, b_size)
    push!(mass_conn_3pt_mean, jack[1])
    push!(mass_conn_3pt_err, jack[2])
end
t_vals_3pt = collect(1:length(mass_conn_3pt_mean)) .- 1 
image_mass_conn_3pt = scatter(
    t_vals_3pt,
    mass_conn_3pt_mean, 
    yerror = mass_conn_3pt_err, 
    xticks = t_vals_3pt,
    title = "3-pt. Masses of conn. plaq. corr.\n β = $β, L = $N_t, n_smear = $n_stout, ρ = $ρ",
    label = latexstring("\$m_{\\textrm{eff}}(t) = \\textrm{acosh}\\,\\left(\\frac{C(t-a) + C(t+a)}{2C(t)} \\right) \$"),
    xlabel = latexstring("\$ t\$"),
    ylim = (0.5, 1.5),
    legend = :topleft
    # legendfontsize = 10
)
    











#=

mean_vals_clover = readdlm(mean_vals_clover_path)

b_sizes = [Int(round(2*auto_corr_time(mean_vals_clover[:,i])+1, RoundNearestTiesAway)) for i = 1:L]
jacks = [jackknife(mean_vals_clover[:,i], b_sizes[i]) for i = 1:L]
clover_means = [jacks[i][1] for i = 1:L]
clover_errs = [jacks[i][2] for i = 1:L]

clov_corr_mean = Vector{Float64}(undef,N_t);
clov_corr_err = Vector{Float64}(undef,N_t);
for t = 1:N_t
    clov_corrs = readdlm(corr_mat_clover_paths[t])[:,1]
    b_size = round(Int, 2*auto_corr_time(clov_corrs)+1)
    clov_corr_mean[t], clov_corr_err[t] = jackknife(clov_corrs, b_size)
end
image_full_clov = plot_corrs(N_t, clov_corr_mean, clov_corr_err, :false)
image_full_clov = plot!(
    title = "Full Clover Corr. \n β = $β, n_meas = $n_meas, n_smear = $n_stout, ρ = $ρ"
)
display(image_full_clov)

clov_corr_con_mean = Vector{Float64}(undef,N_t);
clov_corr_con_err = Vector{Float64}(undef,N_t);
b_sizes = Vector{Float64}(undef,N_t);
for t = 1:N_t
    # t = 1
    clov_corrs = readdlm(corr_mat_clover_paths[t])[:,1]
    clov_means = mean_vals_clover[:,1]
    b_size = maximum([round(Int, 2*auto_corr_time(clov_corrs)+1), round(Int, 2*auto_corr_time(clov_means)+1)])
    b_sizes[t] = b_size
    clov_corr_con_mean[t], clov_corr_con_err[t] = jack_conn_corr_self(clov_corrs, clov_means, b_size)
end
image_con_clov = plot_corrs(N_t, clov_corr_con_mean, clov_corr_con_err, :false)
image_con_clov = plot!(
    title = "Connected Clover Corr. \n β = $β, n_meas = $n_meas, n_smear = $n_stout, ρ = $ρ"
)
display(image_con_clov)

clov_mass_conn_2pt_mean = []
clov_mass_conn_2pt_err = []
for t = 1:N_t
    loop_corrs_t1 = readdlm(corr_mat_clover_paths[t])[:,1]
    loop_corrs_t2 = readdlm(corr_mat_clover_paths[mod1(t+1,N_t)])[:,1]
    loop_means = mean_vals[:,1]
    τ_corrs_t1 = auto_corr_time(loop_corrs_t1)
    τ_corrs_t2 = auto_corr_time(loop_corrs_t2)
    τ_means = auto_corr_time(loop_means)
    τ_int = maximum([τ_corrs_t1, τ_corrs_t2, τ_means])
    b_size = round(Int, 2*τ_int+1)
    jack = jack_mass_conn_corr_self_2pt(loop_corrs_t1, loop_corrs_t2, loop_means, b_size)
    push!(clov_mass_conn_2pt_mean, jack[1])
    push!(clov_mass_conn_2pt_err, jack[2])
end
t_vals_clov_2pt = collect(1:length(clov_mass_conn_2pt_mean)) .- 0.5
image_clov_mass_conn_2pt = scatter(
    t_vals_clov_2pt, 
    clov_mass_conn_2pt_mean, 
    yerror = clov_mass_conn_2pt_err, 
    # xticks = t_vals_clov_2pt,
    title = "2-pt. Masses of conn. plaq. corr.\n β = $β, n_meas = $n_meas, n_smear = $n_stout, ρ = $ρ",
    label = latexstring("\$m_{\\textrm{eff}} = \\log\\left(\\frac{C(t)}{C(t+a)} \\right) \$"),
    xlabel = latexstring("\$ t\$"),
    # legendfontsize = 10
)

clov_mass_conn_2pt_2_mean = []
clov_mass_conn_2pt_2_err = []
for t = 1:N_t
    loop_corrs_t1 = readdlm(corr_mat_clover_paths[t])[:,1]
    loop_corrs_t2 = readdlm(corr_mat_clover_paths[mod1(t+2,N_t)])[:,1]
    loop_means = mean_vals[:,1]
    τ_corrs_t1 = auto_corr_time(loop_corrs_t1)
    τ_corrs_t2 = auto_corr_time(loop_corrs_t2)
    τ_means = auto_corr_time(loop_means)
    τ_int = maximum([τ_corrs_t1, τ_corrs_t2, τ_means])
    b_size = round(Int, 2*τ_int+1)
    jack = jack_mass_conn_corr_self_2pt(loop_corrs_t1, loop_corrs_t2, loop_means, b_size) ./2
    push!(clov_mass_conn_2pt_2_mean, jack[1])
    push!(clov_mass_conn_2pt_2_err, jack[2])
end
t_vals_clov_2pt_2 = collect(1:length(clov_mass_conn_2pt_2_mean))
image_clov_mass_conn_2pt_2 = scatter(
    t_vals_clov_2pt_2,
    clov_mass_conn_2pt_2_mean, 
    yerror = clov_mass_conn_2pt_2_err, 
    xticks = t_vals_clov_2pt,
    title = "2-pt. Masses of conn. plaq. corr.\n β = $β, n_meas = $n_meas, n_smear = $n_stout, ρ = $ρ",
    label = latexstring("\$m_{\\textrm{eff}} = \\frac{1}{2} \\log\\left(\\frac{C(t)}{C(t+2a)} \\right) \$"),
    xlabel = latexstring("\$ t\$"),
    # legendfontsize = 10
)

clov_mass_conn_3pt_mean = [];
clov_mass_conn_3pt_err = [];
for t = 1:N_t
    loop_corrs_t1 = readdlm(corr_mat_clover_paths[t])[:,1]
    loop_corrs_t2 = readdlm(corr_mat_clover_paths[mod1(t+1,N_t)])[:,1]
    loop_corrs_t3 = readdlm(corr_mat_clover_paths[mod1(t+2,N_t)])[:,1]
    loop_means = mean_vals[:,1]
    τ_corrs_t1 = auto_corr_time(loop_corrs_t1)
    τ_corrs_t2 = auto_corr_time(loop_corrs_t2)
    τ_corrs_t3 = auto_corr_time(loop_corrs_t3)
    τ_means = auto_corr_time(loop_means)
    τ_int = maximum([τ_corrs_t1, τ_corrs_t2, τ_means])
    b_size = round(Int, 2*τ_int+1)
    jack = jack_mass_conn_corr_self_3pt(loop_corrs_t1, loop_corrs_t2, loop_corrs_t3, loop_means, b_size)
    push!(clov_mass_conn_3pt_mean, jack[1])
    push!(clov_mass_conn_3pt_err, jack[2])
end
t_vals_clov_3pt = collect(1:length(clov_mass_conn_3pt_mean)) 
image_clov_mass_conn_3pt = scatter(
    t_vals_clov_3pt,
    clov_mass_conn_3pt_mean, 
    yerror = clov_mass_conn_3pt_err, 
    xticks = t_vals_clov_3pt,
    title = "3-pt. Masses of conn. plaq. corr.\n β = $β, n_meas = $n_meas, n_smear = $n_stout, ρ = $ρ",
    label = latexstring("\$m_{\\textrm{eff}} = \\textrm{acosh}\\,\\left(\\frac{C(t) + C(t+2a)}{C(t+a)} \\right) \$"),
    xlabel = latexstring("\$ t\$"),
    # legendfontsize = 10
)

=#








#=
# For each t we want corr_mat_ar[t] to contain the time series of the correlation
# matrix at Euclidean time t
corr_mat_ar = Array{Array{Matrix{Float64}}}(undef, N_t)
for t = 1:N_t
    corr_mats = readdlm(corr_mat_paths[t])
    corr_mat_ar[t] = [corr_mats[Int((i-1)*L+1) : Int(i*L), 1:L] for i = 1:n_meas]
    corr_mat_ar[t] = (corr_mat_ar[t] .+ transpose.(corr_mat_ar[t])) ./ 2
    # for meas = 1:n_meas
    #     corr_mat_ar[t][meas] ./= [i^2*j^2 for i = 1:L, j = 1:L]
    # end
end

plaq_corr_mean = Vector{Float64}(undef,N_t);
plaq_corr_err = Vector{Float64}(undef,N_t);
for t = 1:N_t
    corr_mats = readdlm(corr_mat_paths[t])
    plaq_corrs_t = [corr_mats[Int((i-1)*L+1), 1] for i = 1:n_meas]
    b_size = round(Int, 2*auto_corr_time(plaq_corrs_t)+1)
    plaq_corr_mean[t], plaq_corr_err[t] = jackknife(plaq_corrs_t, b_size)
end
plot_corrs(N_t, plaq_corr_mean, plaq_corr_err, :false)
# plot_corrs(N_t, plaq_corr_mean.-loop_means[1]^2, plaq_corr_err, :false)
=#







e_val_means = Matrix{Float64}(undef, L, N_t)
e_val_errs = Matrix{Float64}(undef, L, N_t)
for t = 1:N_t
    e_val_mat = Matrix{Float64}(undef, n_meas, L)
    for i = 1:n_meas
        # e_val_mat[i,:] = abs.(eigen(corr_mat_ar[t][i]).values)
        # e_val_mat[i,:] = eigen(corr_mat_ar[t][i]).values
        e_val_mat[i,:] = real.(eigen(corr_mat_ar[t][i]).values)    # evals of sym. matr. are real
    end
    b_sizes = [Int(round(2*auto_corr_time(e_val_mat[:,i])+1, RoundNearestTiesAway)) for i = 1:L]
    jacks = [jackknife(e_val_mat[:,i], b_sizes[i]) for i = 1:L]
    e_val_jack_means = [jacks[i][1] for i = 1:L]
    e_val_jack_errs = [jacks[i][2] for i = 1:L]
    e_val_means[:,t] = e_val_jack_means
    e_val_errs[:,t] = e_val_jack_errs
end

gen_e_val_means = Matrix{Float64}(undef, L, N_t)
gen_e_val_errs = Matrix{Float64}(undef, L, N_t)
for t = 1:N_t
    gen_e_val_mat = Matrix{Float64}(undef, n_meas, L)
    for i = 1:n_meas
        gen_e_val_mat[i,:] = real.(eigen(corr_mat_ar[t][i], corr_mat_ar[N_t][i]).values)
    end
    b_sizes = [Int(round(2*auto_corr_time(gen_e_val_mat[:,i])+1, RoundNearestTiesAway)) for i = 1:L]
    jacks = [jackknife(gen_e_val_mat[:,i], b_sizes[i]) for i = 1:L]
    gen_e_val_jack_means = [jacks[i][1] for i = 1:L]
    gen_e_val_jack_errs = [jacks[i][2] for i = 1:L]
    gen_e_val_means[:,t] = gen_e_val_jack_means
    gen_e_val_errs[:,t] = gen_e_val_jack_errs
end

# for t = 1:N_t
#     @assert isapprox(imag.(gen_e_val_means[:,t]), zeros(L))
# end

let ev_nr = 6
    # ev_nr = 3
    display(scatter(1:N_t, real.(gen_e_val_means[ev_nr,:]), yerror = real.(gen_e_val_errs[ev_nr,:]), label = "Gen. EV Nr = $ev_nr, real", legend = :top))
    display(scatter(1:N_t, real.(e_val_means[ev_nr,:]), yerror = real.(e_val_errs[ev_nr,:]), label = "EV Nr = $ev_nr, real", legend = :top))
    # display(scatter(1:N_t, imag.(e_val_means[ev_nr,:]), yerror = imag.(e_val_errs[ev_nr,:]), label = "EV Nr = $ev_nr, imag", legend = :top))
    # display(scatter(1:N_t, imag.(e_val_means[ev_nr,:]) ./ real.(e_val_means[ev_nr,:]), label = "EV Nr = $ev_nr, real/imag"))
end

let
    image_ev = scatter(1:N_t, e_val_means[1,:], yerror = e_val_errs[1,:], label = "EV nr. 1")
    for l = 2:L-1
        image_ev = scatter!(1:N_t, e_val_means[l,:], yerror = e_val_errs[l,:], label = "EV nr. $l")
    end
    image_ev = plot!(
        title = "EV's of corr. matrices of Wilson loops\n 3D U(2), β = $β, N_meas = $n_meas",
        legendfontsize = 7
    )
    display(image_ev)
end

let
    image_gen_ev = scatter(1:N_t, gen_e_val_means[1,:], yerror = gen_e_val_errs[1,:], label = "EV nr. 1")
    for l = 2:L-1
        image_ev = scatter!(1:N_t, gen_e_val_means[l,:], yerror = gen_e_val_errs[l,:], label = "EV nr. $l")
    end
    image_ev = plot!(
        title = "Gen. EV's of corr. matrices of Wilson loops\n 3D U(2), β = $β, N_meas = $n_meas",
        legendfontsize = 7
    )
    display(image_gen_ev)
end

# gen_e_val_means[:,N_t]

# scatter(1:10, e_val_means[16,1:10], yerror = e_val_errs[16,1:10])






small_corr_mat_ar = [] #Array{Array{Matrix{Float64}}}(undef,N_t)
for t = 1:N_t
    push!(small_corr_mat_ar, [])
    corr_mats = readdlm(corr_mat_paths[t])
    # corr_mat_ar[t] = [corr_mats[Int((i-1)*L+1) : Int(i*L), 1:L] for i = 1:n_meas]
    for i = 1:n_meas
        large_corr_mat = corr_mats[Int((i-1)*L+1) : Int(i*L), 1:L]
        small_corr_mat = [large_corr_mat[1+3*j, 1+3*k]/((j+1)*(k+1))^2 for j = 0:5, k = 0:5]
        # small_corr_mat = [large_corr_mat[1+3*j, 1+3*k] for j = 0:5, k = 0:5]
        small_corr_mat = (small_corr_mat + transpose(small_corr_mat))/2
        push!(small_corr_mat_ar[t], small_corr_mat)
    end
end

LLL = 6

e_val_means = Matrix{Float64}(undef, LLL, N_t)
e_val_errs = Matrix{Float64}(undef, LLL, N_t)
for t = 1:N_t
    e_val_mat = Matrix{Float64}(undef, n_meas, LLL)
    for i = 1:n_meas
        e_val_mat[i,:] = real.(eigen(small_corr_mat_ar[t][i]).values)    # evals of sym. matr. are real
    end
    b_sizes = [Int(round(2*auto_corr_time(e_val_mat[:,i])+1, RoundNearestTiesAway)) for i = 1:LLL]
    jacks = [jackknife(e_val_mat[:,i], b_sizes[i]) for i = 1:LLL]
    e_val_jack_means = [jacks[i][1] for i = 1:LLL]
    e_val_jack_errs = [jacks[i][2] for i = 1:LLL]
    e_val_means[:,t] = e_val_jack_means
    e_val_errs[:,t] = e_val_jack_errs
end

let
    image_ev = plot_corrs(N_t, e_val_means[1,:], e_val_errs[1,:], "EV nr. $LLL")
    for l = 2:LLL
        plot_corrs!(N_t, e_val_means[l,:], e_val_errs[l,:], "EV nr. $(LLL+1-l)", image_ev)
    end
    image_ev = plot!(
        title = "EV's of corr. matrices of Wilson squares\n 3D U(2), β = $β, N_meas = $n_meas",
        legendfontsize = 7,
        xlabel = latexstring("Temporal extent \$\\tau\$")
    )
    display(image_ev)
end

# let
#     image_gen_ev = plot_corrs(N_t, gen_e_val_means[1,:], gen_e_val_errs[1,:], "EV nr. $LLL")
#     for l = 2:LLL-1
#         plot_corrs!(N_t, gen_e_val_means[l,:], gen_e_val_errs[l,:], "EV nr. $(LLL+1-l)", image_gen_ev)
#     end
#     image_gen_ev = plot!(
#         title = "Gen. EV's of corr. matrices of Wilson loops\n 3D U(2), β = $β, N_meas = $n_meas",
#         legendfontsize = 7,
#         xlabel = latexstring("Temporal extent \$\\tau\$")
#     )
#     display(image_gen_ev)
# end



let
    global LLL = 6
    global t_start = 4
    global gen_e_val_means = Matrix{Float64}(undef, LLL, N_t)
    global gen_e_val_errs = Matrix{Float64}(undef, LLL, N_t)
    for t = t_start+1:N_t
        gen_e_val_mat = Matrix{Float64}(undef, n_meas, LLL)
        for i = 1:n_meas
            gen_e_val_mat[i,:] = real.(eigen(small_corr_mat_ar[t][i], small_corr_mat_ar[t_start][i]).values)
        end
        if t == 2
            global bla = gen_e_val_mat
        end
        b_sizes = [Int(round(2*auto_corr_time(gen_e_val_mat[:,i])+1, RoundNearestTiesAway)) for i = 1:LLL]
        jacks = [jackknife(gen_e_val_mat[:,i], b_sizes[i]) for i = 1:LLL]
        gen_e_val_jack_means = [jacks[i][1] for i = 1:LLL]
        gen_e_val_jack_errs = [jacks[i][2] for i = 1:LLL]
        gen_e_val_means[:,t] = gen_e_val_jack_means
        gen_e_val_errs[:,t] = gen_e_val_jack_errs
    end
end

let
    interval = t_start+1:N_t
    start_l = 1
    end_l   = LLL
    image_gen_ev = scatter(interval, gen_e_val_means[start_l,interval], yerror = gen_e_val_errs[start_l,interval], label = "EV nr. $end_l")
    for l = start_l+1:end_l
        image_ev = scatter!(interval, gen_e_val_means[l,interval], yerror = gen_e_val_errs[l,interval], label = "EV nr. $(LLL+1-l)")
    end
    image_ev = plot!(
        title = "Gen. EV's of corr. matrices of Wilson loops\n 3D U(2), β = $β, N_meas = $n_meas",
        legendfontsize = 7
    )
    display(image_gen_ev)
end















let
    N_t = N_x   = 8
    ρ           = 0.22
    sim_count   = 2
    n_meas      = 1000
    image_s_wil = plot(title = "Action Density with Timeslice-Smearing \n L = $N_t, N_meas = 1000, ρ = $ρ", xlabel = "N_stout")

    for β in [4.0, 7.0, 10.0]
        s_wil_smeared_means = []
        s_wil_smeared_errs  = []
        for n_stout = 0:2:8
            if n_stout == 0
                sim_count = 1
            else
                sim_count = 2
            end
            base_path           = "D:\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
            params_path         = string(base_path,"\\params.txt") #
            acceptances_path    = string(base_path,"\\acceptances.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
            last_conf_path      = string(base_path,"\\last_config.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\last_config_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
            corr_mat_paths      = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["D:\\Physik Uni\\julia_projects\\SU2\\data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
            mean_vals_path      = string(base_path,"\\mean_vals.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"

            mean_vals = readdlm(mean_vals_path)
            n_meas = size(mean_vals,1)
            L = size(mean_vals,2)

            ses = s_wil.(mean_vals[:,1])
            b_size = round(Int, 2*auto_corr_time(ses)+1)
            s_wil_mean, s_wil_err = jackknife(ses, b_size)
            push!(s_wil_smeared_means, s_wil_mean)
            push!(s_wil_smeared_errs,  s_wil_err)
        end
        image_s_wil = scatter!(
            0:2:8, 
            xticks = 0:2:8,
            s_wil_smeared_means, 
            yerror = s_wil_smeared_errs,
            label = "β = $β"
        )
    end
display(image_s_wil)
end