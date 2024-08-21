include("SU2_analyze_head.jl")

using LaTeXStrings



β           = 7.0
# D           = 3
N_t = N_x   = 32
# ϵ           = 0.2
n_stout     = 0
ρ           = 0.1
sim_count   = 2
loops     = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
# N_metro   = 1
# N_over    = 3

base_path           = "D:\\Physik Uni\\julia_projects\\SU2_data\\3d_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
params_path         = string(base_path,"\\params.txt") #
acceptances_path    = string(base_path,"\\acceptances.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\acceptances_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
last_conf_path      = string(base_path,"\\last_config.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\last_config_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
corr_mat_paths      = [string(base_path,"\\corrs_t_$t.txt") for t = 1:N_t] #["D:\\Physik Uni\\julia_projects\\SU2\\data\\corrs_t_$t._eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt" for t = 1:N_t]
mean_vals_path      = string(base_path,"\\mean_vals.txt") #"D:\\Physik Uni\\julia_projects\\SU2\\data\\mean_vals_eps_$ϵ._beta_$β._L_$N_t._n_stout_$n_stout._rho_$ρ.txt"
# mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")
# last_conf_path      = string(last_base_path,"\\last_config.txt")

mean_vals = readdlm(mean_vals_path)
n_meas = size(mean_vals,1)
# n_meas = 92 #414
L = length(loops)

# acc = readdlm(acceptances_path)
# last(acc) /( (N_metro+N_over) * D * N_t * N_x^(D-1) * n_meas )
# for i = 1:n_meas
#     acc[i] /= (N_metro+N_over) * D * N_t * N_x^(D-1) * i
# end

# last(acc)
# plot(acc)


function plot_corrs(N_t, corr_means, corr_errs, series_label)
    return scatter(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label)
end

function plot_corrs!(N_t, corr_means, corr_errs, series_label, image)
    image = scatter!(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label)
    return nothing
end

function s_from_plaq(β, plaq_mean)
    return 3 * β * (2-plaq_mean) / 2
end

actions = s_from_plaq.(1, mean_vals[:,1]);
plot(actions)
# bla = gaugefield_SU2_cube(8, 8, true);
# pmean = (mean(tr.([plaq_12(bla, x, y, t) for x = 1:8, y = 1:8, t = 1:8])) + mean(tr.([plaq_23(bla, x, y, t) for x = 1:8, y = 1:8, t = 1:8])) + mean(tr.([plaq_13(bla, x, y, t) for x = 1:8, y = 1:8, t = 1:8]))) / 3
# s_from_plaq(1, pmean)
# action_cube(bla, 1)/8^3

# mean_vals = readdlm(mean_vals_path)

b_sizes = [Int(round(2*auto_corr_time(mean_vals[:,i])+1, RoundNearestTiesAway)) for i = 1:L]
jacks = [jackknife(mean_vals[:,i], b_sizes[i]) for i = 1:L]
loop_means = [jacks[i][1] for i = 1:L]
loop_errs = [jacks[i][2] for i = 1:L]

scatter(string.(loops), loop_means, yerror = loop_errs, label = "Wilson loop exp. values")



# For each t we want corr_mat_ar[t] to contain the time series of the correlation
# matrix at Euclidean time t
corr_mat_ar = Array{Array{Matrix{Float64}}}(undef, N_t)
corr_mat_ar_sym = Array{Array{Matrix{Float64}}}(undef, N_t)
for t = 1:N_t
    corr_mats = readdlm(corr_mat_paths[t])
    corr_mat_ar[t] = [corr_mats[Int((i-1)*L+1) : Int(i*L), 1:L] for i = 1:n_meas]
    # corr_mat_ar_sym[t] = [(corr_mats[Int((i-1)*L+1) : Int(i*L), 1:L] + transpose(corr_mats[Int((i-1)*L+1) : Int(i*L), 1:L]))/2 for i = 1:n_meas]
    corr_mat_ar_sym[t] = (corr_mat_ar[t] .+ transpose.(corr_mat_ar[t])) ./ 2
end

e_val_means = Matrix{Float64}(undef, L, N_t)
e_val_errs = Matrix{Float64}(undef, L, N_t)
for t = 1:N_t
    e_val_mat = Matrix{Float64}(undef, n_meas, L)
    for i = 1:n_meas
        # e_val_mat[i,:] = abs.(eigen(corr_mat_ar[t][i]).values)
        # e_val_mat[i,:] = eigen(corr_mat_ar[t][i]).values
        e_val_mat[i,:] = real.(eigen(corr_mat_ar_sym[t][i]).values)    # evals of sym. matr. are real
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
        gen_e_val_mat[i,:] = real.(eigen(corr_mat_ar_sym[t][i], corr_mat_ar_sym[N_t][i]).values)
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

let ev_nr = 16
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






plaq_corr_mean = Vector{Float64}(undef,N_t);
plaq_corr_err = Vector{Float64}(undef,N_t);
for t = 1:N_t
    corr_mats = readdlm(corr_mat_paths[t])
    plaq_corrs_t = [corr_mats[Int((i-1)*L+1), 1] for i = 1:n_meas]
    b_size = round(Int, 2*auto_corr_time(plaq_corrs_t)+1)
    plaq_corr_mean[t], plaq_corr_err[t] = jackknife(plaq_corrs_t, b_size)
end
plot_corrs(N_t, plaq_corr_mean, plaq_corr_err, :false)


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
    for l = 2:LLL-1
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
    global t_start = 1
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
    start_l = 2
    end_l   = LLL-1
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
