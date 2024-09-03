using Plots
using StatsBase
using DelimitedFiles
using LsqFit
using SpecialFunctions
using LaTeXStrings

include("D:\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("D:\\Physik Uni\\julia_projects\\SU2\\observables\\observables_hex.jl")
include("D:\\Physik Uni\\julia_projects\\SU2\\observables\\observables_cube.jl")



enu_endings = ["st", "nd", "rd"]
for i = 1:50
    push!(enu_endings, "th")
end

function bootstrap(obs, b_size, N_boot)
    N_blocks = Int(div(length(obs), b_size, RoundDown))
    # Step 1: blocking
    new_obs = zeros(N_blocks)
    for i = 1:N_blocks
        start = Int(1 + (i-1)*b_size)
        ending = Int(start+b_size-1)
        new_obs[i] = mean(obs[start:ending])
    end
    # Step 2: bootstrap
    means = zeros(N_boot)
    for i = 1:N_boot
        strap = zeros(N_blocks)
        r = rand(1:N_blocks, N_blocks)
        for j = 1:N_blocks
            strap[j] = new_obs[r[j]]
        end
        means[i] = mean(strap)
    end
    # return [mean(means), std(means)]
    return [mean(obs), std(means)]
end

function jackknife(obs, b_size)
    N_blocks = Int(div(length(obs), b_size, RoundDown))
    first_mean = mean(obs[b_size+1:b_size*N_blocks])
    last_mean = mean(obs[1:b_size*N_blocks - b_size])
    jack_means = [first_mean, last_mean]
    for i = 2:N_blocks-1
        push!(jack_means, mean(vcat(obs[1:(i-1)*b_size], obs[i*b_size+1:b_size*N_blocks])))
    end
    obs_mean = mean(obs)
    σ = sqrt((N_blocks-1) * mean((jack_means.-obs_mean).^2 ))

    return [obs_mean, σ]
end

# function to determine the autocorrelation of an observable (stored in an
# array "obs") at a simulation time t
function auto_corr(obs, t)
    M = length(obs)
    C = 0.0
    obs_mean = mean(obs)
    for i = 1:Int(M-t)
        C += (obs[i]-obs_mean)*(obs[Int(i+t)]-obs_mean)
    end
    return C/(M-t)
    # return mean((obs[1:M-t] .- mean(obs[1:M-t])) .* (obs[t+1:M] .- mean(obs[t+1:M])) )
end

function auto_corr_norm(obs, t)
    return auto_corr(obs, t)/auto_corr(obs, 0)
end

# 🐌🐌🐌 Could be faster 🐌🐌🐌
function auto_corr_time(obs)
    τ = 0.5
    for i = 1:length(obs)
        τ_add = auto_corr_norm(obs,i)
        if τ_add <= 0.0
            break
        end
        τ += τ_add
    end
    return τ
end

function mass_2pt(corrs)
    L = length(corrs)
    masses = Array{Float64}(undef,L)
    for i = 1:L
        i_p = mod1(i+1,L)
        masses[i] = log(corrs[i]/corrs[i_p])
    end
    return masses
end

function mass_3pt(corrs)
    L = length(corrs)
    masses = zeros(L)
    for i = 1:L
        i_p = mod1(i+1,L)
        i_m = mod1(i-1,L)
        arg = (corrs[i_p] + corrs[i_m])/(2*corrs[i])
        if arg > 0.0
            masses[i] = acosh(arg)
        end
    end
    return masses
end


#=
function creutz_err(a, b, c, d, aerr, berr, cerr, derr)
    return sqrt((aerr*d/c/b)^2 + (derr*a/c/b)^2 + (cerr*a*b/c^2/b)^2 + (berr*a*b/c/b^2)^2  )
end

function string_err(a, b, c, d, aerr, berr, cerr, derr)
    return sqrt((aerr/a)^2 + (derr/d)^2 + (cerr/c)^2 + (berr/b)^2)
end
=#


function plot_corrs(N_t, corr_means, corr_errs, series_label)
    return scatter(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label)
end

function plot_corrs!(N_t, corr_means, corr_errs, series_label, image)
    image = scatter!(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label)
    return nothing
end

function plot_corrs(N_t, corr_means, corr_errs, series_label, farbe)
    return scatter(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label, color = farbe)
end

function plot_corrs!(N_t, corr_means, corr_errs, series_label, image, farbe)
    image = scatter!(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label, color = farbe)
    return nothing
end

function s_from_plaq(β, plaq_mean)
    return 3 * β * (2-plaq_mean) / 2
end

function s_wil(plaq_mean)
    return 1 - plaq_mean/2
end

function jack_conn_corr_self(corrs, loop_means, b_size)
    N_blocks = Int(div(length(corrs), b_size, RoundDown))

    first_loop = mean(loop_means[b_size+1:b_size*N_blocks])
    first_corr = mean(corrs[b_size+1:b_size*N_blocks])
    first_con = first_corr - first_loop^2

    last_loop = mean(loop_means[1:b_size*N_blocks - b_size])
    last_corr = mean(corrs[1:b_size*N_blocks - b_size])
    last_con = last_corr - last_loop^2

    jack_cons = [first_con, last_con]
    for i = 2:N_blocks-1
        temp_loop = mean(vcat(loop_means[1:(i-1)*b_size], loop_means[i*b_size+1:b_size*N_blocks]))
        temp_corr = mean(vcat(corrs[1:(i-1)*b_size], corrs[i*b_size+1:b_size*N_blocks]))
        push!(jack_cons, temp_corr - temp_loop^2)
    end

    con_mean = mean(corrs) - mean(loop_means)^2
    σ = sqrt((N_blocks-1) * mean((jack_cons .- con_mean).^2 ))
    return con_mean, σ
end

function jack_mass_conn_corr_self_2pt(corrs_t1, corrs_t2, loop_means, b_size)
    N_blocks = Int(div(length(corrs_t1), b_size, RoundDown))
    jack_masses = Vector{Float64}(undef,N_blocks)
    for i = 1:N_blocks
        temp_loop = mean(vcat(loop_means[1:(i-1)*b_size], loop_means[i*b_size+1:b_size*N_blocks]))
        temp_corr_t1 = mean(vcat(corrs_t1[1:(i-1)*b_size], corrs_t1[i*b_size+1:b_size*N_blocks]))
        temp_corr_t2 = mean(vcat(corrs_t2[1:(i-1)*b_size], corrs_t2[i*b_size+1:b_size*N_blocks]))
        jack_masses[i] = log((temp_corr_t1-temp_loop^2) / (temp_corr_t2-temp_loop^2))
    end
    loop = mean(loop_means)
    mass_mean = log((mean(corrs_t1)-loop^2) / (mean(corrs_t2)-loop^2))    
    σ = sqrt((N_blocks-1) * mean((jack_masses .- mass_mean).^2 ))
    return mass_mean, σ
end

function jack_mass_conn_corr_self_3pt(corrs_t1, corrs_t2, corrs_t3, loop_means, b_size)
    N_blocks = Int(div(length(corrs_t1), b_size, RoundDown))
    jack_masses = Vector{Float64}(undef,N_blocks)
    for i = 1:N_blocks
        temp_loop = mean(vcat(loop_means[1:(i-1)*b_size], loop_means[i*b_size+1:b_size*N_blocks]))
        temp_corr_t1 = mean(vcat(corrs_t1[1:(i-1)*b_size], corrs_t1[i*b_size+1:b_size*N_blocks]))
        temp_corr_t2 = mean(vcat(corrs_t2[1:(i-1)*b_size], corrs_t2[i*b_size+1:b_size*N_blocks]))
        temp_corr_t3 = mean(vcat(corrs_t3[1:(i-1)*b_size], corrs_t3[i*b_size+1:b_size*N_blocks]))
        jack_masses[i] = acosh((temp_corr_t1+temp_corr_t3-2*temp_loop^2) / (2*(temp_corr_t2-temp_loop^2)))
    end
    loop = mean(loop_means)
    mass_mean = acosh((mean(corrs_t1)+mean(corrs_t3)-2*loop^2) / (2*(mean(corrs_t2)-loop^2)))    
    σ = sqrt((N_blocks-1) * mean((jack_masses .- mass_mean).^2 ))
    return mass_mean, σ
end

function jack_corr_mat_ev(corr_mats, b_size)
    num_vals = size(corr_mats[1],1)
    N_blocks = Int(div(length(corr_mats), b_size, RoundDown))
    jack_evs = Array{Float64}(undef, N_blocks, num_vals)
    for i = 1:N_blocks 
        temp_corr = mean(vcat(corr_mats[1:(i-1)*b_size], corr_mats[i*b_size+1:b_size*N_blocks]))
        jack_evs[i,:] = eigen(temp_corr).values
    end
    mean_evs = eigen(mean(corr_mats)).values
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end
