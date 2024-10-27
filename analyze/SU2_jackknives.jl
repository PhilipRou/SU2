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

function jack_conn_corr_self(corrs, loop_means, b_size)
    N_blocks = Int(div(length(corrs), b_size, RoundDown))

    # first_loop = mean(loop_means[b_size+1:b_size*N_blocks])
    # first_corr = mean(corrs[b_size+1:b_size*N_blocks])
    # first_con = first_corr - first_loop^2

    # last_loop = mean(loop_means[1:b_size*N_blocks - b_size])
    # last_corr = mean(corrs[1:b_size*N_blocks - b_size])
    # last_con = last_corr - last_loop^2

    # jack_cons = [first_con, last_con]
    jack_cons = []
    for i = 1:N_blocks
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
    N_blocks = Int(div(length(corr_mats), b_size, RoundDown))
    num_vals = size(corr_mats[1],1)
    jack_evs = Array{Float64}(undef, N_blocks, num_vals)
    for i = 1:N_blocks 
        temp_corr = mean(vcat(corr_mats[1:(i-1)*b_size], corr_mats[i*b_size+1:b_size*N_blocks]))
        jack_evs[i,:] = eigen(temp_corr).values
    end
    mean_evs = eigen(mean(corr_mats)).values
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end

function jack_conn_corr_mat_ev(corr_mats, mean_vals, b_size)
    N_blocks = Int(div(length(corr_mats), b_size, RoundDown))
    num_vals = size(corr_mats[1],1)
    mean_vals_foo = [mean_vals[i,:] for i = 1:Int(b_size*N_blocks)] # just to have an array of arrays instead of a [n_meas:num_vals]-matrix
    jack_evs = Array{Float64}(undef, N_blocks, num_vals)
    for i = 1:N_blocks
        temp_means    = mean(vcat(mean_vals_foo[1:(i-1)*b_size], mean_vals_foo[i*b_size+1:b_size*N_blocks]))
        temp_corr     = mean(vcat(corr_mats[1:(i-1)*b_size], corr_mats[i*b_size+1:b_size*N_blocks]))
        jack_evs[i,:] = eigen(temp_corr-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    end
    mean_corr  = mean(corr_mats)
    mean_means = mean(mean_vals_foo)
    mean_evs   = eigen(mean_corr - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end

# bla = collect(1:5);
# [bla[i] * bla[j]  for i = 1:5, j = 1:5]
# bla = [mean_vals[i,:] for i in eachindex(mean_vals[:,1])];
# size(bla)
# vcat(bla, bla)
# mean(bla)
# mean(mean_vals[:,5])
# bla[1]

# function cheat_log(x, alternative_result)
#     if x > 0.0
#         return log(x)
#     else
#         return alternative_result
#     end
# end

function jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1, corr_mats_t2, mean_vals, b_size)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    mean_vals_foo = [mean_vals[i,:] for i = 1:Int(b_size*N_blocks)] # just to have an array of arrays instead of a [n_meas:num_vals]-matrix
    jack_masses = Array{Float64}(undef, N_blocks, num_vals)
    for i = 1:N_blocks
        temp_means    = mean(vcat(mean_vals_foo[1:(i-1)*b_size], mean_vals_foo[i*b_size+1:b_size*N_blocks]))
        temp_corr_t1  = mean(vcat(corr_mats_t1[1:(i-1)*b_size], corr_mats_t1[i*b_size+1:b_size*N_blocks]))
        temp_corr_t2  = mean(vcat(corr_mats_t2[1:(i-1)*b_size], corr_mats_t2[i*b_size+1:b_size*N_blocks]))
        jack_evs_t1   = eigen(temp_corr_t1-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
        jack_evs_t2   = eigen(temp_corr_t2-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
        jack_masses[i,:] = log.(jack_evs_t1 ./ jack_evs_t2)
    end
    mean_means  = mean(mean_vals_foo)
    mean_evs_t1 = eigen(mean(corrs_t1) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    mean_evs_t2 = eigen(mean(corrs_t2) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    mass_means  = log.(mean_evs_t1 ./ mean_evs_t2)
    σ = [sqrt((N_blocks-1) * mean((jack_masses[:,i] .- mass_means[i]).^2 )) for i = 1:num_vals]
    return mass_means, σ
end

function jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1, corr_mats_t2, mean_vals, b_size, ev_nr)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    ev_nr    = num_vals - ev_nr +1  # So that Nr. 1 is the largest
    mean_vals_foo = [mean_vals[i,:] for i = 1:Int(b_size*N_blocks)] # just to have an array of arrays instead of a [n_meas:num_vals]-matrix
    jack_masses = Array{Float64}(undef, N_blocks)
    for i = 1:N_blocks
        temp_means   = mean(vcat(mean_vals_foo[1:(i-1)*b_size], mean_vals_foo[i*b_size+1:b_size*N_blocks]))
        temp_corr_t1 = mean(vcat(corr_mats_t1[1:(i-1)*b_size], corr_mats_t1[i*b_size+1:b_size*N_blocks]))
        temp_corr_t2 = mean(vcat(corr_mats_t2[1:(i-1)*b_size], corr_mats_t2[i*b_size+1:b_size*N_blocks]))
        jack_ev_t1   = eigen(temp_corr_t1-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values[ev_nr]
        jack_ev_t2   = eigen(temp_corr_t2-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values[ev_nr]
        jack_masses[i] = log(jack_ev_t1 / jack_ev_t2)
    end
    mean_means = mean(mean_vals_foo)
    mean_ev_t1 = eigen(mean(corr_mats_t1) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values[ev_nr]
    mean_ev_t2 = eigen(mean(corr_mats_t2) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values[ev_nr]
    mass_mean  = log(mean_ev_t1 / mean_ev_t2)
    σ = sqrt((N_blocks-1) * mean((jack_masses .- mass_mean).^2 ))
    return mass_mean, σ
end

# function jack_corr_mat_ev_mass_2pt(corr_mats_t1, corr_mats_t2, b_size, ev_nr)
#     N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
#     num_vals = size(corr_mats_t1[1],1)
#     ev_nr    = num_vals - ev_nr +1  # So that Nr. 1 is the largest 
#     jack_masses = Array{Float64}(undef, N_blocks)
#     for i = 1:N_blocks
#         temp_corr_t1 = mean(vcat(corr_mats_t1[1:(i-1)*b_size], corr_mats_t1[i*b_size+1:b_size*N_blocks]))
#         temp_corr_t2 = mean(vcat(corr_mats_t2[1:(i-1)*b_size], corr_mats_t2[i*b_size+1:b_size*N_blocks]))
#         jack_ev_t1   = eigen(temp_corr_t1).values[ev_nr]
#         jack_ev_t2   = eigen(temp_corr_t2).values[ev_nr]
#         jack_masses[i] = log(jack_ev_t1 / jack_ev_t2)
#     end
#     mean_ev_t1 = eigen(mean(corr_mats_t1)).values[ev_nr]
#     mean_ev_t2 = eigen(mean(corr_mats_t2)).values[ev_nr]
#     mass_mean  = log(mean_ev_t1 / mean_ev_t2)
#     σ = sqrt((N_blocks-1) * mean((jack_masses .- mass_mean).^2 ))
#     return mass_mean, σ
# end

function jack_conn_corr_mat_GEV(corr_mats_t2, corr_mats_t1, mean_vals, b_size)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    mean_vals_foo = [mean_vals[i,:] for i = 1:Int(b_size*N_blocks)] # just to have an array of arrays instead of a [n_meas:num_vals]-matrix
    jack_evs = Array{Float64}(undef, N_blocks, num_vals)
    for i = 1:N_blocks
        temp_means     = mean(vcat(mean_vals_foo[1:(i-1)*b_size], mean_vals_foo[i*b_size+1:b_size*N_blocks]))
        temp_means_mat = [temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        temp_corr_t2   = mean(vcat(corr_mats_t2[1:(i-1)*b_size], corr_mats_t2[i*b_size+1:b_size*N_blocks]))
        temp_corr_t1   = mean(vcat(corr_mats_t1[1:(i-1)*b_size], corr_mats_t1[i*b_size+1:b_size*N_blocks]))
        jack_evs[i,:]  = eigen(temp_corr_t2-temp_means_mat, temp_corr_t1-temp_means_mat).values
    end
    mean_corr_t2   = mean(corr_mats_t2)
    mean_corr_t1   = mean(corr_mats_t1)
    mean_means     = mean(mean_vals_foo)
    mean_means_mat = [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    mean_evs   = eigen(mean_corr_t2-mean_means_mat, mean_corr_t1-mean_means_mat).values
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end

