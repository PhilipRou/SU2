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
    N_blocks   = Int(div(length(obs), b_size, RoundDown))
    jack_means = Vector{Float64}(undef,N_blocks)
    blocked_means = [mean(obs[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    temp_means    = blocked_means[2:end]
    for i = 1:N_blocks-1
        jack_means[i] = mean(temp_means)
        temp_means[i] = blocked_means[i]
    end
    jack_means[N_blocks] = mean(temp_means)

    obs_mean = mean(obs)
    σ = sqrt((N_blocks-1) * mean((jack_means.-obs_mean).^2 ))
    return [obs_mean, σ]
end
### benchmark on rand(1000), b_size = 5: 12.5 μs

#=
function jackknife1(obs, b_size)
    N_blocks = Int(div(length(obs), b_size, RoundDown))
    jack_means = Vector{Float64}(undef,N_blocks)
    temp_array = obs[b_size+1:end]
    for i = 1:N_blocks-1
        jack_means[i] = mean(temp_array)
        temp_array[(i-1)*b_size+1:i*b_size] = obs[(i-1)*b_size+1:i*b_size]
    end
    jack_means[N_blocks] = mean(temp_array)
    obs_mean = mean(obs)
    σ = sqrt((N_blocks-1) * mean((jack_means.-obs_mean).^2 ))
    return [obs_mean, σ]
end
### benchmark on rand(1000), b_size = 5: 28.1 μs

function jackknife3(obs, b_size)
    N_blocks = Int(div(length(obs), b_size, RoundDown))
    jack_means = Vector{Float64}(undef,N_blocks)
    temp_array = @views obs[b_size+1:end]
    for i = 1:N_blocks-1
        jack_means[i] = mean(temp_array)
        temp_array[(i-1)*b_size+1:i*b_size] = @views obs[(i-1)*b_size+1:i*b_size]
    end
    jack_means[N_blocks] = mean(temp_array)
    obs_mean = mean(obs)
    σ = sqrt((N_blocks-1) * mean((jack_means.-obs_mean).^2 ))
    return [obs_mean, σ]
end
=#
### benchmark on rand(1000), b_size = 5: 18.6 μs

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

function jack_creutz(means::Array, a, b, c, d, b_size)
    N_blocks   = Int(div(length(means[:,1]), b_size, RoundDown))
    jack_means = Array{Float64}(undef,N_blocks)
    blocked_loop_means = [mean(means[(i-1)*b_size+1:i*b_size,j]) for i = 1:N_blocks, j in [a,b,c,d]]
    temp_means    = blocked_loop_means[2:end, :]
    for i = 1:N_blocks-1
        loop_means = [mean(temp_means[:,j]) for j = 1:4]
        jack_means[i] = loop_means[1]*loop_means[4] / (loop_means[2]*loop_means[3])
        temp_means[i,:] = blocked_loop_means[i,:]
    end
    loop_means = [mean(temp_means[:,j]) for j = 1:4]
    jack_means[N_blocks] = loop_means[1]*loop_means[4] / (loop_means[2]*loop_means[3])

    creutz_mean = mean(means[:,a]) * mean(means[:,d]) / (mean(means[:,b])  * mean(means[:,c]))
    σ = sqrt((N_blocks-1) * mean((jack_means.-creutz_mean).^2 ))
    return [creutz_mean, σ]
end

function jack_string_creutz(means::Array, a, b, c, d, b_size)
    N_blocks   = Int(div(length(means[:,1]), b_size, RoundDown))
    jack_means = Array{Float64}(undef,N_blocks)
    blocked_loop_means = [mean(means[(i-1)*b_size+1:i*b_size,j]) for i = 1:N_blocks, j in [a,b,c,d]]
    temp_means    = blocked_loop_means[2:end, :]
    for i = 1:N_blocks-1
        loop_means = [mean(temp_means[:,j]) for j = 1:4]
        jack_means[i] = -log(loop_means[1]*loop_means[4] / (loop_means[2]*loop_means[3]))
        temp_means[i,:] = blocked_loop_means[i,:]
    end
    loop_means = [mean(temp_means[:,j]) for j = 1:4]
    jack_means[N_blocks] = -log(loop_means[1]*loop_means[4] / (loop_means[2]*loop_means[3]))

    string_mean = -log(mean(means[:,a]) * mean(means[:,d]) / (mean(means[:,b])  * mean(means[:,c])))
    σ = sqrt((N_blocks-1) * mean((jack_means.-string_mean).^2 ))
    return [string_mean, σ]
end

function jack_conn_corr_self(corrs, loop_means, b_size)
    N_blocks  = Int(div(length(corrs), b_size, RoundDown))
    jack_cons = Vector{Float64}(undef,N_blocks)
    blocked_loop_means = [mean(loop_means[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    blocked_corr_means = [mean(corrs[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    temp_loop_means    = blocked_loop_means[2:end]
    temp_corr_means    = blocked_corr_means[2:end]
    for i = 1:N_blocks-1
        jack_cons[i]       = mean(temp_corr_means) - mean(temp_loop_means)^2
        temp_loop_means[i] = blocked_loop_means[i]
        temp_corr_means[i] = blocked_corr_means[i]
    end
    jack_cons[N_blocks]       = mean(temp_corr_means) - mean(temp_loop_means)^2

    con_mean = mean(corrs) - mean(loop_means)^2
    σ = sqrt((N_blocks-1) * mean((jack_cons .- con_mean).^2 ))
    return con_mean, σ
end

#=
function jack_conn_corr_self(corrs, loop_means, b_size)
    N_blocks = Int(div(length(corrs), b_size, RoundDown))
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
=#

### 🚧👷 Has to be made faster! 👷🚧
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

### 🚧👷 Has to be made faster! 👷🚧
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
    blocked_corrs = [mean(corr_mats[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    temp_corrs    = blocked_corrs[2:end]
    for i = 1:N_blocks-1
        jack_evs[i,:] = reverse(eigvals(mean(temp_corrs)))
        temp_corrs[i] = blocked_corrs[i]
    end
    jack_evs[N_blocks,:] = reverse(eigvals(mean(temp_corrs)))
    
    mean_evs = reverse(eigvals(mean(corr_mats)))
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end

#=
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
=#

function jack_conn_corr_mat_ev(corr_mats, mean_vals, b_size)
    N_blocks = Int(div(length(corr_mats), b_size, RoundDown))
    num_vals = size(corr_mats[1],1)
    jack_evs = Array{Float64}(undef, N_blocks, num_vals)
    blocked_means = Array{Float64}(undef, N_blocks, num_vals)
    for op = 1:num_vals
        blocked_means[:,op] =  [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    blocked_corrs =  [mean(corr_mats[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    temp_means =  blocked_means[2:end,:]
    temp_corrs =  blocked_corrs[2:end]
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        jack_evs[i,:]   = reverse(eigvals(mean(temp_corrs)-temp_means_mat))
        temp_means[i,:]  =  blocked_means[i,:]
        temp_corrs[i]    =  blocked_corrs[i]
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    jack_evs[N_blocks,:] = reverse(eigvals(mean(temp_corrs)-temp_means_mat))

    mean_corr  =  mean(corr_mats)
    mean_means =  [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_mean_mat =  [mean_means[op1] * mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    mean_evs   = reverse(eigvals(mean_corr - mean_mean_mat))
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end

function jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1, corr_mats_t2, mean_vals, b_size, n_masses)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    jack_evs_t1 = Array{Float64}(undef, N_blocks, num_vals)
    jack_evs_t2 = Array{Float64}(undef, N_blocks, num_vals)
    jack_masses = Array{Float64}(undef, N_blocks, n_masses) # n_masses instead of num_vals
    blocked_means = Array{Float64}(undef, N_blocks, num_vals)
    for op = 1:num_vals
        blocked_means[:,op] =  [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    blocked_corrs_t1 = [mean(corr_mats_t1[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    blocked_corrs_t2 = [mean(corr_mats_t2[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    temp_means =    blocked_means[2:end,:]
    temp_corrs_t1 = blocked_corrs_t1[2:end]
    temp_corrs_t2 = blocked_corrs_t2[2:end]
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        jack_evs_t1[i,:] = reverse(eigvals(mean(temp_corrs_t1)-temp_means_mat))
        jack_evs_t2[i,:] = reverse(eigvals(mean(temp_corrs_t2)-temp_means_mat))
        ### jack_masses[i,:] = log.(jack_evs_t1 ./ jack_evs_t2) ### What if arg of log is negative???
        for j = 1:n_masses
            mass_ind = j # num_vals - j +1
            bla = jack_evs_t1[i,mass_ind] / jack_evs_t2[i,mass_ind]
            if bla > 0.0
                jack_masses[i,j] = log(bla)
            else
                jack_masses[i,j] = NaN
            end
        end
        temp_means[i,:]  = blocked_means[i,:]
        temp_corrs_t1[i] = blocked_corrs_t1[i]
        temp_corrs_t2[i] = blocked_corrs_t2[i]
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    jack_evs_t1[N_blocks,:] = reverse(eigvals(mean(temp_corrs_t1)-temp_means_mat))
    jack_evs_t2[N_blocks,:] = reverse(eigvals(mean(temp_corrs_t2)-temp_means_mat))
    ### jack_masses[N_blocks,:] = log.(jack_evs_t1 ./ jack_evs_t2)
    for j = 1:n_masses
        mass_ind = j # num_vals - j +1
        bla = jack_evs_t1[N_blocks,mass_ind] / jack_evs_t2[N_blocks,mass_ind]
        if bla > 0.0
            jack_masses[N_blocks,j] = log(bla)
        else
            jack_masses[N_blocks,j] = NaN
        end
    end

    mean_means  = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_evs_t1 = reverse(eigvals(mean(corr_mats_t1) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]))
    mean_evs_t2 = reverse(eigvals(mean(corr_mats_t2) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]))
    mass_means = [NaN for j = 1:n_masses]
    for j = 1:n_masses
        mass_ind = j # num_vals - j +1
        bla = mean_evs_t1[mass_ind]/mean_evs_t2[mass_ind]
        if bla > 0.0
            mass_means[j] = log(bla)
        end
    end
    
    σ_m = [sqrt((N_blocks-1) * mean((jack_masses[:,i] .- mass_means[i]).^2 )) for i = 1:n_masses]
    σ_evs = [sqrt((N_blocks-1) * mean((jack_evs_t1[:,i] .- mean_evs_t1[i]).^2 )) for i = 1:num_vals]
    return mass_means, σ_m, mean_evs_t1, σ_evs
end

function jack_conn_corr_mat_ev_mass_2pt_allofem(corr_mats_array, mean_vals, b_size, n_masses)
    N_blocks = Int(div(length(corr_mats_array[1]), b_size, RoundDown))
    num_vals = size(corr_mats_array[1][1],1)
    T_max = length(corr_mats_array)
    jack_evs = Array{Float64}(undef, N_blocks, num_vals, T_max)
    jack_masses = Array{Float64}(undef, N_blocks, n_masses, T_max-1) # n_masses instead of num_vals
    blocked_means = Array{Float64}(undef, N_blocks, num_vals)
    for op = 1:num_vals
        blocked_means[:,op] = [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    temp_means = blocked_means[2:end,:]
    blocked_corrs = Vector{Vector{Matrix{Float64}}}(undef, T_max)
    temp_corrs = Vector{Vector{Matrix{Float64}}}(undef, T_max)
    for t = 1:T_max
        blocked_corrs[t] = [mean(corr_mats_array[t][(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
        temp_corrs[t] = blocked_corrs[t][2:end]
    end
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        for t = 1:T_max
            jack_evs[i,:,t] = reverse(eigvals(mean(temp_corrs[t]) - temp_means_mat))
        end

        for t = 1:T_max-1
            for j = 1:n_masses
                bla = jack_evs[i,j,t] / jack_evs[i,j,t+1]
                if bla > 0.0
                    jack_masses[i,j,t] = log(bla)
                else
                    jack_masses[i,j,t] = NaN
                end
            end
        end
        temp_means[i,:]  = blocked_means[i,:]
        for t = 1:T_max
            temp_corrs[t][i] = blocked_corrs[t][i]
        end
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    for t = 1:T_max
        jack_evs[N_blocks,:,t] = reverse(eigvals(mean(temp_corrs[t]) - temp_means_mat))
    end

    for t = 1:T_max-1
        for j = 1:n_masses
            bla = jack_evs[N_blocks,j,t] / jack_evs[N_blocks,j,t+1]
            if bla > 0.0
                jack_masses[N_blocks,j,t] = log(bla)
            else
                jack_masses[N_blocks,j,t] = NaN
            end
        end
    end

    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_means_mat = [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    mean_evs = Array{Float64}(undef, num_vals, T_max)
    for t = 1:T_max
        mean_evs[:,t] = reverse(eigvals(mean(corr_mats_array[t]) - mean_means_mat))
    end
    mass_means = Array{Float64}(undef, n_masses, T_max-1)
    for t = 1:T_max-1
        for j = 1:n_masses
            bla = mean_evs[j,t]/mean_evs[j,t+1]
            if bla > 0.0
                mass_means[j,t] = log(bla)
            else
                mass_means[j,t] = NaN
            end
        end
    end
    
    σ_m = Array{Float64}(undef, n_masses, T_max-1)
    for t = 1:T_max-1
        for j = 1:n_masses
            σ_m[j,t] = sqrt((N_blocks-1) * mean((jack_masses[:,j,t] .- mass_means[j,t]).^2 ))
        end
    end
    σ_evs = Array{Float64}(undef, num_vals, T_max)
    for t = 1:T_max
        for i = 1:num_vals
            σ_evs[i,t] = sqrt((N_blocks-1) * mean((jack_evs[:,i,t] .- mean_evs[i,t]).^2 )) 
        end
    end
    return transpose(mass_means), transpose(σ_m), transpose(mean_evs), transpose(σ_evs)
end


### NOT TO BE USED, WHY DID I WRITE THIS THING?!
function jack_conn_corr_mat_ev_mass_2pt_allofem_plateau(corr_mats_array, mean_vals, b_size, n_masses, massplateau_range_1, massplateau_range_2)
    N_blocks = Int(div(length(corr_mats_array[1]), b_size, RoundDown))
    num_vals = size(corr_mats_array[1][1],1)
    T_max = length(corr_mats_array)
    jack_evs = Array{Float64}(undef, N_blocks, num_vals, T_max)
    jack_masses = Array{Float64}(undef, N_blocks, n_masses, T_max-1) # n_masses instead of num_vals
    jack_plateau_masses_1 = Array{Float64}(undef, N_blocks)
    jack_plateau_masses_2 = Array{Float64}(undef, N_blocks)
    blocked_means = Array{Float64}(undef, N_blocks, num_vals)
    for op = 1:num_vals
        blocked_means[:,op] = [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    temp_means = blocked_means[2:end,:]
    blocked_corrs = Vector{Vector{Matrix{Float64}}}(undef, T_max)
    temp_corrs = Vector{Vector{Matrix{Float64}}}(undef, T_max)
    for t = 1:T_max
        blocked_corrs[t] = [mean(corr_mats_array[t][(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
        temp_corrs[t] = blocked_corrs[t][2:end]
    end
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        for t = 1:T_max
            jack_evs[i,:,t] = reverse(eigvals(mean(temp_corrs[t]) - temp_means_mat))
        end

        for t = 1:T_max-1
            for j = 1:n_masses
                bla = jack_evs[i,j,t] / jack_evs[i,j,t+1]
                if bla > 0.0
                    jack_masses[i,j,t] = log(bla)
                else
                    jack_masses[i,j,t] = NaN
                end
            end
        end

        jack_plateau_masses_1[i] = mean(jack_masses[i, 1, massplateau_range_1])
        jack_plateau_masses_2[i] = mean(jack_masses[i, 2, massplateau_range_2])

        temp_means[i,:]  = blocked_means[i,:]
        for t = 1:T_max
            temp_corrs[t][i] = blocked_corrs[t][i]
        end
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    for t = 1:T_max
        jack_evs[N_blocks,:,t] = reverse(eigvals(mean(temp_corrs[t]) - temp_means_mat))
    end

    for t = 1:T_max-1
        for j = 1:n_masses
            bla = jack_evs[N_blocks,j,t] / jack_evs[N_blocks,j,t+1]
            if bla > 0.0
                jack_masses[N_blocks,j,t] = log(bla)
            else
                jack_masses[N_blocks,j,t] = NaN
            end
        end
    end

    jack_plateau_masses_1[N_blocks] = mean(jack_masses[N_blocks, 1, massplateau_range_1])
    jack_plateau_masses_2[N_blocks] = mean(jack_masses[N_blocks, 2, massplateau_range_2])

    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_means_mat = [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    mean_evs = Array{Float64}(undef, num_vals, T_max)
    for t = 1:T_max
        mean_evs[:,t] = reverse(eigvals(mean(corr_mats_array[t]) - mean_means_mat))
    end
    mass_means = Array{Float64}(undef, n_masses, T_max-1)
    for t = 1:T_max-1
        for j = 1:n_masses
            bla = mean_evs[j,t]/mean_evs[j,t+1]
            if bla > 0.0
                mass_means[j,t] = log(bla)
            else
                mass_means[j,t] = NaN
            end
        end
    end
    plateaumass_means = [mean(mass_means[1,massplateau_range_1]), mean(mass_means[2,massplateau_range_2])]
    
    σ_evs = Array{Float64}(undef, num_vals, T_max)
    for t = 1:T_max
        for i = 1:num_vals
            σ_evs[i,t] = sqrt((N_blocks-1) * mean((jack_evs[:,i,t] .- mean_evs[i,t]).^2 )) 
        end
    end
    σ_m = Array{Float64}(undef, n_masses, T_max-1)
    for t = 1:T_max-1
        for j = 1:n_masses
            σ_m[j,t] = sqrt((N_blocks-1) * mean((jack_masses[:,j,t] .- mass_means[j,t]).^2 ))
        end
    end
    σ_m_plat = [sqrt((N_blocks-1) * mean((jack_plateau_masses_1 .- plateaumass_means[1]).^2 )) ]
    push!(σ_m_plat, sqrt((N_blocks-1) * mean((jack_plateau_masses_2 .- plateaumass_means[2]).^2 )))
    
    return transpose(mass_means), transpose(σ_m), transpose(mean_evs), transpose(σ_evs), transpose(plateaumass_means), transpose(σ_m_plat)
end



#=
function jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1, corr_mats_t2, mean_vals, b_size, n_masses)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    # n_masses_intern = num_vals - n_masses +1
    jack_evs_t1 = Array{Float64}(undef, N_blocks, num_vals)
    jack_evs_t2 = Array{Float64}(undef, N_blocks, num_vals)
    jack_masses = Array{Float64}(undef, N_blocks, n_masses) # n_masses instead of num_vals
    blocked_means = Array{Float64}(undef, N_blocks, num_vals)
    for op = 1:num_vals
        blocked_means[:,op] =  [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    blocked_corrs_t1 = [mean(corr_mats_t1[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    blocked_corrs_t2 = [mean(corr_mats_t2[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    temp_means =    blocked_means[2:end,:]
    temp_corrs_t1 = blocked_corrs_t1[2:end]
    temp_corrs_t2 = blocked_corrs_t2[2:end]
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        jack_evs_t1[i,:] = eigen(mean(temp_corrs_t1)-temp_means_mat).values
        jack_evs_t2[i,:] = eigen(mean(temp_corrs_t2)-temp_means_mat).values
        ### jack_masses[i,:] = log.(jack_evs_t1 ./ jack_evs_t2) ### What if arg of log is negative???
        for j = 1:n_masses
            mass_ind = num_vals - j +1
            bla = jack_evs_t1[i,mass_ind] / jack_evs_t2[i,mass_ind]
            if bla > 0.0
                jack_masses[i,j] = log(bla)
            else
                jack_masses[i,j] = NaN
            end
        end
        temp_means[i,:]  = blocked_means[i,:]
        temp_corrs_t1[i] = blocked_corrs_t1[i]
        temp_corrs_t2[i] = blocked_corrs_t2[i]
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    jack_evs_t1[N_blocks,:] = eigen(mean(temp_corrs_t1)-temp_means_mat).values
    jack_evs_t2[N_blocks,:] = eigen(mean(temp_corrs_t2)-temp_means_mat).values
    ### jack_masses[N_blocks,:] = log.(jack_evs_t1 ./ jack_evs_t2)
    for j = 1:n_masses
        mass_ind = num_vals - j +1
        bla = jack_evs_t1[N_blocks,mass_ind] / jack_evs_t2[N_blocks,mass_ind]
        if bla > 0.0
            jack_masses[N_blocks,j] = log(bla)
        else
            jack_masses[N_blocks,j] = NaN
        end
    end

    mean_means  = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_evs_t1 = eigen(mean(corr_mats_t1) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    mean_evs_t2 = eigen(mean(corr_mats_t2) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    mass_means = [NaN for j = 1:n_masses]
    for j = 1:n_masses
        mass_ind = num_vals - j +1
        bla = mean_evs_t1[mass_ind]/mean_evs_t2[mass_ind]
        if bla > 0.0
            mass_means[j] = log(bla)
        end
    end    
    
    σ_m = [sqrt((N_blocks-1) * mean((jack_masses[:,i] .- mass_means[i]).^2 )) for i = 1:n_masses]
    σ_evs = [sqrt((N_blocks-1) * mean((jack_evs_t1[:,i] .- mean_evs_t1[i]).^2 )) for i = 1:num_vals]
    return mass_means, σ_m, mean_evs_t1, σ_evs
end

function jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1, corr_mats_t2, mean_vals, b_size, ev_nr)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    ev_nr_intern = num_vals - ev_nr +1
    jack_masses = Array{Float64}(undef, N_blocks, ev_nr_intern)
    blocked_means = Array{Float64}(undef, N_blocks, ev_nr_intern)
    for op = 1:num_vals
        blocked_means[:,op] =  [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    blocked_corrs_t1 = [mean(corr_mats_t1[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    blocked_corrs_t2 = [mean(corr_mats_t2[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    temp_means    = blocked_means[2:end,:]
    temp_corrs_t1 = blocked_corrs_t1[2:end]
    temp_corrs_t2 = blocked_corrs_t2[2:end]
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        jack_evs_t1 = eigen(mean(temp_corrs_t1)-temp_means_mat).values[1:ev_nr_intern]
        jack_evs_t2 = eigen(mean(temp_corrs_t2)-temp_means_mat).values[1:ev_nr_intern]
        jack_masses[i,:] = log.(jack_evs_t1 ./ jack_evs_t2)
        temp_means[i,:]   =  blocked_means[i,:]
        temp_corrs_t1[i]  =  blocked_corrs_t1[i]
        temp_corrs_t2[i]  =  blocked_corrs_t2[i]
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    jack_evs_t1 = eigen(mean(temp_corrs_t1)-temp_means_mat).values[1:ev_nr_intern]
    jack_evs_t2 = eigen(mean(temp_corrs_t2)-temp_means_mat).values[1:ev_nr_intern]
    jack_masses[N_blocks,:] = log.(jack_evs_t1 ./ jack_evs_t2)

    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_evs_t1 = eigen(mean(corrs_t1) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    mean_evs_t2 = eigen(mean(corrs_t2) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    mass_means  = log.(mean_evs_t1 ./ mean_evs_t2)
    σ = [sqrt((N_blocks-1) * mean((jack_masses[:,i] .- mass_means[i]).^2 )) for i = 1:num_vals]
    return mass_means, σ
end

### 🚧👷 Has to be made faster! 👷🚧
### 🪄🔮 Has been  made faster! 🔮🪄
function jack_conn_corr_mat_ev(corr_mats, mean_vals, b_size)
    N_blocks = Int(div(length(corr_mats), b_size, RoundDown))
    num_vals = size(corr_mats[1],1)
    # mean_vals_foo = [mean_vals[i,:] for i = 1:Int(b_size*N_blocks)] # just to have an array of arrays instead of a [n_meas:num_vals]-matrix
    jack_evs = Array{Float64}(undef, N_blocks, num_vals)
    for i = 1:N_blocks
        # temp_means   = mean(vcat(mean_vals_foo[1:(i-1)*b_size], mean_vals_foo[i*b_size+1:b_size*N_blocks]))
        temp_means   = [mean(vcat(mean_vals[1:(i-1)*b_size, op], mean_vals[i*b_size+1:b_size*N_blocks, op])) for op = 1:num_vals]
        temp_corr     = mean(vcat(corr_mats[1:(i-1)*b_size], corr_mats[i*b_size+1:b_size*N_blocks]))
        jack_evs[i,:] = eigen(temp_corr-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    end
    # mean_means = mean(mean_vals_foo)
    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_corr  = mean(corr_mats)
    mean_evs   = eigen(mean_corr - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end

### 🚧👷 Has to be made faster! 👷🚧
### 🪄🔮 Has been  made faster! 🔮🪄
function jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1, corr_mats_t2, mean_vals, b_size)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    # mean_vals_foo = [mean_vals[i,:] for i = 1:Int(b_size*N_blocks)] # just to have an array of arrays instead of a [n_meas:num_vals]-matrix
    jack_masses = Array{Float64}(undef, N_blocks, num_vals)
    for i = 1:N_blocks
        # temp_means   = mean(vcat(mean_vals_foo[1:(i-1)*b_size], mean_vals_foo[i*b_size+1:b_size*N_blocks]))
        temp_means   = [mean(vcat(mean_vals[1:(i-1)*b_size, op], mean_vals[i*b_size+1:b_size*N_blocks, op])) for op = 1:num_vals]
        temp_corr_t1  = mean(vcat(corr_mats_t1[1:(i-1)*b_size], corr_mats_t1[i*b_size+1:b_size*N_blocks]))
        temp_corr_t2  = mean(vcat(corr_mats_t2[1:(i-1)*b_size], corr_mats_t2[i*b_size+1:b_size*N_blocks]))
        jack_evs_t1   = eigen(temp_corr_t1-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
        jack_evs_t2   = eigen(temp_corr_t2-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
        jack_masses[i,:] = log.(jack_evs_t1 ./ jack_evs_t2)
    end
    # mean_means = mean(mean_vals_foo)
    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_evs_t1 = eigen(mean(corrs_t1) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    mean_evs_t2 = eigen(mean(corrs_t2) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values
    mass_means  = log.(mean_evs_t1 ./ mean_evs_t2)
    σ = [sqrt((N_blocks-1) * mean((jack_masses[:,i] .- mass_means[i]).^2 )) for i = 1:num_vals]
    return mass_means, σ
end

### 🚧👷 Has to be made faster! 👷🚧
### 🪄🔮 Has been  made faster! 🔮🪄
function jack_conn_corr_mat_ev_mass_2pt(corr_mats_t1, corr_mats_t2, mean_vals, b_size, ev_nr)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    ev_nr    = num_vals - ev_nr +1  # So that Nr. 1 is the largest
    # mean_vals_foo = [mean_vals[i,:] for i = 1:Int(b_size*N_blocks)] # just to have an array of arrays instead of a [n_meas:num_vals]-matrix
    jack_masses = Array{Float64}(undef, N_blocks)
    for i = 1:N_blocks
        # temp_means   = mean(vcat(mean_vals_foo[1:(i-1)*b_size], mean_vals_foo[i*b_size+1:b_size*N_blocks]))
        temp_means   = [mean(vcat(mean_vals[1:(i-1)*b_size, op], mean_vals[i*b_size+1:b_size*N_blocks, op])) for op = 1:num_vals]
        temp_corr_t1 = mean(vcat(corr_mats_t1[1:(i-1)*b_size], corr_mats_t1[i*b_size+1:b_size*N_blocks]))
        temp_corr_t2 = mean(vcat(corr_mats_t2[1:(i-1)*b_size], corr_mats_t2[i*b_size+1:b_size*N_blocks]))
        jack_ev_t1   = eigen(temp_corr_t1-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values[ev_nr]
        jack_ev_t2   = eigen(temp_corr_t2-[temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values[ev_nr]
        jack_masses[i] = log(jack_ev_t1 / jack_ev_t2)
    end
    # mean_means = mean(mean_vals_foo)
    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_ev_t1 = eigen(mean(corr_mats_t1) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values[ev_nr]
    mean_ev_t2 = eigen(mean(corr_mats_t2) - [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]).values[ev_nr]
    mass_mean  = log(mean_ev_t1 / mean_ev_t2)
    σ = sqrt((N_blocks-1) * mean((jack_masses .- mass_mean).^2 ))
    return mass_mean, σ
end
=#

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

function jack_conn_corr_mat_GEV(corr_mats_t0, corr_mats_t2, mean_vals, b_size)
    N_blocks = Int(div(length(corr_mats_t0), b_size, RoundDown))
    num_vals = size(corr_mats_t0[1],1)
    jack_evs = Array{Float64}(undef, N_blocks, num_vals)
    blocked_means = Array{Float64}(undef, N_blocks, num_vals)
    for op = 1:num_vals
        blocked_means[:,op] = [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    blocked_corrs_t0 =  [mean(corr_mats_t0[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    blocked_corrs_t2 =  [mean(corr_mats_t2[(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
    temp_means    =  blocked_means[2:end,:]
    temp_corrs_t0 =  blocked_corrs_t0[2:end]
    temp_corrs_t2 =  blocked_corrs_t2[2:end]
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        jack_evs[i,:]   = reverse(eigvals(mean(temp_corrs_t2)-temp_means_mat, mean(temp_corrs_t0)-temp_means_mat))
        temp_means[i,:]  = blocked_means[i,:]
        temp_corrs_t0[i] = blocked_corrs_t0[i]
        temp_corrs_t2[i] = blocked_corrs_t2[i]
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    jack_evs[N_blocks,:] = reverse(eigvals(mean(temp_corrs_t2)-temp_means_mat, mean(temp_corrs_t0)-temp_means_mat))

    mean_means =  [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_means_mat = [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    # mean_corr_t2   = mean(corr_mats_t2)
    # mean_corr_t0   = mean(corr_mats_t0)
    mean_evs   = reverse(eigvals(mean(corr_mats_t2)-mean_means_mat, mean(corr_mats_t0)-mean_means_mat))
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end



function jack_conn_corr_mat_GEV_mass_2pt_allofem(t0, corr_mats_array_extern, mean_vals, b_size, n_masses)
    corr_mats_array = corr_mats_array_extern[t0:end]
    N_blocks = Int(div(length(corr_mats_array[1]), b_size, RoundDown))
    num_vals = size(corr_mats_array[1][1],1)
    T_max = length(corr_mats_array)
    jack_GEVs = Array{Float64}(undef, N_blocks, num_vals, T_max-1)
    jack_masses = Array{Float64}(undef, N_blocks, n_masses, T_max-2) # n_masses instead of num_vals
    blocked_means = Array{Float64}(undef, N_blocks, num_vals)
    for op = 1:num_vals
        blocked_means[:,op] = [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    temp_means = blocked_means[2:end,:]
    blocked_corrs = Vector{Vector{Matrix{Float64}}}(undef, T_max)
    temp_corrs = Vector{Vector{Matrix{Float64}}}(undef, T_max)
    for t = 1:T_max
        blocked_corrs[t] = [mean(corr_mats_array[t][(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
        temp_corrs[t] = blocked_corrs[t][2:end]
    end
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        temp_conn_corr_t0 = mean(temp_corrs[1]) - temp_means_mat
        for t = 1:T_max-1
            # jack_evs[i,:,t] = reverse(eigvals(mean(temp_corrs[t]) - temp_means_mat))
            jack_GEVs[i,:,t] = reverse(eigvals(mean(temp_corrs[t+1]) - temp_means_mat, temp_conn_corr_t0))
        end

        for t = 1:T_max-2
            for j = 1:n_masses
                bla = jack_GEVs[i,j,t] / jack_GEVs[i,j,t+1]
                if bla > 0.0
                    jack_masses[i,j,t] = log(bla)
                else
                    jack_masses[i,j,t] = NaN
                end
            end
        end
        temp_means[i,:]  = blocked_means[i,:]
        for t = 1:T_max
            temp_corrs[t][i] = blocked_corrs[t][i]
        end
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    temp_conn_corr_t0 = mean(temp_corrs[1]) - temp_means_mat
    for t = 1:T_max-1
        # jack_evs[N_blocks,:,t] = reverse(eigvals(mean(temp_corrs[t]) - temp_means_mat))
        jack_GEVs[N_blocks,:,t] = reverse(eigvals(mean(temp_corrs[t+1]) - temp_means_mat, temp_conn_corr_t0))
    end

    for t = 1:T_max-2
        for j = 1:n_masses
            bla = jack_GEVs[N_blocks,j,t] / jack_GEVs[N_blocks,j,t+1]
            if bla > 0.0
                jack_masses[N_blocks,j,t] = log(bla)
            else
                jack_masses[N_blocks,j,t] = NaN
            end
        end
    end

    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_means_mat = [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    mean_conn_corr_t0 = mean(corr_mats_array[1]) - mean_means_mat
    mean_GEVs = Array{Float64}(undef, num_vals, T_max-1)
    for t = 1:T_max-1
        # mean_evs[:,t] = reverse(eigvals(mean(corr_mats_array[t]) - mean_means_mat))
        mean_GEVs[:,t] = reverse(eigvals(mean(corr_mats_array[t+1]) - mean_means_mat, mean_conn_corr_t0))
    end
    mass_means = Array{Float64}(undef, n_masses, T_max-2)
    for t = 1:T_max-2
        for j = 1:n_masses
            bla = mean_GEVs[j,t]/mean_GEVs[j,t+1]
            if bla > 0.0
                mass_means[j,t] = log(bla)
            else
                mass_means[j,t] = NaN
            end
        end
    end
    
    σ_m = Array{Float64}(undef, n_masses, T_max-2)
    for t = 1:T_max-2
        for j = 1:n_masses
            σ_m[j,t] = sqrt((N_blocks-1) * mean((jack_masses[:,j,t] .- mass_means[j,t]).^2 ))
        end
    end
    σ_GEVs = Array{Float64}(undef, num_vals, T_max-1)
    for t = 1:T_max-1
        for i = 1:num_vals
            σ_GEVs[i,t] = sqrt((N_blocks-1) * mean((jack_GEVs[:,i,t] .- mean_GEVs[i,t]).^2 ))
        end
    end
    return transpose(mass_means), transpose(σ_m), transpose(mean_GEVs), transpose(σ_GEVs)
end



function jack_conn_corr_mat_GEV_mass_2pt_allofem_plateau(t0, corr_mats_array_extern, mean_vals, b_size, n_masses, massplateau_range_1, massplateau_range_2)
    corr_mats_array = corr_mats_array_extern[t0:end]
    N_blocks = Int(div(length(corr_mats_array[1]), b_size, RoundDown))
    num_vals = size(corr_mats_array[1][1],1)
    T_max = length(corr_mats_array)
    jack_GEVs = Array{Float64}(undef, N_blocks, num_vals, T_max-1)
    jack_masses = Array{Float64}(undef, N_blocks, n_masses, T_max-2) # n_masses instead of num_vals
    jack_plateau_masses_1 = Array{Float64}(undef, N_blocks)
    jack_plateau_masses_2 = Array{Float64}(undef, N_blocks)
    blocked_means = Array{Float64}(undef, N_blocks, num_vals)
    for op = 1:num_vals
        blocked_means[:,op] = [mean(mean_vals[(i-1)*b_size+1:i*b_size, op]) for i = 1:N_blocks]
    end
    temp_means = blocked_means[2:end,:]
    blocked_corrs = Vector{Vector{Matrix{Float64}}}(undef, T_max)
    temp_corrs = Vector{Vector{Matrix{Float64}}}(undef, T_max)
    for t = 1:T_max
        blocked_corrs[t] = [mean(corr_mats_array[t][(i-1)*b_size+1:i*b_size]) for i = 1:N_blocks]
        temp_corrs[t] = blocked_corrs[t][2:end]
    end
    for i = 1:N_blocks-1
        temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
        temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        temp_conn_corr_t0 = mean(temp_corrs[1]) - temp_means_mat
        for t = 1:T_max-1
            # jack_evs[i,:,t] = reverse(eigvals(mean(temp_corrs[t]) - temp_means_mat))
            jack_GEVs[i,:,t] = reverse(eigvals(mean(temp_corrs[t+1]) - temp_means_mat, temp_conn_corr_t0))
        end

        for t = 1:T_max-2
            for j = 1:n_masses
                bla = jack_GEVs[i,j,t] / jack_GEVs[i,j,t+1]
                if bla > 0.0
                    jack_masses[i,j,t] = log(bla)
                else
                    jack_masses[i,j,t] = NaN
                end
            end
        end

        jack_plateau_masses_1[i] = mean(jack_masses[i, 1, massplateau_range_1])
        jack_plateau_masses_2[i] = mean(jack_masses[i, 2, massplateau_range_2])

        temp_means[i,:]  = blocked_means[i,:]
        for t = 1:T_max
            temp_corrs[t][i] = blocked_corrs[t][i]
        end
    end
    temp_mean_means = [mean(temp_means[:,op]) for op = 1:num_vals]
    temp_means_mat  = [temp_mean_means[op1] * temp_mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    temp_conn_corr_t0 = mean(temp_corrs[1]) - temp_means_mat
    for t = 1:T_max-1
        # jack_evs[N_blocks,:,t] = reverse(eigvals(mean(temp_corrs[t]) - temp_means_mat))
        jack_GEVs[N_blocks,:,t] = reverse(eigvals(mean(temp_corrs[t+1]) - temp_means_mat, temp_conn_corr_t0))
    end

    for t = 1:T_max-2
        for j = 1:n_masses
            bla = jack_GEVs[N_blocks,j,t] / jack_GEVs[N_blocks,j,t+1]
            if bla > 0.0
                jack_masses[N_blocks,j,t] = log(bla)
            else
                jack_masses[N_blocks,j,t] = NaN
            end
        end
    end
    
    jack_plateau_masses_1[N_blocks] = mean(jack_masses[N_blocks, 1, massplateau_range_1])
    jack_plateau_masses_2[N_blocks] = mean(jack_masses[N_blocks, 2, massplateau_range_2])

    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_means_mat = [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    mean_conn_corr_t0 = mean(corr_mats_array[1]) - mean_means_mat
    mean_GEVs = Array{Float64}(undef, num_vals, T_max-1)
    for t = 1:T_max-1
        # mean_evs[:,t] = reverse(eigvals(mean(corr_mats_array[t]) - mean_means_mat))
        mean_GEVs[:,t] = reverse(eigvals(mean(corr_mats_array[t+1]) - mean_means_mat, mean_conn_corr_t0))
    end
    mass_means = Array{Float64}(undef, n_masses, T_max-2)
    for t = 1:T_max-2
        for j = 1:n_masses
            bla = mean_GEVs[j,t]/mean_GEVs[j,t+1]
            if bla > 0.0
                mass_means[j,t] = log(bla)
            else
                mass_means[j,t] = NaN
            end
        end
    end
    plateaumass_means = [mean(mass_means[1,massplateau_range_1]), mean(mass_means[2,massplateau_range_2])]
    
    σ_GEVs = Array{Float64}(undef, num_vals, T_max-1)
    for t = 1:T_max-1
        for i = 1:num_vals
            σ_GEVs[i,t] = sqrt((N_blocks-1) * mean((jack_GEVs[:,i,t] .- mean_GEVs[i,t]).^2 ))
        end
    end
    σ_m = Array{Float64}(undef, n_masses, T_max-2)
    for t = 1:T_max-2
        for j = 1:n_masses
            σ_m[j,t] = sqrt((N_blocks-1) * mean((jack_masses[:,j,t] .- mass_means[j,t]).^2 ))
        end
    end
    σ_m_plat = [sqrt((N_blocks-1) * mean((jack_plateau_masses_1 .- plateaumass_means[1]).^2 )) ]
    push!(σ_m_plat, sqrt((N_blocks-1) * mean((jack_plateau_masses_2 .- plateaumass_means[2]).^2 )))
    
    return transpose(mass_means), transpose(σ_m), transpose(mean_GEVs), transpose(σ_GEVs), transpose(plateaumass_means), transpose(σ_m_plat)
end

#=
function jack_conn_corr_mat_GEV(corr_mats_t2, corr_mats_t1, mean_vals, b_size)
    N_blocks = Int(div(length(corr_mats_t1), b_size, RoundDown))
    num_vals = size(corr_mats_t1[1],1)
    jack_evs = Array{Float64}(undef, N_blocks, num_vals)
    # mean_vals_foo = [mean_vals[i,:] for i = 1:Int(b_size*N_blocks)] # just to have an array of arrays instead of a [n_meas:num_vals]-matrix
    for i = 1:N_blocks
        # temp_means   = mean(vcat(mean_vals_foo[1:(i-1)*b_size], mean_vals_foo[i*b_size+1:b_size*N_blocks]))
        temp_means   = [mean(vcat(mean_vals[1:(i-1)*b_size, op], mean_vals[i*b_size+1:b_size*N_blocks, op])) for op = 1:num_vals]
        temp_means_mat = [temp_means[op1]*temp_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
        temp_corr_t2   = mean(vcat(corr_mats_t2[1:(i-1)*b_size], corr_mats_t2[i*b_size+1:b_size*N_blocks]))
        temp_corr_t1   = mean(vcat(corr_mats_t1[1:(i-1)*b_size], corr_mats_t1[i*b_size+1:b_size*N_blocks]))
        jack_evs[i,:] = eigen(temp_corr_t2-temp_means_mat, temp_corr_t1-temp_means_mat).values
    end
    # mean_means = mean(mean_vals_foo)
    mean_means = [mean(mean_vals[:,op]) for op = 1:num_vals]
    mean_means_mat = [mean_means[op1]*mean_means[op2] for op1 = 1:num_vals, op2 = 1:num_vals]
    mean_corr_t2   = mean(corr_mats_t2)
    mean_corr_t1   = mean(corr_mats_t1)
    mean_evs   = eigen(mean_corr_t2-mean_means_mat, mean_corr_t1-mean_means_mat).values
    σ = [sqrt((N_blocks-1) * mean((jack_evs[:,i] .- mean_evs[i]).^2 )) for i = 1:num_vals]
    return mean_evs, σ
end
=#
