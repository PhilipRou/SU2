include("observables_head.jl")



# Following three functions 
function plaq_12(U, x, y, t)
    NX = size(U,2)
    # NT = size(U,4)
    # t_p = mod1(t+1,NT)
    x_p = mod1(x+1,NX)
    y_p = mod1(y+1,NX)
    return U[1,x,y,t] * U[2,x_p,y,t] * adjoint(U[1,x,y_p,t])  * adjoint(U[2,x,y,t])
end

function plaq_13(U, x, y, t)
    NX = size(U,2)
    NT = size(U,4)
    t_p = mod1(t+1,NT)
    x_p = mod1(x+1,NX)
    # y_p = mod1(y+1,NX)
    return U[1,x,y,t] * U[3,x_p,y,t] * adjoint(U[1,x,y,t_p])  * adjoint(U[3,x,y,t])
end

function plaq_23(U, x, y, t)
    NX = size(U,2)
    NT = size(U,4)
    t_p = mod1(t+1,NT)
    # x_p = mod1(x+1,NX)
    y_p = mod1(y+1,NX)
    return U[2,x,y,t] * U[3,x,y_p,t] * adjoint(U[2,x,y,t_p])  * adjoint(U[3,x,y,t])
end

# Compute the Wilson gauge action of a 3-dim. cubic config
function action_cube(U, β)
    NX = size(U,2)
    NT = size(U,4)
    S = 3*2*NX*NX*NT    # 2: because SU(2), N_c = 2    AND    3: because we sum over 3 plaq's per site
    for t = 1:NT
        t_p = mod1(t+1,NT)
        for x = 1:NX
            x_p = mod1(x+1,NX)
            for y = 1:NX
                y_p = mod1(y+1,NX)
                # Don't use the plaq's because then significantly more mod1()-computations
                S -= tr(U[1,x,y,t] * U[3,x_p,y,t] * adjoint(U[1,x,y,t_p])  * adjoint(U[3,x,y,t]))
                S -= tr(U[2,x,y,t] * U[3,x,y_p,t] * adjoint(U[2,x,y,t_p])  * adjoint(U[3,x,y,t]))
                S -= tr(U[1,x,y,t] * U[2,x_p,y,t] * adjoint(U[1,x,y_p,t])  * adjoint(U[2,x,y,t]))
            end
        end
    end
    return S*β/2        # Again 2 because N_c = 2
end
        

# 
function loop_mat_cube(U, l_1, l_2)
    NX = size(U,2)
    NT = size(U,4)
    res_t = [coeffs_SU2(1.0,0.0,0.0,0.0) for x = 1:NX, y = 1:NX,  t = 1:NT]
    # res_x = [coeffs_SU2(1.0,0.0,0.0,0.0) for t = 1:NT, x = 1:NX, y = 1:NX]
    # res_y = [coeffs_SU2(1.0,0.0,0.0,0.0) for t = 1:NT, x = 1:NX, y = 1:NX]
    # t_arr = collect(1:NT)
    x_arr = collect(1:NX)
    y_arr = collect(1:NX)

    # The y-x-loops, stored in res_t:
    # for i = 1:l_1
    #     res_t = res_t .* U[2,x_arr,y_arr,:]
    #     circshift!(y_arr,-1)    
    # end
    # for i = 1:l_2
    #     res_t = res_t .* U[1,x_arr,y_arr,:]
    #     circshift!(x_arr,-1)   
    # end
    # for i = 1:l_1
    #     circshift!(y_arr,1)
    #     res_t = res_t .* adjoint.(U[2,x_arr,y_arr,:])
    # end
    # for i = 1:l_2
    #     circshift!(x_arr,1)
    #     res_t = res_t .* adjoint.(U[1,x_arr,y_arr,:])
    # end
    
    for i = 1:l_1
        res_t = res_t .* U[1,x_arr,y_arr,:]
        circshift!(x_arr,-1)    
    end
    for i = 1:l_2
        res_t = res_t .* U[2,x_arr,y_arr,:]
        circshift!(y_arr,-1)   
    end
    for i = 1:l_1
        circshift!(x_arr,1)
        res_t = res_t .* adjoint.(U[1,x_arr,y_arr,:])
    end
    for i = 1:l_2
        circshift!(y_arr,1)
        res_t = res_t .* adjoint.(U[2,x_arr,y_arr,:])
    end

    # # The x-t-loops, stored in res_y:
    # for i = 1:l_1
    #     res_y = res_y .* U[1,x_arr,:,t_arr]
    #     circshift!(x_arr,-1)    
    # end
    # for i = 1:l_2
    #     res_y = res_y .* U[3,x_arr,:,t_arr]
    #     circshift!(t_arr,-1)   
    # end
    # for i = 1:l_1
    #     circshift!(x_arr,1)
    #     res_y = res_y .* adjoint.(U[1,x_arr,:,t_arr])
    # end
    # for i = 1:l_2
    #     circshift!(t_arr,1)
    #     res_y = res_y .* adjoint.(U[3,x_arr,:,t_arr])
    # end

    # # The y-t-loops, stored in res_x:
    # for i = 1:l_1
    #     res_x = res_x .* U[2,:,y_arr,t_arr]
    #     circshift!(y_arr,-1)    
    # end
    # for i = 1:l_2
    #     res_x = res_x .* U[3,:,y_arr,t_arr]
    #     circshift!(t_arr,-1)   
    # end
    # for i = 1:l_1
    #     circshift!(y_arr,1)
    #     res_x = res_x .* adjoint.(U[2,:,y_arr,t_arr])
    # end
    # for i = 1:l_2
    #     circshift!(t_arr,1)
    #     res_x = res_x .* adjoint.(U[3,:,y_arr,t_arr])
    # end

    return res_t #, res_x, res_y
end

# bla = gaugefield_SU2_cube(8,8,true);
# res = tr.(loop_mat_cube(bla,1,1));
# res_true = [tr(plaq_12(bla,x,y,t)) for x = 1:8, y = 1:8, t = 1:8];
# isapprox(res,res_true)

# bla = gaugefield_SU2_cube(8,8,true);
# bli = deepcopy(bla);
# bli[1:2,:,:,1] = temp_gauge(bli[1:2,:,:,1],"SU2");
# mat_bla = loop_mat_cube(bla,3,4);
# mat_bli = loop_mat_cube(bli,3,4);
# isapprox(tr.(mat_bla), tr.(mat_bli))


# A function which measures everything one can measure (yet) using Wilson loops.
# The loops are specified in the array 'loops' in which tuples of [n_t, n_x] 
# are specified, i.e. loops = [[1,1], [1,2], [1,4], [3,4], ...]
function measure_RT_loops_corrs_cube(U, loops::Array, n_stout, ρ) 
    # NX = size(U,2)
    NT = size(U,4)
    L = length(loops)
    t_arr = collect(1:NT)
    results_t = Vector{Array{Float64}}(undef,L)
    V = stout_cube_timeslice(U,n_stout,ρ)
    for i = 1:L
        # "results_t" conatins matrices corresponding to different resp. loops;
        # these matrices contain the trace of
        # the loop at the resp. space-time point
        # mat_t = loop_mat_cube(stout_midpoint_cube_timeslice(U,n_stout,ρ), loops[i][1],loops[i][2])
        mat_t = loop_mat_cube(V, loops[i][1],loops[i][2])
        results_t[i] = tr.(mat_t)
    end
    # Now for each loop we want to obtain a column in "summed":
    # the row-index in that column is equal to the t-index of the time slice 
    # over which we sum our loop-observable
    summed_t = Matrix{Float64}(undef,NT,L)
    for i = 1:L
        summed_t[:,i] = [mean(results_t[i][:,:,t]) for t = 1:NT]
    end

    # The config-mean of each loop will be stored in "mean_vals_conf":
    # mean_vals_conf = Vector{Float64}(undef,L)
    # for l = 1:L
    #     mean_vals_conf[l] = mean(summed_t[:,l])
    # end
    mean_vals_conf = [mean(summed_t[:,l]) for l = 1:L]

    # A vector containing correlation matrices
    # corrs_t = Array{Float64}(undef,L,L,NT)
    # for τ = 1:NT
    #     corrs_t[:,:,τ] = [mean(summed_t[:,i] .* summed_t[circshift(t_arr,-τ),j]) for i = 1:L, j = 1:L] # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    # end
    corrs_t = [mean(summed_t[:,l1] .* summed_t[circshift(t_arr,-τ),l2]) for l1 = 1:L, l2 = 1:L, τ = 1:NT]

    return corrs_t, mean_vals_conf
end
# So the result of measure_RT_loops is an array containing corrs and mean_vals_conf:
#   corrs[:,:,t] is the correlation matrix at phys. time t
#   mean_vals_conf[i] is the mean value of the i-th loop (average of the lattice)

# same as measure_RT_loops_corrs_cube, but without cross-correlations
# (self-correlations only)
function measure_RT_loops_corrs_cube_selfonly(U, loops::Array, n_stout, ρ) 
    # NX = size(U,2)
    NT = size(U,4)
    L = length(loops)
    t_arr = collect(1:NT)
    results_t = Vector{Array{Float64}}(undef,L)
    V = stout_cube_timeslice(U, n_stout, ρ)
    for l = 1:L
        # "results_t" conatins matrices corresponding to different resp. loops;
        # these matrices contain the trace of
        # the loop at the resp. space-time point
        # mat_t = loop_mat_cube(stout_midpoint_cube_timeslice(U, n_stout, ρ), loops[l][1],loops[l][2])
        mat_t = loop_mat_cube(V, loops[l][1],loops[l][2])
        results_t[l] = tr.(mat_t)
    end
    # Now for each loop we want to obtain a column in "summed":
    # the row-index in that column is equal to the t-index of the time slice 
    # over which we sum our loop-observable
    summed_t = Matrix{Float64}(undef,NT,L)
    for l = 1:L
        summed_t[:,l] = [mean(results_t[l][:,:,t]) for t = 1:NT]
    end

    # The config-mean of each loop will be stored in "mean_vals_conf":
    # mean_vals_conf = Vector{Float64}(undef,L)
    # for l = 1:L
    #     mean_vals_conf[l] = mean(summed_t[:,l])
    # end
    mean_vals_conf = [mean(summed_t[:,l]) for l = 1:L]

    # A vector containing correlation matrices
    # corrs_t = Array{Float64}(undef,L,NT)
    corrs_t = [mean(summed_t[:,l] .* summed_t[circshift(t_arr,-τ),l]) for l = 1:L, τ = 1:NT]

    return corrs_t, mean_vals_conf
end

function s_wil_timeslice(U, x, y, t)
    return 1-tr(plaq_12(U,x,y,t))/2
end

function measure_s_wil(U, n_stout::Int, ρ)
    NX = size(U,2)
    NT = size(U,4)
    t_arr = collect(1:NT)
    V = stout_cube_timeslice(U,n_stout,ρ)
    results = [s_wil_timeslice(V,x,y,t) for x = 1:NX, y = 1:NX, t = 1:NT]
    summed = [mean(results[:,:,t]) for t = 1:NT]
    corrs = [mean(summed[:] .* summed[circshift(t_arr,-τ)]) for τ = 1:NT]
    return corrs, mean(summed)
end

function measure_s_wil(U, smear_nums::Vector, ρ)
    NX = size(U,2)
    NT = size(U,4)
    t_arr = collect(1:NT)
    smear_help = vcat([0],smear_nums)

    V = deepcopy(U)
    smeared_sums = Array{Float64}(undef, NT, length(smear_nums))
    smeared_means = Vector{Float64}(undef, length(smear_nums))
    for smear in eachindex(smear_nums)
        n_smear = smear_nums[smear] - smear_help[smear]
        V = stout_cube_timeslice(V,n_smear,ρ)
        results = [s_wil_timeslice(V,x,y,t) for x = 1:NX, y = 1:NX, t = 1:NT]
        summed = [mean(results[:,:,t]) for t = 1:NT]
        smeared_sums[:,smear] = summed
        smeared_means[smear]  = mean(summed)
    end

    corrs = [mean(smeared_sums[:,s1] .* smeared_sums[circshift(t_arr,-τ),s2]) for s1 in eachindex(smear_nums), s2 in eachindex(smear_nums), τ = 1:NT]
    return corrs, smeared_means
end

# bla = gaugefield_SU2_cube(8,8,true);
# bli = stout_cube_timeslice(bla, 0.24);
# ble = stout_cube_timeslice(bla, 2, 0.24);
# blub = measure_s_wil(bla, [0,1,2], 0.24);

# mean([s_wil_timeslice(ble, x, y, t) for x = 1:8, y = 1:8, t = 1:8])
# blub[2]

function clover_cube_timeslice(U,x,y,t)
    NX = size(U,2)
    # NT = size(U,4)
    x_p = mod1(x+1,NX)
    x_m = mod1(x-1,NX)
    y_p = mod1(y+1,NX)
    y_m = mod1(y-1,NX)
    tmq =  U[1,x,y,t] * U[2,x_p,y,t] * adjoint(U[1,x,y_p,t]) * adjoint(U[2,x,y,t])
    tmq += U[2,x,y,t] * adjoint(U[1,x_m,y_p,t]) * adjoint(U[2,x_m,y,t]) * U[1,x_m,y,t]
    tmq += adjoint(U[1,x_m,y,t]) * adjoint(U[2,x_m,y_m,t]) * U[1,x_m,y_m,t] * U[2,x,y_m,t]
    tmq += adjoint(U[2,x,y_m,t]) * U[1,x,y_m,t] * U[2,x_p,y_m,t] * adjoint(U[1,x,y,t])
    tmq = (tmq-adjoint(tmq)) # /8 /im
    return -tr(tmq*tmq)/128
end

function measure_clover(U, n_stout::Int, ρ)
    NX = size(U,2)
    NT = size(U,4)
    t_arr = collect(1:NT)
    V = stout_cube_timeslice(U,n_stout,ρ)
    results = [clover_cube_timeslice(V,x,y,t) for x = 1:NX, y = 1:NX, t = 1:NT]
    summed = [mean(results[:,:,t]) for t = 1:NT]
    corrs = [mean(summed[:] .* summed[circshift(t_arr,-τ)]) for τ = 1:NT]
    return corrs, mean(summed)
end

function measure_clover(U, smear_nums::Vector, ρ)
    NX = size(U,2)
    NT = size(U,4)
    t_arr = collect(1:NT)
    smear_help = vcat([0],smear_nums)

    V = deepcopy(U)
    smeared_sums = Array{Float64}(undef, NT, length(smear_nums))
    smeared_means = Vector{Float64}(undef, length(smear_nums))
    for smear in eachindex(smear_nums)
        n_smear = smear_nums[smear] - smear_help[smear]
        V = stout_cube_timeslice(V,n_smear,ρ)
        results = [clover_cube_timeslice(V,x,y,t) for x = 1:NX, y = 1:NX, t = 1:NT]
        summed = [mean(results[:,:,t]) for t = 1:NT]
        smeared_sums[:,smear] = summed
        smeared_means[smear]  = mean(summed)
    end

    corrs = [mean(smeared_sums[:,s1] .* smeared_sums[circshift(t_arr,-τ),s2]) for s1 in eachindex(smear_nums), s2 in eachindex(smear_nums), τ = 1:NT]
    return corrs, smeared_means
end

function measure_plaq_12(U, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,4)
    t_arr = collect(1:NT)
    V = stout_cube_timeslice(U,n_stout,ρ)
    results = [tr(plaq_12(V,x,y,t)) for x = 1:NX, y = 1:NX, t = 1:NT]
    summed = [mean(results[:,:,t]) for t = 1:NT]
    corrs = [mean(summed[:] .* summed[circshift(t_arr,-τ)]) for τ = 1:NT]
    return corrs, mean(summed)
end

function measure_s_wil_x_clover(U, smear_nums::Vector, ρ)
    NX = size(U,2)
    NT = size(U,4)
    t_arr = collect(1:NT)
    smear_help = vcat([0],smear_nums)

    V = deepcopy(U)
    L = length(smear_nums)
    smeared_sums = Array{Float64}(undef, NT, 2*L)
    smeared_means = Vector{Float64}(undef, 2*L)
    for smear = 1:L
        n_smear = smear_nums[smear] - smear_help[smear]
        V = stout_cube_timeslice(V,n_smear,ρ)
        results_s_wil  = [s_wil_timeslice(V,x,y,t) for x = 1:NX, y = 1:NX, t = 1:NT]
        results_clover = [clover_cube_timeslice(V,x,y,t) for x = 1:NX, y = 1:NX, t = 1:NT]
        summed_s_wil  = [mean(results_s_wil[:,:,t]) for t = 1:NT]
        summed_clover = [mean(results_clover[:,:,t]) for t = 1:NT]
        smeared_sums[:,smear]   = summed_s_wil
        smeared_sums[:,smear+L] = summed_clover
        smeared_means[smear]   = mean(summed_s_wil)
        smeared_means[smear+L] = mean(summed_clover)
    end

    corrs = [mean(smeared_sums[:,s1] .* smeared_sums[circshift(t_arr,-τ),s2]) for s1 = 1:2*L, s2 = 1:2*L, τ = 1:NT]
    return corrs, smeared_means
end

