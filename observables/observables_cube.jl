include("observables_head.jl")



# Following three functions 
function plaq_12(U, x, y, t)
    NX = size(U,2)
    NT = size(U,4)
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
function action_cube(U)
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
    for i = 1:l_1
        res_t = res_t .* U[2,x_arr,y_arr,:]
        circshift!(y_arr,-1)    
    end
    for i = 1:l_2
        res_t = res_t .* U[1,x_arr,y_arr,:]
        circshift!(x_arr,-1)   
    end
    for i = 1:l_1
        circshift!(y_arr,1)
        res_t = res_t .* adjoint.(U[2,x_arr,y_arr,:])
    end
    for i = 1:l_2
        circshift!(x_arr,1)
        res_t = res_t .* adjoint.(U[1,x_arr,y_arr,:])
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

# A function which measures everything one can measure (yet) using Wilson loops.
# The loops are specified in the array 'loops' in which tuples of [n_t, n_x] 
# are specified, i.e. loops = [[1,1], [1,2], [1,4], [3,4], ...]
function measure_RT_loops_corrs_cube(U, loops::Array) # ⭕ Implement smearing!
    NX = size(U,2)
    NT = size(U,4)
    L = length(loops)
    t_arr = collect(1:NT)
    # x_arr = collect(1:NX)
    # y_arr = collect(1:NX)
    # Transpose loop_mat because Julia is column-major, see below
    results_t = Vector{Array{Float64}}(undef,L)
    # results_x = Vector{Array{Float64}}(undef,L)
    # results_y = Vector{Array{Float64}}(undef,L)
    for i = 1:L
        # "results" conatins matrices corresponding to different resp. loops;
        # these matrices contain the trace (2× first entry of quaternion array) of
        # the loop at the resp. space-time point
        # mat_t, mat_x, mat_y = loop_mat_3d(U,loops[i][1],loops[i][2])
        mat_t = loop_mat_cube(U,loops[i][1],loop[i][2])
        results_t[i] = tr.(mat_t)
        # results_x[i] = tr.(mat_x)
        # results_y[i] = tr.(mat_y)
    end
    # Now for each loop we want to obtain a column in "summed":
    # the position in that column is equal to the t-index of the time slice 
    # over which we sum our loop-observable
    summed_t = Matrix{Float64}(undef,NT,L)
    # summed_x = Matrix{Float64}(undef,NX,L)
    # summed_y = Matrix{Float64}(undef,NX,L)
    for i = 1:L
        summed_t[:,i] = [sum(results_t[i][:,:,t]) for t = 1:NT]
        # summed_x[:,i] = [sum(results_x[i][:,x,:]) for x = 1:NX]
        # summed_y[:,i] = [sum(results_y[i][:,:,y]) for y = 1:NX]
    end

    # A vector containing correlation matrices
    corrs_t = Array{Float64}(undef,L,L,NT)
    # corrs_x = Array{Float64}(undef,L,L,NX)
    # corrs_y = Array{Float64}(undef,L,L,NX)
    for τ = 1:NT
        corrs_t[:,:,τ] = [sum(summed_t[:,i] .* summed_t[circshift(t_arr,-τ),j]) for i = 1:L, j = 1:L] # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    end
    # for τ = 1:NX
    #     corrs_x[:,:,τ] = [sum(summed_x[:,i] .* summed_x[circshift(x_arr,-τ),j]) for i = 1:L, j = 1:L] # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    #     corrs_y[:,:,τ] = [sum(summed_y[:,i] .* summed_y[circshift(y_arr,-τ),j]) for i = 1:L, j = 1:L] # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    # end
    corrs_t ./= NX^4*NT
    # corrs_x ./= NX^3*NT^2
    # corrs_y ./= NX^3*NT^2

    mean_vals_t = Vector{Float64}(undef,L)
    # mean_vals_x = Vector{Float64}(undef,L)
    # mean_vals_y = Vector{Float64}(undef,L)
    for i = 1:L
        mean_vals_t[i] = sum(summed_t[:,i])
        # mean_vals_x[i] = sum(summed_x[:,i])
        # mean_vals_y[i] = sum(summed_y[:,i])
    end
    mean_vals_t ./= NT*NX^2
    # mean_vals_x ./= NT*NX^2
    # mean_vals_y ./= NT*NX^2

    # return corrs_t, corrs_x, corrs_y, mean_vals_t, mean_vals_x, mean_vals_y
    return corrs_t, mean_vals_t
end
# So the result of measure_RT_loops is an array containing corrs and mean_vals:
#   corrs[:,:,t] is the correlation matrix at phys. time t
#   mean_vals[i] is the mean value of the i-th loop (average of the lattice)

