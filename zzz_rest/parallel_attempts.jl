using LinearAlgebra
using Statistics
# using Plots
using DelimitedFiles
using BenchmarkTools
using Base.Threads

# Crazy outdated⭕⭕⭕
# include("SU2_gaugefields.jl")
# include("SU2_observables.jl")
# include("SU2_updates.jl")

# Create an array of tasks
tasks = [
    ()->task_1()  # Define the first task as a function
    ()->task_2()  # Define the second task as a function
    # ... add more tasks as needed
]

# Create an array to store the results
results = Any[nothing for _ in 1:length(tasks)]

# Create an array of threads
threads = [Thread() for _ in 1:length(tasks)]

# Start each thread and assign a task to it
for i in 1:length(tasks)
    @threadcall threads[i] = ()->results[i] = tasks[i]()  # Assign the task to the thread
end

# Wait for all threads to finish
for thread in threads
    wait(thread)
end

# Process the results
for i in 1:length(results)
    process_result(results[i])
end







function task_1()
    # Task 1 code
    return result
end

function task_2()
    # Task 2 code
    return result
end

# Spawn tasks using Threads.@spawn
result1 = @spawn task_1()
result2 = @spawn task_2()

# Alternatively, you can explicitly specify the thread using Threads.@spawnat
result1 = @spawnat 2 task_1()  # Spawn task 1 on thread 2
result2 = @spawnat 4 task_2()  # Spawn task 2 on thread 4

# Retrieve the results
result1 = fetch(result1)
result2 = fetch(result2)




function partial_glue_1x2_corr(U::gaugefield_SU2, start_t, fin_t)
    corrs = zeros(fin_t-start_t+1)
    for Δt = start_t:fin_t
        c = 0.0
        for t = 1:U.N_t
            T = (t+Δt-1)%U.N_t +1   # T = mod1(t+Δt, U.N_t)
            for x1 = 1:U.N_x
                for x2 = 1:U.N_x
                    c += (loop_1x2(U,t,x1)[1] + loop_2x1(U,t,x1)[1]) * (loop_1x2(U,T,x2)[1] + loop_2x1(U,T,x2)[1])
                end
            end
        end
        corrs[Δt+1] = c
    end
    return corrs
end

function partial_glue_1x2_corr(U::gaugefield_SU2, times)
    corrs = Vector{Float64}(undef,length(times))
    i = 1
    for Δt in times
        c = 0.0
        for t = 1:U.N_t
            T = (t+Δt-1)%U.N_t +1   # T = mod1(t+Δt, U.N_t)
            for x1 = 1:U.N_x
                for x2 = 1:U.N_x
                    c += (loop_1x2(U,t,x1)[1] + loop_2x1(U,t,x1)[1]) * (loop_1x2(U,T,x2)[1] + loop_2x1(U,T,x2)[1])
                end
            end
        end
        corrs[i] = c
        i += 1
    end
    return corrs
end

function non_para_glue_1x2(U::gaugefield_SU2)
    chunks = Iterators.partition(0:U.N_t-1, div(U.N_t, nthreads()))
    tasks = map(chunks) do chunk
        partial_glue_1x2_corr(U,chunk)
    end
    return tasks
end

non_para_glue_1x2(test_field)

partial_glue_1x2_corr(test_field, 0:3)


# function sum_multi_good(a)
#     chunks = Iterators.partition(a, div(length(a), nthreads()))
#     tasks = map(chunks) do chunk
#         Threads.@spawn sum_single(chunk)
#     end
#     chunk_sums = fetch.(tasks)
#     return sum_single(chunk_sums)
# end

function parallel_glue_1x2(U::gaugefield_SU2)
    chunks = Iterators.partition(0:U.N_t-1, div(U.N_t, nthreads()))
    tasks = map(chunks) do chunk
        @spawn partial_glue_1x2_corr(U,chunk)
    end
    partial_corrs = fetch.(tasks)
    return partial_corrs
end

test_field = gaugefield_SU2(64,64,true);

@benchmark non_para_glue_1x2(test_field)    # (1.12 ± 0.04)s, 64x64: single result with 17.838 s
@benchmark parallel_glue_1x2(test_field)    # (1.80 ± 0.03)s for 8 threads, (1.60 ± 0.12)s for 4, (1.78 ± 0.33) for 2
#                                           # 64x64: single result for 4 threads: 26.838 s, for 8: 40.79 s, BUT WHYYY


function para_glue_1x2_corr(U::gaugefield_SU2)
    ntr = nthreads()            # number of threads
    @assert U.N_t % ntr == 0
    trlen = Int(U.N_t/ntr)           # thread length
    corrs = zeros(trlen,ntr)
    tasks = Array{Task}(undef,ntr)
    for i = 0:ntr-1
        tasks[i+1] = @spawn for Δt = (i*trlen):((1+i)*trlen-1)
            c = 0.0
            for t = 1:U.N_t
                T = (t+Δt-1)%U.N_t +1   # T = mod1(t+Δt, U.N_t)
                for x1 = 1:U.N_x
                    for x2 = 1:U.N_x
                        c += (loop_1x2(U,t,x1)[1] + loop_2x1(U,t,x1)[1]) * (loop_1x2(U,T,x2)[1] + loop_2x1(U,T,x2)[1])
                    end
                end
            end
            # println(bla)
            corrs[Δt+1,i+1] = c
        end 
    end
    partial_corrs = fetch.(tasks)
    # return reshape(corrs,1,:)[1,:]
    return partial_corrs
end

test_field = gaugefield_SU2(32,32,true);

isapprox(glue_1x2_corr(test_field),  para_glue_1x2_corr(test_field))

glue_1x2_corr(test_field)
test_corrs = para_glue_1x2_corr(test_field)
for i = 1:length(test_corrs)
    if test_corrs[i] != 0.0
        println(i)
    end
end

test_vec1 = Array(1:16)
test_vec = reshape(test_vec, 4, 4)
reshape(test_vec,1,:)[1,:] == Vector(1:16)
print(test_vec)
print(Vector(1:16))



# Several simulations simultaneously



function non_para_sims(L, iters)
    U1 = gaugefield_SU2(L, L, true)
    U2 = gaugefield_SU2(L, L, true)
    U3 = gaugefield_SU2(L, L, true)
    for i = 1:iters
        chess_metro!(U1,0.3,3.0)
    end
    for i = 1:iters
        chess_metro!(U2,0.3,4.0)
    end
    for i = 1:iters
        chess_metro!(U3,0.3,5.0)
    end
    return nothing
end

function para_sims(L, iters)
    U1 = gaugefield_SU2(L, L, true)
    U2 = gaugefield_SU2(L, L, true)
    U2 = gaugefield_SU2(L, L, true)
    task_1 = Threads.@spawn for i = 1:iters
        chess_metro!(U1,0.3,3.0)
    end
    task_2 = Threads.@spawn for i = 1:iters
        chess_metro!(U2,0.3,4.0)
    end
    task_3 = Threads.@spawn for i = 1:iters
        chess_metro!(U2,0.3,5.0)
    end
    fetch(task_1)
    fetch(task_2)
    fetch(task_3)
    return nothing
end

@benchmark non_para_sims(8,20) samples = 100000  #(48 ± 4)ms
@benchmark para_sims(8,20) samples = 100000      #(40 ± 22)ms


test_config = gaugefield_SU2(16,16,true)

function non_para(M)
    results = Vector{Vector}(undef,M)
    for i = 1:M
        results[i] = glue_nxm_corr(test_config, 1, i)
    end
    return results
end

function my_para(M)
    results = Vector{Vector}(undef,M)
    @threads for i = 1:M
        results[i] = glue_nxm_corr(test_config, 1, i)
    end
    return results
end

glue_nxm_corr(test_config, 1, 1)
@benchmark non_para(8)  # μ±σ = (173±9)ms
@benchmark alex_para(4)  # μ±σ = 

@benchmark my_para(8)   # μ±σ = (263±22)ms



# How about parallelize the entire flipping simulation


function test_sim(M)
    for i = 1:M
        N_t = 32
        N_x = 32
        β   = N_t*N_x/128
        hot = true
        ϵ   = 0.2 
        
        N_metro = 3     # N_metro-many Metropolois sweeps followed by...
        N_over = 1      # ...N_over-many overrelaxation sweeps will be performed...
        N_therm = 50    # ... for N_therm times...
        N_meas = 400/(N_metro+N_over) # ...and N_meas times,
        #                               # i.e. (N_therm + N_meas) ⋅ (N_metro + N_over) sweeps in total
        
        counter = 0        
        
        U = gaugefield_SU2(N_t, N_x, hot)
        
        for i = 1:N_therm
            for j = 1:N_metro
                chess_metro!(U,ϵ,β)
            end
            for j = 1:N_over
                chess_overrelax!(U)
            end
        end
        
        # for i = 1:N_meas 
        #     if i%(Int(N_meas/100)) == 1
        #         println(" ")
        #         println("We're already ", counter, "% deep in the simulation with N_t = $N_t, β = $β and ϵ = $ϵ !")
        #         counter += 1
        #     end
        #     for metro = 1:N_metro
        #         chess_metro!(U,ϵ,β)
        #     end
        #     for over = 1:N_over
        #         chess_overrelax!(U)
        #     end
        # end
    end
    return nothing
end

@benchmark test_sim(4) samples = 10      # 
@benchmark test_sim_para(4) samples = 10 # 
function test_sim_para(M)
    @threads for i = 1:M
        N_t = 32
        N_x = 32
        β   = N_t*N_x/128
        hot = true
        ϵ   = 0.2 
        
        N_metro = 3     # N_metro-many Metropolois sweeps followed by...
        N_over = 1      # ...N_over-many overrelaxation sweeps will be performed...
        N_therm = 50    # ... for N_therm times...
        N_meas = 400/(N_metro+N_over) # ...and N_meas times,
        #                               # i.e. (N_therm + N_meas) ⋅ (N_metro + N_over) sweeps in total
        
        counter = 0        
        
        U = gaugefield_SU2(N_t, N_x, hot)
        
        for i = 1:N_therm
            for j = 1:N_metro
                chess_metro!(U,ϵ,β)
            end
            for j = 1:N_over
                chess_overrelax!(U)
            end
        end
        
        # for i = 1:N_meas 
        #     if i%(Int(N_meas/100)) == 1
        #         println(" ")
        #         println("We're already ", counter, "% deep in the simulation with N_t = $N_t, β = $β and ϵ = $ϵ !")
        #         counter += 1
        #     end
        #     for metro = 1:N_metro
        #         chess_metro!(U,ϵ,β)
        #     end
        #     for over = 1:N_over
        #         chess_overrelax!(U)
        #     end
        # end
    end
    return nothing
end