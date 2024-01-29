include("test_SU2_gaugefields.jl")
include("test_hex.jl")
include("test_observables_square.jl")
include("test_observables_hex.jl")
include("test_smearing.jl")
include("test_updates.jl")
using BenchmarkTools

# benchmark two possibilities: alternative is mine, original according to Stephan
# ⭕⭕⭕ Not in use anymore! (Or soon to be)
function delta_S_gauge(U::gaugefield_SU2, t::Int64, x::Int64, dir::Int64, old_coeffs::Array, new_coeffs::Array, β::Float64)
    return β*0.5*real(tr((coeffs2grp(old_coeffs) - coeffs2grp(new_coeffs)) * staple_dag(U,t,x,dir)))
end

function local_metro!_alt(U::gaugefield_SU2, t::Int64, x::Int64, dir::Int64, step::Float64, β::Float64)
    old_coeffs = deepcopy(U.U[dir,t,x])
    new_coeffs = mult_SU2(ran_SU2(step), old_coeffs)
    # staple_d = staple_dag(U,t,x,dir)
    # S_old = β*0.5*real(tr(coeffs2grp(old_coeffs) * staple_d))
    # S_new = β*0.5*real(tr(coeffs2grp(new_coeffs) * staple_d))
    if rand() < exp(-delta_S_gauge(U, t, x, dir, old_coeffs, new_coeffs, β))
        U.U[dir,t,x] = new_coeffs
        U.acc_count += 1
    end
    return nothing
end

test_field = gaugefield_SU2(8, 8, true)
t = rand(1:8)
x = rand(1:8)
dir = rand(1:2)
@benchmark local_metro!(test_field, t, x, dir, 0.2, 1.0)        # (5.6 ± 15)μs
@benchmark local_metro!_alt(test_field, t, x, dir, 0.2, 1.0)    # (5.9 ± 16)μs





# Benchmark a multiplication on a struct for the 4 quaternionic coefficients of an SU2 matrix:

mutable struct coeffs_SU2
    coeffs::Vector{Float64}
    function coeffs_SU2(
        a0::Float64,
        a1::Float64,
        a2::Float64,
        a3::Float64
    )
    return new([a0,a1,a2,a3])
    end
end

function mult_SU2(a::coeffs_SU2, b::coeffs_SU2)
    c0 = a.coeffs[1]*b.coeffs[1] - a.coeffs[2]*b.coeffs[2] - a.coeffs[3]*b.coeffs[3] - a.coeffs[4]*b.coeffs[4] 
    c1 = a.coeffs[1]*b.coeffs[2] + a.coeffs[2]*b.coeffs[1] - a.coeffs[3]*b.coeffs[4] + a.coeffs[4]*b.coeffs[3]
    c2 = a.coeffs[1]*b.coeffs[3] + a.coeffs[2]*b.coeffs[4] + a.coeffs[3]*b.coeffs[1] - a.coeffs[4]*b.coeffs[2]
    c3 = a.coeffs[1]*b.coeffs[4] - a.coeffs[2]*b.coeffs[3] + a.coeffs[3]*b.coeffs[2] + a.coeffs[4]*b.coeffs[1]
    return coeffs_SU2(c0,c1,c2,c3)
end

function Base.:*(a::coeffs_SU2, b::coeffs_SU2)
    c0 = a.coeffs[1]*b.coeffs[1] - a.coeffs[2]*b.coeffs[2] - a.coeffs[3]*b.coeffs[3] - a.coeffs[4]*b.coeffs[4] 
    c1 = a.coeffs[1]*b.coeffs[2] + a.coeffs[2]*b.coeffs[1] - a.coeffs[3]*b.coeffs[4] + a.coeffs[4]*b.coeffs[3]
    c2 = a.coeffs[1]*b.coeffs[3] + a.coeffs[2]*b.coeffs[4] + a.coeffs[3]*b.coeffs[1] - a.coeffs[4]*b.coeffs[2]
    c3 = a.coeffs[1]*b.coeffs[4] - a.coeffs[2]*b.coeffs[3] + a.coeffs[3]*b.coeffs[2] + a.coeffs[4]*b.coeffs[1]
    return coeffs_SU2(c0,c1,c2,c3)
end


test_coeffs_1 = ran_SU2(0.5)
test_coeffs_2 = ran_SU2(0.5)
alt_coeffs_1 = coeffs_SU2(test_coeffs_1[1],test_coeffs_1[2],test_coeffs_1[3],test_coeffs_1[4])
alt_coeffs_2 = coeffs_SU2(test_coeffs_2[1],test_coeffs_2[2],test_coeffs_2[3],test_coeffs_2[4])
@benchmark mult_SU2(test_coeffs_1, test_coeffs_2)   # μ±σ = (66 ± 144)ns
@benchmark mult_SU2(alt_coeffs_1, alt_coeffs_2)     # μ±σ = (66 ± 158)ns
@benchmark alt_coeffs_1*alt_coeffs_2                # μ±σ = (79 ± 154)ns

# So the real advantage, given by the overload of the * operator, is not worth it

function Base.:*(a::Array{Float64}, b::Array{Float64})
    c0 = a[1]*b[1] - a[2]*b[2] - a[3]*b[3] - a[4]*b[4] 
    c1 = a[1]*b[2] + a[2]*b[1] - a[3]*b[4] + a[4]*b[3]
    c2 = a[1]*b[3] + a[2]*b[4] + a[3]*b[1] - a[4]*b[2]
    c3 = a[1]*b[4] - a[2]*b[3] + a[3]*b[2] + a[4]*b[1]
    return [c0,c1,c2,c3]
end

test_coeffs_1 = ran_SU2(0.5)
test_coeffs_2 = ran_SU2(0.5)
@benchmark mult_SU2(test_coeffs_1, test_coeffs_2)   # μ±σ = (150±344)ns
@benchmark test_coeffs_1 * test_coeffs_2            # μ±σ = (180±440)ns

# Even without struct the mult_SU2 is just faster. For which ever reason

test_coeffs_1 = ran_SU2(0.5)
test_coeffs_2 = ran_SU2(0.5)
test_coeffs_3 = ran_SU2(0.5)
@benchmark mult_SU2(test_coeffs_1, mult_SU2(test_coeffs_2, test_coeffs_3))  # μ±σ = (130±239)ns
@benchmark test_coeffs_1 * test_coeffs_2 * test_coeffs_3                    # μ±σ = (127±257)ns

test_coeffs_1 = ran_SU2(0.5)
test_coeffs_2 = ran_SU2(0.5)
test_coeffs_3 = ran_SU2(0.5)
test_coeffs_4 = ran_SU2(0.5)
@benchmark mult_SU2(test_coeffs_1, mult_SU2(test_coeffs_2, mult_SU2(test_coeffs_3, test_coeffs_4))) # μ±σ = (127±257)ns
@benchmark test_coeffs_1 * test_coeffs_2 * test_coeffs_3 * test_coeffs_4                            # μ±σ = (186±313)ns




# 

N_t = 8
N_x = 8
acc_count = [0]

function staple_dag_alt(U::Array, t::Int64, x::Int64, dir::Int64)
    a = zeros(4)
    b = zeros(4)
    t_plu = t%N_t +1                 
    x_plu = x%N_x +1                 
    t_min = (t + N_t -2)%N_t +1   
    x_min = (x + N_x -2)%N_x +1   

    if dir == 1
        a = mult_SU2(U[2,t_plu,x], mult_SU2(adjoint(U[1,t,x_plu]), adjoint(U[2,t,x])))
        b = mult_SU2(adjoint(U[2,t_plu,x_min]), mult_SU2(adjoint(U[1,t,x_min]), U[2,t,x_min]))
    elseif dir == 2
        a = mult_SU2(U[1,t,x_plu], mult_SU2(adjoint(U[2,t_plu,x]), adjoint(U[1,t,x])))
        b = mult_SU2(adjoint(U[1,t_min,x_plu]), mult_SU2(adjoint(U[2,t_min,x]), U[1,t_min,x]))
    end
    return coeffs2grp(a) + coeffs2grp(b)
end

function local_metro_alt_alt!(U::Array{Vector{Float64}, 3}, t::Int64, x::Int64, dir::Int64, step::Float64, β::Float64)
    old_coeffs = deepcopy(U[dir,t,x])
    new_coeffs = mult_SU2(ran_SU2(step), old_coeffs)
    staple_d = staple_dag_alt(U,t,x,dir)
    S_old = β*0.5*real(tr(coeffs2grp(old_coeffs) * staple_d))
    S_new = β*0.5*real(tr(coeffs2grp(new_coeffs) * staple_d))
    if rand() < exp(S_new-S_old)
        U[dir,t,x] = new_coeffs
        acc_count[1] += 1
    end
    return nothing
end

U = gaugefield_SU2(N_t, N_x, true)
U_alt = deepcopy(U.U)

t = rand(1:N_t)
x = rand(1:N_x)
dir = rand(1:2)
β = 1.0
step = 0.5*rand()

@benchmark local_metro!(U,t,x,dir,step,β)               # μ±σ = (17±44)μs, t ∈ [8μs,...,3.1ms]
@benchmark local_metro_alt_alt!(U_alt,t,x,dir,step,β)   # μ±σ = (19±85)μs, t ∈ [10μs,...8.5ms]




@benchmark mywrite("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\test_data\\test_mywrite.txt", 5.0) # (247 ± 177)μs
test_array = []
@benchmark push!(test_array,5.0)                                                                                             # (54 ± 316)ns



# mutable struct SU2_coeffs
#     a::Float64
#     b::Float64
#     c::Float64
#     d::Float64
#     function SU2_coeffs(
#         a::Float64,
#         b::Float64,
#         c::Float64,
#         d::Float64
#         )
#         return new(a,b,c,d)
#     end
# end

# function Base.:*(X::SU2_coeffs, Y::SU2_coeffs)
#     a = X.a*Y.a - X.b*Y.b - X.c*Y.c - X.d*Y.d 
#     b = X.a*Y.b + X.b*Y.a - X.c*Y.d + X.d*Y.c
#     c = X.a*Y.c + X.b*Y.d + X.c*Y.a - X.d*Y.b
#     d = X.a*Y.d - X.b*Y.c + X.c*Y.b + X.d*Y.a
#     return SU2_coeffs(a,b,c,d)
# end

# function ran_SU2(ϵ)
#     # r0 = rand([-1,1])
#     # r0 = 1.0
#     r1 = 2 * (rand()-0.5)
#     r2 = 2 * (rand()-0.5)
#     r3 = 2 * (rand()-0.5)
#     absr = sqrt(r1^2 + r2^2 + r3^2)
#     # return [r0*sqrt(1-ϵ^2), ϵ*r1/absr, ϵ*r2/absr, ϵ*r3/absr]
#     return SU2_coeffs(sqrt(1-ϵ^2), ϵ*r1/absr, ϵ*r2/absr, ϵ*r3/absr)
# end

coeffs1 = ran_SU2(rand())
coeffs2 = ran_SU2(rand())

coeffs1 * coeffs2

@benchmark coeffs1 * coeffs2        # μ±σ = (34±38)ns

coeffs11 = [coeffs1.a,coeffs1.b,coeffs1.c,coeffs1.d]
coeffs22 = [coeffs2.a,coeffs2.b,coeffs2.c,coeffs2.d]

@benchmark coeffs1 * coeffs2                # μ±σ = (34±38)ns
@benchmark mult_SU2(coeffs11, coeffs22)     # μ±σ = (50±30)ns


N_t = N_x = 64
β = 6.0
acc = [0]
test_field = gaugefield_SU2(64,64,true)
@benchmark chess_metro!(test_field, 0.1)    # (20±2)ms

function test_measure_func()
    N_t = N_x = 64
    β = 6.0
    acc = [0]
    test_field = gaugefield_SU2(64,64,true)
    measure_loops(test_field, [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4]], 0, 0)
end

@benchmark test_measure_func()      # (20±2)ms for [[1,1]]
#                                   # (213±4)ms for [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4]]

test_field = gaugefield_SU2(64,64,true)
@benchmark measure_loops(test_field, [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4]], 0, 0) # (202 ± 7)ms

function ran_SU2(ϵ)
    r1 = 2 * (rand()-0.5)
    r2 = 2 * (rand()-0.5)
    r3 = 2 * (rand()-0.5)
    absr = sqrt(r1^2 + r2^2 + r3^2)
    return coeffs_SU2(sqrt(1-ϵ^2), ϵ*r1/absr, ϵ*r2/absr, ϵ*r3/absr)
    # r = 2 .* (rand(3) .- 0.5)
    # absr = sqrt(sum(r.^2))
    # r .*= ϵ/absr
    # return coeffs_SU2(sqrt(1-ϵ^2), r[1], r[2], r[3])
end

function ran_SU2_proposed(ϵ)
    # r1 = 2 * (rand()-0.5)
    # r2 = 2 * (rand()-0.5)
    # r3 = 2 * (rand()-0.5)
    # absr = sqrt(r1^2 + r2^2 + r3^2)
    # return coeffs_SU2(sqrt(1-ϵ^2), ϵ*r1/absr, ϵ*r2/absr, ϵ*r3/absr)
    r = 2 .* (rand(3) .- 0.5)
    absr = sqrt(sum(r.^2))
    r .*= ϵ/absr
    return coeffs_SU2(sqrt(1-ϵ^2), r[1], r[2], r[3])
end

@benchmark ran_SU2(0.5)             # (25±48)ns, median 21ns
@benchmark ran_SU2_proposed(0.5)    # (127±189)ns, median 103ns


########    why does measure_loops() take so long compared to chess_metro()    ########


function cum_loop_mat(U, loops)
    for loop in loops
        loop_mat(U, loop[1], loop[2])
    end
end

function measure_loops_bench(U, loops::Array)
    NT = size(U,2)
    NX = size(U,3)
    L = length(loops) #     ⬇ no stout(U), that's the trick!
    results = [tr.(loop_mat(U, loop[1], loop[2])) for loop in loops]
    mean_vals = [sum(results[i]) for i = 1:L] ./(NT*NX)
    return mean_vals
end

function measure_loops_bench(U, loops::Array, n_stout, ρ)
    NT = size(U,2)
    NX = size(U,3)
    L = length(loops)
    results = [tr.(loop_mat(stout(U,n_stout,ρ), loop[1], loop[2])) for loop in loops]
    mean_vals = [sum(results[i]) for i = 1:L] ./(NT*NX)
    return mean_vals
end

test_field = gaugefield_SU2(64,64,true)
loops = [[1,1], [1,2], [2,1]]
β = 6.0
acc = [0]

@benchmark chess_metro!(test_field,0.2) # (21 ± 2)ms
@benchmark cum_loop_mat(test_field,loops) #(1.9 ± 1.3) ms
@benchmark measure_loops_bench(test_field, loops, 0, 0.0)    # (63 ± 6)ms
@benchmark measure_loops_bench(test_field, loops)   # (1.8 ± 1.0)ms


#### ⇒ It's the hecking Stout smearing! Even if n_stout = 0

function stout_bench(U, n_stout, ρ)
    V = deepcopy(U)
    for i = 1:n_stout
        for μ = 1:2
            for trip = 1:2  # For the chequer board pattern
                for t = 1:N_t
                    for x = (1+mod(t+trip,2)):2:N_x
                        stap_link = staple(V,μ,t,x) * adjoint(V[μ,t,x])
                        stap_link = 0.5*ρ*stap_link
                        V[μ,t,x] = exp_traceless(stap_link) * V[μ,t,x]
                        # stap_link = coeffs2grp(staple(V,μ,t,x) * adjoint(V[μ,t,x]))
                        # V[μ,t,x] = grp2coeffs(exp(α*0.5*(stap_link - adjoint(stap_link)))) * V[μ,t,x]
                    end
                end
            end
        end
    end
    return V
end

@benchmark stout_bench(test_field, 0, 0.0) # (20 ± 2)ms

function deep_bench(U)
    return deepcopy(U)
end

@benchmark deep_bench(test_field)   # (20 ± 2)ms    ⇒    fair enough...
# That makes perfect sense! We had 3 loops with no smearing, i.e. we called 
# deepcopy() three times, hence the ∼ +60ms

function similar_bench(U)
    return similar(U)
end

@benchmark similar_bench(test_field) # (16 ± 39) μs

# It's not only slow, it's wrong! Let's correct that:

function stout_single_bench(U, ρ)
    NT = size(U,2)
    NX = size(U,3)
    V = similar(U)
    # for i = 1:n_stout
        for μ = 1:2
            for trip = 1:2  # For the chequer board pattern
                for t = 1:NT
                    for x = (1+mod(t+trip,2)):2:NX
                        stap_link = staple(U,μ,t,x) * adjoint(U[μ,t,x])
                        stap_link = 0.5*ρ*stap_link
                        V[μ,t,x] = exp_traceless(stap_link) * U[μ,t,x]
                        # stap_link = coeffs2grp(staple(V,μ,t,x) * adjoint(V[μ,t,x]))
                        # V[μ,t,x] = grp2coeffs(exp(α*0.5*(stap_link - adjoint(stap_link)))) * V[μ,t,x]
                    end
                end
            end
        end
    # end
    return V
end

@benchmark stout_single_bench(test_field, 0.1)  # (1.4 ± 1.1)ms

function stout_single_2_bench(U, ρ)
    NT = size(U,2)
    NX = size(U,3)
    V = similar(U)
    for μ = 1:2
        for t = 1:NT
            for x = 1:NX
                stap_link = staple(U,μ,t,x) * adjoint(U[μ,t,x])
                stap_link = 0.5*ρ*stap_link
                V[μ,t,x] = exp_traceless(stap_link) * U[μ,t,x]
                # stap_link = coeffs2grp(staple(V,μ,t,x) * adjoint(V[μ,t,x]))
                # V[μ,t,x] = grp2coeffs(exp(α*0.5*(stap_link - adjoint(stap_link)))) * V[μ,t,x]
            end
        end
    end
    return V
end

false ∈ isapprox.(stout_single_bench(test_field,0.1), stout_single_2_bench(test_field, 0.1))

# Corrected stout, let's benchmark again
function measure_loops_bench(U, loops::Array, n_stout, ρ)
    NT = size(U,2)
    NX = size(U,3)
    L = length(loops)
    results = [tr.(loop_mat(stout(U,n_stout,ρ), loop[1], loop[2])) for loop in loops]
    mean_vals = [sum(results[i]) for i = 1:L] ./(NT*NX)
    return mean_vals
end

measure_loops_bench(test_field,loops,0,0.0)

@benchmark measure_loops_bench(test_field, loops, 0, 0.0)   # (1.5 ± 0.8)ms

test_hex = hexfield_SU2(64, 64, true)
@benchmark loop_mat_hex(test_hex, 1, 2) # (0.9 ± 1.0)ms






function local_metro_bench!(U, μ, t, x, step)
    old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * old_coeffs
    staple_d = staple_dag(U,μ,t,x)
    S_old = β*0.5*tr(old_coeffs * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function local_metro_bench_2!(U, μ, t, x, step)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,t,x]
    staple_d = staple_dag(U,μ,t,x)
    S_old = β*0.5*tr(U[μ,t,x] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

@benchmark local_metro_bench!(test_field, 1, 1, 1, 0.2) # (3±4) μs
@benchmark local_metro_bench_2!(test_field, 1, 1, 1, 0.2) # (0.4±0.5) ns   # Brooo!


@benchmark chess_metro!(test_field,0.2) # (5±3)ms


# using Plots
arvis(D,n,r) = sqrt(r^2 + 2*pi*(n- (D-2)/24))
arvis_taylor(D,n,r) = r + pi*(n-(D-2)/24)/r

function plot_arvis(r_start, inter, resol, D, n_max)
    r_vals = collect(r_start:1/resol:inter)
    y_vals = [[arvis(D,n,r) for r in r_vals] for n in 0:n_max]
    image = plot(title = "Arvis Potentials Vₙ(R) in D = $D dim.", xlabel = "Distance R")
    for n = 0:n_max
        image = plot!(r_vals, y_vals[n+1], label = "V_$n")
    end    
    return image
end
function plot_arvis_taylor(r_start, inter, resol, D, n_max)
    r_vals = collect(r_start:1/resol:inter)
    y_vals = [[arvis_taylor(D,n,r) for r in r_vals] for n in 0:n_max]
    image = plot(title = "Arvis Potentials Vₙ(R) in D = $D dim., Taylor exp.", xlabel = "Distance R")
    for n = 0:n_max
        image = plot!(r_vals, y_vals[n+1], label = "V_$n")
    end    
    return image
end

# display(plot_arvis_taylor(0.1,5,100,26,2))
# display(plot_arvis(0,5,100,2,2))

function my_mod1(t, N_t)
    return (t-1)%N_t + 1
end

bla = rand(1:2^10)
@benchmark mod1(bla,64) # (8 ± 5) ns
@benchmark my_mod1(bla,64) # (23±70)ns  ❗ That has changed ❗




coeffs1 = ran_SU2(rand())
coeffs2 = ran_SU2(rand())
mat1 = coeffs2grp(coeffs1)
mat2 = coeffs2grp(coeffs2)
@benchmark coeffs1 * coeffs2    # (30 ± 40)ns
@benchmark mat1 * mat2          # (80 ± 40)ns


function metro_bench!(U, μ, t, x, step, β, acc)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,t,x]
    staple_d = staple_dag(U,μ,t,x)
    S_old = β*0.5*tr(U[μ,t,x] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function metro_bench_2!(U, μ, t, x, step)
    # old_coeffs = deepcopy(U[μ,t,x])
    new_coeffs = ran_SU2(step) * U[μ,t,x]
    staple_d = staple_dag(U,μ,t,x)
    S_old = β*0.5*tr(U[μ,t,x] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,t,x] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

β = 6.0
acc = [0]
test_field = gaugefield_SU2(32,32,true)
μ = rand(1:2)
t = rand(1:32)
x = rand(1:32)
@benchmark metro_bench!(test_field,μ,t,x,0.2,β,acc) #(630 ± 540)ns
@benchmark metro_bench_2!(test_field,μ,t,x,0.2) #(680 ± 620)ns


function auto_corr_bench_1(obs, t)
    M = length(obs)
    C = 0.0
    obs_mean = mean(obs)
    for i = 1:Int(M-t)
        C += (obs[i]-obs_mean)*(obs[Int(i+t)]-obs_mean)
    end
    return C/(M-t)
    # obs_m = obs .- mean(obs)
    # return mean(obs_m .* circshift(obs_m, t))
end

function auto_corr_bench_2(obs, t)
    # M = length(obs)
    # C = 0.0
    # obs_mean = mean(obs)
    # for i = 1:Int(M-t)
    #     C += (obs[i]-obs_mean)*(obs[Int(i+t)]-obs_mean)
    # end
    # return C/(M-t)
    obs_m = obs .- mean(obs)
    return mean(obs_m .* circshift(obs_m, t))
end

obs = rand(1000)
t = rand(1:1000)

@benchmark auto_corr_bench_1(obs,t) # (370 ± 70) ns
@benchmark auto_corr_bench_2(obs,t) # (12 ± 0.04) μs

H = hexfield_SU2(64,64,true)
@benchmark temp_gauge_hex(H)    # (960 ± 950) µs
@benchmark temp_gauge(H)        # (270 ± 460) µs




test_field = gaugefield_SU2(128, 128, true);
@benchmark loop_mat(test_field, 4, 4)   # for 32²: (0.9 ± 0.6)ms    for 128²: (14±2)ms
@benchmark [RT_loop(test_field, 4, 4, x, t) for x in 1:128, t in 1:128]   #for 32²: (0.8±0.5)   for 128²: (14±2)ms



function ran_mat_U2(ϵ)
    ϕ = pi*rand()
    return exp(im*ϕ/2) * coeffs2grp(ran_SU2(ϵ))
end

mat1 = ran_U2(rand())
mat2 = ran_U2(rand())
abs(det(mat1))
mat1 * adjoint(mat1)

c1 = ran_SU2(rand())
c2 = ran_SU2(rand())

@benchmark c1 * c2      # (33±50)ns
@benchmark mat1 * mat2  # (83±35)ns


mutable struct coeffs_U2{T <: Real}
    a::T
    b::T
    c::T
    d::T
    ϕ::T
    function coeffs_U2(
        a::T,
        b::T,
        c::T,
        d::T,
        ϕ::T
        ) where {T <: Real}
        return new{T}(a,b,c,d,ϕ)
    end
end


# Multiply two U(2)-matrices whose coefficients are given in X and Y
function Base.:*(X::coeffs_U2, Y::coeffs_U2)
    a = X.a*Y.a - X.b*Y.b - X.c*Y.c - X.d*Y.d 
    b = X.a*Y.b + X.b*Y.a - X.c*Y.d + X.d*Y.c
    c = X.a*Y.c + X.b*Y.d + X.c*Y.a - X.d*Y.b
    d = X.a*Y.d - X.b*Y.c + X.c*Y.b + X.d*Y.a
    ϕ = mod(X.ϕ + Y.ϕ, 2*π)     # This is supposed to store the phase of the determinant
    return coeffs_U2(a,b,c,d,ϕ)
end

function ran_U2(ϵ)
    # r1 = 2 * (rand()-0.5)
    # r2 = 2 * (rand()-0.5)
    # r3 = 2 * (rand()-0.5)
    r1, r2, r3 = 2 .* (rand(3) .- 0.5)
    ϕ = ϵ * 2 * π * rand() 
    absr = sqrt(r1^2 + r2^2 + r3^2)
    return coeffs_U2(sqrt(1-ϵ^2), ϵ*r1/absr, ϵ*r2/absr, ϵ*r3/absr, ϕ)
end

c3 = ran_U2(rand())
c4 = ran_U2(rand())


@benchmark c3*c4                # (37±50) ns
@benchmark ran_U2(0.1*pi)       # (105±184) ns
@benchmark ran_mat_U2(0.1*pi)   # (660±863) ns



mutable struct coeffs_U2{T <: Number}
    a::T
    b::T
    c::T
    d::T
    function coeffs_U2(
        a::T,
        b::T,
        c::T,
        d::T
        ) where {T <: Number}
        return new{T}(a,b,c,d)
    end
end


# Multiply two U(2)-matrices whose coefficients are given in X and Y
function Base.:*(X::coeffs_U2, Y::coeffs_U2)
    a = X.a*Y.a - X.b*Y.b - X.c*Y.c - X.d*Y.d 
    b = X.a*Y.b + X.b*Y.a - X.c*Y.d + X.d*Y.c
    c = X.a*Y.c + X.b*Y.d + X.c*Y.a - X.d*Y.b
    d = X.a*Y.d - X.b*Y.c + X.c*Y.b + X.d*Y.a
    # ϕ = mod(X.ϕ + Y.ϕ, 2*π)     # This is supposed to store the phase of the determinant
    return coeffs_U2(a,b,c,d)#,ϕ)
end

function ran_U2(ϵ)
    # r1 = 2 * (rand()-0.5)
    # r2 = 2 * (rand()-0.5)
    # r3 = 2 * (rand()-0.5)
    r1, r2, r3 = 2 .* (rand(3) .- 0.5)
    phase = exp(ϵ*im*π*rand()) 
    absr = sqrt(r1^2 + r2^2 + r3^2)
    return coeffs_U2(phase*sqrt(1-ϵ^2), phase*ϵ*r1/absr, phase*ϵ*r2/absr, phase*ϵ*r3/absr)#, ϕ)
end

c1 = ran_U2(rand())
c2 = ran_U2(rand())


@benchmark c1*c2               # (45±60) ns