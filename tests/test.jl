################################################################################
# using Plots
using StatsBase
using BenchmarkTools
using LinearAlgebra
using DelimitedFiles

# Crazy outdated⭕⭕⭕
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("SU2_observables.jl")
include("SU2_observables_hex.jl")
include("SU2_updates.jl")
include("SU2_smearing.jl")
include("SU2_d_dim.jl")


#=
If any changes in the code have been made, one can run this entire file. All
of the contents of the files included above have a test function below, which
will be executed right below its definition.

Exceptions include:
    SU2_gaugefields.jl:
        gaugefield_SU2 (what is there to test?)
        mult_field_SU2 (mult_SU2 works though)

    SU2_observables.jl:
        action (will have to be checked more elaborately)
=#



########    Preliminaries        ########



# A lot of the functions use N_t, N_x . To avoid having
# to give these as arguments we need to specify them at some point.
# For testing purposes they shall be:
N_t = 64
N_x = 64
β   = 6.0
acc = [0]

# An SU(2) matrix has to be unitary and have determinant = 1.0
function is_SU2(mat::Matrix)
    return isapprox(det(mat), 1.0) * isapprox(mat*adjoint(mat), [1 0; 0 1])
end

function is_SU2(X::coeffs_SU2)
    return isapprox(X.a^2 + X.b^2 + X.c^2 + X.d^2, 1.0)
end

# Already defined elsewhere
# function Base.isapprox(X::coeffs_SU2, Y::coeffs_SU2)
#     return isapprox(get_array(X), get_array(Y))
# end



########    Gauge Fields        ########



function alg2grp_SU2_test()
    mat = alg2grp_SU2(rand(3).-0.5)
    @assert is_SU2(mat) "Something went wrong while testing alg2grp_SU2"
    return true
end

alg2grp_SU2_test()

function ran_SU2_test()
    coeffs = ran_SU2(rand()-0.5)
    @assert isapprox(sum(get_array(coeffs).^2), 1.0) "Something went wrong while testing ran_SU2"
    return true
end
# It was also tested that for ϵ → 0 no -id can appear anymore

ran_SU2_test()

function mult_SU2_test()
    coeffs_1 = ran_SU2(rand())
    coeffs_2 = ran_SU2(rand())
    mat_1 = coeffs2grp(coeffs_1)
    mat_2 = coeffs2grp(coeffs_2)
    prod_coeffs = coeffs_1 * coeffs_2
    prod_mat = coeffs2grp(prod_coeffs)
    @assert isapprox(prod_mat, mat_1*mat_2) "Something went wrong while testing mult_SU2"
    return true
end

mult_SU2_test()

function tr_test()
    coeffs = ran_SU2(rand())
    mat = coeffs2grp(coeffs)
    @assert tr(coeffs) == mat[1,1] + mat[2,2]
    return true
end

tr_test()

function det_test()
    coeffs = rand() * ran_SU2(rand())
    mat = coeffs2grp(coeffs)
    @assert isapprox(det(mat), det(coeffs))
    return true
end

det_test()

function proj2man!_test()
    coeffs = rand() * ran_SU2(rand())
    @assert is_SU2(proj2man!(coeffs))
    return true
end

proj2man!_test()

function coeffs2grp_test()
    coeffs = ran_SU2(rand())
    @assert is_SU2(coeffs2grp(coeffs)) "Something went wrong while testing coeffs2grp"
    return true
end

coeffs2grp_test()

function grp2coeffs_test()
    coeffs = rand() * ran_SU2(rand())
    mat = coeffs2grp(coeffs)
    new_coeffs = grp2coeffs(mat)
    # @assert isapprox(coeffs.a, new_coeffs.a) * isapprox(coeffs.b, new_coeffs.b) * isapprox(coeffs.c, new_coeffs.c) * isapprox(coeffs.d, new_coeffs.d)
    @assert isapprox(coeffs, new_coeffs)
    return true
end

grp2coeffs_test()

function adjoint_test()
    coeffs = ran_SU2(rand())
    adj_coeffs = adjoint(coeffs)
    mat = coeffs2grp(coeffs)
    adj_mat = coeffs2grp(adj_coeffs)
    @assert isapprox(adj_mat, adjoint(mat)) "Something went wrong while testing adjoint"
    return true
end

adjoint_test()

# function read_last_config_test()
#     test_config = gaugefield_SU2(4,4,true)
#     writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\test_data\\test_config.txt", test_config.U)
#     recon = read_last_config("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\test_data\\test_config.txt")
#     # ⭕ For some reason (test_config == recon) = false, hence check all attributes individually: ⭕
#     @assert test_config.U           == recon.U
#     @assert test_config.N_t         == recon.N_t
#     @assert test_config.N_x         == recon.N_x
#     @assert test_config.V           == recon.V
#     @assert test_config.acc_count   == recon.acc_count
#     # @assert test_config.hot         == recon.hot      # irrelevant, but true when function written that way
#     return true
# end

# read_last_config_test()



########        Observables        ########



function mod1_homebrew(t,N_t)
    return (t-1)%N_t +1
end

function mod1_homebrew_test()
    N_t = 4
    for i = 1:12
        @assert mod1(i,N_t) == mod1_homebrew(i,N_t) "Homebrew mod1 does not seem to work"
    end
    return true
end

mod1_homebrew_test()

# function test_mywrite()
#     rm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\test_data\\test_mywrite.txt")
#     for i = 1:10
#         mywrite("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\test_data\\test_mywrite.txt",i)
#     end
#     bla = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\test_data\\test_mywrite.txt")
#     @assert Int.(bla[:,1]) == collect(1:10)
#     return true
# end

# test_mywrite()

#=
function mywrite_last_conf_test()
    V = gaugefield_SU2(N_t,N_x,true)
    mywrite("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\test_data\\mywrite_config_test.txt", V)
    W = read_last_conf("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\test_data\\mywrite_config_test.txt")
    for μ = 1:2
        for t = 1:N_t
            for x = 1:N_x
                @assert isapprox(V[μ,t,x], W[μ,t,x])
            end
        end
    end
    return true
end


mywrite_last_conf_test()
=#

function plaq_test()
    test_field = gaugefield_SU2(N_t,N_x,true)
    plaquette = coeffs2grp(plaq(test_field, rand(1:N_t), rand(1:N_x)))
    @assert is_SU2(plaquette) "Something went wrong while testing plaq"
    return true
end

plaq_test()

function gen_rect_test()
    U = gaugefield_SU2(N_t,N_x,true)
    @assert isapprox(plaq(U,1,1) , gen_rect(U,1,1,1,1)) "Something wrong with gen_rect"
    return true
end

gen_rect_test()

function loop_mat_test()
    l_t = rand(1:N_t)
    l_x = rand(1:N_x)
    U = gaugefield_SU2(N_t,N_x,true)
    to_be_tested = loop_mat(U,l_t,l_x)
    mat = Matrix{coeffs_SU2}(undef,N_t,N_x)
    for i = 1:8
        for j = 1:8
            mat[i,j] = gen_rect(U,i,j,l_t,l_x)
        end
    end
    for i = 1:8
        for j = 1:8
            @assert isapprox(mat[i,j], to_be_tested[i,j]) "Something went wrong while testing loop_mat"
        end
    end
    return true
end

loop_mat_test()

# Still need to test measure_RT_loops_corrs...

function measure_RT_loops_test()
    N_t = N_x = 32
    test_field = gaugefield_SU2(N_t, N_x, true)
    T1 = rand(1:N_t)
    T2 = rand(1:N_t)
    R1 = rand(1:N_t)
    R2 = rand(1:N_t)
    loops = [[T1,R1], [T2,R2]]
    @assert isapprox(measure_RT_loops_corrs(test_field, loops, 0, 0.0)[2], measure_RT_loops(test_field, loops, 0, 0.0))
    return true
end

measure_RT_loops_test()



# auto_corr, auto_corr_norm, auto_corr_time



########        Updates        ########



function delta_S_gauge_test()
    test_field_1 = gaugefield_SU2(N_t,N_x,true)
    test_field_2 = deepcopy(test_field_1)
    μ = rand(1:2)
    t = rand(1:N_t)
    x = rand(1:N_x)
    test_field_2[μ,t,x] = ran_SU2(rand())
    true_dif =  action(test_field_2)-action(test_field_1)
    test_dif = delta_S_gauge(test_field_1, μ,t,x, test_field_1[μ,t,x], test_field_2[μ,t,x], β)
    @assert isapprox(true_dif, test_dif) "Something went wrong while testing delta_S_gauge"
    return true
end

delta_S_gauge_test()

function staple_dag_D_test()
    test_field = gaugefield_SU2(64, 64, true)
    μ = rand(1:2)
    t = rand(1:64)
    x = rand(1:64)
    @assert isapprox(staple_dag(test_field, μ,t,x), staple_dag_D(test_field, μ, [t,x]))
    return true
end

staple_dag_D_test()

acc = [0]
function local_metro!_test()
    test_field = gaugefield_SU2(N_t,N_x,true)
    test_test_field = deepcopy(test_field)
    for i = 1:50
        metro!(test_field, 1,1,1,0.1,β,acc)
    end
    @assert test_test_field != test_field "local_metro! didn't change the config it was given"
    @assert acc[1] != 0 "local_metro! didn't increment the acceptance"
    return true
end

local_metro!_test()

# function gaugefield_SU2_acc_count_test()
#     test_field = gaugefield_SU2(8, 8, true)
#     for i = 1:50
#         local_metro!(test_field, 1, 1, 1, 0.1)
#     end
#     @assert test_field.acc_count != 0
#     return true
# end

# gaugefield_SU2_acc_count_test()


# Only the indices of chess_metro! are tested, as the only remaining thing,
# local_metro!, has already been tested above.
function chess_indices_test()
    indices = []
    for μ = 1:2
        for t = 1:N_t
            for x = 1:N_x
                push!(indices, [μ,t,x])
            end
        end
    end

    test_indices = []
    for dir = 1:2
        for pass = 1:2
            for t = 1:N_t
                for x = (1+mod(t+pass,2)):2:N_x
                    push!(test_indices, [dir,t,x])
                end
            end
        end
    end

    @assert issetequal(indices, test_indices) "The indices of chess_metro! seem weird"
    return true
end

chess_indices_test()

function overrelax!_test()
    test_field = gaugefield_SU2(N_t,N_x,true)
    old_field = deepcopy(test_field)
    old_action = action(test_field)
    for i = 1:10
        t = rand(1:8)
        x = rand(1:8)
        μ = rand(1:2)
        overrelax!(test_field,μ,t,x)
    end
    @assert old_field != test_field "There was no update in overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action(test_field)) "overrelax! didn't leave the action invariant"
    return true
end

overrelax!_test()

acc = [0]
function lexico_overrelax!_test()
    test_field = gaugefield_SU2(N_t,N_x,true)
    old_field = deepcopy(test_field)
    old_action = action(test_field)
    lexico_overrelax!(test_field,acc)
    @assert old_field != test_field "There was no update in lexico_overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action(test_field)) "lexico_overrelax! didn't leave the action invariant"
    return true
end

lexico_overrelax!_test()

function chess_overrelax!_test()
    test_field = gaugefield_SU2(N_t,N_x,true)
    old_field = deepcopy(test_field)
    old_action = action(test_field)
    chess_overrelax!(test_field,acc)
    @assert old_field != test_field "There was no update in overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action(test_field)) "lexico_overrelax! didn't leave the action invariant"
    return true
end

lexico_overrelax!_test()



########        Smearing        ########



function staple_test()
    U = gaugefield_SU2(N_t,N_x,true)
    μ = rand(1:2)
    t = rand(1:N_t)
    x = rand(1:N_x)
    @assert isapprox(staple(U,μ,t,x), adjoint(staple_dag(U,μ,t,x)))
    return true
end

staple_test()

function traceless_test()
    coeffs = ran_SU2(rand())
    @assert tr(coeffs - adjoint(coeffs)) == 0.0
    return true
end

traceless_test()

function exp_traceless_test()
    coeffs = ran_SU2(rand())
    coeffs_tr_less = coeffs - adjoint(coeffs)
    @assert isapprox(exp_traceless(coeffs), grp2coeffs(exp(coeffs2grp(coeffs_tr_less))))
    return true
end

exp_traceless_test()

function old_VS_new_stout_test()
    function old_stout(U, n_stout, α)
        V = deepcopy(U)
        for i = 1:n_stout
            for μ = 1:2
                for trip = 1:2  # For the chequer board pattern
                    for t = 1:N_t
                        for x = (1+mod(t+trip,2)):2:N_x
                            # stap_link = staple(V,μ,t,x) * adjoint(V[μ,t,x])
                            # stap_link *= α*0.5
                            # V[μ,t,x] = exp_traceless(stap_link) * V[μ,t,x]
                            stap_link = coeffs2grp(staple(V,μ,t,x) * adjoint(V[μ,t,x]))
                            V[μ,t,x] = grp2coeffs(exp(α*0.5*(stap_link - adjoint(stap_link)))) * V[μ,t,x]
                        end
                    end
                end
            end
        end
        return V
    end
    function new_stout(U, n_stout, α)
        V = deepcopy(U)
        for i = 1:n_stout
            for μ = 1:2
                for trip = 1:2  # For the chequer board pattern
                    for t = 1:N_t
                        for x = (1+mod(t+trip,2)):2:N_x
                            stap_link = staple(V,μ,t,x) * adjoint(V[μ,t,x])
                            stap_link = α*0.5*stap_link
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
    test_field = gaugefield_SU2(N_t,N_x,true)
    old = old_stout(test_field, 1, 0.1)
    new = new_stout(test_field, 1, 0.1)
    for x = 1:N_x
        for t = 1:N_t
            for μ = 1:2
                @assert isapprox(old[μ,t,x], new[μ,t,x])
            end
        end
    end
    return true
end

old_VS_new_stout_test()


# ⭕⭕⭕ Not really a test function
function stout_test()
    U = gaugefield_SU2(N_t,N_x,true)
    V = stout(U,1,0.12)
    return true 
end

stout_test()



####    Hexagonal shenanigans    ####



function delta_S_gauge_hex_test()
    # NX = N_x>>1
    test_field_1 = hexfield_SU2(N_t, N_x, true)
    test_field_2 = deepcopy(test_field_1)
    μ = rand(1:3)
    t = rand(1:N_t)
    x = rand(1:N_x)
    test_field_2[μ,t,x] = ran_SU2(rand())
    true_dif =  action_hex(test_field_2)-action_hex(test_field_1)
    test_dif = delta_S_gauge_hex(test_field_1, μ,t,x, test_field_1[μ,t,x], test_field_2[μ,t,x],β)
    @assert isapprox(true_dif, test_dif) "Something went wrong while testing delta_S_gauge_hex"
    # print("μ = ", μ,", t = ",t,", x = ",x,",    true_dif = ", true_dif, ",    test_dif = ", test_dif)
    return true
end

delta_S_gauge_hex_test()

#=
for i = 1:100
    push!(actions, action_hex(test_hex))
    chess_metro_hex!(test_hex,0.2)
end

plot(actions)
=#

function overrelax_hex!_test()
    test_field = hexfield_SU2(N_t,N_x,true)
    old_field = deepcopy(test_field)
    old_action = action_hex(test_field)
    for i = 1:10
        t = rand(1:8)
        x = rand(1:8)
        μ = rand(1:3)
        overrelax_hex!(test_field,μ,t,x)
    end
    @assert old_field != test_field "There was no update in overrelax!, might wanna redo that test"
    @assert isapprox(old_action, action_hex(test_field)) "overrelax! didn't leave the action invariant"
    return true
end

overrelax_hex!_test()





################################################################################

println("Test file executed successfully!")

################################################################################






# test_mat = coeffs2grp(ran_SU2(2*(rand()-0.5)))
# is_SU2(test_mat)
# test_mat_exp = test_mat^rand()
# is_SU2(test_mat_exp)

# tr(test_mat)
# tr(test_mat_exp)

#=
function chess_metro!_time_series()
    N_t = 24
    N_x = 24
    β   = 5.0
    hot = true
    ϵ   = 0.1
    # acc = 0 
    N_therm = 1000
    actions = []

    U = gaugefield_SU2(N_t, N_x, hot)

    for i = 1:N_therm
        chess_metro!(U,ϵ,β)
        push!(actions, action(U, β))
    end

    acc_rate = U.acc_count/(N_t*N_x*2*N_therm)

    image = plot(1:length(actions), actions)
    image = plot!(
        size=(750,600),
        title = "Testing of action and chess Metropolis:
        times series of action, acc. rate = $acc_rate",
        legend = :false
        )
    return image
end

chess_metro!_time_series()


function chess_metro!_print_ind_stephan(N_t, N_x)
    for dir = 1:2
        for pass = 1:2
            for t = 1:N_t
                for x = (1+mod(t+pass,2)):2:N_x
                    println("[μ,t,x] = [", dir, ",", t, ",", x, "]")
                end
            end
        end
    end
    return nothing
end
=#

# function overrelax!_test()
#     test_field_1 = gaugefield_SU2(8, 8, true)
#     test_field_2 = deepcopy(test_field_1)
#     dir = rand(1:2)
#     t = rand(1:8)
#     x = rand(1:8)
#     staple_d = staple_dag(test_field_1,t,x,dir)
#     test_field_1.U[dir,t,x] = 
#     return true
# end

# bla = open("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\test_data\\test3.txt", "a")
# write(bla, "Surely this will work \n")
# close(bla)
# bla = open("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\test_data\\test3.txt", "a")
# write(bla, "Surely this will work ALSO \n")
# close(bla)

# function mywrite(sth, path)
#     bla = open(path, "a")
#     write(bla, sth)
#     close(bla)
# end


