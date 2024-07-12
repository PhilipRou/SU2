include("observables_head.jl")

# function mywrite(path, obs)
#     bla = open(path, "a")
#     write(bla, "$obs\n")
#     close(bla)
#     return nothing
# end

# function mywrite(path, obs::Vector)
#     bla = open(path, "a")
#     for i = 1:length(obs)
#         blu = obs[i]
#         write(bla, "$blu\t")
#     end
#     write(bla, "\n")
#     close(bla)
#     return nothing
# end

# function mywrite(path, obs)
#     bla = open(path,"a")
#     writedlm(bla, obs)
#     # writedlm(bla, "\n")
#     close(bla)
#     return nothing
# end

# function mywrite(path, obs::Array)
#     bla = open(path,"a")
#     writedlm(bla, transpose(obs))
#     # writedlm(bla, "\n")
#     close(bla)
#     return nothing
# end

function mywrite(path, obs::Array{coeffs_SU2}, N_x, N_t)
    mat = [get_array(obs[μ,x,t]) for μ = 1:2, x = 1:N_x, t = 1:N_t ]
    bla = open(path,"w")
    writedlm(bla, mat)
    close(bla)
    return nothing
end

function read_last_conf(path, N_x, N_t)
    V1 = readdlm(path)
    V1 = V1[1:size(V1,1)-1, 1:size(V1,2)-1]
    return reshape([coeffs_SU2(V1[i,1], V1[i,2], V1[i,3], V1[i,4]) for i = 1:size(V1,1)], (2,N_x,N_t))
end

#
function plaq(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    x_p = mod1(x+1, NX) # x%NX + 1
    t_p = mod1(t+1, NT) # t%NT + 1
    # return mult_SU2(U.U[2,t,x], mult_SU2(U.U[1,t,x_p], mult_SU2(adjoint(U.U[2,t_p,x]), adjoint(U.U[1,t,x]))))
    return U[1,x,t] * U[2,x_p,t] * adjoint(U[1,x,t_p]) * adjoint(U[2,x,t])
end

function action(U, β)
    NX = size(U,2)
    NT = size(U,3)
    S = 2*NX*NT   # later generalization: N_colour * NT * (NX)^d_s
    for t = 1:NT
        for x = 1:NX
            S -= real(tr(plaq(U,x,t)))
        end
    end
    return β*S/2    # later generalization: β*S/N_colour
end

# A method to calculate the action using the clover instead
# of the plaquette
function action_clover(U,β)
    NX = size(U,2)
    NT = size(U,3)
    S = 0
    for t = 1:NT
        t_p = mod1(t+1,NT)
        t_m = mod1(t-1,NT)
        for x = 1:NX
            x_p = mod1(x+1,NX)
            x_m = mod1(x-1,NX)
            tmq =  U[1,x,t] * U[2,x_p,t] * adjoint(U[1,x,t_p]) * adjoint(U[2,x,t])
            tmq += U[2,x,t] * adjoint(U[1,x_m,t_p]) * adjoint(U[2,x_m,t]) * U[1,x_m,t]
            tmq += adjoint(U[1,x_m,t]) * adjoint(U[2,x_m,t_m]) * U[1,x_m,t_m] * U[2,x,t_m]
            tmq += adjoint(U[2,x,t_m]) * U[1,x,t_m] * U[2,x_p,t_m] * adjoint(U[1,x,t])
            tmq = (tmq-adjoint(tmq)) / (8*im)
            S += β * real(tr(tmq*tmq)) / 4
        end
    end
    return S
end

# A (2×3)-Wilson loop written by hand for debugging purposes
function loop_2x3_square(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp1 = mod1(x+1,NX)
    xp2 = mod1(x+2,NX)
    tp1 = mod1(t+1,NT)
    tp2 = mod1(t+2,NT)
    tp3 = mod1(t+3,NT)
    
    # n   n   n   n   n   |  a   a   a   a   a 
    # 1   1   2   2   2   |  1   1   2   2   2 
    # x   xp1 xp2 xp2 xp2 |  xp1 x   x   x   x 
    # t   t   t   tp1 tp2 |  tp3 tp3 tp2 tp1 t 

    res = U[1,x,t] * U[1,xp1,t] * U[2,xp2,t] * U[2,xp2,tp1] * U[2,xp2,tp2]
    res = res * adjoint(U[1,xp1,tp3]) * adjoint(U[1,x,tp3]) * adjoint(U[2,x,tp2]) * adjoint(U[2,x,tp1]) * adjoint(U[2,x,t])
    return res
end

# Compute an (R×T)- Wilson loop at the coordinate (x,t)
function RT_loop(U, R, T, x, t)
    NX = size(U,2)
    NT = size(U,3)
    # loop = coeffs_Id_SU2()
    loop = U[1,x,t]
    x = mod1(x+1,NX)
    for i = 2:R
        loop *= U[1,x,t]
        x = mod1(x+1, NX)
    end
    for i = 1:T
        loop *= U[2,x,t]
        t = mod1(t+1,NT)
    end
    for i = 1:R
        x = mod1(x-1, NX)
        loop *= adjoint(U[1,x,t])
    end
    for i = 1:T
        t = mod1(t-1,NT)
        loop *= adjoint(U[2,x,t])
    end
    return loop
end

# Same as RT_loop(), but uses 
# C.Michael, NPB 259, 58, eq.(6)
function RT_loop_mike(U, R, T, x, t, avg_U)
    NX = size(U,2)
    NT = size(U,3)

    res = U[1,x,t]
    x = mod1(x+1,NX)
    for i = 2:R-1
        res *= avg_U[1,x,t]
        x = mod1(x+1, NX)
    end
    res *= U[1,x,t]
    x = mod1(x+1,NX)
    
    res *= U[2,x,t]
    t = mod1(t+1,NT)
    for i = 2:T-1
        res *= avg_U[2,x,t]
        t = mod1(t+1,NT)
    end
    res *= U[2,x,t]
    t = mod1(t+1,NT)

    x = mod1(x-1,NX)
    res *= adjoint(U[1,x,t])
    for i = 2:R-1
        x = mod1(x-1, NX)
        res *= adjoint(avg_U[1,x,t])
    end
    x = mod1(x-1,NX)
    res *= adjoint(U[1,x,t])

    t = mod1(t-1,NT)
    res *= adjoint(U[2,x,t])
    for i = 2:T-1
        t = mod1(t-1,NT)
        res *= adjoint(avg_U[2,x,t])
    end
    t = mod1(t-1,NT)
    res *= adjoint(U[2,x,t])

    return res
end

#=
# Returns an (Nₓ × Nₜ)-matrix whose entries carry the rectangular (R×T)-loop 
# at the respective (x,t)-points of the lattice. 
function loop_mat(U, R, T)
    NX = size(U,2)
    NT = size(U,3)

    res = [coeffs_Id_SU2() for x = 1:NX, t = 1:NT] # coeffs of the identity in (NT×NX)-matrix
    x_arr = collect(1:NX)
    t_arr = collect(1:NT)
    for i = 1:R
        res = res .* U[1,x_arr,t_arr]
        circshift!(x_arr,-1)    # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡
    end
    for i = 1:T
        res = res .* U[2,x_arr,t_arr]
        circshift!(t_arr,-1)   
    end
    for i = 1:R
        circshift!(x_arr,1)
        res = res .* adjoint.(U[1,x_arr,t_arr])
    end
    for i = 1:T
        circshift!(t_arr,1)
        res = res .* adjoint.(U[2,x_arr,t_arr])
    end
    return res
end
=#

# Returns an (Nₓ × Nₜ)-matrix whose entries carry the rectangular (R×T)-loop 
# at the respective (x,t)-points of the lattice. 
function loop_mat(U, R, T)
    NX = size(U,2)
    NT = size(U,3)
    return [RT_loop(U,R,T,x,t) for x = 1:NX, t = 1:NT]
end

# Returns an (Nₜ × Nₓ)-matrix just like loop_mat(), but uses 
# C.Michael, NPB 259, 58, eq.(6)
# for noise reduction. This can only be done for those links
# of the loop which are not part of the corners. Hence we will
# have to get an if-statement for each Wilson-line (i.e. the 
# four "if T/x > 2"-blocks). If the Wilson line is long enough,
# we can drag the first and last link of the Wilson line out of
# the inner for-loop and replace the links inside with avg_U[...].
# Regarding tests: only possible by comparing simulation results.
function loop_mat_mike(U, R, T, β)
    NX = size(U,2)
    NT = size(U,3)

    if R<=2 && T<=2
        return loop_mat(U,R,T)
    end
    
    stap_field = [staple(U,μ,x,t) for μ = 1:2, x = 1:NX, t = 1:NT]
    # d_field contains the variable d of the paper, evaluated at each coordinate 
    d_field = sqrt.(det.(stap_field))
    avg_U = [(besseli(2,β*d_field[μ,x,t]) / (besseli(1,β*d_field[μ,x,t]) * d_field[μ,x,t])) * stap_field[μ,x,t] for μ = 1:2, x = 1:NX, t = 1:NT]

    return [RT_loop_mike(U,R,T,x,t,avg_U) for x = 1:NX, t = 1:NT]
end

# A function just to get the mean values 'cause it's faster
function measure_RT_loops(U, loops::Array, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,3)
    # L = length(loops)
    results = [real.(tr.(loop_mat(stout(U,n_stout,ρ), loop[1], loop[2]))) for loop in loops]
    mean_vals = [sum(results[i]) for i = 1:length(loops)] ./(NX*NT)
    return mean_vals
end

# Same as measure_RT_loops(), but using loop_mat_mike() instead of loop_mat()
function measure_RT_loops_mike(U, loops::Array, β)
    NT = size(U,3)
    NX = size(U,2)
    # L = length(loops)
    results = [real.(tr.(loop_mat_mike(U, loops[i][1], loops[i][2], β))) for i in loops]
    mean_vals = [sum(results[i]) for i = 1:length(loops)] ./(NX*NT)
    return mean_vals
end

# A function which measures everything one can measure (yet) using Wilson loops.
# The loops are specified in the array 'loops' in which tuples of [n_t, n_x] 
# are specified, i.e. loops = [[1,1], [1,2], [1,4], [3,4], ...]
function measure_RT_loops_corrs(U, loops::Array, n_stout, ρ)
    NX = size(U,2)
    NT = size(U,3)
    L = length(loops)
    
    # "results" conatins matrices corresponding to different resp. loops;
    # these matrices contain the trace of the loop at the resp. 
    # space-time point
    results = [real.(tr.(loop_mat(stout(U,n_stout,ρ), loop[1], loop[2]))) for loop in loops]
    
    # Now for each loop we want to obtain a column in "summed":
    # the position in that column is equal to the t-index of the time slice 
    # over which we sum our loop-observable
    summed = [sum(results[i][:,t]) for t = 1:NT, i = 1:L]
    
    # Construct 'corrs', a vector containing correlation matrices
    t_arr = collect(1:NT)
    corrs = [sum(summed[:,i] .* summed[circshift(t_arr,-τ),j]) for i = 1:L, j = 1:L, τ = 1:NT] ./ (NX^2*NT) # 😡 circshift and circshift! DO NOT shift in opposite ways ANYMORE 😡

    mean_vals = [sum(summed[:,i]) for i = 1:L] ./ (NX*NT)

    return corrs, mean_vals
end

# So the result of measure_RT_loops is an array containing corrs and mean_vals:
#   corrs[:,:,t] is the correlation matrix at phys. time t
#   mean_vals[i] is the mean value of the i-th loop (average of the lattice)

# By edge we mean Wilson loops of the shape below (uncomment the shape for more clearness)
#   _
#  | |_
#  |_._|
function edge_loop_1(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    # n   n   n   a   n   a   a   a  
    # 1   1   2   1   2   1   2   2
    # x   xp  xpp xp  xp  x   x   x
    # t   t   t   tp  tp  tpp tp  t
    xp  = mod1(x+1, NX) # x%NX + 1
    tp  = mod1(t+1, NT) # t%NT + 1
    xpp = mod1(x+2, NX) 
    tpp = mod1(t+2, NT) 
    return U[1,x,t] * U[1,xp,t] * U[2,xpp,t] * adjoint(U[1,xp,tp]) * U[2,xp,tp] * adjoint(U[1,x,tpp]) * adjoint(U[2,x,tp]) * adjoint(U[2,x,t])
end

# test_field = gaugefield_SU2(32,32,true)
# X = rand(1:32)
# T = rand(1:32)
# edge_loop_1(test_field,X,T)

# By "L" we mean Wilson loops of the shape below (uncomment the shape for more clearness)
#   _
#  | |
#  | |_
#  |_._|
function L_loop_1(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    # n   n   n   a   n   n   |  a   a   a   a
    # 1   1   2   1   2   2   |  1   2   2   2
    # x   xp  xpp xp  xp  xp  |  x   x   x   x 
    # t   t   t   tp  tp  tpp |  tp3 tpp tp  t
    xp  = mod1(x+1, NX) # x%NX + 1
    xpp = mod1(x+2, NX) 
    tp  = mod1(t+1, NT) # t%NT + 1
    tpp = mod1(t+2, NT) 
    tp3 = mod1(t+3, NT) 
    return U[1,x,t] * U[1,xp,t] * U[2,xpp,t] * adjoint(U[1,xp,tp]) * U[2,xp,tp] * U[2,xp,tpp] * adjoint(U[1,x,tp3]) * adjoint(U[2,x,tpp]) * adjoint(U[2,x,tp]) * adjoint(U[2,x,t])
end

#
function rhomb_half_loop_square(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp  = mod1(x+1, NX) # x%NX + 1
    tp  = mod1(t+1, NT) # t%NT + 1
    # xpp = mod1(x+2, NX) 
    tpp = mod1(t+2, NT) 
    a = U[1,x,t]*U[2,xp,t] + U[2,x,t]*U[1,x,tp]
    b = adjoint(U[1,x,tp])*U[2,x,tp] + U[2,xp,tp]*adjoint(U[1,x,tpp])
    # c = a*b/det(a*b)
    # return c*adjoint(U[2,x,tp])*adjoint(U[2,x,t])
    return proj2man(a*b*adjoint(U[2,x,tp])*adjoint(U[2,x,t]))
end

#
function rhomb_loop_square(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    xp  = mod1(x+1, NX) # x%NX + 1
    tp  = mod1(t+1, NT) # t%NT + 1
    xpp = mod1(x+2, NX) 
    tpp = mod1(t+2, NT) 
    a = adjoint(U[2,x,t])*U[1,x,t] + U[1,x,tp]*adjoint(U[2,xp,t])
    b = U[1,xp,t]*U[2,xpp,t] + U[2,xp,t]*U[1,xp,tp]
    c = adjoint(U[1,xp,tp])*U[2,xp,tp] + U[2,xpp,tp]*adjoint(U[1,xp,tpp])
    d = adjoint(U[2,xp,tp])*adjoint(U[1,x,tp]) + adjoint(U[1,x,tpp])*adjoint(U[2,x,tp])
    # return a*b*c*d/det(a*b*c*d)
    return proj2man(a*b*c*d)
end



# clover action
# topological charge (even though non-existent in 2D SU(2))
# Polyakov loop???





########    U(2) Shenanigans    ########





# A neat definition of the topological charge utilizing the
# geometrical features of 2D QFT, yielding an integer-valued
# topological charge (which is not necessarily more 
# correct than the field theoretic definition below)
function top_charge_U2(U)
    NX = size(U,2)
    NT = size(U,3)
    return sum([imag(log(det(plaq(U, x, t)))) for x = 1:NX, t = 1:NT]) / 2 / π
end

# The field theoretic definition of the topological charge
function top_charge_U2_wil(U)
    NX = size(U,2)
    NT = size(U,3)
    q = 0
    for t = 1:NT
        t_p = mod1(t+1,NT)
        # t_m = mod1(t-1,NT)
        for x = 1:NX
            x_p = mod1(x+1,NX)
            # x_m = mod1(x-1,NX)
            plaq =  U[1,x,t] * U[2,x_p,t] * adjoint(U[1,x,t_p]) * adjoint(U[2,x,t])
            q += imag(tr(plaq))
        end
    end
    return q/2/π
end


#=
function one_link_integral_U2(A, β, n_sum_max, r_sum_max)
    Z = 0.0
    M = β^2 * A * adjoint(A)
    for n = 1:n_sum_max
        for r = 1:r_sum_max
            Z += tr(M)^n * det(M)^r / (factorial(big(n+2*r+1)) * factorial(big(n)) * factorial(big(r))^2)
        end
    end
    @assert imag(Z) < 10.0^(-10) "Z had non-negligible imaginary part for input A = $A"
    return real(Z)
end

let
    A = ran_U2(rand()) + ran_U2(rand())
    n_max = 11
    r_max = 11
    start_val_for_heatmap = 6
    Z_eval = [one_link_integral_U2(A,1.0,N,R) for N = 1:n_max, R = 1:r_max]
    display(heatmap(Z_eval, xticks = 1:2:n_max, yticks = 1:2:r_max))
    display(heatmap(
        start_val_for_heatmap:n_max,
        start_val_for_heatmap:r_max,
        Z_eval[start_val_for_heatmap:n_max,start_val_for_heatmap:r_max]
    ))
end


function one_link_expect_value_U2(U, μ, x, t, β, n_sum_max, r_sum_max, p_sum_max, s_sum_max)
    A = staple(U,μ,x,t)
    A_dagmod = coeffs_U2(conj(A.a), conj(A.b), conj(A.c), conj(A.d)) #  (-1)ⁱ⁺ʲ A†[̸j,̸i]
    # M = A*adjoint(A)
    # A_dagmod = conj.([A[2,2] -A[2,1]; -A[1,2] A[1,1]])
    trM = tr(A*adjoint(A))
    detA = det(A)
    detA_dag = adjoint(detA)
    detM = detA*detA_dag # det(A*adjoint(A))
    sum_normal = 0.0
    sum_dagmod = 0.0
    for n = 1:n_sum_max
        for r = 1:r_sum_max
            for p = 1:p_sum_max
                for s = 1:s_sum_max
                    # Z += tr(M)^n * det(M)^r / (factorial(big(n+2*r+1)) * factorial(big(n)) * factorial(big(r))^2)
                    denom = factorial(big(n+2*r+1))*factorial(big(p+2*s+1))*factorial(big(n))*factorial(big(p))*factorial(big(r))^2*factorial(big(s))^2
                    sum_normal += p * trM^(p+n-1) * detM^(r+s) / denom
                    sum_dagmod += s * trM^(p+n) * detA^(r+s) * detA_dag^(r+s-1) / denom
                end
            end
        end
    end
    return convert(ComplexF64,sum_normal)*A + convert(ComplexF64,sum_dagmod)*A_dagmod
end

N_x = N_t = L = 16
β = 1
U = gaugefield_U2(N_x, N_t, true);
for i = 1:200 chess_metro!(U,0.1,β,[0.0],"U2") end
μ = rand([1,2])
x, t = rand(Array(1:L), 2)
Us = []
V = deepcopy(U);
for i = 1:100000
    push!(Us, V[μ,x,t])
    metro!(V,μ,x,t,0.1,β,[0.0],"U2")
end
Us = real.(tr.(Us))
b_size = round(Int, 2*auto_corr_time(Us) + 1)
jackknife(Us,b_size)
U_exp = one_link_expect_value_U2(U,μ,x,t,β,6,6,6,6)
real(tr(U_exp))

# aaa = coeffs_U2(1.522347318101908 + 0.541280686115554im, 0.07674475624372301 + 0.14067257572225875im, 0.3587265810130666 + 0.15678099746753316im, 0.08346597802398936 - 0.08003358792801255im)
# zzz = 0.0
# mmm = aaa*adjoint(aaa)
# for n = 1:10
#     for r = 1:10
#         zzz += tr(mmm)^n * det(mmm)^r / factorial(big(n+2*r+1)) / factorial(big(n)) / factorial(big(r))^2
#     end
# end
# zzz


# function one_link_expect_value_SU2(U, μ, x, t, β)
#     stap = staple(U,μ,x,t)
#     d = sqrt(det(stap))
#     return 1/d*besseli(2,β*d)/besseli(1,β*d)*stap
# end

# N_x = N_t = L = 16
# β = 1
# U = gaugefield_SU2(N_x, N_t, true);
# for i = 1:200 chess_metro!(U,0.1,β,[0.0],"SU2") end
# μ = rand([1,2])
# x, t = rand(Array(1:L), 2)
# Us = []
# V = deepcopy(U);
# for i = 1:100000
#     push!(Us, V[μ,x,t])
#     metro!(V,μ,x,t,0.1,β,[0.0],"SU2")
# end
# Us = real.(tr.(Us))
# b_size = round(Int, 2*auto_corr_time(Us) + 1)
# jackknife(Us,b_size)
# U_exp = one_link_expect_value_SU2(U,μ,x,t,β)
# real(tr(U_exp))

=#