using Statistics
using DelimitedFiles
using LinearAlgebra
using Plots

λ1 = Complex.([0 1 0; 1 0 0; 0 0 0])
λ2 = Complex.([0 -im 0; im 0 0; 0 0 0])
λ3 = Complex.([1 0 0; 0 -1 0; 0 0 0])
λ4 = Complex.([0 0 1; 0 0 0; 1 0 0])
λ5 = Complex.([0 0 -im; 0 0 0; im 0 0])
λ6 = Complex.([0 0 0; 0 0 1; 0 1 0])
λ7 = Complex.([0 0 0; 0 0 -im; 0 im 0])
λ8 = Complex.([1 0 0; 0 1 0; 0 0 -2]) / sqrt(3)

λ0 = Complex.([1 0 0; 0 1 0; 0 0 1])
λ_zero = Complex.([0 0 0; 0 0 0; 0 0 0])

Λ = [λ0, λ1, λ2, λ3, λ4, λ5, λ6, λ7, λ8]

function ran_U3(ϵ)
    return exp(ϵ * im * sum((2*rand(9) .- 1.0).*Λ))
end

function gaugefield_U3(N_x, N_t, hot)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    if hot
        U = [ran_U3(2*rand()-1) for μ = 1:2, x = 1:N_x, t = 1:N_t]
    else
        U = [λ0 for μ = 1:2, x = 1:N_x, t = 1:N_t]
    end
    return U
end

function insta_U3(N_x, N_t, Q)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    w = 0.0
    if mod(Q%3,3) == 1
        w = -2*π/sqrt(3)
    elseif mod(Q%3,3) == 2
        w = 2*π/sqrt(3)
    end
    U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((-im*t*w)/(N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp((im*x*w)/(N_x) * λ8) for x = 1:N_x]
    return U
end

function plaq(U, x, t)
    NX = size(U,2)
    NT = size(U,3)
    x_p = mod1(x+1, NX) # x%NX + 1
    t_p = mod1(t+1, NT) # t%NT + 1
    # return mult_SU2(U.U[2,t,x], mult_SU2(U.U[1,t,x_p], mult_SU2(adjoint(U.U[2,t_p,x]), adjoint(U.U[1,t,x]))))
    return U[1,x,t] * U[2,x_p,t] * adjoint(U[1,x,t_p]) * adjoint(U[2,x,t])
end

function action_U3(U, β)
    NX = size(U,2)
    NT = size(U,3)
    S = 3*NX*NT   # later generalization: N_colour * NT * (NX)^d_s
    for t = 1:NT
        for x = 1:NX
            S -= real(tr(plaq(U,x,t)))
        end
    end
    return β*S/3    # later generalization: β*S/N_colour
end

function top_charge(U)
    NX = size(U,2)
    NT = size(U,3)
    return sum([imag(log(det(plaq(U, x, t)))) for x = 1:NX, t = 1:NT]) / 2 / π
end


################################################################################


function staple_dag(U, μ, x, t)
    NX = size(U,2)
    NT = size(U,3)
    a = λ0
    b = λ0
    x_p = mod1(x+1, NX) # x%NX +1                 
    t_p = mod1(t+1, NT) # t%NT +1                 
    x_m = mod1(x-1, NX) # (x + NX -2)%NX +1   
    t_m = mod1(t-1, NT) # (t + NT -2)%NT +1   

    # 🐌 More efficient: only use adjoint once 🐌 (but less human-readable, no?)
    if μ == 1
        a = U[2,x_p,t] * adjoint(U[1,x,t_p]) * adjoint(U[2,x,t])
        b = adjoint(U[2,x_p,t_m]) * adjoint(U[1,x,t_m]) * U[2,x,t_m]
    else #if μ == 2
        a = U[1,x,t_p] * adjoint(U[2,x_p,t]) * adjoint(U[1,x,t])
        b = adjoint(U[1,x_m,t_p]) * adjoint(U[2,x_m,t]) * U[1,x_m,t]
    end
    return a + b 
end

function metro_U3!(U, μ, x, t, step, β, acc)
    new_coeffs = ran_U3(step) * U[μ,x,t]
    staple_d = staple_dag(U,μ,x,t)
    S_old = β*0.5*real(tr(U[μ,x,t] * staple_d))
    S_new = β*0.5*real(tr(new_coeffs * staple_d))
    if rand() < exp(S_new-S_old)
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function chess_metro_U3!(U, step, β, acc)
    NX = size(U,2)
    NT = size(U,3)
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    metro_U3!(U,μ,x,t,step, β, acc)
                end
            end
        end
    end
    return nothing
end

function insta_update_U3!(U, β, acc, ΔQ)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U3(NX,NT,rand([-ΔQ,ΔQ])) .* U
    if rand() < exp(action_U3(U,β) - action_U3(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end

function insta_U3_tryout(N_x, N_t, up)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    w = 4*π/sqrt(3)
    Q = -1
    if up
        w = -4*π/sqrt(3)
        Q = 1
    end
    U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((-im*t*w)/(N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp((im*x*w)/(N_x) * λ8) for x = 1:N_x]
    return U
end

function insta_update_U3_tryout!(U, β, acc)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U3_tryout(NX,NT,rand(Bool)) .* U
    if rand() < exp(action_U3(U,β) - action_U3(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end

function cool_U3!(U, μ, x, t, step, acc)
    new_coeffs = ran_U3(step) * U[μ,x,t]
    staple_d = staple_dag(U,μ,x,t)
    S_old = real(tr(U[μ,x,t] * staple_d)) # ❗ technically missing global factors like β,
    S_new = real(tr(new_coeffs * staple_d)) # but we don't need them here
    if S_old < S_new
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function chess_cool_U3!(U, step, acc)
    NX = size(U,2)
    NT = size(U,3)
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    cool_U3!(U,μ,x,t,step,acc)
                end
            end
        end
    end
    return nothing
end



################################################################################



N_x     = 32
N_t     = N_x
β       = 16.0
hot     = true

N_therm = 500
N_up    = 1000
N_metro = 1
N_insta = 1
ϵ       = 0.05
acc_wish = 0.8
ΔQ      = 2


actions = []
charges = []


acc_metro_therm = [0]
acc_metro = [0]
acc_insta = [0]
count_therm = [0]
count = [0]

U = gaugefield_U3(N_x, N_t, hot)

for i = 1:N_therm
    acc_copy = acc_metro_therm[1]
    chess_metro_U3!(U, ϵ, β, acc_metro_therm)
    push!(actions,action_U3(U, β))
    push!(charges, top_charge(U))
    if mod(i,Int(N_therm/20)) == 0
        println("Progress Thermalization: ", count_therm[1], "%")
        count_therm[1] += 5
    end
    ϵ *= sqrt((acc_metro_therm[1]-acc_copy)  /2/N_t/N_x/N_metro / acc_wish) # only update ϵ acc. to Metropolis
    # global epsilon = deepcopy(ϵ)
    # println(epsilon)
end

for i = 1:N_up
    for j = 1:N_metro
        chess_metro_U3!(U, ϵ, β, acc_metro)
    end
    for i = 1:N_insta
        # insta_update_U3_tryout!(U, β, acc_insta)
        insta_update_U3!(U, β, acc_insta, ΔQ)
    end
    if mod(i,Int(N_up/20)) == 0
        println("Progress Simulation: ", count[1], "%")
        count[1] += 5
    end
    push!(actions, action_U3(U,β))
    push!(charges, top_charge(U))
end





accrate_metro = acc_metro[1]/N_x/N_t/2/N_metro/N_up
accrate_insta = acc_insta[1]/N_x/N_t/2/N_metro/N_up
accrate_metro_print = round(100*accrate_metro, digits = 1)
accrate_insta_print = round(100*accrate_insta, digits = 1)
Q_mean = mean(charges[N_therm+1:end])
Q_mean_print = round(Q_mean, digits = 2)

let
    image_actions = plot(
        actions,
        title = "Actions in 2D U(3): β = $β, N_x = N_t = $N_x \n N_metro = $N_metro, N_insta = $N_insta, ΔQ = $ΔQ",
        label = :false, #"acc. rate Metrop. after therm: $accrate_metro_print %",
        xlabel = "Monte Carlo time"
    )
    image_actions = vline!([N_therm], label = "Thermalization (only Metrop.)", color = :red)
    display(image_actions)
end
let
    image_charges = plot(
        charges,
        title = "Top. charges in 2D U(3): β = $β, N_x = N_t = $N_x \n N_metro = $N_metro, N_insta = $N_insta, ΔQ = $ΔQ",
        label = "acc. rate insta. after therm: $accrate_insta_print %",
        xlabel = "Monte Carlo time"
    )
    image_charges = vline!([N_therm], label = "Thermalization (only Metrop.)", color = :red)
    display(image_charges)
end
let
    hist_charges = histogram(
        round.(Int, charges[N_therm+1:end]),
        title = "Top. charges in 2D U(3): β = $β, N_x = N_t = $N_x \n N_metro = $N_metro, N_insta = $N_insta, ΔQ = $ΔQ",
        label = "Charges after therm. \n⟨Q⟩ = $Q_mean_print"
    )
    display(hist_charges)
end


#=
top_charge(U)
action_U3(U,β)
action_U3(insta_U3_tryout(N_x,N_t,false).*U,β)


for i = 1:100000
    chess_cool_U3!(U, ϵ/50, [0])
end
action_U3(U,β)


function write_conf_U3(U, path)
    NX = size(U,2)
    NT = size(U,3)
    bla = []
    for μ = 1:2
        for x = 1:NX
            for t = 1:NT
                push!(bla, U[μ,x,t])
            end
        end
    end
    writedlm(path, bla, ',')
    return nothing
end

function read_config_U3(path)
    blub = readdlm(path, ',', ComplexF64)
    L = Int(sqrt(size(blub,1)/2))
    N = Int(sqrt(size(blub,2)))
    return [reshape(blub[t+(x-1)*L+(μ-1)*L^2,:], N, N) for μ = 1:2, x = 1:L, t = 1:L] 
end

# write_conf_U3(U,"C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\conf3.txt")
# bla = read_config_U3("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\conf3.txt");
# # U == bla # true


top_charge(bla)
action_U3(bla, 1)
action_U3(insta_U3(N_x,N_t,3),1)
mean([plaq(bla,x,t) for x = 1:N_x, t = 1:N_t])
std([imag.(plaq(bla,x,t)) for x = 1:N_x, t = 1:N_t])

mean([real(tr(plaq(bla,x,t))) for x = 1:N_x, t = 1:N_t])
std([real(tr(plaq(bla,x,t))) for x = 1:N_x, t = 1:N_t])
heatmap([real(tr(plaq(bla,x,t))) for x = 1:N_x, t = 1:N_t])

maximum([real(tr(plaq(bla,x,t))) for x = 1:N_x, t = 1:N_t])
minimum([real(tr(plaq(bla,x,t))) for x = 1:N_x, t = 1:N_t])

# action_U3(U,β)
# insta_action_min_U3(β, N_x, N_t, 3, 3)
# insta_action_min(β, 3, N_x, N_t, 3, 3)


function insta_U3_tryout(N_x, N_t, up)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    w = 4*π/sqrt(3)
    Q = -1
    if up
        w = -4*π/sqrt(3)
        Q = 1
    end
    U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((-im*t*w)/(N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp((im*x*w)/(N_x) * λ8) for x = 1:N_x]
    return U
end

bla1 = insta_U3(N_x,N_t,1);
bla2 = insta_U3(N_x,N_t,2);
bla3 = insta_U3(N_x,N_t,3);
top_charge(bla3)
top_charge(bla1 .* insta_U3_tryout(N_x,N_t,true))
action_U3(bla2,1)
action_U3(insta_U3_tryout(N_x,N_t,true) .*bla1 , 1)

plaq(insta_U3(N_x,N_t,3),rand(1:N_x),rand(1:N_t))
plaq(bla,rand(1:N_x),rand(1:N_t))
=#