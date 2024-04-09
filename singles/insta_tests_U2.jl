include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")

using Plots
using LsqFit
using LaTeXStrings
using LinearAlgebra


# Create a minimum of the gauge action in the topological
# sector of charge Q. Unfortunately the link values of the 
# odd-Q-instantons are not elements of the center of U(2)
# (i.e. ∉ U(1)), so that 
function insta_U2(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    if iseven(Q)
        U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * coeffs_Id_U2() for x = 1:N_x]
    else
        U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*π/N_x/N_t)*coeffs_Id_U2() - sin(t*π/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*π/N_x)*coeffs_Id_U2() + sin(x*π/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    end
    return U
end

# Create a naive (N_x × N_t)-Q-instanton configuration, in
# the sense that it is not an actual minimum of the gauge
# action for odd Q (see insta_U2() below)
function insta_U2_naive(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * coeffs_Id_U2() for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * coeffs_Id_U2() for x = 1:N_x]
    return U
end

function def_not_insta_U2(N_x, N_t, Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*π/N_x/N_t)*coeffs_Id_U2() - sin(t*π/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*π/N_x)*coeffs_Id_U2() + sin(x*π/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    return U
end

function insta_U2_split(N_x, N_t, Q, v)
    # M = grp2coeffs_U2(exp(sum(im .* v .* [σ1,σ2,σ3])))
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * grp2coeffs_U2(exp(sum(-(im*t*π/N_x/N_t) .* v .* [σ1,σ2,σ3]))) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * grp2coeffs_U2(exp(sum((im*x*π/N_x) .* v .* [σ1,σ2,σ3]))) for x = 1:N_x]
    return U
end

# insta_update() doesn't really work, but cooling the proposal
# only a couple of times (N_insta_cool) seems to do the trick
function insta_cool_update_U2!(U, step, β, N_insta_cool, acc)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U2(NX,NT,rand([-2,2])) .* U
    for i = 1:N_insta_cool
        chess_cool!(U_prop,step,β,[0],"U2")
    end
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end

# insta_update() doesn't really work, but cooling the proposal
# only a couple of times (N_insta_cool) seems to do the trick
# (which does not mean that the results obtained are correct
# in any way).
function insta_ran_cool_update_U2!(U, step, β, N_insta_cool, cool_degree_in_percent, acc)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U2_naive(NX,NT,rand([-2,2])) .* U
    for i = 1:N_insta_cool
        ran_cool!(U_prop,cool_degree_in_percent,step,β,[0],"U2")
    end
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT*cool_degree_in_percent/100
    end
    return nothing
end

function insta_flow_update_U2!(U, N_stout, acc, Q)
    NX = size(U,2)
    NT = size(U,3)
    ρ  = 1/N_stout
    U_prop = stout(U,ρ)
    for i = 1:N_stout-1
        U_prop = stout(U_prop,ρ)
    end
    U_prop = insta_U2_naive(NX,NT,rand([-Q,Q])) .* U_prop
    for i = 1:N_stout
        U_prop = stout(U_prop,-ρ)
    end
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end



# β = 1.0
# N_t = N_x = 32

# config = insta_U2(N_t,N_x,1);
# S = action(config,β)
# println("Action of the 1-instanton: ", S)
# for i = 1:10000 chess_cool!(config,0.01,β,[0],"U2") end
# println("Difference after 10'000 cooling steps: ", -S + action(config,β))
# μ = rand(1:2)
# t = rand(1:N_t)
# x = rand(1:N_x)
# smol = gaugefield_U2(N_t,N_x,false);
# smol[μ,t,x] = ran_U2(0.01)
# config = smol.*config
# println("Difference after small perturbation: ", -S + action(config,β))
# for i = 1:10000 chess_cool!(config,0.01,β,[0],"U2") end
# println("Difference after another 10'000 cooling steps: ", -S + action(config,β))



β = 12.0
ϵ = 0.03
N_x = 32
N_t = N_x
Q_bound = 4

# bla = gaugefield(N_x, N_t,true,"U2","square");
bla = insta_U2(N_x, N_t, 0);
for i = 1:300 chess_metro!(bla, ϵ, β, [0], "U2") end
# bla = insta_U2(N_x, N_t, -round(Int, top_charge_U2(bla))) .* bla
# for i = 1:1000 chess_metro!(bla, ϵ, β, [0], "U2") end
# action(bla,12)
# top_charge_U2(bla)

# bla = insta_U2(N_x, N_t, 0);

actions = []
actions_naive = []
Qs = Vector(-Q_bound:0.01:Q_bound)
charges = []
charges_naive = []
for q in Qs
    config = def_not_insta_U2(N_x, N_t, q) .* bla
    config_naive = insta_U2_naive(N_x, N_t, q) .* bla
    push!(actions, action(config, 1.0 ))
    push!(actions_naive, action(config_naive, 1.0 ))
    push!(charges, top_charge_U2(config) )
    push!(charges_naive, top_charge_U2(config_naive) )
end

actions = actions .- minimum(actions_naive)
actions_naive = actions_naive .- minimum(actions_naive)

image_actions = plot(
    Qs,
    actions_naive,
    label = "abelian instantons",
    legend = :right,
    xlabel = "Argument for top. charge Q of the instanton",
    title = "S/β after Mult. with \"even\" and \"odd\" Instantons\n 2D U(2), N_x = N_t = $N_x "
)
image_actions = plot!(
    Qs,
    actions,
    label = "non-ab. instantons"
)

display(image_actions)

# display(histogram2d(round.(Int,charges), round.(Int,charges_naive), normalize = :true))

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\insta_update_actions_2_U2_N_t_$N_t.pdf")


#=
N_x = N_t = 32
even_actions = []
odd_actions = []
insta_actions = []
even_Qs = []
odd_Qs = []
Q_bound = 10
for q in -Q_bound:Q_bound
    push!(insta_actions, action(insta_U2(N_x, N_t, q), 1))
    if iseven(q)
        push!(even_actions, action(insta_U2(N_x, N_t, q), 1))
        push!(even_Qs,q)
    elseif isodd(q)
        push!(odd_actions, action(insta_U2(N_x, N_t, q), 1))
        push!(odd_Qs,q)
    end
end

# model(x, p) = p[1] .* x.^2 .+ p[2]
model(x, p) = p[1] .* (1 .- cos.(p[2] .*x)) .+ p[3]
p0 = [insta_actions[6], π/(N_x*N_t), 0.1]
# fit = curve_fit(model, Vector(-Q_bound:Q_bound), insta_actions, p0)
# a = round(fit.param[1], sigdigits = 2)
# c = round(fit.param[2], sigdigits = 2)
# a_err = round(stderror(fit)[1], sigdigits = 2, RoundUp)
# c_err = round(stderror(fit)[2], sigdigits = 2, RoundUp)

fit_even = curve_fit(model, even_Qs, even_actions, p0)
a_even = round(fit_even.param[1], sigdigits = 2)
c_even = round(fit_even.param[2], sigdigits = 2)
a_err_even = round(stderror(fit_even)[1], sigdigits = 2, RoundUp)
c_err_even = round(stderror(fit_even)[2], sigdigits = 2, RoundUp)

fit_odd = curve_fit(model, odd_Qs, odd_actions, p0)
a_odd = round(fit_odd.param[1], sigdigits = 2)
c_odd = round(fit_odd.param[2], sigdigits = 2)
a_err_odd = round(stderror(fit_odd)[1], sigdigits = 2, RoundUp)
c_err_odd = round(stderror(fit_odd)[2], sigdigits = 2, RoundUp)

# image_insta_actions =  plot(
#     Vector(-Q_bound:0.01:Q_bound),
#     [model(x, fit.param) for x in Vector(-Q_bound:0.01:Q_bound)],
#     color = :grey25,
#     linestyle = :dashdot,
#     # label = "Overall fit(x) = c + a⋅x²,\nc = $c ± $c_err,   a = $a ± $a_err",
#     label = L"Overall fit $f(x) = a\cdot \cos(bQ) + c$ ",
#     legend = :top,
#     xlabel = "Top. charge Q of the instanton",
#     # ylabel = "Action S",
#     title = "S/β for Instanton Configurations\n 2D U(2), N_x = N_t = $N_x",
#     xticks = -Q_bound:2:Q_bound
# )
image_insta_actions =  plot(
    Vector(-Q_bound:0.01:Q_bound),
    [model(x, fit_even.param) for x in Vector(-Q_bound:0.01:Q_bound)],
    color = palette(:default)[1],
    # label = "Even fit \nc = $c_even ± $c_err_even,   a = $a_even ± $a_err_even",
    label = L"Even fit $f(x) = a\cdot \cos(bQ) + c$ ",
    legend = :top,
    xlabel = "Top. charge Q of the instanton",
    title = "S/β for Instanton Configurations\n 2D U(2), N_x = N_t = $N_x",
    xticks = -Q_bound:2:Q_bound
)
image_insta_actions =  plot!(
    Vector(-Q_bound:0.01:Q_bound),
    [model(x, fit_odd.param) for x in Vector(-Q_bound:0.01:Q_bound)],
    linestyle = :dash,
    color = palette(:default)[2],
    # label = "Odd fit \nc = $c_odd ± $c_err_odd,   a = $a_odd ± $a_err_odd"
    label = L"Odd fit $f(x) = a\cdot \cos(bQ) + c$ ",
)
image_insta_actions =  scatter!(
    even_Qs,
    even_actions,
    color = palette(:default)[1],
    label = "Even Q"
)
image_insta_actions =  scatter!(
    odd_Qs,
    odd_actions,
    color = palette(:default)[2],
    markershape = :diamond,
    label = "Odd Q"
)
# image_insta_actions = plot!(
#     Vector(-Q_bound:0.01:Q_bound),
#     [1024*(1-cos(q*pi/1024)) for q = -Q_bound:0.01:Q_bound]
# )
display(image_insta_actions)
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\insta_actions_U2_N_t_$N_t.Q_bound_$Q_bound.pdf")



N_x = N_t = 32
Q_bound = 4

actions_naive = []
actions_non_abel = []
Qs = Vector(-Q_bound:0.01:Q_bound)
for q in Qs
    config_naive = insta_U2_naive(N_x, N_t, q) 
    config_non_abel = def_not_insta_U2(N_x, N_t, q)
    push!(actions_naive, action(config_naive, 1.0 )/(N_x*N_t))
    push!(actions_non_abel, action(config_non_abel, 1.0 )/(N_x*N_t))
end

spacing = 0.1
Qs_1 = Vector(-Q_bound:spacing:Q_bound)
Qs_2 = Vector(-Q_bound:spacing:Q_bound) .+ spacing/4
Qs_3 = Vector(-Q_bound:spacing:Q_bound) .+ (2*spacing/4)
Qs_rand = Vector(-Q_bound:spacing:Q_bound) .+ (3*spacing/4)
actions_3 = []
actions_2 = []
actions_1 = []
actions_rand = []
v = rand(3)
v = v./(sqrt(sum(v.^2)))
for q in Qs_1
    config_1 = insta_U2_split(N_x, N_t, q, [1.0, 0.0, 0.0])
    config_2 = insta_U2_split(N_x, N_t, q + spacing/4, [0.0, 1.0, 0.0])
    config_3 = insta_U2_split(N_x, N_t, q + 2*spacing/4, [0.0, 0.0, 1.0])
    config_rand = insta_U2_split(N_x, N_t, q + 3*spacing/4, v)
    push!(actions_1, action(config_1, 1.0 )/(N_x*N_t))
    push!(actions_2, action(config_2, 1.0 )/(N_x*N_t))
    push!(actions_3, action(config_3, 1.0 )/(N_x*N_t))
    push!(actions_rand, action(config_rand, 1.0 )/(N_x*N_t))
end

image_actions = plot(
    Qs,
    actions_naive,
    label = "abelian instantons",
    color = :dodgerblue2,
    legend = :right,
    xlabel = "Argument for top. charge Q of the instanton",
    title = "S/βV of \"even\" and \"odd\" Instantons\n 2D U(2), N_x = N_t = $N_x "
)
# image_actions = plot!(
#     Qs,
#     actions_non_abel,
#     label = :false,
#     color = :red
# )
image_actions = scatter!(
    Qs_1,
    actions_1,
    markershape = :circle,
    # color = :red,
    # markersize = 5,
    label = L"$\sigma_1$-instantons"
)
image_actions = scatter!(
    Qs_2,
    actions_2,
    markershape = :diamond,
    # color = :red,
    # markersize = 5,
    label = L"$\sigma_2$-instantons"
)
image_actions = scatter!(
    Qs_3,
    actions_3,
    markershape = :utriangle,
    # color = :red,
    # markersize = 5,
    label = L"$\sigma_3$-instantons"
)
image_actions = scatter!(
    Qs_rand,
    actions_rand,
    markershape = :star4,
    # color = :red,
    markersize = 5,
    label = "mixed instantons"
)

display(image_actions)

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\insta_actions_U2_eo.pdf")




# function mod_pi_half(X::coeffs_U2)
#     ϕ = imag(log(X.a/abs(X.a)))
#     # return exp(im*(sign(ϕ) * mod(abs(ϕ),π/2) - ϕ)) * X
#     return exp(im*(mod(ϕ+π/2,π) - π/2 - ϕ)) * X
# end

# function mod_pi_half(X::Number)
#     ϕ = imag(log(X/abs(X)))
#     # return exp(im*(sign(ϕ) * mod(abs(ϕ),π/2) - ϕ)) * X
#     return exp(im*(mod(ϕ+π/2,π) - π/2 - ϕ)) * X
# end

# function top_charge_U2_alt(U)
#     NX = size(U,2)
#     NT = size(U,3)
#     return sum([imag(log(det(mod_pi_half(plaq(U, x, t))))) for x = 1:NX, t = 1:NT]) / 2 / π
# end

=#
#=
bla = gaugefield(32,32,true,"U2","square");
for i = 1:300 chess_metro!(bla, 0.03, 12, [0], "U2") end
action(bla,12)
top_charge_U2(bla)
Q_bla = round(Int, top_charge_U2(bla))
bla = insta_U2(32,32,-Q_bla) .* bla;
action(bla,12)
for i = 1:300 chess_metro!(bla, 0.03, 12, [0], "U2") end


# insta_ran_cool_update_U2!(bla, 0.03, 12, 1, 50, acc_insta_cool)

cool_deg = 30
acc_wish_insta = 0.45
acc_insta_cool = [0.0]
cool_degs = []
for i = 1:1000
    acc_copy = deepcopy(acc_insta_cool[1])
    # for j = 1:20
        for j = 1:3
            chess_metro!(bla, 0.03, 12, [0],"U2")
        end
        insta_ran_cool_update_U2!(bla, 0.03, 12, 1, cool_deg, acc_insta_cool)
    # end
    # acc_rate_insta = acc_insta_cool[1]/5/32/32/2/cool_deg*100
    accepted = (acc_insta_cool[1] - acc_copy)/32/32/2/cool_deg*100
    if accepted < acc_wish_insta
        cool_deg += 0.3
    else
        cool_deg -= 0.3
    end
    # cool_deg *= sqrt(acc_wish_insta/accepted)
    # cool_deg += 0.5 * (acc_wish_insta - accepted)
    # println(acc_insta_cool)
    # if mod(i,10) == 1
    #     println(accepted)
    #     println(cool_deg)
    # end
    # println(cool_deg)
    push!(cool_degs,cool_deg)
    # ϵ *= sqrt((acc_insta[1]-acc_copy)  /2/N_t/N_x/N_metro / acc_wish) # only update ϵ acc. to Metropolis
end

plot(cool_degs)
sum(cool_degs[200:1000])/802


# bla_prop = insta_U2(32,32,1) .* bla;
# top_charge_U2(bla_prop)
# action(bla_prop,12)
# acc = [0]
# # chess_metro!(bla_prop,0.03,12,acc,"U2")
# ran_cool!(bla_prop,50,0.03,12,acc,"U2")

# acc[1]/32/32/2/50*100
action(bla,12) - action(bla_prop,12)

=#


λ1 = Complex.([0 1 0; 1 0 0; 0 0 0])
λ2 = Complex.([0 -im 0; im 0 0; 0 0 0])
λ3 = Complex.([1 0 0; 0 -1 0; 0 0 0])
λ4 = Complex.([0 0 1; 0 0 0; 1 0 0])
λ5 = Complex.([0 0 -im; 0 0 0; im 0 0])
λ6 = Complex.([0 0 0; 0 0 1; 0 1 0])
λ7 = Complex.([0 0 0; 0 0 -im; 0 im 0])
λ8 = Complex.([1 0 0; 0 1 0; 0 0 -2]) / sqrt(3)

Λ = [λ1, λ2, λ3, λ4, λ5, λ6, λ7, λ8]

λ0 = Complex.([1 0 0; 0 1 0; 0 0 1])
λ_zero = Complex.([0 0 0; 0 0 0; 0 0 0])

N_x = N_t = 32

function insta_U3(N_x, N_t, Q)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    if Q%3 == 0
        U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * λ0 for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * λ0 for x = 1:N_x]
    elseif (Q-1)%3 == 0
        U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp(-(im*t*2*π)/(3*N_x*N_t) * λ3) for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp(im*x*2*π/(3*N_x) * λ3) for x = 1:N_x]
    elseif (Q-2)%3 == 0
        U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp(-(im*t*2*π)/(3*N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
        U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
        U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp(im*x*2*π/(3*N_x) * λ8) for x = 1:N_x]
    end
    return U
end

function ran_U3(N_x, N_t, ϵ)
    return [exp(ϵ * im * rand()) * exp(ϵ * im * sum(rand(8).*Λ)) for μ = 1:2, x = 1:N_x, t = 1:N_t]
end
function ran_U3(ϵ)
    return exp(ϵ * im * sum(rand(8).*Λ))
end



# function top_charge_U3(U)
#     NX = size(U,2)
#     NT = size(U,3)
#     return sum([imag(log(det(plaq(U, x, t)))) for x = 1:NX, t = 1:NT]) / 2 / π
# end

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

# Q_bound = 6
# scatter(
#     -Q_bound:Q_bound, 
#     [action_U3(insta_U3(N_x, N_t, q), 1) for q = -Q_bound:Q_bound],
#     xticks = -Q_bound:Q_bound
# )

# top_charge_U2(insta_U3(N_x, N_t, 0))



# action_U3(insta_U3(N_x, N_t, 1),1)
# action_U3(ran_U3(N_x, N_t, 0.01) .* insta_U3(N_x, N_t, 1), 1)

# bla = insta_U3(N_x, N_t, 2);
# action_U3(bla,1)
# bla[rand(1:2), rand(1:N_x), rand(1:N_t)] *= ran_U3(0.01)
# action_U3(bla, 1)

function cool!(U, μ, x, t, step, β, acc, group)
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    elseif group == "U3"
        new_coeffs = ran_U3(step) * new_coeffs
    end
    staple_d = staple_dag(U,μ,x,t)
    S_old = real(tr(U[μ,x,t] * staple_d))
    S_new = real(tr(new_coeffs * staple_d))
    if S_old < S_new
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

bla1 = insta_U3(N_x, N_t, 1)
bla2 = insta_U3(N_x, N_t, 2)
acc1 = [0]
acc2 = [0]

N_cool = 1000
count = [0]
for i = 1:N_cool
    chess_cool!(bla1, 0.01, 1.0, acc1, "U3")
    chess_cool!(bla2, 0.01, 1.0, acc2, "U3")
    if Int(i%(N_cool/100)) == 1
        println("Cooling progress: ", count[1], "%" )
        count[1] += 1
    end
end

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

heatmap([real(tr(plaq(insta_U2(N_x, N_t, 2), x, t)))/2 for t = 1:N_t, x = 1:N_x])
heatmap([real(tr(plaq(insta_U3(N_x, N_t, 0), x, t)))/3 for t = 1:N_t, x = 1:N_x])
heatmap([real(tr(plaq(bla2, x, t)))/3 for t = 1:N_t, x = 1:N_x])

