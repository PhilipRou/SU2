include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")

using Plots
using LsqFit


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



#=
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



N_x = N_t = 32
even_actions = []
odd_actions = []
insta_actions = []
even_Qs = []
odd_Qs = []
Q_bound = 6
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

model(x, p) = p[1] .* x.^2
p0 = [insta_actions[5]]
fit = curve_fit(model, Vector(-Q_bound:Q_bound), insta_actions, p0)
c = round(1000*fit.param[1], digits = 1)
c_err = round(1000* stderror(fit)[1], digits = 1, RoundUp)

image_insta_actions =  plot(
    Vector(-Q_bound:0.01:Q_bound),
    [model(x, fit.param) for x in Vector(-Q_bound:0.01:Q_bound)],
    color = :grey25,
    linestyle = :dash,
    label = "fit(x) = c⋅x², c = ($c ± $c_err)e-3",
    legend = :top,
    xlabel = "Top. charge Q of the instanton",
    # ylabel = "Action S",
    title = "Norm. Action S/β for Instanton Configurations\n 2D U(2), N_x = N_t = $N_x",
    xticks = -Q_bound:2:Q_bound
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
    label = "Odd Q"
)

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\insta_actions_U2_N_t_$N_t.pdf")

=#


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

