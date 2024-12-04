include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")

using StatsBase

function action_2x2(U,β)
    NX = size(U,2)
    NT = size(U,3)
    S = 2*NX*NT
    for t = 1:NT
        for x = 1:NX
            S -= real(tr(RT_loop(U,2,2,x,t)))
        end
    end
    return β*S/32
end



L         = 32
N_x = N_t = L
β         = 8.0
N_therm   = 500
N_meas    = 5e3
N_metro   = 4
N_over    = 3
group     = "SU2"
ϵ         = 0.1
acc_wish  = 0.8

acc_metro = [0.0]
acc_over  = [0.0]
acc_insta = [0.0]

actions     = []
actions_2x2 = []
sm_actions     = []
sm_actions_2x2 = []

U = gaugefield_SU2(L, L, true)

for therm = 1:N_therm
    chess_metro!(U,ϵ,β,acc_metro,group)
    ϵ *= sqrt(acc_metro[1] / acc_wish) # only update ϵ acc. to Metropolis
    global epsilon = deepcopy(ϵ)
    if therm%Int(N_therm/10) == 0
        println("\t\t Acceptance: $(round(acc_metro[1],digits = 3)), ϵ: $epsilon ")
    end
end

count = 0

for meas = 1:N_meas
    for metro = 1:N_metro
        chess_metro!(U,ϵ,β,acc_metro,group)
    end
    for over  = 1:N_over
        chess_overrelax!(U,acc_over)
    end
    V = stout(U,15,0.2)
    push!(actions, action(U,β))
    push!(actions_2x2, action_2x2(U,β))
    push!(sm_actions, action(V,β))
    push!(sm_actions_2x2, action_2x2(V,β))
    # push!(actions, mean([real(tr(plaq(U,x,t))) for x = 1:N_x, t = 1:N_t]))
    # push!(actions_2x2, mean([real(tr(RT_loop(U,2,2,x,t))) for x = 1:N_x, t = 1:N_t]))
    # push!(actions, mean([real(tr(plaq(V,x,t))) for x = 1:N_x, t = 1:N_t]))
    # push!(actions_2x2, mean([real(tr(RT_loop(V,2,2,x,t))) for x = 1:N_x, t = 1:N_t]))
    if meas%Int(N_meas/20) == 0
        count += 5
        println("\t\t We're already $count% deep in the sim")
    end
end

# auto_corr_time(actions_2x2)
# auto_corr_time(actions) 
# auto_corr_time(sm_actions_2x2)
# auto_corr_time(sm_actions) 

jackknife(sm_actions_2x2,5)
jackknife(sm_actions,5)



# mean([real(tr(plaq(U,x,t))) for x = 1:N_x, t = 1:N_t])
# mean([real(tr(RT_loop(U,2,2,x,t))) for x = 1:N_x, t = 1:N_t])