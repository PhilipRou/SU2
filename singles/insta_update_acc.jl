include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
# include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\analyze\\SU2_analyze_head.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\analyze\\SU2_jackknives.jl")

using DelimitedFiles
using StatsBase



for L  = 32:8:64
    println("Started L = $L")
    for Δq = 0:6
        println("\t Started Δq = $Δq")
        # L        = 32
        # Δq       = 1
        global β        = 5.0 * L^2 / 32^2
        global N_therm  = 500
        global N_meas   = 5e3
        global N_metro  = 4
        global N_over   = 3
        global N_insta  = 1
        global group    = "U2"
        global ϵ        = 0.1
        global acc_wish = 0.8

        global acc_metro = [0.0]
        global acc_over  = [0.0]
        global acc_insta = [0.0]

        global metro_acceptances = []
        global over_acceptances  = []
        global insta_acceptances = []

        global actions = []
        global charges = []

        global U = gaugefield_U2(L, L, true)

        for therm = 1:N_therm
            chess_metro!(U,ϵ,β,acc_metro,group)
            ϵ *= sqrt(acc_metro[1] / acc_wish) # only update ϵ acc. to Metropolis
            global epsilon = deepcopy(ϵ)
            println("\t\t Acceptance: $(round(acc_metro[1],digits = 3)), ϵ: $epsilon ")
        end

        println("\t L = $L, Δq = $Δq")

        count = 0

        for meas = 1:N_meas
            for metro = 1:N_metro
                chess_metro!(U,ϵ,β,acc_metro,group)
                push!(metro_acceptances,acc_metro[1])
            end
            for over  = 1:N_over
                chess_overrelax!(U,acc_over)
                push!(over_acceptances,acc_over[1])
            end
            for insta = 1:N_insta
                insta_update_U2!(U,β,acc_insta,Δq)
                push!(insta_acceptances,acc_insta[1])
            end
            push!(actions, action(U,β))
            push!(charges, top_charge_U2(U))
            if meas%Int(N_meas/20) == 0
                count += 5
                println("\t\t We're already $count% deep in the sim")
            end
        end

        metro_acc_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\acc_metro_L_$(L)_deltaq_$Δq.txt"
        over_acc_path  = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\acc_over_L_$(L)_deltaq_$Δq.txt"
        insta_acc_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\acc_insta_L_$(L)_deltaq_$Δq.txt"
        actions_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\actions_L_$(L)_deltaq_$Δq.txt"
        charges_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\charges_L_$(L)_deltaq_$Δq.txt"
        writedlm(metro_acc_path, metro_acceptances)
        writedlm(over_acc_path, over_acceptances)
        writedlm(insta_acc_path, insta_acceptances)
        writedlm(actions_path, actions)
        writedlm(charges_path, charges)
    end
end


Ls = Vector(32:8:64)
qs = Vector(0:6)

insta_acc_all     = Array{Float64}(undef,length(qs),length(Ls))
insta_acc_err_all = Array{Float64}(undef,length(qs),length(Ls))

susc_all     = Array{Float64}(undef,length(qs),length(Ls))
susc_err_all = Array{Float64}(undef,length(qs),length(Ls))

for L_ind = 1:length(Ls)
    for q_ind = 1:length(qs)
        L  = Ls[L_ind]
        Δq = qs[q_ind]
        metro_acc_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\acc_metro_L_$(L)_deltaq_$Δq.txt"
        over_acc_path  = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\acc_over_L_$(L)_deltaq_$Δq.txt"
        insta_acc_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\acc_insta_L_$(L)_deltaq_$Δq.txt"
        actions_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\actions_L_$(L)_deltaq_$Δq.txt"
        charges_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\acc\\charges_L_$(L)_deltaq_$Δq.txt"
        metro_acc = readdlm(metro_acc_path)
        over_acc = readdlm(over_acc_path)
        insta_acc = readdlm(insta_acc_path)
        actions = readdlm(actions_path)
        charges = readdlm(charges_path)

        b_size_insta = round(Int, 2*auto_corr_time(insta_acc)+1, RoundUp)
        insta_acc_all[q_ind,L_ind], insta_acc_err_all[q_ind,L_ind] = jackknife(insta_acc, b_size_insta)

        b_size_q = round(Int, 2*auto_corr_time(charges)+1, RoundUp)
        susc_all[q_ind,L_ind], susc_err_all[q_ind,L_ind] = jackknife(charges.^2 ./L^2, b_size_q)
    end
end