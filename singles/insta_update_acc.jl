include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
# include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\observables\\smearing.jl")

using DelimitedFiles



#=
for L  = 32:8:64
    println("Started L = $L")
    # for Δq = 0:6
    Δq = 2
        println("\t Started Δq = $Δq")
        # L        = 32
        # Δq       = 1
        global β        = 5.0 * L^2 / 32^2
        global N_therm  = 500
        global N_meas   = 5e4
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
            if therm%Int(N_therm/10) == 0
                println("\t\t Acceptance: $(round(acc_metro[1],digits = 3)), ϵ: $epsilon ")
            end
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

        base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\U2_data_new\\acc_3"
        metro_acc_path = string(base_path,"\\acc_metro_L_$(L)_deltaq_$Δq.txt")
        over_acc_path  = string(base_path,"\\acc_over_L_$(L)_deltaq_$Δq.txt")
        insta_acc_path = string(base_path,"\\acc_insta_L_$(L)_deltaq_$Δq.txt")
        actions_path = string(base_path,"\\actions_L_$(L)_deltaq_$Δq.txt")
        charges_path = string(base_path,"\\charges_L_$(L)_deltaq_$Δq.txt")
        writedlm(metro_acc_path, metro_acceptances)
        writedlm(over_acc_path, over_acceptances)
        writedlm(insta_acc_path, insta_acceptances)
        writedlm(actions_path, actions)
        writedlm(charges_path, charges)
    # end
end
=#





include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\analyze\\SU2_analyze_head.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\SU2\\analyze\\SU2_jackknives.jl")



Ls = Vector(32:8:64)
# qs = Vector(0:6)
qs = [2]
betas = 5.0 .* Ls.^2 ./ 32^2

insta_acc_all     = Array{Float64}(undef,length(qs),length(Ls))
insta_acc_err_all = Array{Float64}(undef,length(qs),length(Ls))

susc_all     = Array{Float64}(undef,length(qs),length(Ls))
susc_err_all = Array{Float64}(undef,length(qs),length(Ls))

even_susc_all     = Array{Float64}(undef,length(qs),length(Ls))
even_susc_err_all = Array{Float64}(undef,length(qs),length(Ls))

odd_susc_all     = Array{Float64}(undef,length(qs),length(Ls))
odd_susc_err_all = Array{Float64}(undef,length(qs),length(Ls))

swil_all     = Array{Float64}(undef,length(qs),length(Ls))
swil_err_all = Array{Float64}(undef,length(qs),length(Ls))

for L_ind = 1:length(Ls)
    L  = Ls[L_ind]
    β  = betas[L_ind]
    for q_ind = 1:length(qs)
        Δq = qs[q_ind]        
        base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\julia_projects\\U2_data_new\\acc_2"
        metro_acc_path = string(base_path,"\\acc_metro_L_$(L)_deltaq_$Δq.txt")
        over_acc_path  = string(base_path,"\\acc_over_L_$(L)_deltaq_$Δq.txt")
        insta_acc_path = string(base_path,"\\acc_insta_L_$(L)_deltaq_$Δq.txt")
        actions_path = string(base_path,"\\actions_L_$(L)_deltaq_$Δq.txt")
        charges_path = string(base_path,"\\charges_L_$(L)_deltaq_$Δq.txt")
        metro_acc = readdlm(metro_acc_path)
        over_acc = readdlm(over_acc_path)
        insta_acc = readdlm(insta_acc_path)
        charges      = readdlm(charges_path)
        susc         = charges.^2 ./L^2 .* β ./4
        even_charges = charges[(x->iseven(x)).(round.(Int,charges))]
        odd_charges  = charges[(x->isodd(x)).(round.(Int,charges))]
        even_susc = even_charges.^2 ./L^2 .* β ./4
        odd_susc  = odd_charges.^2 ./L^2 .* β ./4
        # swil = readdlm(actions_path) ./ (β*L^2)

        if isnan(auto_corr_time(insta_acc))
            insta_acc_all[q_ind,L_ind]     = insta_acc[1]
            insta_acc_err_all[q_ind,L_ind] = 0
        else
            b_size_insta = round(Int, 2*auto_corr_time(insta_acc)+1, RoundUp)
            insta_acc_all[q_ind,L_ind], insta_acc_err_all[q_ind,L_ind] = jackknife(insta_acc, b_size_insta)
        end
        
        b_size_q = round(Int, 2*auto_corr_time(charges)+1, RoundUp) ### auto time of charges on purpose cause it's longer
        @show b_size_q
        susc_all[q_ind,L_ind], susc_err_all[q_ind,L_ind] = jackknife(susc, b_size_q)
        if length(even_susc)>0
            even_susc_all[q_ind,L_ind], even_susc_err_all[q_ind,L_ind] = jackknife(even_susc, b_size_q)
        else
            even_susc_all[q_ind,L_ind], even_susc_err_all[q_ind,L_ind] = [NaN,NaN]
        end
        if length(odd_susc)>0
            odd_susc_all[q_ind,L_ind], odd_susc_err_all[q_ind,L_ind] = jackknife(odd_susc, b_size_q)
        else
            odd_susc_all[q_ind,L_ind], odd_susc_err_all[q_ind,L_ind] = [NaN,NaN]
        end
        
        # b_size_swil = round(Int, 2*auto_corr_time(swil)+1, RoundUp)
        # swil_all[q_ind,L_ind], swil_err_all[q_ind,L_ind] = jackknife(swil, b_size_swil)

        # display(plot(charges, title = "top. charge time series \n L = $L, Δq = $Δq, β = $β", label = :false))
    end
end






let
    shapes = [:circle, :utriangle, :diamond, :rtriangle, :star5, :dtriangle, :star8]
    image_acc = plot(
        xlabel = latexstring("inverse coupling \$\\beta\$"),
        ylabel = "inst. acceptance rate",
        # legend = :bottomright,
        background_color_legend = nothing,
        tickfontsize = 10,
        labelfontsize = 15,
        legendfontsize = 11,
    )
    for L_ind = 1:length(Ls)
        L = Ls[L_ind]
        β = betas[L_ind]
        if L_ind == 1
            for q_ind = 1:length(qs)
                Δq = qs[q_ind]
                image_acc = scatter!(
                    [β],
                    [insta_acc_all[q_ind,L_ind]],
                    yerror = [insta_acc_err_all[q_ind,L_ind]],
                    label  = latexstring("\$\\Delta q = $Δq\$"),
                    markersize = 6,
                    color  = cb_colors[q_ind],
                    markerstrokecolor = cb_colors[q_ind],
                    markershape = shapes[q_ind]
                )
            end
        else
            for q_ind = 1:length(qs)
                Δq = qs[q_ind]
                image_acc = scatter!(
                    [β],
                    [insta_acc_all[q_ind,L_ind]],
                    yerror = [insta_acc_err_all[q_ind,L_ind]],
                    # label  = latexstring("\$\\Delta q = $Δq\$"),
                    label = :false,
                    markersize = 6,
                    color  = cb_colors[q_ind],
                    markerstrokecolor = cb_colors[q_ind],
                    markershape = shapes[q_ind]
                )
            end
        end
    end
    display(image_acc)
end

fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\insta_updates\\acc_insta_update_2.pdf"
# savefig(fig_path)

tab_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\tabellen\\insta_updates\\acc_insta_update_2.txt"
bla = open(tab_path, "w")
write(bla, "\\begin{table}[H]\n\t\\centering\n\t\\hline\n\t\\begin{tabular}{|c||c|c|c|c|c|}\n")
write(bla, "\t\t & \t \$L=$(Ls[1])\$ & \t \$L=$(Ls[2])\$ & \t \$L=$(Ls[3])\$ & \t \$L=$(Ls[4])\$ & \t \$L=$(Ls[5])\$ \t \\\\\n")
write(bla, "\t\t \$\\Delta q\$ & \t \$L=$(betas[1])\$ & \t \$L=$(betas[2])\$ & \t \$L=$(betas[3])\$ & \t \$L=$(betas[4])\$ & \t \$L=$(betas[5])\$ \t \\\\\\hline\\hline\n")
error_digs     = 2
anal_digs_pre  = 1
anal_digs_post = 8
insta_acc_str = format_x_err.(insta_acc_all, insta_acc_err_all, error_digs)
for q_ind = 1:length(qs)
    write(bla, "\t\t $(qs[q_ind]) & \t $(insta_acc_str[q_ind,1]) & \t $(insta_acc_str[q_ind,2]) & \t $(insta_acc_str[q_ind,3]) & \t $(insta_acc_str[q_ind,4]) & \t $(insta_acc_str[q_ind,5]) & \t \\\\\\hline\n ")
end
write(bla,"\t\\end{tabular}\n\t\\caption{Caption}\n\t\\label{tab:my_label}\n\\end{table}")
close(bla)





#=
function analytic_susc_U2(β)
    nasty(α)   = besseli(1,β*cos(α))/cos(α)
    nastier(α) = α^2 * besseli(1,β*cos(α))/cos(α)
    return quadgk(nastier,-π/2,π/2)[1] / quadgk(nasty,-π/2,π/2)[1] / π^2
end

betas_anal = first(betas)-0.5:0.1:last(betas)+0.5
susc_anal = [analytic_susc_U2(β) for β in betas_anal]

let
    # greys = reverse([:grey84, :grey72, :grey60, :grey48,:grey36])
    greys  = reverse([:grey70, :grey60, :grey50, :grey40,:grey30])
    shapes = [:circle, :utriangle, :diamond, :rtriangle, :star5, :dtriangle, :star8]
    image_susc = plot(
        xlabel = latexstring("inverse coupling \$\\beta\$"),
        ylabel = latexstring("\$ \\chi_\\mathrm{top} a^2 \$"),
        legend = :outerright,
        # foreground_color_legend = nothing,
        # background_color_legend = nothing,
        tickfontsize = 10,
        labelfontsize = 15,
        legendfontsize = 12,
        size = (750,400),
        leftmargin = 4mm,
        bottommargin = 2mm
    )
   image_susc = plot!(
        betas_anal,
        susc_anal,
        label = :false,# "analytic",
        color = :black
   )
    
    for L_ind = 1:length(Ls)
        L = Ls[L_ind]
        β = betas[L_ind]
        if L_ind == 1
            for q_ind = 1:length(qs)
                Δq = qs[q_ind]
                image_susc = scatter!(
                    [β],
                    [susc_all[q_ind,L_ind]],
                    yerror = [susc_err_all[q_ind,L_ind]],
                    label  = latexstring("\$\\Delta q = $Δq\$"),
                    markersize = 6,
                    color  = cb_colors[q_ind],
                    markerstrokecolor = cb_colors[q_ind],
                    markershape = shapes[q_ind]
                )
            end
        else
            for q_ind = 1:length(qs)
                Δq = qs[q_ind]
                image_susc = scatter!(
                    [β],
                    [susc_all[q_ind,L_ind]],
                    yerror = [susc_err_all[q_ind,L_ind]],
                    # label  = latexstring("\$\\Delta q = $Δq\$"),
                    label  = :false,
                    markersize = 6,
                    color  = cb_colors[q_ind],
                    markerstrokecolor = cb_colors[q_ind],
                    markershape = shapes[q_ind]
                )
            end
        end
    end
    image_susc = plot!(
         [betas_anal[1]],
         [susc_anal[1]],
         label = "analytic",
         color = :black
    )
    display(image_susc)
end

fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\insta_updates\\susc_insta_update.pdf"
# savefig(fig_path)
=#

tab_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\tabellen\\insta_updates\\susc_insta_update_2.txt"
bla = open(tab_path, "w")
write(bla, "\\begin{table}[H]\n\t\\centering\n\t\\hline\n\t\\begin{tabular}{|c||c|c|c|c|c|}\n")
write(bla, "\t\t & \t \$L=$(Ls[1])\$ & \t \$L=$(Ls[2])\$ & \t \$L=$(Ls[3])\$ & \t \$L=$(Ls[4])\$ & \t \$L=$(Ls[5])\$ \t \\\\\n")
write(bla, "\t\t \$\\Delta q\$ & \t \$L=$(betas[1])\$ & \t \$L=$(betas[2])\$ & \t \$L=$(betas[3])\$ & \t \$L=$(betas[4])\$ & \t \$L=$(betas[5])\$ \t \\\\\\hline\\hline\n")
error_digs     = 2
anal_digs_pre  = 1
anal_digs_post = 8
susc_str      = format_x_err.(susc_all, susc_err_all, error_digs)
even_susc_str = format_x_err.(even_susc_all, even_susc_err_all, error_digs)
odd_susc_str  = format_x_err.(odd_susc_all, odd_susc_err_all, error_digs)
for q_ind = 1:length(qs)-1
    write(bla, "\t\t $(qs[q_ind]) & \t $(susc_str[q_ind,1]) & \t $(susc_str[q_ind,2]) & \t $(susc_str[q_ind,3]) & \t $(susc_str[q_ind,4]) & \t $(susc_str[q_ind,5]) & \t \\\\\\hline\n ")
end
write(bla, "\t\t $(qs[end]) & \t $(susc_str[end,1]) & \t $(susc_str[end,2]) & \t $(susc_str[end,3]) & \t $(susc_str[end,4]) & \t $(susc_str[end,5]) & \t \\\\\\hline\\hline\n ")
write(bla, "\t\t $(qs[end]), even & \t $(even_susc_str[end,1]) & \t $(even_susc_str[end,2]) & \t $(even_susc_str[end,3]) & \t $(even_susc_str[end,4]) & \t $(even_susc_str[end,5]) & \t \\\\\\hline\\hline\n ")
write(bla, "\t\t $(qs[end]), odd & \t $(odd_susc_str[end,1]) & \t $(odd_susc_str[end,2]) & \t $(odd_susc_str[end,3]) & \t $(odd_susc_str[end,4]) & \t $(odd_susc_str[end,5]) & \t \\\\\\hline\\hline\n ")
susc_anal_str = [@sprintf("%*.*f", anal_digs_pre, anal_digs_post, analytic_susc_U2(β)) for β in betas]
write(bla, "\t\t anal. & \t $(susc_anal_str[1]) & \t $(susc_anal_str[2]) & \t $(susc_anal_str[3]) & \t $(susc_anal_str[4]) & \t $(susc_anal_str[5]) & \t \\\\\\hline\n")
write(bla,"\t\\end{tabular}\n\t\\caption{Caption}\n\t\\label{tab:my_label}\n\\end{table}")
close(bla)



function analytic_susc_U2(β)
    nasty(α)   = besseli(1,β*cos(α))/cos(α)
    nastier(α) = α^2 * besseli(1,β*cos(α))/cos(α)
    return quadgk(nastier,-π/2,π/2)[1] / quadgk(nasty,-π/2,π/2)[1] / π^2
end


let
    betas_anal = 1:0.1:700
    susc_anal = [analytic_susc_U2(β)*β/4 for β in betas_anal]

    image_susc = plot(
        xlabel = latexstring("\$1/\\beta = (ag)^2/4\$"),
        ylabel = latexstring("\$ \\chi_\\mathrm{top}/g^2 \$"),
        tickfontsize = 10,
        labelfontsize = 15,
        legendfontsize = 11,
        leftmargin = 2mm
    )
    image_susc = plot!(
        1 ./ betas_anal,
        susc_anal,
        label = :false,
        color = cb_orange,
        linewidth = 1.5
    )
    image_susc = scatter!(
        1 ./ betas,
        susc_all[3,:],
        yerror = susc_err_all[3,:],
        label = "Simulation data",
        markershape = :diamond,
        markersize = 5,
        color = cb_blue
    )
    image_susc = plot!(
        [1/betas_anal[1]],
        [susc_anal[1]],
        label = "Analytic",
        color = cb_orange
    )
end

let
    versatz = 0.002
    betas_anal = 3:0.1:700
    susc_anal = [analytic_susc_U2(β)*β/4 for β in betas_anal]

    image_susc = plot(
        xlabel = latexstring("\$1/\\beta = (ag)^2/4\$"),
        ylabel = latexstring("\$ \\chi_\\mathrm{top}/g^2 \$"),
        tickfontsize = 10,
        labelfontsize = 15,
        legendfontsize = 11,
        leftmargin = 2mm
    )
    image_susc = plot!(
        1 ./ betas_anal,
        susc_anal,
        label = :false,
        color = :black,
        linewidth = 1.5
    )
    # image_susc = scatter!(
    #     1 ./ betas,
    #     susc_all[:],
    #     yerror = susc_err_all[:],
    #     label = latexstring("all \$q\$"),
    #     markershape = :diamond,
    #     markersize = 5,
    #     color = cb_blue,
    #     # markerstrokecolor = cb_blue
    # )
    image_susc = scatter!(
        1 ./ betas .+ versatz,
        even_susc_all[:],
        yerror = even_susc_err_all[:],
        label = latexstring("even \$q\$ only"),
        markershape = :rtriangle,
        markersize = 8,
        color = cb_orange,
        # markerstrokecolor = cb_orange
    )
    image_susc = scatter!(
        1 ./ betas .- versatz,
        odd_susc_all[:],
        yerror = odd_susc_err_all[:],
        label = latexstring("odd \$q\$ only"),
        markershape = :ltriangle,
        markersize = 8,
        color = cb_green,
        # markerstrokecolor = cb_green
    )
    image_susc = scatter!(
        1 ./ betas,
        susc_all[:],
        yerror = susc_err_all[:],
        label = latexstring("all \$q\$"),
        markershape = :diamond,
        markersize = 6,
        color = cb_blue,
        # markerstrokecolor = cb_blue
    )
    image_susc = plot!(
        [1/betas_anal[1]],
        [susc_anal[1]],
        label = "Analytic",
        color = :black
    )
end

fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\insta_updates\\susc_insta_update_2_eo.pdf"
# savefig(fig_path)