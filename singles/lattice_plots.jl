include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\analyze\\SU2_analyze_head.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\analyze\\SU2_jackknives.jl")


using Plots
# using LsqFit
using LaTeXStrings
using LinearAlgebra
# using Roots
using DelimitedFiles
using Statistics
# using Optim
using QuadGK
using Measures

data_path ="C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Lattice_projects\\Lattice2024\\data" 
# fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Lattice_projects\\Lattice2024\\plots"
fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Lattice_projects\\Lattice2024\\PoS"

function insta_action_U2_min(β, N_x, N_t, q)
    return β*N_x*N_t*(1-cos(q*π/N_x/N_t))
end

function insta_action_U2(β, N_x, N_t, q, z)
    N_c = 2
    ReTrP = (N_c-1)*cos(2*π*z/N_x/N_t) + cos(2*π*(q - (N_c-1)*z)/N_x/N_t ) 
    return β*N_x*N_t*(1 - ReTrP/N_c)
end

function insta_action_Nc(β, N_c, N_x, N_t, Q, z)
    ReTrP = (N_c-1)*cos(2*π*z/N_x/N_t) + cos(2*π*(Q - (N_c-1)*z)/N_x/N_t ) 
    return β*N_x*N_t*(1 - ReTrP/N_c)
end

function analytic_susc_U1(β)
    kern(ϕ) = ϕ^2 * exp(β*cos(ϕ))
    return quadgk(kern, -π, π)[1] / besseli(0,β) / (2*π)^3
end

function analytic_susc_U2(β)
    nasty(α)   = besseli(1,β*cos(α))/cos(α)
    nastier(α) = α^2 * besseli(1,β*cos(α))/cos(α)
    return quadgk(nastier,0,π/2)[1] / quadgk(nasty,0,π/2)[1] / π^2
end

function analytic_plaq_U2(β)
    numer(α) = besseli(0,β*cos(α)) + besseli(2,β*cos(α))
    denom(α) = 2*besseli(1,β*cos(α))/cos(α)
    return quadgk(numer,0,π)[1]/quadgk(denom,0,π)[1] - 1/β
end

function insta_U2_z(N_x, N_t, q, z)
    # w = π*(2*z-q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*q*t*π/N_x/N_t) * exp_u2(-t*2*π/N_x/N_t * (z-q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*q*x*π/N_x) * exp_u2(x*2*π/N_x * (z-q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x]
    return U
end

cb_colors = parse.(Colorant, ["#377eb8", "#ff7f00", "#4daf4a", "#f781bf", "#984ea3", "#999999", "#e41a1c"]);
cb_blue, cb_orange, cb_green, cb_pink, cb_grey, cb_purple, cb_red = cb_colors;





######## Top. freezing ########





# let
    N_meas = 1000
    L_fine =   64
    L_medium = 48
    L_coarse = 32
    U_fine =   gaugefield_U2(L_fine, L_fine, true); β_fine   = 20.0;
    U_medium = gaugefield_U2(L_medium, L_medium, true); β_medium = β_fine*L_medium^2/L_fine^2;
    U_coarse = gaugefield_U2(L_coarse, L_coarse, true); β_coarse = β_fine*L_coarse^2/L_fine^2;

    actions_therm_fine =   []
    actions_therm_medium = []
    actions_therm_coarse = []
    acc_therm_fine =   [0.0]
    acc_therm_medium = [0.0]
    acc_therm_coarse = [0.0]
    ϵ_fine =   0.1
    ϵ_medium = 0.1
    ϵ_coarse = 0.1
    acc_wish =  0.8
    for therm = 1:300
        chess_metro!(U_fine, ϵ_fine, β_fine, acc_therm_fine, "U2")
        chess_metro!(U_medium, ϵ_medium, β_medium, acc_therm_medium, "U2")
        chess_metro!(U_coarse, ϵ_coarse, β_coarse, acc_therm_coarse, "U2")
        ϵ_fine *=   sqrt(acc_therm_fine[1] / acc_wish)
        ϵ_medium *= sqrt(acc_therm_medium[1] / acc_wish)
        ϵ_coarse *= sqrt(acc_therm_coarse[1] / acc_wish)
        push!(actions_therm_fine,   action(U_fine, β_fine))
        push!(actions_therm_medium, action(U_medium, β_medium))
        push!(actions_therm_coarse, action(U_coarse, β_coarse))
    end
    U_fine = insta_U2(L_fine, L_fine, -round(Int,top_charge_U2(U_fine))) .* U_fine
    U_medium = insta_U2(L_medium, L_medium, -round(Int,top_charge_U2(U_medium))) .* U_medium
    U_coarse = insta_U2(L_coarse, L_coarse, -round(Int,top_charge_U2(U_coarse))) .* U_coarse
    for therm = 1:100
        chess_metro!(U_fine, ϵ_fine, β_fine, acc_therm_fine, "U2")
        chess_metro!(U_medium, ϵ_medium, β_medium, acc_therm_medium, "U2")
        chess_metro!(U_coarse, ϵ_coarse, β_coarse, acc_therm_coarse, "U2")
        ϵ_fine *=   sqrt(acc_therm_fine[1] / acc_wish)
        ϵ_medium *= sqrt(acc_therm_medium[1] / acc_wish)
        ϵ_coarse *= sqrt(acc_therm_coarse[1] / acc_wish)
        push!(actions_therm_fine,   action(U_fine, β_fine))
        push!(actions_therm_medium, action(U_medium, β_medium))
        push!(actions_therm_coarse, action(U_coarse, β_coarse))
    end
    image_therm = plot(actions_therm_fine/L_fine^2, label = "fine")
    image_therm = plot!(actions_therm_medium/L_medium^2, label = "medium")
    image_therm = plot!(actions_therm_coarse/L_coarse^2, label = "coarse")
    display(image_therm)
    # println("ϵ_fine = $(round(ϵ_fine, digits = 3)), ϵ_medium = $(round(ϵ_medium, digits = 3)), ϵ_coarse = $(round(ϵ_coarse, digits = 3)), ")


    charges_coarse = []
    charges_medium = []
    charges_fine =   []
    for i = 1:N_meas
        chess_metro!(U_fine, ϵ_fine, β_fine, acc_therm_fine, "U2")
        chess_metro!(U_medium, ϵ_medium, β_medium, acc_therm_medium, "U2")
        chess_metro!(U_coarse, ϵ_coarse, β_coarse, acc_therm_coarse, "U2")
        push!(charges_fine,   top_charge_U2(U_fine))
        push!(charges_medium, top_charge_U2(U_medium))
        push!(charges_coarse, top_charge_U2(U_coarse))
    end
    charges_fine_path = string(data_path, "\\charges_fine.txt")
    charges_medium_path = string(data_path, "\\charges_medium.txt")
    charges_coarse_path = string(data_path, "\\charges_coarse.txt")
    # # writedlm(charges_fine_path, charges_fine)
    # # writedlm(charges_medium_path, charges_medium)
    # # writedlm(charges_coarse_path, charges_coarse)
    charges_fine = readdlm(charges_fine_path)
    charges_medium = readdlm(charges_medium_path)
    charges_coarse = readdlm(charges_coarse_path)

    plot_range = 1:N_meas

    image_charge_series = plot(
        charges_coarse[plot_range], 
        label = latexstring("\$L = $L_coarse,\\, \\beta = $β_coarse\$"),
        # color = cb_grey,
        # color = :salmon,
        color = :skyblue1,
        alpha = 0.8,
        linewidth = 1.5,
    )
    image_charge_series = plot!(
        charges_medium[plot_range],
        label = latexstring("\$L = $L_medium,\\, \\beta = $β_medium\$"),
        # color = cb_blue,
        # color = :firebrick1,
        # color = :dodgerblue2,
        color = :royalblue1,
        # alpha = 0.9,
        linewidth = 1.5,
    )
    image_charge_series = plot!(
        charges_fine[plot_range],
        label = latexstring("\$L = $L_fine,\\, \\beta = $β_fine\$"),
        # color = cb_orange,
        # color = :firebrick,
        # color = :navyblue,
        color = :mediumblue,
        # alpha = 0.9,
        linewidth = 1.5,
    )
    image_charge_series = plot!(
        grid = :false,
        legend = :topleft,
        background_color_legend = RGBA(1.0,1.0,1.0,0.85),
        yticks = -16:4:16,
        ylabel = latexstring("Top. Charge \$\\, q\$"),
        xlabel = "MC Time",
        size = (700,200),
        left_margin = 3mm,
        bottom_margin = 3mm,
    )
    display(image_charge_series)
    image_charge_series_path = string(fig_path, "\\charge_series_6.pdf")
    # savefig(image_charge_series_path)
# end





######## Verify formula for susc. ########





beta_16 = 1.0
N_meas_16 = 50000
N_sepa = 10
Ls = Array(16:8:48)
betas = beta_16/16^2 .* Ls.^2
N_meas = round.(Int, N_meas_16 * 16^2 .* (Ls.^(-2)) )
# N_meas = round.(Int, N_meas_16 * 16^2 .* (Ls.^(-2)) .* betas.^2 )
# queues = Array{Float64}(undef, N_meas, length(Ls))
queues  = []
actions = []

for i in eachindex(Ls)
    # i = 1
    N_meas = 5000 # N_meas[i]
    push!(queues, [])
    push!(actions, [])
    L = Ls[i]
    U = gaugefield_U2(L,L,true)
    ϵ = 0.1
    acc_therm = [0.0]
    acc_wish = 0.8
    for therm = 1:300
        chess_metro!(U, ϵ, betas[i], acc_therm, "U2")
        ϵ *=   sqrt(acc_therm[1] / acc_wish)
    end
    count = 0
    for meas = 1:N_meas
        for sepa = 1:N_sepa
            chess_metro!(U, ϵ, betas[i], [0.0], "U2")
        end
        # push!(queues[i], top_charge_U2(U))
        push!(actions[i], action(U, 1.0))
        if meas%Int(N_meas/100) == 0
            count += 1
            println("L = $(Ls[i]), Progress: $count%")
        end
        # if meas%Int(N_meas/2) == 0
        #     count += 50
        #     println("L = $(Ls[i]), Progress: $count%")
        # end
    end
end

# for i in eachindex(Ls)
#     queues_path = string(fig_path, "\\queues_L_$(Ls[i]).txt")
#     writedlm(queues_path, queues[i])
# end

# queues = []
# for i in eachindex(Ls)
#     push!(queues, readdlm(string(fig_path, "\\queues_L_$(Ls[i]).txt"))[:,1])
# end

susc_vals = []
susc_errs = []
for i in eachindex(Ls)
    b_size = round(Int, 2*auto_corr_time(queues[i])+1)
    println("Block size of q at L = $(Ls[i]): ", b_size)
    jack = jackknife(queues[i].^2 ./ Ls[i]^2 .* betas[i] ./2, b_size)
    push!(susc_vals, jack[1])
    push!(susc_errs, jack[2])
end

# beta_range = Array(minimum(betas):0.01:maximum(betas))
# beta_range = Array(0.8:0.01:9.5)
beta_range = Array(1/1.05:0.01:700) # previous upper limit: 10
susc_anal = [analytic_susc_U2(β)*β/2 for β in beta_range]
beta_range_inv = 1 ./ beta_range

let
    image_susc = plot(
        beta_range_inv,
        susc_anal,
        label = "Analytical result",
        linecolor = palette(:default)[2], # cb_grey # cb_orange
        linewidth = 1.5,
    )
    image_susc = scatter!(
        1 ./ betas,
        susc_vals,
        yerror = susc_errs,
        # label = latexstring("Num. result for \$L=32\$"),
        label = "Num. results",
        markershape = :diamond,
        markersize = 6,
        # markercolor = palette(:default)[1], # cb_red, #cb_blue
        # markerstrokecolor = :black # palette(:default)[1], # cb_red, #cb_blue
        color = palette(:default)[1],
    )
    image_susc = plot!(
        xlabel = latexstring("\$1 / \\beta\$"),
        # ylabel = latexstring("\$\\frac{\\chi_\\mathrm{top}}{g^2}\$"),
        ylabel = L"\chi_{\textbf{top}} \, / g^2",
        # xticks = (Array(1:2:9), string.(Array(1.0:2:9.0))),
        legend = :topright,
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        left_margin = 2mm,
        # xlim = (0.0,1.1),
    )
    display(image_susc)
end

image_susc_path = string(fig_path, "\\susc_over_beta_inv.pdf")
# savefig(image_susc_path)


s_wil_vals = []
s_wil_errs = []
for i in eachindex(Ls)
    b_size = round(Int, 2*auto_corr_time(actions[i])+1)
    println("Block size of s_wil at L = $(Ls[i]): ", b_size)
    jack = jackknife(actions[i] ./ Ls[i]^2 .* betas[i] ./ 2, b_size)
    push!(s_wil_vals, jack[1])
    push!(s_wil_errs, jack[2])
end

beta_range = Array(1/1.05:0.01:700) # previous upper limit: 10
s_wil_anal = [(1-analytic_plaq_U2(β))*β/2 for β in beta_range]
beta_range_inv = 1 ./ beta_range

let
    image_s_wil = plot(
        beta_range_inv,
        s_wil_anal,
        label = "Analytical result",
        linecolor = palette(:default)[2], # cb_grey # cb_orange
        linewidth = 1.5,
    )
    image_s_wil = scatter!(
        1 ./ betas,
        s_wil_vals,
        yerror = s_wil_errs,
        # label = latexstring("Num. result for \$L=32\$"),
        label = "Num. results",
        markershape = :square,
        markersize = 4,
        # markercolor = palette(:default)[1], # cb_red, #cb_blue
        # markerstrokecolor = :black # palette(:default)[1], # cb_red, #cb_blue
        color = palette(:default)[1],
    )
    image_s_wil = plot!(
        xlabel = latexstring("\$1 / \\beta\$"),
        # ylabel = latexstring("\$\\frac{\\chi_\\mathrm{top}}{g^2}\$"),
        ylabel = L"s_{\mathrm{wil}} \, / g^2",
        # xticks = (Array(1:2:9), string.(Array(1.0:2:9.0))),
        legend = :topright,
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        left_margin = 2mm,
        # xlim = (0.0,1.1),
    )
    display(image_s_wil)
end

image_s_wil_path = string(fig_path, "\\s_wil_over_beta_inv.pdf")
# savefig(image_s_wil_path)




######## Global Bogomolny ########





let
    L = 32
    q_bound = 5
    resol = 0.05
    margin = 0.1
    insta_actions = [insta_action_U2_min(1.0, L, L, q) for q = -q_bound:q_bound]
    insta_curve = [insta_action_U2_min(1.0, L, L, q) for q = -(q_bound+margin):resol:q_bound+margin]
    image_insta = plot(
        -(q_bound+margin):resol:q_bound+margin,
        insta_curve,
        # label = L"N_{x} N_{t} \left(1 - \cos\left\{\frac{q\pi}{N_x N_t} \right\}\right)",
        # label = latexstring("Eq.\$\\,\\,\$(7)"),
        label = "Lower bound",
        color = palette(:default)[2],
        linewidth = 1.5
    )
    image_insta = scatter!(
        -q_bound:q_bound,
        insta_actions,
        label = "Insta Actions",
        color = palette(:default)[1],
        markershape = :diamond,
        markersize = 6,
    )
    image_insta = plot!(
        legend = :top,
        xlabel = latexstring("Top. Charge \$q\$"),
        ylabel = L"S/\beta",
        xticks = -q_bound:q_bound,
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        # left_margin = 2mm
    )
end
image_insta_path = string(fig_path,"\\insta.pdf")
# savefig(image_insta_path)






######## Smearing Configs of q = 3 (for multiple pairs ρ*N_smear = τ) ########





# L = 32
# β = 5.0
# N_therm = 300
# N_meas  = 10
# N_sepa  = 10
# U = gaugefield_U2(L, L, true)
# ϵ = 0.1
# acc_therm = [0.0]
# acc_wish = 0.8
# for therm = 1:N_therm
#     chess_metro!(U, ϵ, β, acc_therm, "U2")
#     ϵ *= sqrt(acc_therm[1] / acc_wish)
# end

# τ = 50
# rhos   = [0.02, 0.05, 0.1]
# smears = Int.(τ ./rhos)
# meta_measures = []
# for i in eachindex(rhos)
#     push!(meta_measures, [])
# end
# for meas = 1:N_meas
#     for sepa = 1:N_sepa
#         chess_metro!(U, ϵ, β, [0.0], "U2")
#     end
#     if round(Int,abs(top_charge_U2(U))) == 1
#         for i in eachindex(rhos)
#             V = stout_midpoint(U, rhos[i])
#             sms = [action(V,β)/L^2]
#             for smear in smears[i]-1
#                 V = stout_midpoint(U, rhos[i])
#                 push!(sms, action(V,β)/L^2)
#             end
#             push!(meta_measures[i], sms)
#         end
#     end
# end
# meta_measures

let
    # @assert 1==0 "Do we really want to start a smearing run???"
    β   = 6.0
    L   = 32
    N_x = L
    N_t = L
    τ   = 10000
    rhos = [0.1, 0.05]
    N_smears = Int.(τ ./ rhos)
    N_therm   = 500
    N_meas    = 50
    N_sepa    = 100
    acc_wish  = 0.8
    ϵ         = 0.1
    base_path = string(data_path, "\\smearing") # "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_data_18"

    acc_therm = [0.0]
    U = gaugefield(N_x, N_t, true, "U2", "square")
    for therm = 1:N_therm
        chess_metro!(U,ϵ,β,acc_therm,"U2")
        ϵ *= sqrt(acc_therm[1] / acc_wish) # only update ϵ acc. to Metropolis
    end
    
    meas_count = 3

    for meas = 1:N_meas

        for sepa = 1:N_sepa
            chess_metro!(U,ϵ,β,[0.0],"U2")
        end

        ρ = rhos[1]
        N_smear = N_smears[1]
        q = top_charge_U2(U)
        smeared_actions = [action(U,β)]
        smeared_charges = [q]
        V = stout_midpoint(U,ρ)
        count = 0
        for smear = 1:1000
            q = top_charge_U2(V)
            push!(smeared_actions,action(V,β))
            push!(smeared_charges,q)
            V = stout_midpoint(V,ρ)
            if smear%Int(N_smear/100) == 0
                count += 1
                println("Measurement Nr.: $meas, rho_nr.: 1, Smearing Progress: $count%")
            end
        end

        if round(Int, abs(last(smeared_charges))) == 1
            for smear = 1:N_smear-1000
                q = top_charge_U2(V)
                push!(smeared_actions,action(V,β))
                push!(smeared_charges,q)
                V = stout_midpoint(V,ρ)
                if smear%Int(N_smear/100) == 0
                    count += 1
                    println("Measurement Nr.: $meas, rho_nr.: 1, Smearing Progress: $count%")
                end
            end
            meas_count += 1
            S_path = string(base_path,"\\sms_$meas_count.rho_nr_1.txt")
            Q_path = string(base_path,"\\smq_$meas_count.rho_nr_1.txt")
            writedlm(S_path,smeared_actions)
            writedlm(Q_path,smeared_charges)
            for i = 2:length(rhos)
                ρ = rhos[i]
                N_smear = N_smears[i]
                q = top_charge_U2(U)
                smeared_actions = [action(U,β)]
                smeared_charges = [q]
                V = stout_midpoint(U,ρ)
                count = 0
                for smear = 1:N_smear
                    q = top_charge_U2(V)
                    push!(smeared_actions,action(V,β))
                    push!(smeared_charges,q)
                    V = stout_midpoint(V,ρ)
                    if smear%Int(N_smear/100) == 0
                        count += 1
                        println("Measurement Nr.: $meas, rho_nr.: $i, Smearing Progress: $count%")
                    end
                end
                S_path = string(base_path,"\\sms_$meas_count.rho_nr_$i.txt")
                Q_path = string(base_path,"\\smq_$meas_count.rho_nr_$i.txt")
                writedlm(S_path,smeared_actions)
                writedlm(Q_path,smeared_charges)
            end
        end
        if meas_count == 6
            break
        end
    end
    println("We're done here!")
end



let
    L   = 32
    β   = 6.0
    N_x = L
    N_t = L
    τ   = 10000 # 5000 is enough bro...
    τ0  = 500
    resol = 5
    # rhos = [0.2, 0.1, 0.05]
    rhos = [0.1, 0.05]
    N_smears =   Int.(τ ./ rhos)
    start_inds = Int.(τ0./ rhos)
    N_measurements = 3

    greys = [:grey30, :grey50, :grey69]
    base_path = string(data_path, "\\smearing")
    # base_path = data_path
    # let
    image_smear = plot(
        # [insta_action_U2_min(1, N_x, N_t, 1)],
        # color = cb_red,
        # xticks = 0:resol:τ,
        xlabel = latexstring("Flow Time \$\\tau\$"),
        ylabel = L"S/\beta",
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        # left_margin = 2mm
    )
    for i in eachindex(rhos)
        x_vals = rhos[i] .* Array(0:N_smears[i])
        sms = readdlm(string(base_path,"\\sms_1.rho_nr_$i.txt")) ./β
        image_smear = plot!(
            x_vals[start_inds[i]:end],
            sms[start_inds[i]:end],
            # color = greys[i],
            color = palette(:default)[1],
            label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$"),
            alpha = 1/i,
        )
        for meas = 2:N_measurements
            sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt")) ./β
            image_smear = plot!(
                x_vals[start_inds[i]:end],
                sms[start_inds[i]:end],
                # color = greys[i],
                color = palette(:default)[meas],
                label = :false,
                alpha = 1/i
            )
        end
        # println("$(greys[i])")
    end
    image_smear = hline!(
        [insta_action_U2_min(1, N_x, N_t, 1)],
        # label = latexstring("\$S_{\\textrm{Insta}}\$ with \$|q| = 1\$"),
        label = latexstring("Eq.\$\\,\\,\$(7)"),
        color = cb_red,
        linestyle = :dash,
    )
    display(image_smear)
    # end
end

image_smear_path = string(fig_path,"\\smeared_actions.pdf")
# savefig(image_smear_path)


let
    L   = 32
    β   = 6.0
    N_x = L
    N_t = L
    τ   = 10000
    τ0  = 500
    τ1  = 2000
    rhos = [0.1, 0.05]
    N_smears =   Int.(τ ./ rhos)
    start_inds = Int.(τ0./ rhos)
    end_inds = Int.(τ1./ rhos)
    N_measurements = 3

    base_path = string(data_path, "\\smearing")

    image_smear = plot(
        # [insta_action_U2_min(1, N_x, N_t, 1)],
        # color = cb_red,
        # xticks = 0:resol:τ,
        xlabel = latexstring("Flow Time \$\\tau\$"),
        # ylabel = L"(S-S_{\textrm{insta}})\,/\,\beta",
        ylabel = latexstring("\$ s-s_{\\mathrm{lower bound}} \$"), # L"s-\mathrm{Eq.}(7)/V",
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        yaxis = :log,
        left_margin = 2mm
    )

    col_ind = 0
    for meas in [1,2,6]
        col_ind +=1
        for i in eachindex(rhos)
            N_meas = N_smears[i]
            x_vals = rhos[i] .* Array(0:N_smears[i])
            sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt"))
            sms = (sms.-insta_action_U2_min(β,L,L,1)) ./ (L^2)
            image_smear = plot!(
                x_vals[start_inds[i]:end_inds[i]],
                sms[start_inds[i]:end_inds[i]],
                # color = greys[i],
                # color = palette(:default)[meas],
                color = cb_colors[col_ind],
                # label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$"),
                label = latexstring("Conf. Nr. $col_ind, \$\\rho = $(rhos[i])\$"),
                alpha = 1/i,
                yticks = [10.0^(-i) for i = 4:12],
                linewidth = 1.5,
            )
        end
    end
    display(image_smear)
end

image_smear_path = string(fig_path,"\\smeared_actions.pdf")
# savefig(image_smear_path)


        # i = 1
        # x_vals = rhos[i] .* Array(0:N_smears[i])
        # sms = readdlm(string(base_path,"\\sms_1.rho_nr_$i.txt"))
        # image_smear = plot!(
        #     x_vals[start_inds[i]:end],
        #     sms[start_inds[i]:end],
        #     color = :grey30,
        #     label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$")
        # )
        # for meas = 2:N_measurements
        #     sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt"))
        #     image_smear = plot!(
        #         x_vals[start_inds[i]:end],
        #         sms[start_inds[i]:end],
        #         color = :grey30,
        #         label = :false
        #     )
        # end
        
        # i = 2
        # x_vals = rhos[i] .* Array(0:N_smears[i])
        # sms = readdlm(string(base_path,"\\sms_1.rho_nr_$i.txt"))
        # image_smear = plot!(
        #     x_vals[start_inds[i]:end],
        #     sms[start_inds[i]:end],
        #     color = :grey50,
        #     label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$")
        # )
        # for meas = 2:N_measurements
        #     sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt"))
        #     image_smear = plot!(
        #         x_vals[start_inds[i]:end],
        #         sms[start_inds[i]:end],
        #         color = :grey50,
        #         label = :false
        #     )
        # end
        
        # i = 3
        # x_vals = rhos[i] .* Array(0:N_smears[i])
        # sms = readdlm(string(base_path,"\\sms_1.rho_nr_$i.txt"))
        # image_smear = plot!(
        #     x_vals[start_inds[i]:end],
        #     sms[start_inds[i]:end],
        #     color = :grey69,
        #     label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$")
        # )
        # for meas = 2:N_measurements
        #     sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt"))
        #     image_smear = plot!(
        #         x_vals[start_inds[i]:end],
        #         sms[start_inds[i]:end],
        #         color = :grey69,
        #         label = :false
        #     )
        # end





#=
i = 1
    x_vals = rhos[i] .* Array(0:N_smears[i])
    sms = readdlm(string(base_path,"\\sms_1.rho_nr_$i.txt"))
    image_smear = plot!(
        x_vals[start_inds[i]:end],
        sms[start_inds[i]:end],
        color = :grey30,
        label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$")
    )
    for meas = 2:N_measurements
        sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt"))
        image_smear = plot!(
            x_vals[start_inds[i]:end],
            sms[start_inds[i]:end],
            color = :grey30,
            label = :false
        )
    end
    
    i = 2
    x_vals = rhos[i] .* Array(0:N_smears[i])
    sms = readdlm(string(base_path,"\\sms_1.rho_nr_$i.txt"))
    image_smear = plot!(
        x_vals[start_inds[i]:end],
        sms[start_inds[i]:end],
        color = :grey50,
        label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$")
    )
    for meas = 2:N_measurements
        sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt"))
        image_smear = plot!(
            x_vals[start_inds[i]:end],
            sms[start_inds[i]:end],
            color = :grey50,
            label = :false
        )
    end
    
    i = 3
    x_vals = rhos[i] .* Array(0:N_smears[i])
    sms = readdlm(string(base_path,"\\sms_1.rho_nr_$i.txt"))
    image_smear = plot!(
        x_vals[start_inds[i]:end],
        sms[start_inds[i]:end],
        color = :grey69,
        label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$")
    )
    for meas = 2:N_measurements
        sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt"))
        image_smear = plot!(
            x_vals[start_inds[i]:end],
            sms[start_inds[i]:end],
            color = :grey69,
            label = :false
        )
    end
=#





######## Local Bogomolny ########





let
    L = 32
    z_bound = 3
    q_bound = 5
    q_over  = 0.2
    q_resol = 0.05
    y_bound_lower = -0.01
    y_bound_upper = 0.3
    x_bound_lower = -5.3
    x_bound_upper = 5.3
    q_vals = Array(-(q_bound+q_over):q_resol:(q_bound+q_over))
    z_vals = Array(-z_bound:z_bound)
    q_vals_insta = Array(-q_bound:q_bound)
    actions     = [insta_action_U2(1, L, L, q, z) for q in q_vals, z in z_vals ]
    actions_min = [insta_action_U2_min(1, L, L, q) for q in q_vals]
    insta_actions = [insta_action_U2(1, L, L, q, z) for q in q_vals_insta, z in z_vals ]
    markers = [:circle, :rect, :diamond, :star4]
    image_local = plot(
            q_vals,
            actions_min,
            color = :black,
            # linestyle = :dash,
            linewidth = 1.2,
            # label = latexstring("Eq.\$\\,\\,\$(7)"),
            label = "Lower bound",
            xticks = -q_bound:q_bound,
            xlabel = latexstring("Top. Charge \$q\$"),
            ylabel = L"S/\beta",
            size = (700,300),
            # tickfontsize = 9,
            labelfontsize = 13,
            legendfontsize = 9,
            # ylim = (y_bound_lower, y_bound_upper),
            # xlim = (x_bound_lower, x_bound_upper),
            grid = :false,
            legend = :outerright,
            # foreground_color_legend = :false,
            # background_color_legend = :false,
            # background_color_legend = RGBA(1.0,1.0,1.0,0.8),
            left_margin = 3mm,
            bottom_margin = 3mm,
    )
    for z in z_vals
        image_local = plot!(
            q_vals,
            actions[:,z+z_bound+1],
            label = :false, # latexstring("\$z = $z\$"),
            color = palette(:default)[z+z_bound+1],
            linestyle = :dash,
            linewidth = 1.2
        )
        image_local = scatter!(
            q_vals_insta,
            insta_actions[:,z+z_bound+1],
            label = latexstring("\$z = $z\$"), # :false,
            markershape = markers[mod(z,length(markers))+1],
            markersize = 5,
            color = palette(:default)[z+z_bound+1],
            alpha = 0.6,
        )
    end
    image_local = lens!(
        [0.8,1.201],
        [0.0,0.015],
        inset = (1, bbox(0.16,0.0,0.15,0.15)),
    )
    image_local = lens!(
        [1.8,2.2],
        [0.01,0.025],
        inset = (1, bbox(0.44,0.0,0.15,0.15)),
    )
    display(image_local)
end

image_local_path = string(fig_path,"\\local_insta.pdf")
# savefig(image_local_path)





######## Disturb the locals 🤙 ######## 





L = 32
q = 1
z = 5
N_smear = 10^4
ρ = 0.1
epsilons = round.(vcat([0.0], [10.0^(-i) for i = 2:7]), digits = 8)
dist_actions = Array{Float64}(undef,N_smear+1,length(epsilons))

#=
for i in eachindex(epsilons)
    ϵ = epsilons[i]
    U = [ran_U2(ϵ*rand()) for μ = 1:2, x = 1:L, t = 1:L ] .* insta_U2_z(L, L, q, z)
    smeared_actions = [action(U,1)]
    V = stout_midpoint(U,ρ)
    for smear = 1:N_smear
        push!(smeared_actions, action(V,1))
        V = stout_midpoint(V,ρ)
    end
    dist_actions[:,i] = smeared_actions
end
=#

dist_actions_path = string(data_path,"\\dist_actions.txt")
# # writedlm(dist_actions_path, dist_actions)
dist_actions = readdlm(dist_actions_path)

let
    τ0 = 2
    start_ind = Int(τ0/ρ)
    image_dist = plot(
        xlabel = latexstring("Flow time \$\\tau\$"),
        ylabel = L"S/\beta",
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 12,
        legend = :top,
        # foreground_color_legend = :false,
        background_color_legend = RGBA(1.0,1.0,1.0,0.8),
        # left_margin = 2mm,
    )
    # image_dist = plot!(
    #     ρ.*Array(0:N_smear)[start_ind:end],
    #     dist_actions[start_ind:end,1],
    #     label = latexstring("\$\\epsilon = 0.0\$"),
    #     color = cb_colors[7],
    #     linewidth = 2,
    #     alpha = 0.8
    # )
    # greys = (:grey77, :grey70, :grey63, :grey56, :grey49, :grey42, :grey84)
    # greys = (:grey45, :grey40, :grey35, :grey30, :grey25, :grey20, :grey50)
    greys = (:grey84, :grey77, :grey70, :grey63, :grey56, :grey49, :grey42)
    greys = (:grey77, :grey70, :grey63, :grey56, :grey49, :grey42, :grey30)
    # greys = (:grey75, :grey70, :grey65, :grey60, :grey55, :grey50, :grey45)
    for i = 1:length(epsilons)-1 # in eachindex(epsilons)
        j = length(epsilons)-i+1
        # j = i
        image_dist = plot!(
            ρ.*Array(0:N_smear)[start_ind:end],
            dist_actions[start_ind:end,j],
            label = latexstring("\$\\epsilon = $(epsilons[j])\$"),
            color = greys[i],
            linewidth = 2.5,
            alpha = 1
        )
    end
    image_dist = plot!(
        ρ.*Array(0:N_smear)[start_ind:end],
        dist_actions[start_ind:end,1],
        label = latexstring("\$\\epsilon = 0.0\$"),
        color = :black,
        linewidth = 2,
        alpha = 0.8
    )
    image_dist = hline!(
        [insta_action_U2_min(1,L,L,q)],
        label = latexstring("Lower bound for \$q=1\$"),
        color = cb_red,
        linestyle = :dot,
        legend = :right,
        linewidth = 2.5,
        alpha = 0.8,
        # grid = :false
    )
    display(image_dist)
end

image_dist_path = string(fig_path,"\\disturbed_locals.pdf")
# savefig(image_dist_path)





######## An unplanned endeavor 👀 ########




# let
    N_meas = 10000
    L = 32
    q = 1
    z = 5
    insta = insta_U2_z(L, L, q, z);
    epsilons = round.([10.0^(-i) for i = 1:7], digits = 8)
    prop_actions = []
    prop_errs = []
    for ϵ in epsilons
        actions = []
        for meas = 1:N_meas
            push!(actions, action([ran_U2(ϵ*rand()) for μ = 1:2, x = 1:L, t = 1:L ].*insta,1))
        end
        push!(prop_actions, mean(actions))
        push!(prop_errs, std(actions))
    end
    image_props = scatter(
        epsilons,
        prop_actions .- insta_action_U2(1, L, L, q, z),
        yerror = prop_errs,
        xaxis = :log,
        yaxis = :log,
        xticks = [10.0^(-i) for i = 1:9],
        yticks = [10.0^(-i) for i = 0:2:12],
        legend = :false,
        # size = (700,200),
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        xlabel = latexstring("Step Size \$\\epsilon\$"),
        ylabel = L"\Delta S/\beta",
        left_margin = 1mm,
        bottom_margin = 1mm,
        # ylim = (0.5*10.0^(-13), 0.0),
        markershape = :rtriangle,
        markersize = 8,
        # markercolor = cb_blue,
        # markerstrokecolor = cb_blue,
        )
# end
# last(prop_actions)
# last(prop_errs)
image_props_path = string(fig_path,"\\props.pdf")
# savefig(image_props_path)