include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\analyze\\SU2_analyze_head.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\analyze\\SU2_jackknives.jl")


using Plots
using LsqFit
using LaTeXStrings
using LinearAlgebra
# using Roots
using DelimitedFiles
using Statistics
# using Optim
using QuadGK
using Measures

data_path ="C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Lattice_projects\\Lattice2024\\data"
# data_path ="C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Lattice_projects\\Lattice2024\\PoS" 
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

function analytic_plaq_U2(β)
    # b = big(β)
    b = β
    numer(α) = besseli(0,b*cos(α)) + besseli(2,b*cos(α))
    denom(α) = 2*besseli(1,b*cos(α))/cos(α)
    return quadgk(numer,0,π/2)[1]/quadgk(denom,0,π/2)[1] - 1/b
end

function analytic_plaq_U2_approx(β)
    b = big(β)
    numer(α) = exp(b*cos(α)) / sqrt(cos(α)) * 2
    denom(α) = 2 * exp(b*cos(α)) / sqrt(cos(α))^3
    return quadgk(numer,0,π/2)[1]/quadgk(denom,0,π/2)[1] - 1/b
end

function analytic_plaq_U2_approx_2(β)
    b = big(β)
    numer(α) = exp(b*cos(α)) / sqrt(cos(α)) * (2 - 14/(8*b*cos(α)))
    denom(α) = 2 * exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)))
    return quadgk(numer,0,π/2)[1]/quadgk(denom,0,π/2)[1] - 1/b
end

function analytic_plaq_U2_approx_3(β)
    b = big(β)
    numer(α) = exp(b*cos(α)) / sqrt(cos(α)) * (2 - 14/(8*b*cos(α)) + 21/(32*(b*cos(α))^2))
    denom(α) = 2 * exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)) -15/(128*(b*cos(α))^2))
    return quadgk(numer,0,π/2)[1]/quadgk(denom,0,π/2)[1] - 1/b
end

function analytic_susc_U2(β)
    b = β
    denom(α)   = besseli(1,b*cos(α))/cos(α)
    numer(α) = α^2 * besseli(1,b*cos(α))/cos(α)
    return quadgk(numer,0,π/2)[1] / quadgk(denom,0,π/2)[1] / π^2
end

function analytic_susc_U2_approx(β)
    b = big(β)
    denom(α)   = exp(b*cos(α)) / sqrt(cos(α))^3
    numer(α) = α^2 * exp(b*cos(α)) / sqrt(cos(α))^3
    return quadgk(numer,0,π/2)[1] / quadgk(denom,0,π/2)[1] / π^2
end

function analytic_susc_U2_approx_2(β)
    b = big(β)
    denom(α)   = exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)))
    numer(α) = α^2 * exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)))
    return quadgk(numer,0,π/2)[1] / quadgk(denom,0,π/2)[1] / π^2
end

function analytic_susc_U2_approx_3(β)
    b = big(β)
    denom(α)   = exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)) - 15/(128*(b*cos(α))^2))
    numer(α) = α^2 * exp(b*cos(α)) / sqrt(cos(α))^3 * (1 - 3/(8*b*cos(α)) - 15/(128*(b*cos(α))^2))
    return quadgk(numer,0,π/2)[1] / quadgk(denom,0,π/2)[1] / π^2
end

function insta_U2_z(N_x, N_t, q, z)
    # w = π*(2*z-q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*q*t*π/N_x/N_t) * exp_u2(-t*2*π/N_x/N_t * (z-q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*q*x*π/N_x) * exp_u2(x*2*π/N_x * (z-q/2) * coeffs_U2(0.0im, 0.0im, 0.0im, complex(1.0))) for x = 1:N_x]
    return U
end

cb_colors = parse.(Colorant, ["#377eb8", "#ff7f00", "#4daf4a", "#e41a1c", "#999999", "#984ea3", "#f781bf"]);
cb_blue, cb_orange, cb_green, cb_red, cb_grey, cb_purple, cb_pink  = cb_colors;





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
#     actions_path = string(fig_path, "\\actions_L_$(Ls[i]).txt")
#     # # writedlm(queues_path, queues[i])
#     # # writedlm(actions_path, actions[i])
# end

# queues = []
# actions = []
# for i in eachindex(Ls)
#     push!(queues, readdlm(string(fig_path, "\\queues_L_$(Ls[i]).txt"))[:,1])
#     push!(actions, readdlm(string(fig_path, "\\actions_L_$(Ls[i]).txt"))[:,1])
# end

susc_vals = []
susc_errs = []
for i in eachindex(Ls)
    b_size = round(Int, 2*auto_corr_time(queues[i])+1)
    println("Block size of q at L = $(Ls[i]): ", b_size)
    jack = jackknife(queues[i].^2 ./ Ls[i]^2 .* betas[i] ./4, b_size)
    push!(susc_vals, jack[1])
    push!(susc_errs, jack[2])
end

# beta_range = Array(minimum(betas):0.01:maximum(betas))
# beta_range = Array(0.8:0.01:9.5)
beta_range = Array(1/1.05:0.01:700) # previous upper limit: 10
susc_anal = [analytic_susc_U2(β)*β/4 for β in beta_range]
beta_range_inv = 1 ./ beta_range

fits_executed = true

let
    lw = 1.8
    b_inv_fit_range = 0.0:0.05:0.25
    image_susc = plot(
        beta_range_inv,
        susc_anal,
        label = "Analytical result",
        linecolor = cb_orange, # palette(:default)[2], # cb_grey # cb_orange
        linewidth = lw,
    )
    image_susc = scatter!(
        1 ./ betas,
        susc_vals,
        yerror = susc_errs,
        # label = latexstring("Num. result for \$L=32\$"),
        label = "Simulation data",
        markershape = :diamond,
        markersize = 6,
        # markercolor = palette(:default)[1], # cb_red, #cb_blue
        # markerstrokecolor = :black # palette(:default)[1], # cb_red, #cb_blue
        color = cb_blue, # palette(:default)[1],
    )
    image_susc = plot!(
        xlabel = latexstring("\$1 / \\beta = (ag)^2/4\$"),
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
    if fits_executed
        image_susc = plot!(
            b_inv_fit_range,
            [model_lin(b_inv,fit_susc.param) for b_inv in b_inv_fit_range], 
            color = cb_green,
            label = "Slope at continuum value",
            linewidth = lw
        )
        image_susc = scatter!(
            [0.0],
            [1/4/π^2],
            label = latexstring("Continuum value \$1/(2\\pi)^2\$"),
            markersize = 8,
            markershape = :star5,
            color = cb_grey
        )
    end
    display(image_susc)
end

image_susc_path = string(fig_path, "\\susc_over_beta_inv.pdf")
# savefig(image_susc_path)


s_wil_vals = []
s_wil_errs = []
for i in eachindex(Ls)
    b_size = round(Int, 2*auto_corr_time(actions[i])+1)
    println("Block size of s_wil at L = $(Ls[i]): ", b_size)
    jack = jackknife(actions[i] ./ Ls[i]^2 .* betas[i] ./ 4, b_size)
    push!(s_wil_vals, jack[1])
    push!(s_wil_errs, jack[2])
end

beta_range = Array(1/1.05:0.01:700) # previous upper limit: 10
s_wil_anal = [(1-analytic_plaq_U2(β))*β/4 for β in beta_range]
beta_range_inv = 1 ./ beta_range

let
    lw = 1.8
    b_inv_fit_range = 0.0:0.05:0.2
    image_s_wil = plot(
        beta_range_inv,
        s_wil_anal,
        label = "Analytical result",
        linecolor = cb_orange, # palette(:default)[2], # cb_grey # cb_orange
        linewidth = lw,
    )
    image_s_wil = scatter!(
        1 ./ betas,
        s_wil_vals,
        yerror = s_wil_errs,
        # label = latexstring("Num. result for \$L=32\$"),
        label = "Simulation data",
        markershape = :rtriangle,
        markersize = 7,
        # markercolor = palette(:default)[1], # cb_red, #cb_blue
        # markerstrokecolor = :black # palette(:default)[1], # cb_red, #cb_blue
        color = cb_blue, # palette(:default)[1],
    )
    image_s_wil = plot!(
        xlabel = latexstring("\$1 / \\beta = (ag)^2/4 \$"),
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
    if fits_executed
        image_s_wil = plot!(
            b_inv_fit_range,
            [model_lin(b_inv,fit_s_wil.param) for b_inv in b_inv_fit_range], 
            color = cb_green,
            label = "Slope at continuum value",
            linewidth = lw
        )
        image_s_wil = scatter!(
            [0.0],
            [1/2],
            label = "Continuum value \$1/2\$",
            markersize = 7,
            markershape = :star5,
            color = cb_grey,
        )
    end
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
        label = "lower bound",
        color = cb_orange,# palette(:default)[2],
        linewidth = 1.5
    )
    image_insta = scatter!(
        -q_bound:q_bound,
        insta_actions,
        label = "instanton action",
        color = cb_blue,# palette(:default)[1],
        markershape = :diamond,
        markersize = 6,
    )
    image_insta = plot!(
        legend = :top,
        xlabel = latexstring("top. charge \$q\$"),
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
    τ   = 5000
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



# let
#     L   = 32
#     β   = 6.0
#     N_x = L
#     N_t = L
#     τ   = 10000 # 5000 is enough bro...
#     τ0  = 500
#     resol = 5
#     # rhos = [0.2, 0.1, 0.05]
#     rhos = [0.1, 0.05]
#     N_smears =   Int.(τ ./ rhos)
#     start_inds = Int.(τ0./ rhos)
#     N_measurements = 3

#     greys = [:grey30, :grey50, :grey69]
#     base_path = string(data_path, "\\smearing")
#     # base_path = data_path
#     # let
#     image_smear = plot(
#         # [insta_action_U2_min(1, N_x, N_t, 1)],
#         # color = cb_red,
#         # xticks = 0:resol:τ,
#         xlabel = latexstring("flow time \$\\tau/a^2\$"),
#         ylabel = L"S/\beta",
#         tickfontsize = 10,
#         labelfontsize = 17,
#         legendfontsize = 11,
#         # left_margin = 2mm
#     )
#     for i in eachindex(rhos)
#         x_vals = rhos[i] .* Array(0:N_smears[i])
#         sms = readdlm(string(base_path,"\\sms_1.rho_nr_$i.txt")) ./β
#         image_smear = plot!(
#             x_vals[start_inds[i]:end],
#             sms[start_inds[i]:end],
#             # color = greys[i],
#             color = palette(:default)[1],
#             label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$"),
#             alpha = 1/i,
#         )
#         for meas = 2:N_measurements
#             sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt")) ./β
#             image_smear = plot!(
#                 x_vals[start_inds[i]:end],
#                 sms[start_inds[i]:end],
#                 # color = greys[i],
#                 color = palette(:default)[meas],
#                 label = :false,
#                 alpha = 1/i
#             )
#         end
#         # println("$(greys[i])")
#     end
#     image_smear = hline!(
#         [insta_action_U2_min(1, N_x, N_t, 1)],
#         # label = latexstring("\$S_{\\textrm{Insta}}\$ with \$|q| = 1\$"),
#         label = latexstring("Eq.\$\\,\\,\$(7)"),
#         color = cb_red,
#         linestyle = :dash,
#     )
#     display(image_smear)
#     # end
# end

# image_smear_path = string(fig_path,"\\smeared_actions.pdf")
# # savefig(image_smear_path)


let
    L   = 32
    β   = 6.0
    N_x = L
    N_t = L
    τ   = 10000
    τ0  = 1
    τ1  = 1500
    rhos = [0.1, 0.05]
    alpha_for_rhos = [1.0, 0.5]
    linestyle_for_rhos = [:dash, :solid]
    N_smears =   Int.(τ ./ rhos)
    start_inds = Int.(τ0./ rhos)
    end_inds = Int.(τ1./ rhos)
    N_measurements = 3

    base_path = string(data_path, "\\smearing")

    image_smear = plot(
        # [insta_action_U2_min(1, N_x, N_t, 1)],
        # color = cb_red,
        # xticks = 0:resol:τ,
        xlabel = latexstring("flow time \$\\tau/a^2\$"),
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
                alpha = alpha_for_rhos[i],
                linestyle = linestyle_for_rhos[i],
                yticks = [10.0^(-i) for i = 0:2:12],
                linewidth = 1.5,
            )
        end
    end
    display(image_smear)
end

image_smear_path = string(fig_path,"\\smeared_actions.pdf")
# savefig(image_smear_path)





let
    # @assert 1==0 "Do we really want to start a smearing run???"
    β   = 6.0
    L   = 32
    N_x = L
    N_t = L
    τ   = 5000
    rhos = [0.1, 0.05]
    N_smears = Int.(τ ./ rhos)
    N_therm   = 500
    N_meas    = 10
    N_sepa    = 100
    acc_wish  = 0.8
    ϵ         = 0.1
    base_path = string(data_path, "\\smearing_new") # "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\sms\\sms_data_18"

    acc_therm = [0.0]
    U = gaugefield(N_x, N_t, true, "U2", "square")
    for therm = 1:N_therm
        chess_metro!(U,ϵ,β,acc_therm,"U2")
        ϵ *= sqrt(acc_therm[1] / acc_wish) # only update ϵ acc. to Metropolis
    end

    for meas = 1:N_meas
        for sepa = 1:N_sepa
            chess_metro!(U,ϵ,β,[0.0],"U2")
        end
        for rho_ind = 1:length(rhos)
            ρ = rhos[rho_ind]
            N_smear = N_smears[rho_ind]
            S_path = string(base_path,"\\sms_$meas.rho_nr_$rho_ind.txt")
            Q_path = string(base_path,"\\smq_$meas.rho_nr_$rho_ind.txt")
            smeared_actions = [action(U,β)]
            smeared_charges = [top_charge_U2(U)]
            V = stout_midpoint(U,ρ)
            count = 0
            for smear = 1:N_smear
                push!(smeared_actions,action(V,β))
                push!(smeared_charges,top_charge_U2(V))
                V = stout_midpoint(V,ρ)
                if smear%Int(N_smear/100) == 0
                    count += 1
                    println("Measurement Nr.: $meas, rho_nr.: $rho_ind, Smearing Progress: $count%")
                end
            end # smear
            writedlm(S_path,smeared_actions)
            writedlm(Q_path,smeared_charges)
        end # rho_ind
    end
    println("We're done here!")
end

let
    L   = 32
    β   = 6.0
    N_x = L
    N_t = L
    τ   = 5000
    τ0  = 1
    τ1  = 1200
    rhos = [0.1, 0.05]
    alpha_for_rhos = [1.0, 0.5]
    linestyle_for_rhos = [:dash, :solid]
    N_smears =   Int.(τ ./ rhos)
    start_inds = Int.(τ0./ rhos)
    end_inds = Int.(τ1./ rhos)

    base_path = string(data_path, "\\smearing_new")

    image_smear = plot(
        # [insta_action_U2_min(1, N_x, N_t, 1)],
        # color = cb_red,
        # xticks = 0:resol:τ,
        xlabel = latexstring("flow time \$\\tau/a^2\$"),
        # ylabel = L"(S-S_{\textrm{insta}})\,/\,\beta",
        ylabel = latexstring("\$ s-s_{\\mathrm{inst}}(q) \$"), # L"s-\mathrm{Eq.}(7)/V",
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        yaxis = :log,
        xticks = 0:200:1200,
        left_margin = 2mm
    )

    col_ind = 0
    for meas in [1,3,5,6,9]
        col_ind +=1
        for i in eachindex(rhos)
            N_meas = N_smears[i]
            x_vals = rhos[i] .* Array(0:N_smears[i])
            sms = readdlm(string(base_path,"\\sms_$meas.rho_nr_$i.txt")) 
            smq = readdlm(string(base_path,"\\smq_$meas.rho_nr_$i.txt"))
            q = round(Int,last(smq))
            sms = (sms.-insta_action_U2_min(β,L,L,q)) ./ (L^2) .+ 1e-12
            image_smear = plot!(
                x_vals[start_inds[i]:end_inds[i]],
                sms[start_inds[i]:end_inds[i]],
                # color = greys[i],
                # color = palette(:default)[meas],
                color = cb_colors[col_ind],
                # label = latexstring("\$\\rho = $(rhos[i]), \\, N_{\\textrm{smear}} = $(N_smears[i]) \$"),
                label = latexstring("Nr. $meas, \$q = $q\$, \$\\rho = $(rhos[i])\$"),
                alpha = alpha_for_rhos[i],
                linestyle = linestyle_for_rhos[i],
                yticks = [10.0^(-i) for i = 0:2:12],
                linewidth = 1.5,
            )
        end
    end
    display(image_smear)
end

image_smear_path = string(fig_path,"\\smeared_actions.pdf")
# savefig(image_smear_path)





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
            label = "lower bound",
            xticks = -q_bound:q_bound,
            xlabel = latexstring("top. charge \$q\$"),
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
            color = cb_colors[z+z_bound+1],# palette(:default)[z+z_bound+1],
            linestyle = :dash,
            linewidth = 1.2
        )
        image_local = scatter!(
            q_vals_insta,
            insta_actions[:,z+z_bound+1],
            label = latexstring("\$z = $z\$"), # :false,
            markershape = markers[mod(z,length(markers))+1],
            markersize = 5,
            color = cb_colors[z+z_bound+1],# palette(:default)[z+z_bound+1],
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
        xlabel = latexstring("flow time \$\\tau/a^2\$"),
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
        label = latexstring("lower bound for \$q=1\$"),
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





let
    N_meas = 10000
    L = 32
    q = 1
    z = 5
    insta = insta_U2_z(L, L, q, z);
    epsilons = round.([10.0^(-i) for i = 1:7], digits = 8)
    global prop_actions = []
    global prop_errs = []
    for ϵ in epsilons
        actions = []
        for meas = 1:N_meas
            push!(actions, action([ran_U2(ϵ) for μ = 1:2, x = 1:L, t = 1:L ].*insta,1))
        end
        push!(prop_actions, mean(actions))
        push!(prop_errs, std(actions)/sqrt(N_meas))
    end
    image_props = scatter(
        epsilons,
        prop_actions .- insta_action_U2(1, L, L, q, z),
        yerror = prop_errs,
        xaxis = :log,
        yaxis = :log,
        xticks = [10.0^(-i) for i = 1:9],
        yticks = [10.0^(-i) for i = -2:2:12],
        legend = :false,
        # size = (700,200),
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        xlabel = latexstring("step size \$\\epsilon\$"),
        ylabel = L"\Delta S/\beta",
        left_margin = 1mm,
        bottom_margin = 1mm,
        # ylim = (0.5*10.0^(-13), 0.0),
        markershape = :rtriangle,
        markersize = 8,
        # markercolor = cb_blue,
        # markerstrokecolor = cb_blue,
        )
end
# last(prop_actions)
# last(prop_errs)
image_props_path = string(fig_path,"\\props.pdf")
# savefig(image_props_path)





######## Yet another, less planned endeavor 🎩 ########





# let
    N_meas = 100
    L = 16
    q = 1
    z = 5
    insta = insta_U2_z(L, L, q, z);
    epsilons = [0.01] # round.([10.0^(-i) for i = 1:7], digits = 8)
    prop_actions = Array{Float64}(undef, 2, L, L)
    prop_errs = Array{Float64}(undef, 2, L, L)
    for ϵ in epsilons
        for μ = 1:2
        for x = 1:L
        for t = 1:L
        stap = staple_dag(insta, μ, x, t)
        link_old = insta[μ,x,t]
        s_old = -0.5*real(tr(link_old * stap))
        actions = []
            for meas = 1:N_meas
                link_new = ran_U2(ϵ) * link_old
                s_new = -0.5*real((tr(link_new * stap)))
                push!(actions, s_new - s_old)
            end
        prop_actions[μ,x,t] = mean(actions)
        prop_errs[μ,x,t] = std(actions)/sqrt(N_meas)
        end
        end
        end
    end
    image_props_2 = heatmap(prop_actions[2,:,:].*100, bins = 16^2, right_margin = 10mm)
# end
display(image_props_2)

image_props_2_path = string(fig_path,"\\props_2.pdf")
# savefig(image_props_2_path)






#=
function analytic_susc_U2_integrand(α, β)
    denom(α)   = besseli(1,β*cos(α))/cos(α)
    numer(α) = α^2 * besseli(1,β*cos(α))/cos(α)
    return numer(α)
end

betas = [1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0]

image_int = plot(
    0.0:0.01:π/2,
    [analytic_susc_U2_integrand(α,betas[1]) for α = 0.0:0.01:π/2],
    label = "β = $(betas[1])",
    title = "Denominator of χ",
)
for beta in betas[2:end]
    image_int = plot!(
        0.0:0.01:π/2,
        [analytic_susc_U2_integrand(α,beta) for α = 0.0:0.01:π/2],
        label = "β = $beta",
        
    )
end
display(image_int)

for beta in betas
    # beta = 1
    alpha_range = 0.0:0.01:π/2
    numerator_vals = [analytic_susc_U2_integrand(α,beta) for α in alpha_range]
    last_numerator = round(last(numerator_vals), sigdigits = 3)
    image_int = plot(
        alpha_range,
        numerator_vals,
        xticks = (Vector(0:π/8:π/2), ["0.0", "π/8", "π/4", "3π/8", "π/2"]),
        label = "β = $(beta)\nlast val: $last_numerator",
        title = "Numerator of χ",
    )
    display(image_int)
end

for beta in betas
    # beta = 1
    alpha_range = 0.0:0.01:π/2
    denominator_vals = [analytic_susc_U2_integrand(α,beta)/α^2 for α in alpha_range]
    last_denominator = round(last(denominator_vals), sigdigits = 3)
    image_int = plot(
        alpha_range,
        denominator_vals,
        xticks = (Vector(0:π/8:π/2), ["0.0", "π/8", "π/4", "3π/8", "π/2"]),
        label = "β = $(beta)\nlast val: $last_denominator",
        title = "Numerator of χ",
    )
    display(image_int)
end

analytic_susc_U2_integrand(π/2,500)/pi^2*4^2

=#

beta_range_anal = Array(300:1:700) # previous upper limit: 10
susc_anal = [analytic_susc_U2(β)*β/4 for β in beta_range_anal]
beta_range_approx = Array(300:1:1e4)
susc_approx = [analytic_susc_U2_approx(β)*β/4 for β in beta_range_approx]
beta_range_approx_2 = Array(300:1:1e4)
susc_approx_2 = [analytic_susc_U2_approx_2(β)*β/4 for β in beta_range_approx_2]
beta_range_approx_3 = Array(300:1:1e4)
susc_approx_3 = [analytic_susc_U2_approx_3(β)*β/4 for β in beta_range_approx_3]

let
    image_susc_approx = plot(
        1 ./beta_range_anal,
        susc_anal,
        label = "full anal. result",
        linecolor = palette(:default)[1], # cb_grey # cb_orange
        linewidth = 1.5,
        xlabel = latexstring("\$1 / \\beta\$"),
        ylabel = L"\chi_{\textbf{top}} \, / g^2",
        legend = :topleft,
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 10,
        left_margin = 2mm,
        right_margin = 2mm,
        xlim = (-0.0001, 0.0036),
        xticks = 0.0:0.0005:0.0035,
    )
    image_susc_approx = plot!(
        1 ./ beta_range_approx,
        susc_approx,
        label = "large arg. expan. k = 0",
        linecolor = palette(:default)[2], # cb_grey # cb_orange
        linewidth = 1.5,
        linestyle = :dash,
    )
    image_susc_approx = plot!(
        1 ./ beta_range_approx_2,
        susc_approx_2,
        label = "large arg. expan. k = 0,1",
        linecolor = palette(:default)[3], # cb_grey # cb_orange
        linewidth = 1.5,
        linestyle = :dashdot,
    )
    image_susc_approx = plot!(
        1 ./ beta_range_approx_3,
        susc_approx_3,
        label = "large arg. expan. k = 0,1,2",
        linecolor = palette(:default)[4], # cb_grey # cb_orange
        linewidth = 1.5,
        linestyle = :dashdotdot,
    )
    image_susc_approx = plot!(
        1 ./beta_range_anal,
        susc_anal,
        label = :false,# "full analyt. result",
        linecolor = palette(:default)[1], # cb_grey # cb_orange
        linewidth = 1.5,
    )
    display(image_susc_approx)
end





# beta_range_anal = Array(300:1:700) # previous upper limit: 10
s_wil_anal = [(1-analytic_plaq_U2(β))*β/4 for β in beta_range_anal]
beta_range_approx = Array(300:1:1e4)
s_wil_approx = [(1-analytic_plaq_U2_approx(β))*β/4 for β in beta_range_approx]
beta_range_approx_2 = Array(300:1:1e4)
s_wil_approx_2 = [(1-analytic_plaq_U2_approx_2(β))*β/4 for β in beta_range_approx_2]
beta_range_approx_3 = Array(300:1:1e4)
s_wil_approx_3 = [(1-analytic_plaq_U2_approx_3(β))*β/4 for β in beta_range_approx_3]

let
    image_s_wil_approx = plot(
        1 ./beta_range_anal,
        s_wil_anal,
        label = "full anal. result",
        linecolor = palette(:default)[1], # cb_grey # cb_orange
        linewidth = 1.5,
        xlabel = latexstring("\$1 / \\beta\$"),
        ylabel = L"s_{\mathrm{wil}} \, / g^2",
        legend = :right,
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 11,
        left_margin = 2mm,
        xlim = (-0.0001, 0.0036),
        xticks = 0.0:0.0005:0.0035,
        yticks = 0.38:0.02:0.5,
        # xlim = (0.0,1.1),
    )
    image_s_wil_approx = plot!(
        1 ./beta_range_approx,
        s_wil_approx,
        label = "large arg. expan. k = 0",
        linecolor = palette(:default)[2], # cb_grey # cb_orange
        linewidth = 1.5,
        linestyle = :dash
    )
    image_s_wil_approx = plot!(
        1 ./beta_range_approx_2,
        s_wil_approx_2,
        label = "large arg. expan. k = 0,1",
        linecolor = palette(:default)[3], # cb_grey # cb_orange
        linewidth = 1.5,
        linestyle = :dashdot
    )
    image_s_wil_approx = plot!(
        1 ./beta_range_approx_3,
        s_wil_approx_3,
        label = "large arg. expan. k = 0,1,2",
        linecolor = palette(:default)[4], # cb_grey # cb_orange
        linewidth = 1.5,
        linestyle = :dashdotdot
    )
    image_s_wil_approx = plot!(
        1 ./beta_range_anal,
        s_wil_anal,
        label = :false,# "full anal. result",
        linecolor = palette(:default)[1], # cb_grey # cb_orange
        linewidth = 1.5,
    )
    display(image_s_wil_approx)
end

model_lin(x, p) = p[1] .+ p[2]*x
p_susc  = [(2*π)^2, 0.1]
p_s_wil = [0.5, 0.0]
fit_susc = curve_fit(model_lin, 1 ./beta_range_anal, susc_anal, p_susc)
fit_s_wil = curve_fit(model_lin, 1 ./beta_range_anal, s_wil_anal, p_s_wil)
params_s_wil = round.(fit_s_wil.param, sigdigits = 3)
params_susc = round.(fit_susc.param, sigdigits = 3)

inv_beta_range_fit = Vector(0.0:1e-4:0.0025)

let
    image_susc_fit = plot(
        1 ./beta_range_anal,
        susc_anal,
        label = "full anal. result",
        linecolor = palette(:default)[1], # cb_grey # cb_orange
        linewidth = 1.5,
        xlabel = latexstring("\$1 / \\beta\$"),
        ylabel = L"\chi_{\textbf{top}} \, / g^2",
        legend = :topleft,
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 10,
        left_margin = 2mm,
        right_margin = 2mm,
        xlim = (-0.0001, 0.0036),
        xticks = 0.0:0.0005:0.0035,
    )
    image_susc_fit = plot!(
        inv_beta_range_fit,
        [model_lin(b_inv,fit_susc.param) for b_inv in inv_beta_range_fit],
        label = latexstring("\$\\mathrm{fit}(\\beta^{-1}) = a + b\\cdot \\beta^{-1} \$\n\$ a = $(params_susc[1]), b = $(params_susc[2])\$\n\$a\\cdot(2\\pi)^2 = $(fit_susc.param[1]*(2*pi)^2)\$ ")
    )
    display(image_susc_fit)
end

let
    image_s_wil_fit = plot(
        1 ./beta_range_anal,
        s_wil_anal,
        label = "full anal. result",
        linecolor = palette(:default)[1], # cb_grey # cb_orange
        linewidth = 1.5,
        xlabel = latexstring("\$1 / \\beta\$"),
        ylabel = L"s_{\mathrm{wil}} \, / g^2",
        legend = :topleft,
        tickfontsize = 10,
        labelfontsize = 17,
        legendfontsize = 10,
        left_margin = 2mm,
        xlim = (-0.0001, 0.0036),
        xticks = 0.0:0.0005:0.0035,
        # yticks = 0.38:0.02:0.5,
        # xlim = (0.0,1.1),
    )
    image_s_wil_fit = plot!(
        inv_beta_range_fit,
        [model_lin(b_inv,fit_s_wil.param) for b_inv in inv_beta_range_fit],
        label = latexstring("\$\\mathrm{fit}(\\beta^{-1}) = a + b\\cdot \\beta^{-1} \$\n\$ a = $(params_s_wil[1]), b = $(params_s_wil[2])\$")
    )
    display(image_s_wil_fit)
end





########  *How* special are our special configs?  ########





#=
function ran_U2_lie_direction(lie_dir, ϵ)
    v = [0.0, 0.0, 0.0]
    if lie_dir == 0
        return coeffs_U2(exp(im*ϵ/2), complex(0.0), complex(0.0), complex(0.0))
    elseif lie_dir == 1
        return coeffs_U2(complex(cos(ϵ/2)), complex(sin(ϵ/2)), complex(0.0),    complex(0.0))
    elseif lie_dir == 2
        return coeffs_U2(complex(cos(ϵ/2)), complex(0.0),    complex(sin(ϵ/2)), complex(0.0))
    elseif lie_dir == 3
        return coeffs_U2(complex(cos(ϵ/2)), complex(0.0),    complex(0.0),    complex(sin(ϵ/2)))
    else
        error("Group directions 'lie_dir' are labeled from 0 to 3")
    end
end

# abs(det(ran_U2_lie_direction(rand(0:3),rand())))
# bla = ran_U2_lie_direction(rand(0:3),rand())
# coeffs2grp(bla * adjoint(bla))
# for r = 0:3
    # reps = rand()
    # @assert isapprox(coeffs2grp(ran_U2_lie_direction(r,reps)), exp(im*reps*Σ[r+1]))
# end

function vary_config_U2(U, ϵ)
    NX = size(U,2)
    NT = size(U,3)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                stap_dag = staple_dag(U, μ, x, t)
                old_link = U[μ,x,t]
                s_old = -real(tr(old_link*stap_dag))
                for lie_dir in [0,1,2,3]
                    new_link_p = ran_U2_lie_direction(lie_dir,+ϵ) * old_link
                    new_link_m = ran_U2_lie_direction(lie_dir,-ϵ) * old_link
                    delta_s_p = -real(tr(new_link_p*stap_dag)) - s_old
                    delta_s_m = -real(tr(new_link_m*stap_dag)) - s_old
                    if delta_s_p < 0
                        error("Action decrease in μ = $μ, x = $x, t = $t, ϵ = $ϵ in pos. group direction $lie_dir")
                    elseif delta_s_m < 0
                        error("Action decrease in μ = $μ, x = $x, t = $t, ϵ = $ϵ in neg. group direction $lie_dir")
                    end
                end
            end
        end
    end
    return nothing
end

function vary_config_ran_U2(U, ϵ, N_prop)
    NX = size(U,2)
    NT = size(U,3)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                stap_dag = staple_dag(U, μ, x, t)
                old_link = U[μ,x,t]
                s_old = -real(tr(old_link*stap_dag))
                for N = 1:N_prop
                    # r = ran_U2(ϵ*(2*rand()-1))
                    r = ran_U2(ϵ)
                    new_link = r * old_link
                    s_new = -real(tr(new_link*stap_dag))
                    if s_new-s_old < 0 
                        error("Action decrease in μ = $μ, x = $x, t = $t, ϵ = $ϵ for random element $(coeffs2grp(r))")
                    end
                end
            end
        end
    end
end


N_x = N_t = 32
for q = -10:10
    for z = -10:10
        insta = insta_U2_z(N_x, N_t, q, z)
        for pow = 8:8
            vary_config_U2(insta,10.0^(-pow))
            # vary_config_ran_U2(insta,10.0^(-pow),10)
        end
    end
end

let
    N_x = N_t = 32
    q = 2
    z = 1
    ϵ = 1e-7
    insta = insta_U2_z(N_x, N_t, q, z)
    vary_config_U2(insta,ϵ)
end

# bla = ran_U2_lie_direction(rand(0:3),1e-3); coeffs2grp(bla)[1,1]

function vary_special_config_U2(tol, L, q_vals, z_vals, ϵ_powers, file_name)
    mac_prec = eps()
    for q in q_vals
    for z in z_vals; insta = insta_U2_z(L, L, q, z)
    for pow in ϵ_powers; ϵ = 10.0^(-pow)
        for t = 1:L
        for x = 1:L
        for μ = 1:2
            stap_dag = staple_dag(insta, μ, x, t)
            old_link = insta[μ,x,t]
            s_old    = -real(tr(old_link*stap_dag))
            for lie_dir in [0,1,2,3]
                new_link_p = ran_U2_lie_direction(lie_dir,+ϵ) * old_link
                new_link_m = ran_U2_lie_direction(lie_dir,-ϵ) * old_link
                delta_s_p = -real(tr(new_link_p*stap_dag)) - s_old
                delta_s_m = -real(tr(new_link_m*stap_dag)) - s_old
                if delta_s_p < tol * mac_prec
                    open(file_name, "a") do io
                        write(io, "q = $q, z = $z, μ = $μ, x = $x, t = $t, pow = -$pow, lie_dir = +$lie_dir, ΔS/eps = $(delta_s_p/mac_prec)\n")
                    end
                end
                if delta_s_p < tol * mac_prec
                    open(file_name, "a") do io
                        write(io, "q = $q, z = $z, μ = $μ, x = $x, t = $t, pow = -$pow, lie_dir = -$lie_dir, ΔS/eps = $(delta_s_m/mac_prec) \n")
                    end
                end
            end
        end # t
        end # x
        end # μ
    end # q
    end # z
    end # pow
    return nothing
end

for hitsize_pow = 1:7
    # hitsize = 11
    vary_file_name = string(data_path, "\\vary_hitsize_$hitsize_pow.txt")
    vary_special_config_U2(10, 32, 0:5, 0:5, hitsize_pow:hitsize_pow, vary_file_name)
end
=#


function ran_U2_lie_direction(lie_dir, ϵ)
    v = [0.0, 0.0, 0.0]
    if lie_dir == 0
        return coeffs_U2(exp(im*ϵ/2), complex(0.0), complex(0.0), complex(0.0))
    elseif lie_dir == 1
        return coeffs_U2(complex(cos(ϵ/2)), complex(sin(ϵ/2)), complex(0.0),    complex(0.0))
    elseif lie_dir == 2
        return coeffs_U2(complex(cos(ϵ/2)), complex(0.0),    complex(sin(ϵ/2)), complex(0.0))
    elseif lie_dir == 3
        return coeffs_U2(complex(cos(ϵ/2)), complex(0.0),    complex(0.0),    complex(sin(ϵ/2)))
    else
        error("Group directions 'lie_dir' are labeled from 0 to 3")
    end
end

function sdiff_min_of_special_config(L, q, z, μ, x, t, hitsize_range, lie_dir)
    insta    = insta_U2_z(L, L, q, z)
    old_link = insta[μ,x,t]
    stap_dag = staple_dag(insta,μ,x,t)
    s_old    = -real(tr(old_link*stap_dag))
    s_diffs  = []
    for hitsize in hitsize_range
        new_link = ran_U2_lie_direction(lie_dir, hitsize) * old_link
        s_new    = -real(tr(new_link*stap_dag))
        push!(s_diffs, s_new-s_old)
    end
    return s_diffs
end

model_sq(x, p) = p[1] .*x.^2
p_sdiff = [1.0]
all_params = []

let
    L = 32
    # q = 0 #  rand(1:5)
    # z = 0 #  rand(1:5)
    for q = 0:3
    for z = 0:3
    μ_vals = rand(1:2, 4)
    x_vals = [rand(1:L-1), rand(1:L-1), L,           L]
    t_vals = [rand(1:L-1), L,           rand(1:L-1), L]
    hitsize_range  = Vector(-3.5e-7:1e-7:3.5e-7)
    fit_plot_xvals = Vector(-4e-7:1e-8:4e-7)
    for site_ind = 1:4
        μ = μ_vals[site_ind]
        x = x_vals[site_ind]
        t = t_vals[site_ind]
        for lie_dir = 0:3
            sdiff        = sdiff_min_of_special_config(L, q, z, μ, x, t, hitsize_range, lie_dir) ./ eps()
            fit_sdiff    = curve_fit(model_sq, 1e7 .* hitsize_range, sdiff, p_sdiff)
            params_sdiff = round.(fit_sdiff.param .* 1e14, sigdigits = 3)
            image_sdiff = plot(
                title = latexstring("\$ \\Delta S \$ of the special config at \$(q,z) = ($q,$z)\$ \n \$(\\mu, x, t) = ($μ,$x,$t)\$, Lie-dir.: \$ $lie_dir\$"),
                xlabel = "hitsize",
                ylabel = latexstring("\$\\Delta S / \\left(\\beta \\cdot \\epsilon_\\textrm{machine} \\right)\$"),
                rightmargin = 5mm,
                labelfontsize = 15,
                tickfontsize = 10,
                legend = :top
            )
            image_sdiff = plot!(
                fit_plot_xvals,
                model_sq(fit_plot_xvals, 1e14.*fit_sdiff.param),
                label = latexstring("\$\\mathrm{fit}(x) = a\\cdot x^2\$ \n \$a = $(params_sdiff[1])\$"),
                color = cb_orange
            )
            image_sdiff = scatter!(
                hitsize_range,
                sdiff,
                label = "Measurements",
                color = cb_blue
            )
            # display(image_sdiff)
            push!(all_params, fit_sdiff.param[1])
        end # lie_dir
    end # site_ind
    end # z
    end # q
end # let

# mean(all_params)
# std(all_params)/16



let
    L = 32
    q = rand(1:5)
    z = rand(1:5)
    μ = rand(1:2)
    x = rand(1:L)
    t = rand(1:L)
    hitsize_range =  Vector(1e-8:1e-9:5e-7)
    for lie_dir = 0:3
        sdiff_plu = sdiff_min_of_special_config(L, q, z, μ, x, t, hitsize_range, lie_dir) ./ eps()
        sdiff_min = sdiff_min_of_special_config(L, q, z, μ, x, t, -hitsize_range, lie_dir) ./ eps()
        image_sdiff = plot(
            hitsize_range,
            sdiff_plu,
            title = latexstring("\$ \\Delta S \$ of the special config at \$(q,z) = ($q,$z)\$ \n \$(\\mu, x, t) = ($μ,$x,$t)\$, Lie-dir.: \$ $lie_dir\$"),
            label = :false,
            xlabel = "hitsize",
            ylabel = latexstring("\$\\Delta S / \\left(\\beta \\cdot \\epsilon_\\textrm{machine} \\right)\$"),
            rightmargin = 5mm,
            labelfontsize = 15,
            tickfontsize = 10,
            color = palette(:default)[1]
        )
        image_sdiff = plot!(
            -hitsize_range,
            sdiff_min,
            color = palette(:default)[1],
            label = :false
        )
        display(image_sdiff)
    end
end