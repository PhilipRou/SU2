include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")
using LsqFit
using QuadGK

function analytic_plaq_U2(β)
    numer(α) = besseli(0,β*cos(α)) + besseli(2,β*cos(α))
    denom(α) = 2*besseli(1,β*cos(α))/cos(α)
    return quadgk(numer,0,π)[1]/quadgk(denom,0,π)[1] - 1/β
end

function analytic_string_U2(β)
    return -log(analytic_plaq_U2(β))
end

betas = [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0]


all_ev_plat_ms      = Array{Float64}(undef, length(betas), 1)
all_ev_plat_m_errs  = Array{Float64}(undef, length(betas), 1)
all_GEV_plat_ms     = Array{Float64}(undef, length(betas), 2)
all_GEV_plat_m_errs = Array{Float64}(undef, length(betas), 2)

for beta_ind in eachindex(betas)
    # beta_ind = 8
    beta = betas[beta_ind]
    # beta = 20.0
    Nx=Ny=Nz  = round(Int, beta*4)
    beta      = string(beta)*"0"
    n_stout   = 7
    smearlist = [0,1,3,7,15]
    rho_2D    = 0.24
    rho_3D    = 0.16

    
    base_path               = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\beta_$(beta)_Nz_$(Nz)_Nx_$(Nx)_smearlist_rho_2D_$(rho_2D)"
    timeseries_path         = base_path * "_hist_3D.txt"
    timeseries_2Dsmear_path = base_path * "_hist_2D.txt"
    corr_path               = base_path * "_crosscorr.txt"
    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\beta_$(beta)"
    # mkpath(base_path)
    base_path = base_path * "\\beta_$(beta)_Nz_$(Nz)_Nx_$(Nx)_smearlist_rho_2D_$(rho_2D)"


    plat_ranges_1     = [1:5, 1:4, 3:7, 2:5, 1:7, 4:8, 4:8, 3:6]
    plat_ranges_2     = [1:1, 1:1, 1:1, 1:1, 1:1, 1:1, 1:1, 1:1]
    
    plat_ranges_GEV_1 = [1:2, 1:2, 1:3, 1:4, 1:3, 1:4, 2:8, 3:5]
    plat_ranges_GEV_2 = [1:1, 1:2, 1:2, 1:3, 1:3, 1:2, 3:5, 5:12]
        
    plat_range_1 = plat_ranges_1[beta_ind]
    plat_range_2 = plat_ranges_2[beta_ind]
    plot_range_m_ev_dr = minimum([plat_range_1[1], plat_range_2[1]]):maximum([plat_range_1[end],plat_range_2[end]])+1
    plat_range_1_GEV_dr = plat_ranges_GEV_1[beta_ind]
    plat_range_2_GEV_dr = plat_ranges_GEV_2[beta_ind]
    plot_range_m_GEV_dr = minimum([plat_range_1_GEV_dr[1], plat_range_2_GEV_dr[1]]):maximum([plat_range_1_GEV_dr[end],plat_range_2_GEV_dr[end]])+1



    


    evs_dr             = readdlm(base_path * "_evs_8x8_dr.txt")
    evs_errs_dr        = readdlm(base_path * "_evs_errs_8x8_dr.txt")
    m_evs_dr           = readdlm(base_path * "_m_evs_8x8_dr.txt")
    m_evs_errs_dr      = readdlm(base_path * "_m_evs_errs_8x8_dr.txt")
    m_evs_plat_dr      = readdlm(base_path * "_m_evs_plat_8x8_dr.txt")
    m_evs_plat_errs_dr = readdlm(base_path * "_m_evs_plat_errs_8x8_dr.txt")

    all_ev_plat_ms[beta_ind, 1]     = m_evs_plat_dr[1]
    all_ev_plat_m_errs[beta_ind, 1] = m_evs_plat_errs_dr[1]



    let
        image_conn_evs = plot(
            # title  = latexstring("EV's of conn. corr. 8² mats. (bottom right) \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:4 # size(evs,2)
            image_conn_evs = scatter!(
                Vector(0:Nz>>1) .+ (i-1)*0.1,
                evs_dr[:,i],
                yerror = evs_errs_dr[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-16, 10*maximum(evs_dr)),
                yaxis = :log,
                legend = :topright,
                tickfontsize = 10,
                labelfontsize = 15,
                legendfontsize = 11,
                background_color_legend = nothing
            )
        end
        display(image_conn_evs)
        fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\m_plateaus\\evs_beta_$beta.pdf"
        # savefig(fig_path)
    end
    let
        image_conn_evs_masses = plot(
            # title  = latexstring("2pt Masses of EV's of conn. corr. 8² mats. (bottom right) \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2#size(m_evs_dr,2)
            image_conn_evs_masses = scatter!(
                (Vector(0:Nz>>1) .+ 0.5 .+ (i-1)*0.05)[plot_range_m_ev_dr],
                m_evs_dr[plot_range_m_ev_dr,i],
                yerror = m_evs_errs_dr[plot_range_m_ev_dr,i],
                label = latexstring("EV-mass nr. \$$i\$"),
                color = cb_colors[i],
                markerstrokecolor = cb_colors[i],
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :bottomleft,
                markershape = [:circ, :diamond][i],
                # rightmargin = 5mm,
                tickfontsize = 10,
                labelfontsize = 15,
                legendfontsize = 11,
                background_color_legend = nothing
            )
        end
        plat_range_1_plot = vcat(Vector(plat_range_1), [last(plat_range_1)+1])
        image_conn_evs_masses = plot!(
            (Vector(0:Nz>>1) .+ 0.5)[plat_range_1_plot],
            m_evs_plat_dr[1] .* ones(length(plat_range_1_plot)),
            ribbon = m_evs_plat_errs_dr[1] .* ones(length(plat_range_1_plot)),
            color = cb_colors[1],
            lw = 2,
            label = "plateau nr. 1",
        )
        # image_conn_evs_masses = plot!(
        #     (Vector(0:Nz>>1) .+ 0.55)[plat_range_2],
        #     m_evs_plat_dr[2] .* ones(length(plat_range_2)),
        #     ribbon = m_evs_plat_errs_dr[2] .* ones(length(plat_range_2)),
        #     color = cb_colors[2],
        #     lw = 2,
        #     label = "plateau nr. 2",
        #     linestyle = :dash,
        # )
        display(image_conn_evs_masses)
        # fig_path = base_path * "_m_evs_dr_8x8.pdf"
        fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\m_plateaus\\m_plateau_ev_beta_$beta.pdf"
        # savefig(fig_path)
    end





    t0 = 1
    GEVs_dr             = readdlm(base_path * "_GEVs_t0_$(t0)_8x8_dr.txt")
    GEVs_errs_dr        = readdlm(base_path * "_GEVs_errs_t0_$(t0)_8x8_dr.txt")
    m_GEVs_dr           = readdlm(base_path * "_m_GEVs_t0_$(t0)_8x8_dr.txt")
    m_GEVs_errs_dr      = readdlm(base_path * "_m_GEVs_errs_t0_$(t0)_8x8_dr.txt")
    m_GEVs_plat_dr      = readdlm(base_path * "_m_GEVs_plat_t0_$(t0)_8x8_dr.txt")
    m_GEVs_plat_errs_dr = readdlm(base_path * "_m_GEVs_plat_errs_t0_$(t0)_8x8_dr.txt")

    all_GEV_plat_ms[beta_ind,:]     = m_GEVs_plat_dr
    all_GEV_plat_m_errs[beta_ind,:] = m_GEVs_plat_errs_dr

    let
        image_conn_GEVs = plot(
            # title  = latexstring("GEV's of conn. corr. 8² mats. (bottom right) \n \$\\beta = $beta\$, \$L = $Nz\$, \$t_0 = $(t0-1)\$, \$n_\\mathrm{smear} \\in\$ $(smearlist[2:end])"),
            xlabel = latexstring("\$ t\$"),
            )
        for i = 1:4 # size(GEVs,2)
            image_conn_GEVs = scatter!(
                Vector(t0+1:Nz>>1+1) .+ (i-1)*0.02,
                GEVs_dr[:,i],
                yerror = GEVs_errs_dr[:,i],
                label = latexstring("GEV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-4, 10*maximum(GEVs_dr)),
                yaxis = :log,
                # legend = :topright,
                tickfontsize = 10,
                labelfontsize = 15,
                legendfontsize = 11,
                # background_color_legend = nothing
            )
        end
        display(image_conn_GEVs)
        # fig_path = base_path * "_GEVs_dr_8x8.pdf"
        fig_path = fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\m_plateaus\\GEVs_beta_$beta.pdf"
        # savefig(fig_path)
    end

    let
        image_conn_GEVs_masses = plot(
            # title  = latexstring("2pt Masses of GEV's of conn. corr. 8² mats. (bottom right)\n \$\\beta = $beta, L = $Nz, t_0 = $(t0-1), n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2 # size(m_evs_dr,2)
            image_conn_GEVs_masses = scatter!(
                (Vector(t0:Nz>>1) .+ 0.5 .+ (i-1)*0.05)[plot_range_m_GEV_dr],
                m_GEVs_dr[plot_range_m_GEV_dr,i],
                yerror = m_GEVs_errs_dr[plot_range_m_GEV_dr,i],
                label = latexstring("GEV-mass nr. \$$i\$"),
                color = cb_colors[i],
                markerstrokecolor = cb_colors[i],
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                # legend = :topright,
                markershape = [:circ, :diamond][i],
                tickfontsize = 10,
                labelfontsize = 15,
                legendfontsize = 11,
                # background_color_legend = nothing
            )
        end
        image_conn_evs_masses = plot!(
            (Vector(t0:Nz>>1) .+ 0.5)[plat_range_1_GEV_dr],
            m_GEVs_plat_dr[1] .* ones(length(plat_range_1_GEV_dr)),
            ribbon = m_GEVs_plat_errs_dr[1] .* ones(length(plat_range_1_GEV_dr)),
            color = cb_colors[1],
            lw = 2,
            label = "plateau nr. 1",
        )
        image_conn_evs_masses = plot!(
            (Vector(t0:Nz>>1) .+ 0.55)[plat_range_2_GEV_dr],
            m_GEVs_plat_dr[2] .* ones(length(plat_range_2_GEV_dr)),
            ribbon = m_GEVs_plat_errs_dr[2] .* ones(length(plat_range_2_GEV_dr)),
            color = cb_colors[2],
            lw = 2,
            label = "plateau nr. 2",
            linestyle = :dash,
        )
        display(image_conn_GEVs_masses)
        # fig_path = base_path * "_m_plateau_GEVs_dr_8x8.pdf"
        fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\m_plateaus\\m_plateau_GEV_beta_$beta.pdf"
        # savefig(fig_path)
    end
end

# writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\all_ev_plat_ms.txt", all_ev_plat_ms)
# writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\all_ev_plat_m_errs.txt", all_ev_plat_m_errs)
# writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\all_GEV_plat_ms.txt", all_GEV_plat_ms)
# writedlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\all_GEV_plat_m_errs.txt", all_GEV_plat_m_errs)












include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")
using LsqFit
using QuadGK

betas = [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0]
strings = [NaN, 0.3129, 0.2562, NaN, NaN, 0.1179, NaN, NaN]

all_ev_plat_ms = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\all_ev_plat_ms.txt")
all_ev_plat_m_errs = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\all_ev_plat_m_errs.txt")
all_GEV_plat_ms = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\all_GEV_plat_ms.txt")
all_GEV_plat_m_errs = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\prel_16x16_data\\all_GEV_plat_m_errs.txt")

all_GEV_plat_ms[6,2] = 1.0214751284143868
all_GEV_plat_m_errs[6,2] = 0.036703477787555214

all_ev_plat_ms = all_ev_plat_ms .* betas
all_ev_plat_m_errs = all_ev_plat_m_errs .* betas
all_GEV_plat_ms = all_GEV_plat_ms .* betas
all_GEV_plat_m_errs = all_GEV_plat_m_errs .* betas

# all_ev_plat_ms = all_ev_plat_ms ./ strings
# all_ev_plat_m_errs = all_ev_plat_m_errs ./ strings
# all_GEV_plat_ms = all_GEV_plat_ms ./ strings
# all_GEV_plat_m_errs = all_GEV_plat_m_errs ./ strings

# scatter(betas_inv, all_GEV_plat_ms[:,2], yerror = all_GEV_plat_m_errs[:,2])

model_const(x,p)     = p[1] .+ 0.0 .* x
model_lin(x,p)       = p[1] .+ p[2] .* x
model_para_only(x,p) = p[1] .+ p[2] .* x.^2 
model_para(x,p)      = p[1] .+ p[2] .* x .+ p[3] .* x.^2
model_pol(x,p)       = p[1] .+ p[2] .* x .+ p[3] .* x.^2 .+ p[4] .* x.^3

start_p_const     = [1.0]
start_p_lin       = ones(Float64, 2)
start_p_para_only = ones(Float64, 2)
start_p_para      = ones(Float64, 3)
start_p_pol       = ones(Float64, 4)


betas_inv = 1 ./ betas.^2
betas_inv_plot = 0.0:0.001:0.065
wt_ev = all_ev_plat_m_errs[:].^(-2)

ev_fit_lin = curve_fit(model_lin, betas_inv[1:7], all_ev_plat_ms[1:7], wt_ev[1:7], start_p_lin)
lin_params = ev_fit_lin.param
lin_params_err = stderror(ev_fit_lin)

fit_values_plot_lin = model_lin(betas_inv_plot, lin_params) 
upper_band_lin = model_lin(betas_inv_plot, lin_params.+lin_params_err) .- fit_values_plot_lin
lower_band_lin = fit_values_plot_lin .- model_lin(betas_inv_plot, lin_params.-lin_params_err)
cont_value_lin = model_lin([0.0], lin_params)
cont_value_lin_err_upper = model_lin([0.0], lin_params.+lin_params_err) - cont_value_lin
cont_value_lin_err_lower = cont_value_lin - model_lin([0.0], lin_params.-lin_params_err)
chisqdof_lin = round(sum(ev_fit_lin.resid.^2)/dof(ev_fit_lin), sigdigits = 3)

# ev_fit_para_only = curve_fit(model_para_only, betas_inv, all_ev_plat_ms[:], wt_ev, start_p_para_only)
# para_only_params = ev_fit_para_only.param
# para_only_params_err = stderror(ev_fit_para_only)

# fit_values_plot_para_only = model_para_only(betas_inv_plot, para_only_params) 
# upper_band_para_only = model_para_only(betas_inv_plot, para_only_params.+para_only_params_err) .- fit_values_plot_para_only
# lower_band_para_only = fit_values_plot_para_only .- model_para_only(betas_inv_plot, para_only_params.-para_only_params_err)
# cont_value_para_only = model_para_only([0.0], para_only_params)
# cont_value_para_only_err_upper = model_para_only([0.0], para_only_params.+para_only_params_err)
# cont_value_para_only_err_lower = model_para_only([0.0], para_only_params.-para_only_params_err)
# chisqdof_para_only = round(sum(ev_fit_para_only.resid.^2)/dof(ev_fit_para_only), sigdigits = 3)


# ev_fit_para = curve_fit(model_para, betas_inv, all_ev_plat_ms[:], wt_ev, start_p_para)
# para_params = ev_fit_para.param
# para_params_err = stderror(ev_fit_para)

# fit_values_plot_para = model_para(betas_inv_plot, para_params) 
# upper_band_para = model_para(betas_inv_plot, para_params.+para_params_err) .- fit_values_plot_para
# lower_band_para = fit_values_plot_para .- model_para(betas_inv_plot, para_params.-para_params_err)
# cont_value_para = model_para([0.0], para_params)
# cont_value_para_err_upper = model_para([0.0], para_params.+para_params_err)
# cont_value_para_err_lower = model_para([0.0], para_params.-para_params_err)
# chisqdof_para = round(sum(ev_fit_para.resid.^2)/dof(ev_fit_para), sigdigits = 3)






let
    image_ev_plat_ms = scatter(
        betas_inv[1:7],
        all_ev_plat_ms[1:7],
        yerror = all_ev_plat_m_errs[1:7],
        color = cb_blue,
        markerstrokecolor = cb_blue,
        # label = latexstring("\$m\\beta\$ from EVs"),
        label = latexstring("EVs nr. 1: \$m\\beta\$ data"),
        # ylim = (0.0,4.0)
    )
    image_ev_plat_ms = scatter!(
        [betas_inv[8]],
        [all_ev_plat_ms[8]],
        yerror = [all_ev_plat_m_errs[8]],
        color = :white,
        # markerstrokecolor = cb_blue,
        label = :false, # latexstring("\$m\\beta\$ from EVs"),
        # ylim = (0.0,4.0)
    )
    image_ev_plat_ms = plot!(
        betas_inv_plot,
        fit_values_plot_lin,
        lw = 2,
        color = cb_blue,
        ribbon = (upper_band_lin, lower_band_lin),
        tickfontsize = 10,
        labelfontsize = 15,
        legendfontsize = 11,
        markersize = 5,
        xlabel = latexstring("\$1/\\beta^2 = a^2g^4/16 \$"),
        # label = latexstring("fit, \$\\chi^2/\\mathrm{dof} = $chisqdof_lin\$"),
        label = "EVs nr. 1: lin. fit",
    )
    display(image_ev_plat_ms)
    println("EV nr. 1: m_eff = $cont_value_lin + $cont_value_lin_err_upper - $cont_value_lin_err_lower")
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\m_EV_cont.pdf")  

# let
#     image_ev_plat_ms = plot(
#         betas_inv_plot,
#         fit_values_plot_para_only,
#         lw = 2,
#         color = cb_orange,
#         ribbon = (upper_band_para_only, lower_band_para_only),
#         tickfontsize = 10,
#         labelfontsize = 15,
#         legendfontsize = 11,
#         markersize = 5,
#         xlabel = latexstring("\$1/\\beta^2 = a^2g^4/16 \$"),
#         label = latexstring("fit, \$\\chi^2/\\mathrm{dof} = $chisqdof_para_only\$")
#     )
#     image_ev_plat_ms = scatter!(
#         betas_inv,
#         all_ev_plat_ms,
#         yerror = all_ev_plat_m_errs,
#         color = cb_blue,
#         label = latexstring("\$m\\beta\$ from EVs"),
#         # ylim = (0.0,4.0)
#     )
# end

# let
#     image_ev_plat_ms = plot(
#         betas_inv_plot,
#         fit_values_plot_para,
#         lw = 2,
#         color = cb_orange,
#         ribbon = (upper_band_para, lower_band_para),
#         tickfontsize = 10,
#         labelfontsize = 15,
#         legendfontsize = 11,
#         markersize = 5,
#         xlabel = latexstring("\$1/\\beta^2 = a^2g^4/16 \$"),
#         label = latexstring("fit, \$\\chi^2/\\mathrm{dof} = $chisqdof_para\$")
#     )
#     image_ev_plat_ms = scatter!(
#         betas_inv,
#         all_ev_plat_ms,
#         yerror = all_ev_plat_m_errs,
#         color = cb_blue,
#         label = latexstring("\$m\\beta\$ from EVs"),
#         # ylim = (0.0,4.0)
#     )
# end



mask_1 = 1:8
not_mask_1 = Vector(1:8)[(x-> !(x in (Vector(mask_1)))).(Vector(1:8))]

wt_GEV_1 = all_GEV_plat_m_errs[mask_1,1].^(-2)

# GEV_fit_lin_1 = curve_fit(model_lin, betas_inv[mask_1], all_GEV_plat_ms[mask_1,1], wt_GEV_1, start_p_lin)
# lin_params_1 = GEV_fit_lin_1.param
# lin_params_err_1 = stderror(GEV_fit_lin_1)

# fit_values_plot_lin_1 = model_lin(betas_inv_plot, lin_params_1) 
# upper_band_lin_1 = model_lin(betas_inv_plot, lin_params_1.+lin_params_err_1) .- fit_values_plot_lin_1
# lower_band_lin_1 = fit_values_plot_lin_1 .- model_lin(betas_inv_plot, lin_params_1.-lin_params_err_1)
# cont_value_lin_1 = model_lin([0.0], lin_params_1)
# cont_value_lin_err_upper_1 = model_lin([0.0], lin_params_1.+lin_params_err_1) - cont_value_lin_1
# cont_value_lin_err_lower_1 = cont_value_lin_1 - model_lin([0.0], lin_params_1.-lin_params_err_1)
# chisqdof_lin_1 = round(sum(GEV_fit_lin_1.resid.^2)/dof(GEV_fit_lin_1), sigdigits = 3)

wt_GEV_1 = all_GEV_plat_m_errs[mask_1,1].^(-2)

GEV_fit_const_1 = curve_fit(model_const, betas_inv[mask_1], all_GEV_plat_ms[mask_1,1], wt_GEV_1, start_p_const)
const_params_1 = GEV_fit_const_1.param
const_params_err_1 = stderror(GEV_fit_const_1)

fit_values_plot_const_1 = model_const(betas_inv_plot, const_params_1) 
upper_band_const_1 = model_const(betas_inv_plot, const_params_1.+const_params_err_1) .- fit_values_plot_const_1
lower_band_const_1 = fit_values_plot_const_1 .- model_const(betas_inv_plot, const_params_1.-const_params_err_1)
cont_value_const_1 = model_const([0.0], const_params_1)
cont_value_const_err_upper_1 = model_const([0.0], const_params_1.+const_params_err_1) - cont_value_const_1
cont_value_const_err_lower_1 = cont_value_const_1 - model_const([0.0], const_params_1.-const_params_err_1)
chisqdof_const_1 = round(sum(GEV_fit_const_1.resid.^2)/dof(GEV_fit_const_1), sigdigits = 3)



mask_2 = 2:7
not_mask_2 = Vector(1:8)[(x-> !(x in (Vector(mask_2)))).(Vector(1:8))]

wt_GEV_2 = all_GEV_plat_m_errs[mask_2,2].^(-2)

GEV_fit_lin_2 = curve_fit(model_lin, betas_inv[mask_2], all_GEV_plat_ms[mask_2,2], wt_GEV_2, start_p_lin)
lin_params_2 = GEV_fit_lin_2.param
lin_params_err_2 = stderror(GEV_fit_lin_2)

fit_values_plot_lin_2 = model_lin(betas_inv_plot, lin_params_2) 
upper_band_lin_2 = model_lin(betas_inv_plot, lin_params_2.+lin_params_err_2) .- fit_values_plot_lin_2
lower_band_lin_2 = fit_values_plot_lin_2 .- model_lin(betas_inv_plot, lin_params_2.-lin_params_err_2)
cont_value_lin_2 = model_lin([0.0], lin_params_2)
cont_value_lin_err_upper_2 = model_lin([0.0], lin_params_2.+lin_params_err_2) - cont_value_lin_2
cont_value_lin_err_lower_2 = cont_value_lin_2 - model_lin([0.0], lin_params_2.-lin_params_err_2)
chisqdof_lin_2 = round(sum(GEV_fit_lin_2.resid.^2)/dof(GEV_fit_lin_2), sigdigits = 3)

let
    image_GEV_plat_ms = scatter(
        betas_inv[mask_1],
        all_GEV_plat_ms[mask_1,1],
        yerror = all_GEV_plat_m_errs[mask_1,1],
        color = cb_blue,
        markerstrokecolor = cb_blue,
        label = latexstring("GEVs nr. 1: \$m\\beta\$ data"),
        tickfontsize = 10,
        labelfontsize = 15,
        legendfontsize = 11,
        markersize = 5,
        xlabel = latexstring("\$1/\\beta^2 = a^2g^4/16 \$"),
        # ylabel = latexstring("\$m\\beta\$"),
        markershape = :circ,
        # ylim = (0.0,4.0)
    )
    # image_GEV_plat_ms = scatter!(
    #     label = :false,
    #     betas_inv[not_mask_1],
    #     all_GEV_plat_ms[not_mask_1,1],
    #     yerror = all_GEV_plat_m_errs[not_mask_1,1],
    #     markershape = :circ,
    #     color = :white,
    #     markersize = 5,
    #     # alpha = 0.7
    # )
    image_GEV_plat_ms = plot!(
        betas_inv_plot,
        fit_values_plot_const_1,
        lw = 2,
        color = cb_blue,
        ribbon = (lower_band_const_1, upper_band_const_1),
        label = "GEVs nr. 1: const. fit",
    )
    image_GEV_plat_ms = scatter!(
        betas_inv[mask_2],
        all_GEV_plat_ms[mask_2,2],
        yerror = all_GEV_plat_m_errs[mask_2,2],
        color = cb_orange,
        markerstrokecolor = cb_orange,
        label = latexstring("GEVs nr. 2: \$m\\beta\$ data"),
        markershape = :diamond,
        markersize = 5,
        # ylim = (0.0,4.0)
    )
    image_GEV_plat_ms = plot!(
        betas_inv_plot,
        fit_values_plot_lin_2,
        lw = 2,
        color = cb_orange,
        ribbon = (lower_band_lin_2, upper_band_lin_2),
        label = "GEVs nr. 2: lin. fit",
        linestyle = :dash
    )
    not_mask_2 = [7]
    image_GEV_plat_ms = scatter!(
        label = :false,
        betas_inv[not_mask_2],
        all_GEV_plat_ms[not_mask_2,2],
        yerror = all_GEV_plat_m_errs[not_mask_2,2],
        markershape = :diamond,
        # color = :white,
        color = cb_orange,
        markerstrokecolor = cb_orange,
        markersize = 5,
        # alpha = 0.7
    )
    display(image_GEV_plat_ms)
    println("GEV nr. 1: m_eff = $cont_value_const_1 + $cont_value_const_err_upper_1 - $cont_value_const_err_lower_1")
    println("GEV nr. 2: m_eff = $cont_value_lin_2 + $cont_value_lin_err_upper_2 - $cont_value_lin_err_lower_2")
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\plots\\m_GEV_cont.pdf")  