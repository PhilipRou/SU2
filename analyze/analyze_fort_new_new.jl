include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")



# for beta in [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0] #6.0:1.5:15.0
    ### beta = 8.0
    beta = 20.0
    println("\n\n\n")
    println("\t\t\t Started analysis at beta = $beta")
    println("\n\n\n")
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
    mkpath(base_path)
    base_path = base_path * "\\beta_$(beta)_Nz_$(Nz)_Nx_$(Nx)_smearlist_rho_2D_$(rho_2D)"



    @assert 1==0 "plot range 1:... and TITLE weg"

    



    
    m_evs_dr      = blg[1]
    m_evs_errs_dr = blg[2]
    evs_dr        = blg[3]
    evs_errs_dr   = blg[4]
    m_plat_dr     = blg[5]
    m_plat_err_dr = blg[6]

    plot_range_m_ev_dr = 1:5

    let
        image_conn_evs = plot(
            title  = latexstring("EV's of conn. corr. 8² mats. (bottom right) \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:4 # size(evs,2)
            image_conn_evs = scatter!(
                Vector(0:Nz>>1) .+ (i-1)*0.05,
                evs_dr[:,i],
                yerror = evs_errs_dr[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-16, 10*maximum(evs_dr)),
                yaxis = :log,
                legend = :topright
            )
        end
        display(image_conn_evs)
        # fig_path = base_path * "_evs_dr_8x8.pdf"
        # savefig(fig_path)
    end
    let
        image_conn_evs_masses = plot(
            title  = latexstring("2pt Masses of EV's of conn. corr. 8² mats. (bottom right) \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2#size(m_evs_dr,2)
            image_conn_evs_masses = scatter!(
                (Vector(0:Nz>>1) .+ 0.5 .+ (i-1)*0.05)[plot_range_m_GEV_dr],
                m_evs_dr[plot_range_m_GEV_dr,i],
                yerror = m_evs_errs_dr[plot_range_m_GEV_dr,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :topright,
                # rightmargin = 5mm
            )
        end
        for i = 1:2
            image_conn_evs_masses = plot!(
                (Vector(0:Nz>>1) .+ 0.5 .+ (i-1)*0.05)[plot_range_m_GEV_dr],
                m_GEVs_plat_dr[i] .* ones(length(plot_range_m_GEV_dr)),
                ribbon = m_GEVs_plat_errs_dr[i] .* ones(length(plot_range_m_GEV_dr)),
                color = palette(:default)[i],
                lw = 2
            )
        end
        display(image_conn_evs_masses)
        # fig_path = base_path * "_m_evs_dr_8x8.pdf"
        # savefig(fig_path)
    end






    m_GEVs_dr           = blh[1]
    m_GEVs_errs_dr      = blh[2]
    GEVs_dr             = blh[3]
    GEVs_errs_dr        = blh[4]
    m_GEVs_plat_dr      = blh[5]
    m_GEVs_plat_errs_dr = blh[6]

    let
        image_conn_GEVs = plot(
            title  = latexstring("GEV's of conn. corr. 8² mats. (bottom right) \n \$\\beta = $beta\$, \$L = $Nz\$, \$t_0 = $(t0-1)\$, \$n_\\mathrm{smear} \\in\$ $(smearlist[2:end])"),
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
                legend = :topright
            )
        end
        display(image_conn_GEVs)
        # fig_path = base_path * "_GEVs_dr_8x8.pdf"
        # savefig(fig_path)
    end

    let
        image_conn_GEVs_masses = plot(
            title  = latexstring("2pt Masses of GEV's of conn. corr. 8² mats. (bottom right)\n \$\\beta = $beta, L = $Nz, t_0 = $(t0-1), n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2 # size(m_evs_dr,2)
            image_conn_GEVs_masses = scatter!(
                (Vector(t0:Nz>>1) .+ 0.5 .+ (i-1)*0.05)[1:6],
                m_GEVs_dr[1:6,i],
                yerror = m_GEVs_errs_dr[1:6,i],
                label = latexstring("GEV-mass nr. \$$i\$"),
                color = cb_colors[i],
                markerstrokecolor = cb_colors[i],
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :bottomleft,
                markershape = [:circ, :diamond][i],
            )
        end
        image_conn_evs_masses = plot!(
            (Vector(t0:Nz>>1) .+ 0.5)[3:6],
            m_GEVs_plat_dr[1] .* ones(length(3:6)),
            ribbon = m_GEVs_plat_errs_dr[1] .* ones(length(3:6)),
            color = cb_colors[1],
            lw = 2,
            label = "plateau nr. 1",
        )
        image_conn_evs_masses = plot!(
            (Vector(t0:Nz>>1) .+ 0.55)[1:2],
            m_GEVs_plat_dr[2] .* ones(length(1:2)),
            ribbon = m_GEVs_plat_errs_dr[2] .* ones(length(1:2)),
            color = cb_colors[2],
            lw = 2,
            label = "plateau nr. 2",
            linestyle = :dash,
        )
        display(image_conn_GEVs_masses)
        # fig_path = base_path * "_m_GEVs_dr_8x8.pdf"
        # savefig(fig_path)
    end
