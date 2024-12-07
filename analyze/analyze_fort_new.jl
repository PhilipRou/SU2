include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")



for beta in [4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 16.0] #6.0:1.5:15.0
    ### beta = 8.0
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

    
    ### s_wil    s_opt    smeared_s_wil  smeared_clover  accrate_metro  accrate_ovlax  time
    ops_smeared = readdlm(timeseries_2Dsmear_path, skipstart=1)
    n_op = size(ops_smeared,2)
    bsize_ops_smeared = round(Int, 2*maximum([auto_corr_time(ops_smeared[:,j]) for j = 1:n_op])) + 1

    raw_corr  = readdlm(corr_path, skipstart=1)
    n_meas    = minimum([size(raw_corr,1), size(readdlm(timeseries_path, skipstart=1),1)])
    corrmats  = []
    for t = 1:Nz
        corrmats_t = []
        for meas = 1:n_meas
            push!(corrmats_t, reshape(raw_corr[Nz*(meas-1)+t,:], n_op, n_op))
        end
        push!(corrmats, corrmats_t)
    end
    autocorr_corrmats = []
    for t = 1:Nz
    for op1 = 1:n_op
    for op2 = op1:n_op
        push!(autocorr_corrmats, auto_corr_time([corrmats[t][meas][op1,op2] for meas = 1:500]) )
    end
    end
    end
    bsize_corrmats = round(Int, 2*maximum(autocorr_corrmats)) + 1

    ### Symmetrize corr. matrices in t space
    corrmats_symm = []
    push!(corrmats_symm, corrmats[Nz][:]);
    for t = 1:Nz>>1-1
        push!(corrmats_symm, (corrmats[t] .+ corrmats[Nz-t]) ./ 2) 
        # push!(corrmats_symm, corrmats[t][:])  ### for debugging purposes 🚧
    end
    push!(corrmats_symm, corrmats[Nz>>1][:]);
    
    ### Symmetrize every corr. matrix (in operator space)
    println("Symmetrize corr. matrices")
    @time for t = 1:Nz>>1+1
        @show t
        corrmats_symm[t] = (corrmats_symm[t] .+ transpose.(corrmats_symm[t])) ./ 2
        @assert !(false in ishermitian.(corrmats_symm[t])) "There are non-hermitian correlation matrices for t = $t"
    end

    corrmats_symm_8x8_upleft = []
    corrmats_symm_8x8_downright = []
    if n_op == 16
        for t = 1:Nz>>1+1
            push!(corrmats_symm_8x8_upleft, [corrmats_symm[t][meas][1:8,1:8] for meas = 1:n_meas])
            push!(corrmats_symm_8x8_downright, [corrmats_symm[t][meas][9:16,9:16] for meas = 1:n_meas])
        end
    end



####################################################################################################################################



    ### Connected correlators of individual smeared s_wil series
    corr_con_mean = Array{Float64}(undef, length(smearlist)-1, Nz>>1+1);
    corr_con_err = Array{Float64}(undef, length(smearlist)-1, Nz>>1+1);
    for smear_ind = 1:length(smearlist)-1
        n_smear = smearlist[smear_ind+1]
        for t = 1:Nz>>1+1
            temp_corrs = [corrmats_symm[t][meas][smear_ind,smear_ind] for meas = 1:n_meas]
            temp_means = ops_smeared[:,smear_ind]
            corr_con_mean[smear_ind,t], corr_con_err[smear_ind,t] = jack_conn_corr_self(temp_corrs, temp_means, bsize_corrmats)
        end # t
    end # n_smear
    let
        image_con = plot(
            title  = latexstring("Connected \$s_\\mathrm{wil}\$ corr. \n \$β = $beta, L = $Nz, ρ = $rho_2D\$"),
            ylims  = (1e-12, 1e-5),
            yaxis  = :log,
            legend = :topright,
            xlabel = latexstring("\$t\$")
        )
        for smear_ind = 1:length(smearlist)-1
            image_con = scatter!(
                Vector(0:Nz>>1) .+ ((smear_ind-1)*0.1), 
                corr_con_mean[smear_ind,:], 
                yerror = corr_con_err[smear_ind,:], 
                label = latexstring("\$n_\\mathrm{smear} = $(smearlist[smear_ind+1])\$"),
                markerstrokecolor = :auto
            )
        end
        display(image_con)
        fig_path = base_path * "_swil_correlator.pdf"
        savefig(fig_path)
    end
# end # beta



    ### EVs and masses from connected correlator: 8×8, upper left
    
    println("Start jackknife of 8x8 upper left evs")
    @time bla = jack_conn_corr_mat_ev_mass_2pt_allofem(corrmats_symm_8x8_upleft, ops_smeared[:,1:8], 2, 3)
    
    m_evs_ul      = bla[1]
    m_evs_errs_ul = bla[2]
    evs_ul        = bla[3]
    evs_errs_ul   = bla[4]

    let
        image_conn_evs = plot(
            title  = latexstring("EV's of conn. corr. 8² mats. (upper left) \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:4 # size(evs,2)
            image_conn_evs = scatter!(
                Vector(0:Nz>>1) .+ (i-1)*0.05,
                evs_ul[:,i],
                yerror = evs_errs_ul[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-16, 10*maximum(evs_ul)),
                yaxis = :log,
                legend = :topright
            )
        end
        display(image_conn_evs)
        fig_path = base_path * "_evs_ul_8x8.pdf"
        savefig(fig_path)
    end
    let
        image_conn_evs_masses = plot(
            title  = latexstring("2pt Masses of EV's of conn. corr. 8² mats. (upper left) \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2#size(m_evs_ul,2)
            image_conn_evs_masses = scatter!(
                Vector(0:Nz>>1) .+ 0.5 .+ (i-1)*0.05,
                m_evs_ul[:,i],
                yerror = m_evs_errs_ul[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :topright,
                rightmargin = 5mm
            )
        end
        display(image_conn_evs_masses)
        fig_path = base_path * "_m_evs_ul_8x8.pdf"
        savefig(fig_path)
    end

    let
    #     evs_path = base_path * "_evs_8x8_ul.txt"
    #     blu = open(evs_path, "w")
    #     writedlm(blu, evs_ul)
    #     close(blu)
    #     evs_errs_path = base_path * "_evs_errs_8x8_ul.txt"
    #     bli = open(evs_errs_path, "w")
    #     writedlm(bli, evs_errs_ul)
    #     close(bli)
    # end
    # let
    #     m_evs_path = base_path * "_m_evs_8x8_ul.txt"
    #     blu = open(m_evs_path, "w")
    #     writedlm(blu, m_evs_ul)
    #     close(blu)
    #     m_evs_errs_path = base_path * "_m_evs_errs_8x8_ul.txt"
    #     bli = open(m_evs_errs_path, "w")
    #     writedlm(bli, m_evs_errs_ul)
    #     close(bli)
    end



    ### EVs and masses from connected correlator: 8×8, down right
    
    println("Start jackknife of 8x8 bottom right evs")
    @time blb = jack_conn_corr_mat_ev_mass_2pt_allofem(corrmats_symm_8x8_downright, ops_smeared[:,9:16], 2, 3)
    
    m_evs_dr      = blb[1]
    m_evs_errs_dr = blb[2]
    evs_dr        = blb[3]
    evs_errs_dr   = blb[4]

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
        fig_path = base_path * "_evs_dr_8x8.pdf"
        savefig(fig_path)
    end
    let
        image_conn_evs_masses = plot(
            title  = latexstring("2pt Masses of EV's of conn. corr. 8² mats. (bottom right) \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2#size(m_evs_dr,2)
            image_conn_evs_masses = scatter!(
                Vector(0:Nz>>1) .+ 0.5 .+ (i-1)*0.05,
                m_evs_dr[:,i],
                yerror = m_evs_errs_dr[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :topright,
                # rightmargin = 5mm
            )
        end
        display(image_conn_evs_masses)
        fig_path = base_path * "_m_evs_dr_8x8.pdf"
        savefig(fig_path)
    end

    let
        evs_path = base_path * "_evs_8x8_dr.txt"
        blu = open(evs_path, "w")
        writedlm(blu, evs_dr)
        close(blu)
        evs_errs_path = base_path * "_evs_errs_8x8_dr.txt"
        bli = open(evs_errs_path, "w")
        writedlm(bli, evs_errs_dr)
        close(bli)
    end
    let
        m_evs_path = base_path * "_m_evs_8x8_dr.txt"
        blu = open(m_evs_path, "w")
        writedlm(blu, m_evs_dr)
        close(blu)
        m_evs_errs_path = base_path * "_m_evs_errs_8x8_dr.txt"
        bli = open(m_evs_errs_path, "w")
        writedlm(bli, m_evs_errs_dr)
        close(bli)
    end



    ### EVs and masses from connected correlator: 16×16
    
    println("Start jackknife of 16x16 evs")
    @time blc = jack_conn_corr_mat_ev_mass_2pt_allofem(corrmats_symm, ops_smeared, 2, 3)
    
    m_evs      = blc[1]
    m_evs_errs = blc[2]
    evs        = blc[3]
    evs_errs   = blc[4]

    let
        image_conn_evs = plot(
            title  = latexstring("EV's of conn. corr. 16² mats. \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:4 # size(evs,2)
            image_conn_evs = scatter!(
                Vector(0:Nz>>1) .+ (i-1)*0.05,
                evs[:,i],
                yerror = evs_errs[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-16, 10*maximum(evs)),
                yaxis = :log,
                legend = :topright,
                rightmargin = 5mm
            )
        end
        display(image_conn_evs)
        fig_path = base_path * "_evs.pdf"
        savefig(fig_path)
    end
    let
        image_conn_evs_masses = plot(
            title  = latexstring("2pt Masses of EV's of conn. corr. 16² mats. \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2#size(m_evs,2)
            image_conn_evs_masses = scatter!(
                Vector(0:Nz>>1) .+ 0.5 .+ (i-1)*0.05,
                m_evs[:,i],
                yerror = m_evs_errs[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :topright,
                rightmargin = 5mm
            )
        end
        display(image_conn_evs_masses)
        fig_path = base_path * "_m_evs.pdf"
        savefig(fig_path)
    end

    let
        evs_path = base_path * "_evs.txt"
        blu = open(evs_path, "w")
        writedlm(blu, evs)
        close(blu)
        evs_errs_path = base_path * "_evs_errs.txt"
        bli = open(evs_errs_path, "w")
        writedlm(bli, evs_errs)
        close(bli)
    end
    let
        m_evs_path = base_path * "_m_evs.txt"
        blu = open(m_evs_path, "w")
        writedlm(blu, m_evs)
        close(blu)
        m_evs_errs_path = base_path * "_m_evs_errs.txt"
        bli = open(m_evs_errs_path, "w")
        writedlm(bli, m_evs_errs)
        close(bli)
    end


    
####################################################################################################################################



    ### GEVs of the correlation matrices
    ### t0 is the index of the lower bound for the time-index of the GEVP, starting at 1!
    #=
    t0 = 1
    GEVs      = Array{Float64}(undef, Nz>>1+1-t0, n_op)
    GEVs_errs = Array{Float64}(undef, Nz>>1+1-t0, n_op)
        # GEVs      = Array{Float64}(undef, Nz>>1-t0, 4)   ### for debugging purposes 🚧
        # GEVs_errs = Array{Float64}(undef, Nz>>1-t0, 4)   ### for debugging purposes 🚧
        # GEVs      = Array{Float64}(undef, Nz>>1-t0, 8)   ### for debugging purposes 🚧
        # GEVs_errs = Array{Float64}(undef, Nz>>1-t0, 8)   ### for debugging purposes 🚧
    @time for t = t0+1:Nz>>1+1
        @show t
        # t = 3
        jack = jack_conn_corr_mat_GEV(corrmats_symm[t0], corrmats_symm[t], ops_smeared, maximum([bsize_corrmats,bsize_ops_smeared]))
        ### for debugging purposes 🚧:
            # jack = jack_conn_corr_mat_GEV([corrmats_symm[t][meas][1:4,1:4] for meas = 1:n_meas], [corrmats_symm[t0][meas][1:4,1:4] for meas = 1:n_meas], ops_smeared[:,1:4], maximum([bsize_corrmats,bsize_ops_smeared]))
            # jack = jack_conn_corr_mat_GEV([corrmats_symm[t][meas][1:8,1:8] for meas = 1:n_meas], [corrmats_symm[t0][meas][1:8,1:8] for meas = 1:n_meas], ops_smeared[:,1:8], maximum([bsize_corrmats,bsize_ops_smeared]))
            # jack = jack_conn_corr_mat_GEV([corrmats_symm[t][meas][9:16,9:16] for meas = 1:n_meas], [corrmats_symm[t0][meas][9:16,9:16] for meas = 1:n_meas], ops_smeared[:,9:16], maximum([bsize_corrmats,bsize_ops_smeared]))
        GEVs[t-t0,:]  = jack[1]
        GEVs_errs[t-t0,:] = jack[2]
    end
    =#



    ### GEVs and masses of connected correlation matrices: 8×8, upper left

    t0 = 1
    println("Start jackknife of 8x8 upper left GEVs")
    @time bld = jack_conn_corr_mat_GEV_mass_2pt_allofem(t0, corrmats_symm_8x8_upleft, ops_smeared[:,1:8], maximum([bsize_corrmats,bsize_ops_smeared]), 3)

    m_GEVs_ul      = bld[1]
    m_GEVs_errs_ul = bld[2]
    GEVs_ul        = bld[3]
    GEVs_errs_ul   = bld[4]

    let
        image_conn_GEVs = plot(
            title  = latexstring("GEV's of conn. corr. 8² mats. (upper left) \n \$\\beta = $beta\$, \$L = $Nz\$, \$t_0 = $(t0-1)\$, \$n_\\mathrm{smear} \\in\$ $(smearlist[2:end])"),
            xlabel = latexstring("\$ t\$"),
            )
        for i = 1:4 # size(GEVs,2)
            image_conn_GEVs = scatter!(
                Vector(t0+1:Nz>>1+1) .+ (i-1)*0.02,
                GEVs_ul[:,i],
                yerror = GEVs_errs_ul[:,i],
                label = latexstring("GEV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-4, 10*maximum(GEVs_ul)),
                yaxis = :log,
                legend = :topright
            )
        end
        display(image_conn_GEVs)
        fig_path = base_path * "_GEVs_ul_8x8.pdf"
        savefig(fig_path)
    end
    let
        image_conn_GEVs_masses = plot(
            title  = latexstring("2pt Masses of GEV's of conn. corr. 8² mats. (upper left)\n \$\\beta = $beta, L = $Nz, t_0 = $(t0-1), n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2#size(m_evs_ul,2)
            image_conn_GEVs_masses = scatter!(
                Vector(t0:Nz>>1) .+ 0.5 .+ (i-1)*0.05,
                m_GEVs_ul[:,i],
                yerror = m_GEVs_errs_ul[:,i],
                label = latexstring("GEV nr. \$$i\$"),
                markerstrokecolor = :auto,
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :topright
            )
        end
        display(image_conn_GEVs_masses)
        fig_path = base_path * "_m_GEVs_ul_8x8.pdf"
        savefig(fig_path)
    end

    let
        GEVs_path = base_path*"_GEVs_t0_$(t0)_8x8_ul.txt"
        blu = open(GEVs_path, "w")
        writedlm(blu, GEVs_ul)
        close(blu)
        GEVs_errs_path = base_path*"_GEVs_errs_t0_$(t0)_8x8_ul.txt"
        bli = open(GEVs_errs_path, "w")
        writedlm(bli, GEVs_errs_ul)
        close(bli)
    end
    let
        m_GEVs_path = base_path * "_m_GEVs_t0_$(t0)_8x8_ul.txt"
        blu = open(m_GEVs_path, "w")
        writedlm(blu, m_GEVs_ul)
        close(blu)
        m_GEVs_errs_path = base_path * "_m_GEVs_errs_t0_$(t0)_8x8_ul.txt"
        bli = open(m_GEVs_errs_path, "w")
        writedlm(bli, m_GEVs_errs_ul)
        close(bli)
    end



    ### GEVs and masses of connected correlation matrices: 8×8, down right

    println("Start jackknife of 8x8 bottom right GEVs")
    @time ble = jack_conn_corr_mat_GEV_mass_2pt_allofem(1, corrmats_symm_8x8_downright, ops_smeared[:,9:16], maximum([bsize_corrmats,bsize_ops_smeared]), 3)

    m_GEVs_dr      = ble[1]
    m_GEVs_errs_dr = ble[2]
    GEVs_dr        = ble[3]
    GEVs_errs_dr   = ble[4]

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
        fig_path = base_path * "_GEVs_dr_8x8.pdf"
        savefig(fig_path)
    end
    let
        image_conn_GEVs_masses = plot(
            title  = latexstring("2pt Masses of GEV's of conn. corr. 8² mats. (bottom right)\n \$\\beta = $beta, L = $Nz, t_0 = $(t0-1), n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2 # size(m_evs_dr,2)
            image_conn_GEVs_masses = scatter!(
                Vector(t0:Nz>>1) .+ 0.5 .+ (i-1)*0.05,
                m_GEVs_dr[:,i],
                yerror = m_GEVs_errs_dr[:,i],
                label = latexstring("GEV nr. \$$i\$"),
                markerstrokecolor = :auto,
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :topright
            )
        end
        display(image_conn_GEVs_masses)
        fig_path = base_path * "_m_GEVs_dr_8x8.pdf"
        savefig(fig_path)
    end

    let
        GEVs_path = base_path*"_GEVs_t0_$(t0)_8x8_dr.txt"
        blu = open(GEVs_path, "w")
        writedlm(blu, GEVs_dr)
        close(blu)
        GEVs_errs_path = base_path*"_GEVs_errs_t0_$(t0)_8x8_dr.txt"
        bli = open(GEVs_errs_path, "w")
        writedlm(bli, GEVs_errs_dr)
        close(bli)
    end
    let
        m_GEVs_path = base_path * "_m_GEVs_t0_$(t0)_8x8_dr.txt"
        blu = open(m_GEVs_path, "w")
        writedlm(blu, m_GEVs_dr)
        close(blu)
        m_GEVs_errs_path = base_path * "_m_GEVs_errs_t0_$(t0)_8x8_dr.txt"
        bli = open(m_GEVs_errs_path, "w")
        writedlm(bli, m_GEVs_errs_dr)
        close(bli)
    end



    ### GEVs and masses of full 16x16 connected correlation matrices

    println("Start jackknife of 16x16 GEVs")
    @time blf = jack_conn_corr_mat_GEV_mass_2pt_allofem(1, corrmats_symm, ops_smeared, maximum([bsize_corrmats,bsize_ops_smeared]), 3)

    m_GEVs      = blf[1]
    m_GEVs_errs = blf[2]
    GEVs        = blf[3]
    GEVs_errs   = blf[4]

    let
        image_conn_GEVs = plot(
            title  = latexstring("GEV's of conn. corr. 16² mats. \n \$\\beta = $beta\$, \$L = $Nz\$, \$t_0 = $(t0-1)\$, \$n_\\mathrm{smear} \\in\$ $(smearlist[2:end])"),
            xlabel = latexstring("\$ t\$"),
            )
        for i = 1:4 # size(GEVs,2)
            image_conn_GEVs = scatter!(
                Vector(t0+1:Nz>>1+1) .+ (i-1)*0.02,
                GEVs[:,i],
                yerror = GEVs_errs[:,i],
                label = latexstring("GEV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-4, 10*maximum(GEVs)),
                yaxis = :log,
                legend = :topright
            )
        end
        display(image_conn_GEVs)
        fig_path = base_path * "_GEVs.pdf"
        savefig(fig_path)
    end
    let
        image_conn_GEVs_masses = plot(
            title  = latexstring("2pt Masses of GEV's of conn. corr. 16² mats. \n \$\\beta = $beta, L = $Nz, t_0 = $(t0-1), n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:2#size(m_evs,2)
            image_conn_GEVs_masses = scatter!(
                Vector(t0:Nz>>1) .+ 0.5 .+ (i-1)*0.05,
                m_GEVs[:,i],
                yerror = m_GEVs_errs[:,i],
                label = latexstring("GEV nr. \$$i\$"),
                markerstrokecolor = :auto,
                # ylim = (1e-16, 1e-3),
                # yaxis = :log,
                legend = :topright
            )
        end
        display(image_conn_GEVs_masses)
        fig_path = base_path * "_m_GEVs.pdf"
        savefig(fig_path)
    end

    let
        GEVs_path = base_path*"_GEVs_t0_$(t0).txt"
        blu = open(GEVs_path, "w")
        writedlm(blu, GEVs)
        close(blu)
        GEVs_errs_path = base_path*"_GEVs_errs_t0_$(t0).txt"
        bli = open(GEVs_errs_path, "w")
        writedlm(bli, GEVs_errs)
        close(bli)
    end
    let
        m_GEVs_path = base_path * "_m_GEVs_t0_$(t0).txt"
        blu = open(m_GEVs_path, "w")
        writedlm(blu, m_GEVs)
        close(blu)
        m_GEVs_errs_path = base_path * "_m_GEVs_errs_t0_$(t0).txt"
        bli = open(m_GEVs_errs_path, "w")
        writedlm(bli, m_GEVs_errs)
        close(bli)
    end



end # beta