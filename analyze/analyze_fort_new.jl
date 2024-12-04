include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")



# for beta in [4.0] # [4.0,5.0,6.0,8.0,10.0,12.0,16.0] #6.0:1.5:15.0
    beta = 8.0
    beta       = string(beta)*"0"
    # beta       = "6.00"
    n_stout    = 7
    Nx=Ny=Nz   = 32
    smearlist  = [0,1,3,7,15]
    rho_2D     = 0.24
    rho_3D     = 0.16

    n_op = 16 # length(smearlist)

    base_path               = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\fortran_projects\\SU2_3D\\beta_$(beta)_Nz_$(Nz)_Nx_$(Nx)_smearlist_rho_2D_$(rho_2D)"
    timeseries_path         = base_path * "_hist_3D.txt"
    timeseries_2Dsmear_path = base_path * "_hist_2D.txt"
    corr_path               = base_path * "_crosscorr.txt"

    ### s_wil    s_opt    smeared_s_wil  smeared_clover  accrate_metro  accrate_ovlax  time
    ops_smeared = readdlm(timeseries_2Dsmear_path, skipstart=1)
    bsize_ops_smeared = round(Int, 2*maximum([auto_corr_time(ops_smeared[:,j]) for j = 1:n_op])) + 1

    raw_corr  = readdlm(corr_path, skipstart=1);
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
    for t = 1:Nz>>1+1
        corrmats_symm[t] = (corrmats_symm[t] .+ transpose.(corrmats_symm[t])) ./ 2
        @assert !(0 in ishermitian.(corrmats_symm[t])) "There are non-hermitian correlation matrices for t = $t"
    end

    round.(mean(corrmats_symm[1]), sigdigits = 4)
    # eigen(mean(corrmats_symm[1])[1:4,1:4]).values
    # bla = readdlm("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\fortran_projects\\SU2_3D_data\\stephans_mat.txt")
    # eigen(bla[1:4,1:4]).values

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
            legend = :bottomleft,
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
    end
# end # beta


    
    ### EVs of correlation matrices
    # evs      = Array{Float64}(undef, Nz>>1+1, n_op)
    # evs_errs = Array{Float64}(undef, Nz>>1+1, n_op)
        # evs      = Array{Float64}(undef, Nz>>1+1, 4)    ### for debugging purposes 🚧
        # evs_errs = Array{Float64}(undef, Nz>>1+1, 4)    ### for debugging purposes 🚧
        evs      = Array{Float64}(undef, Nz>>1+1, 8)    ### for debugging purposes 🚧
        evs_errs = Array{Float64}(undef, Nz>>1+1, 8)    ### for debugging purposes 🚧
    @time for t = 1:Nz>>1+1
        @show t
        # jack = jack_conn_corr_mat_ev(corrmats_symm[t], ops_smeared, maximum([bsize_corrmats,bsize_ops_smeared]))
        ### for debugging purposes 🚧:
            # jack = jack_conn_corr_mat_ev([corrmats_symm[t][meas][1:4,1:4] for meas = 1:n_meas], ops_smeared[:,1:4], maximum([bsize_corrmats,bsize_ops_smeared]))
            # jack = jack_conn_corr_mat_ev([corrmats_symm[t][meas][1:8,1:8] for meas = 1:n_meas], ops_smeared[:,1:8], maximum([bsize_corrmats,bsize_ops_smeared]))
            jack = jack_conn_corr_mat_ev([corrmats_symm[t][meas][9:16,9:16] for meas = 1:n_meas], ops_smeared[:,9:16], maximum([bsize_corrmats,bsize_ops_smeared]))
        evs[t,:]      = reverse(jack[1])
        evs_errs[t,:] = reverse(jack[2])
    end

    let
        image_conn_evs = plot(
            title  = latexstring("EV's of conn. corr. 16² mats. \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:size(evs,2)
            image_conn_evs = scatter!(
                Vector(0:Nz>>1) .+ (i-1)*0.05,
                evs[:,i],
                yerror = evs_errs[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-16, 1e-3),
                yaxis = :log,
                legend = :outerright
            )
        end
        display(image_conn_evs)
    end

    let
        # evs_path = base_path*"_evs.txt"
        # bla = open(evs_path, "w")
        # writedlm(bla, evs)
        # close(bla)
        # evs_errs_path = base_path*"_evs_errs.txt"
        # bli = open(evs_errs_path, "w")
        # writedlm(bli, evs_errs)
        # close(bli)
    end
# end # beta



    ### GEVs of the correlation matrices
    ### t0 is the index of the lower bound for the time-index of the GEVP, starting at 1!
    t0 = 1
    # GEVs      = Array{Float64}(undef, Nz>>1-t0, n_op)
    # GEVs_errs = Array{Float64}(undef, Nz>>1-t0, n_op)
        # GEVs      = Array{Float64}(undef, Nz>>1-t0, 4)   ### for debugging purposes 🚧
        # GEVs_errs = Array{Float64}(undef, Nz>>1-t0, 4)   ### for debugging purposes 🚧
        GEVs      = Array{Float64}(undef, Nz>>1-t0, 8)   ### for debugging purposes 🚧
        GEVs_errs = Array{Float64}(undef, Nz>>1-t0, 8)   ### for debugging purposes 🚧
    @time for t = t0+1:Nz>>1
        @show t
        # t = 3
        # jack = jack_conn_corr_mat_GEV(corrmats_symm[t], corrmats_symm[t0], ops_smeared, maximum([bsize_corrmats,bsize_ops_smeared]))
        ### for debugging purposes 🚧:
            # jack = jack_conn_corr_mat_GEV([corrmats_symm[t][meas][1:4,1:4] for meas = 1:n_meas], [corrmats_symm[t0][meas][1:4,1:4] for meas = 1:n_meas], ops_smeared[:,1:4], maximum([bsize_corrmats,bsize_ops_smeared]))
            # jack = jack_conn_corr_mat_GEV([corrmats_symm[t][meas][1:8,1:8] for meas = 1:n_meas], [corrmats_symm[t0][meas][1:8,1:8] for meas = 1:n_meas], ops_smeared[:,1:8], maximum([bsize_corrmats,bsize_ops_smeared]))
            jack = jack_conn_corr_mat_GEV([corrmats_symm[t][meas][9:16,9:16] for meas = 1:n_meas], [corrmats_symm[t0][meas][9:16,9:16] for meas = 1:n_meas], ops_smeared[:,9:16], maximum([bsize_corrmats,bsize_ops_smeared]))
            GEVs[t-t0,:]  = reverse(jack[1])
        GEVs_errs[t-t0,:] = reverse(jack[2])
    end

    let
        image_conn_GEVs = plot(
            title  = latexstring("Gen. EV's of conn. corr. 16² mats. \n \$\\beta = $beta\$, \$L = $Nz\$, \$t_0 = $t0\$, \$n_\\mathrm{smear} \\in\$ $(smearlist[2:end])"),
            xlabel = latexstring("\$ t\$"),
            )
        for i = 1:8 # size(GEVs,2)
            image_conn_GEVs = scatter!(
                Vector(t0+1:Nz>>1) .+ (i-1)*0.02,
                GEVs[:,i],
                yerror = GEVs_errs[:,i],
                label = latexstring("GEV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-4, 1.0),
                yaxis = :log,
                legend = :outerright
            )
        end
        display(image_conn_GEVs)
    end

    let
        # GEVs_path = base_path*"_GEVs.txt"
        # bla = open(GEVs_path, "w")
        # writedlm(bla, GEVs)
        # close(bla)
        # GEVs_errs_path = base_path*"_GEVs_errs.txt"
        # bli = open(GEVs_errs_path, "w")
        # writedlm(bli, GEVs_errs)
        # close(bli)
    end

# end # beta

# [format_x_err(evs[t,i], evs_errs[t,i], 2) for t = 1:7, i = 1:16]
