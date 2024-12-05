include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")



# for beta in 6.0:1.5:15.0
    beta = 6.0
    beta       = string(beta)*"0"
    # beta       = "6.00"
    n_stout    = 7
    Nz = Nx    = 32
    smearlist  = [0,1,3,7,15]
    rho_2D     = 0.24

    n_op = length(smearlist)
    small_inds = 2:length(smearlist)

    base_path         = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\fortran_projects\\SU2_3D_data\\good_old_data\\beta_$(beta)_Nz_$(Nz)_Nx_$(Nx)_smearlist_rho_2D_$(rho_2D)"
    timeseries_path   = base_path * ".txt"
    corr_path         = base_path * "_crosscorr.txt"
    swil_smeared_path = base_path * "_swil_smeared.txt"

    ### s_wil    s_opt    smeared_s_wil  smeared_clover  accrate_metro  accrate_ovlax  time
    raw_tseries = readdlm(timeseries_path, skipstart=1)
    swil_smeared = readdlm(swil_smeared_path, skipstart=1)
    swil_smeared_small = swil_smeared[:,small_inds]
    bsize_swil_smeared = round(Int, 2*maximum([auto_corr_time(swil_smeared[:,j]) for j = 1:n_op])) + 1

    raw_corr  = readdlm(corr_path, skipstart=1);
    n_meas    = minimum([size(raw_corr,1), size(raw_tseries,1)])
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
    corrmats_symm_small = []
    push!(corrmats_symm, corrmats[Nz][:]);
    push!(corrmats_symm_small, [corrmats[Nz][meas][small_inds, small_inds] for meas = 1:n_meas]);
    for t = 1:Nz>>1-1
        push!(corrmats_symm, (corrmats[t][:] .+ corrmats[Nz-t][:]) ./ 2) 
        push!(corrmats_symm_small, [(corrmats[t][meas][small_inds, small_inds] .+ corrmats[Nz-t][meas][small_inds, small_inds]) ./ 2 for meas = 1:n_meas ])
    end
    push!(corrmats_symm, corrmats[Nz>>1][:]);
    push!(corrmats_symm_small, [corrmats[Nz>>1][meas][small_inds, small_inds] for meas = 1:n_meas]);
    
    ### Symmetrize every corr. matrix (in operator space)
    for t = 1:Nz>>1+1
        ## push!(corrmats_symm, [corrmats[t][meas] + transpose(corrmats[t][meas])/ 2 for meas = 1:n_meas])
        ## push!(corrmats_symm_small, [corrmats[t][meas][small_inds, small_inds] + transpose(corrmats[t][meas][small_inds, small_inds])/ 2 for meas = 1:n_meas])
        corrmats_symm[t][:] = (corrmats_symm[t][:] .+ transpose.(corrmats_symm[t][:])) ./ 2
        corrmats_symm_small[t][:] = (corrmats_symm_small[t][:] .+ transpose.(corrmats_symm_small[t][:])) ./ 2
    end

    # mean(corrmats_symm_small[1])[1:4,1:4]
    

    ### Connected correlators of individual smeared s_wil series
    swil_corr_con_mean = Array{Float64}(undef, length(smearlist), Nz>>1+1);
    swil_corr_con_err = Array{Float64}(undef, length(smearlist), Nz>>1+1);
    for n_smear in smearlist
        # n_smear = 1
        smear_ind = findall(==(n_smear), smearlist)[1]
        for t = 1:Nz>>1+1
        # t = 17
            temp_corrs = [corrmats_symm[t][meas][smear_ind,smear_ind] for meas = 1:n_meas]
            temp_means = swil_smeared[:,smear_ind]
            swil_corr_con_mean[smear_ind,t], swil_corr_con_err[smear_ind,t] = jack_conn_corr_self(temp_corrs, temp_means, bsize_corrmats)
        end # t
    end # n_smear
    let
        image_con = plot(
            title  = latexstring("Connected \$s_\\mathrm{wil}\$ corr. \n \$β = $beta, L = $Nz, ρ = $rho_2D\$"),
            ylims  = (1e-14, 1e-4),
            yaxis  = :log,
            legend = :bottomleft,
            xlabel = latexstring("\$t\$")
        )
        for smear_ind = 1:length(smearlist)
            image_con = scatter!(
                0:Nz>>1, 
                swil_corr_con_mean[smear_ind,:], 
                yerror = swil_corr_con_err[smear_ind,:], 
                label = latexstring("\$n_\\mathrm{smear} = $(smearlist[smear_ind])\$"),
                markerstrokecolor = :auto
            )
        end
        display(image_con)
    end
# end # beta


    
    ### EVs of correlation matrices
    small    = true
    evs      = Array{Float64}(undef, Nz>>1+1, n_op)
    evs_errs = Array{Float64}(undef, Nz>>1+1, n_op)
    if small
        evs      = Array{Float64}(undef, Nz>>1+1, length(small_inds))
        evs_errs = Array{Float64}(undef, Nz>>1+1, length(small_inds))
    end
    @time for t = 1:Nz>>1+1
        @show t
        jack = [zeros(n_op), zeros(n_op)]
        if small
            jack = jack_conn_corr_mat_ev(corrmats_symm_small[t], swil_smeared_small, maximum([bsize_corrmats,bsize_swil_smeared]))
        else
            jack = jack_conn_corr_mat_ev(corrmats_symm[t], swil_smeared, maximum([bsize_corrmats,bsize_swil_smeared]))
        end
        evs[t,:]      = reverse(jack[1])
        evs_errs[t,:] = reverse(jack[2])
    end

    let
        image_conn_evs = plot(
            title  = latexstring("EV's of conn. corr. mats. of smeared \$s_\\mathrm{wil}\$ \n \$\\beta = $beta, L = $Nz, n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$t\$"),
            )
        for i = 1:size(evs,2)
            image_conn_evs = scatter!(
                Vector(0:Nz>>1) .+ (i-1)*0.1,
                evs[:,i],
                yerror = evs_errs[:,i],
                label = latexstring("EV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-16, 5e-6),
                yaxis = :log,
                legend = :bottomright
            )
        end
        display(image_conn_evs)
    end
# end # beta



    ### GEVs of the correlation matrices
    ### t0 is the index of the lower bound for the time-index of the GEVP, starting at 1!
    t0 = 1
    small     = true
    GEVs      = Array{Float64}(undef, Nz>>1-t0, n_op)
    GEVs_errs = Array{Float64}(undef, Nz>>1-t0, n_op)
    if small
        GEVs      = Array{Float64}(undef, Nz>>1-t0, length(small_inds))
        GEVs_errs = Array{Float64}(undef, Nz>>1-t0, length(small_inds))
    end
    @time for t = t0+1:Nz>>1
        @show t
        jack = [zeros(n_op), zeros(n_op)]
        if small
            jack = jack_conn_corr_mat_GEV(corrmats_symm_small[t], corrmats_symm_small[t0], swil_smeared_small, maximum([bsize_corrmats,bsize_swil_smeared]))
        else
            jack = jack_conn_corr_mat_GEV(corrmats_symm[t], corrmats_symm[t0], swil_smeared, maximum([bsize_corrmats,bsize_swil_smeared]))
        end
        GEVs[t-t0,:]      = reverse(jack[1])
        GEVs_errs[t-t0,:] = reverse(jack[2])
    end

    let
        image_conn_GEVs = plot(
            title  = latexstring("Gen. EV's of conn. corr. mats. of smeared \$s_\\mathrm{wil}\$ \n \$\\beta = $beta\$, \$L = $Nz\$, \$t_0 = $t0\$, \$n_\\mathrm{smear} \\in\$, $(smearlist[2:end])"),
            xlabel = latexstring("\$ t\$"),
            )
        for i = 1:size(GEVs,2)
            image_conn_GEVs = scatter!(
                Vector(t0+1:Nz>>1) .+ (i-1)*0.1,
                GEVs[:,i],
                yerror = GEVs_errs[:,i],
                label = latexstring("GEV nr. \$$i\$"),
                markerstrokecolor = :auto,
                ylim = (1e-5, 10.0),
                yaxis = :log,
                legend = :bottomright
            )
        end
        display(image_conn_GEVs)
    end

# end # beta
