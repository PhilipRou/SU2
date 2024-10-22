include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")



for beta in 6.0:1.5:15.0
    beta        = string(beta)*"0"
    # beta        = "15.00"
    n_stout     = 7
    Nz = Nx     = 32
    smear_nums  = [0,1,3,7,15]
    rho_2D      = 0.24

    base_path       = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\fortran_projects\\SU2\\data\\beta_$(beta)_Nz_$(Nz)_Nx_$(Nx)_n_stout_$(n_stout)_rho_2D_$(rho_2D)"
    timeseries_path = base_path * ".txt"
    corr_path       = base_path * "_corr.txt"

    # s_wil    s_opt    smeared_s_wil  smeared_clover  accrate_metro  accrate_ovlax  time
    raw_tseries    = readdlm(timeseries_path, skipstart=1)
    swil_tseries, sopt_tseries, swil_sm_tseries, sopt_sm_tseries = [raw_tseries[:,j] for j = 1:4]
    tau_swil, tau_sopt, tau_swil_sm, tau_sopt_sm = [auto_corr_time(raw_tseries[:,j]) for  j = 1:4]
    bsize_swil, bsize_sopt, bsize_swil_sm, bsize_sopt_sm = [round(Int, 2*tau+1) for tau in [tau_swil, tau_sopt, tau_swil_sm, tau_sopt_sm]]
    swil,    swil_err    = jackknife(swil_tseries,    bsize_swil)
    swil_sm, swil_sm_err = jackknife(swil_sm_tseries, bsize_swil_sm)


    raw_corr            = readdlm(corr_path, skipstart=1);
    n_meas              = minimum([size(raw_corr,1), size(raw_tseries,1)])
    taus_corr_swil      = [auto_corr_time(raw_corr[1:n_meas,j]) for j = 1:Nz]
    taus_corr_swil_sm   = [auto_corr_time(raw_corr[1:n_meas,j]) for j = Nz+1:2*Nz]
    bsizes_corr_swil    = round.(Int, 2 .* taus_corr_swil .+1)
    bsizes_corr_swil_sm = round.(Int, 2 .* taus_corr_swil_sm .+1)
    # n_op         = Int(sqrt(size(raw_corr, 2)))
    # corr_mats = [reshape(raw_corr[meas,:], n_op, n_op) for meas = 1:n_meas];
    # corrs     = 


    swil_corr_con_mean    = Vector{Float64}(undef,Nz);
    swil_corr_con_err     = Vector{Float64}(undef,Nz);
    swil_sm_corr_con_mean = Vector{Float64}(undef,Nz);
    swil_sm_corr_con_err  = Vector{Float64}(undef,Nz);
    for z = 1:Nz
        # t = 1
        swil_corrs    = raw_corr[:,z]
        swil_sm_corrs = raw_corr[:,z+Nz]
        swil_bsize    = maximum([bsize_swil,    bsizes_corr_swil[z]])
        swil_sm_bsize = maximum([bsize_swil_sm, bsizes_corr_swil_sm[z]])
        swil_corr_con_mean[z], swil_corr_con_err[z]       = jack_conn_corr_self(swil_corrs[1:n_meas],    swil_tseries[1:n_meas], swil_bsize)
        swil_sm_corr_con_mean[z], swil_sm_corr_con_err[z] = jack_conn_corr_self(swil_sm_corrs[1:n_meas], swil_sm_tseries[1:n_meas], swil_sm_bsize)
    end

    x_vals = Vector(0:Nz-1)
    x_inds = 1:Nz
    image_swil_corr = scatter(x_vals[x_inds], 
        circshift(swil_corr_con_mean,1)[x_inds], 
        yerror = circshift(swil_corr_con_err,1)[x_inds],
        title = "Conn. s_wil correlator \n β = $beta, L/a = $Nz",
        label = :false,
        yaxis = :log,
        ylim = (1e-18, 1e-4)
    )
    display(image_swil_corr)

    x_vals = Vector(0:Nz-1)
    x_inds = 1:Nz
    image_swil_sm_corr = scatter(x_vals[x_inds], 
        circshift(swil_sm_corr_con_mean,1)[x_inds], 
        yerror = circshift(swil_sm_corr_con_err,1)[x_inds],
        title = "Conn. smeared s_wil correlator \n β = $beta, L/a = $Nz",
        label = :false, 
        yaxis = :log,
        ylim = (1e-18, 1e-4)
    )
    display(image_swil_sm_corr)
end

# println("Everything ok")