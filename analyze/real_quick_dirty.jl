include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")

beta = 6.0
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
# ops_smeared = readdlm(timeseries_2Dsmear_path, skipstart=1)
# n_op = size(ops_smeared,2)
# bsize_ops_smeared = round(Int, 2*maximum([auto_corr_time(ops_smeared[:,j]) for j = 1:n_op])) + 1

ops = readdlm(timeseries_2Dsmear_path, skipstart=1)
n_op = size(ops,2)
bsize_ops = round(Int, 2*maximum([auto_corr_time(ops[:,j]) for j = 1:n_op])) + 1

swil  = ops[:,1]
plaqs = 1 .- swil
jackknife(plaqs,bsize_ops)
