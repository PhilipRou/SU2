include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")

using Plots
# using LsqFit
using LaTeXStrings
using LinearAlgebra
# using Roots
using DelimitedFiles
using Statistics
# using Optim
using QuadGK

fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Lattice_projects\\Lattice2024\\plots"

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
    return quadgk(nastier,-π/2,π/2)[1] / quadgk(nasty,-π/2,π/2)[1] / π^2
end