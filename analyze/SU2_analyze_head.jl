using Plots
using StatsBase
using DelimitedFiles
using LsqFit
using SpecialFunctions
using LaTeXStrings
using Measures

include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_hex.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_cube.jl")



enu_endings = ["st", "nd", "rd"]
for i = 1:50
    push!(enu_endings, "th")
end


function plot_corrs(N_t, corr_means, corr_errs, series_label)
    return scatter(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label, markerstrokecolor = :auto)
end

function plot_corrs!(N_t, corr_means, corr_errs, series_label, image)
    image = scatter!(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label, markerstrokecolor = :auto)
    return nothing
end

function plot_corrs(N_t, corr_means, corr_errs, series_label, farbe)
    return scatter(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label, color = farbe, markerstrokecolor = farbe)
end

function plot_corrs!(N_t, corr_means, corr_errs, series_label, image, farbe)
    image = scatter!(0:N_t-1, circshift(corr_means,1), yerror = circshift(corr_errs,1), label = series_label, color = farbe, markerstrokecolor = farbe)
    return nothing
end

function s_from_plaq(β, plaq_mean)
    return 3 * β * (2-plaq_mean) / 2
end

function s_wil(plaq_mean)
    return 1 - plaq_mean/2
end