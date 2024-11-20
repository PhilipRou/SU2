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

cb_hexes  = [
    "#377eb8", # blue
    "#ff7f00", # orange
    "#4daf4a", # green
    "#984ea3", # purple
    "#e41a1c", # red
    "#f781bf", # pink
    "#999999", # grey
    "#a65628", # brown
];
cb_colors = parse.(Colorant, cb_hexes);
cb_blue, cb_orange, cb_green, cb_purple, cb_red, cb_pink, cb_grey, cb_brown  = cb_colors;

function show_cb_colors()
    display(cb_colors)
    return nothing
end

function format_x_err(x,x_err,error_digs)
    sig_digs = round(Int, -log10(x_err), RoundUp) + (error_digs-1) 
    rounded_x     = round(x, sigdigits = sig_digs)
    rounded_x_err = round(x_err, sigdigits = sig_digs)
    x_out     = @sprintf("%.*f", sig_digs, rounded_x)
    x_err_out = round(Int, rounded_x_err * 10^sig_digs)
    return "$x_out($x_err_out)"
end

# format_x_err(5.310,0.002,2)