using LaTeXStrings
using Plots
using SpecialFunctions
using QuadGK
using Roots
using Polylogarithms

# testar = Vector(1:100)
# plot!(testar,analytic_susc_U2.(testar))

function analytic_plaq_U1(β)
    return besseli(1,β) / besseli(0,β)
end

function analytic_plaq_SU2(β)
    return (besseli(0,β) + besseli(2,β)) / (2*besseli(1,β)) - 1/β
end

function analytic_plaq_U2(β)
    numer(α) = besseli(0,β*cos(α)) + besseli(2,β*cos(α))
    denom(α) = 2*besseli(1,β*cos(α))/cos(α)
    return quadgk(numer,0,π)[1]/quadgk(denom,0,π)[1] - 1/β
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

function a(k,ν)
    return prod([4*ν^2-(2*i-1)^2 for i = 1:k])/(factorial(big(k))*8^k)
end

function analytic_plaq_U1_largebeta(β)
    # return polylog(1,a(1,0)/β) - polylog(1,a(1,1)/β)
    return (1-a(1,1)/β) / (1-a(1,0)/β)
end

function analytic_plaq_SU2_largebeta(β)
    return (2 - (a(1,0)+a(1,2))/β )/(1-a(1,1)/β)/2 #- 2/β
end
function log_plaq_SU2_largebeta(β)
    # return (2 - (a(1,0)+a(1,2))/β )/(1-a(1,1)/β)/2 #- 2/β
    return (a(1,0)-a(1,2))/(2*β) #- (1-a(1,1)/β)/(2-(a(1,0)+a(1,2))/β)*2/β
end

# plot(1:500, [log(analytic_plaq_SU2(b)) - log_plaq_SU2_largebeta(b) for b = 1:500])

#=
function LCP_beta(β_1, V_1, V_2, group)
    if group == "SU(2)"
        # function F(β_2)
        #     return log(analytic_plaq_SU2(β_2))*V_2 - log(analytic_plaq_SU2(β_1))*V_1
        # end
        return find_zero(β_2 -> log(analytic_plaq_SU2(β_2))*V_2 - log(analytic_plaq_SU2(β_1))*V_1, [β_1, 1.1*β_1/V_1*V_2])
    elseif group == "U(2)"
        # function F(β_2)
        #     return log(analytic_plaq_U2(β_2))*V_2 - log(analytic_plaq_U2(β_1))*V_1
        # end
        return find_zero(β_2 -> log(analytic_plaq_U2(β_2))*V_2 - log(analytic_plaq_U2(β_1))*V_1, [β_1, 1.1*β_1/V_1*V_2])
    elseif group == "U(1)"
        # function F(β_2)
        #     return log(analytic_plaq_U1(β_2))*V_2 - log(analytic_plaq_U1(β_1))*V_1
        # end
        return find_zero(β_2 -> log(analytic_plaq_U1(β_2))*V_2 - log(analytic_plaq_U1(β_1))*V_1, [β_1, 1.1*β_1/V_1*V_2])
    end
end

function test_F(β_1, β_2, V_1, V_2, group)
    if group == "SU(2)"
        return log(analytic_plaq_SU2(β_2))*V_2 - log(analytic_plaq_SU2(β_1))*V_1
    elseif group == "U(2)"
        return log(analytic_plaq_U2(β_2))*V_2 - log(analytic_plaq_U2(β_1))*V_1
    elseif group == "U(1)"
        return log(analytic_plaq_U1(β_2))*V_2 - log(analytic_plaq_U1(β_1))*V_1
    end
end

function LCP_a(a_1, β_1, β_2, group)
    if group == "SU(2)"
        return a_1 * sqrt(log(analytic_plaq_SU2(β_2)) / log(analytic_plaq_SU2(β_1)))
    elseif group == "U(2)"
        return a_1 * sqrt(log(analytic_plaq_U2(β_2)) / log(analytic_plaq_U2(β_1)))
    elseif group == "U(1)"
        return a_1 * sqrt(log(analytic_plaq_U1(β_2)) / log(analytic_plaq_U1(β_1)))
    end
end


# Possible arguments for the group G: 
#   "U(1)"
#   "U(2)"
#   "SU(2)"
function plot_LCP_beta(N1,N2,G)
    V1 = N1^2
    V2 = N2^2
    image_LCP = plot(
        sqrt.(Vector(V1:16:V2)), 
        [2.0*V/V1 for V in V1:16:V2],
        color = :grey25,
        linestyle = :dash,
        label = L"$\beta_2 = \beta_1 \cdot V_2^\Lambda / V_1^\Lambda $" ,
        # title = "Exemplary Lines of Constant Physics \n in 2D $G",
        xlabel = "N_x = N_t",
        ylabel = "β",
        xticks = N1:16:N2
        # foreground_color_legend = nothing
    )
    image_LCP = plot!(
        sqrt.(Vector(V1:16:V2)),
        [4.0*V/V1 for V in V1:16:V2],
        color = :grey25,
        linestyle = :dash,
        label = :false
    )
    image_LCP = plot!(
        sqrt.(Vector(V1:16:V2)),
        [6.0*V/V1 for V in V1:16:V2],
        color = :grey25,
        linestyle = :dash,
        label = :false
    )
    image_LCP = plot!(
        sqrt.(Vector(V1:16:V2)), 
        [LCP_beta(2.0, V1, V, G) for V in V1:16:V2],
        # label = "β = 2.0 for N_x = N_t = 32",
        label = L"$\beta_2$ for $\beta_1 = 2.0$ at $N_x = N_t = 32$",
        color = palette(:default)[1]
    )
    image_LCP = plot!(
        sqrt.(Vector(V1:16:V2)), 
        [LCP_beta(4.0, V1, V, G) for V in V1:16:V2],
        # label = "β = 2.0 for N_x = N_t = 32",
        label = L"$\beta_2$ for $\beta_1 = 4.0$ at $N_x = N_t = 32$",
        color = palette(:default)[2]
    )
    image_LCP = plot!(sqrt.(Vector(V1:16:V2)), 
        [LCP_beta(6.0, V1, V, G) for V in V1:16:V2],
        # label = "β = 2.0 for N_x = N_t = 32",
        label = L"$\beta_2$ for $\beta_1 = 6.0$ at $N_x = N_t = 32$",
        color = palette(:default)[3]
    )
    return image_LCP
end
=#

function LCP_a_U1(a_1, β_1, β_2)
    return a_1 * sqrt(log(analytic_plaq_U1(β_2)) / log(analytic_plaq_U1(β_1)))
end

function LCP_a_SU2(a_1, β_1, β_2)
    return a_1 * sqrt(log(analytic_plaq_SU2(β_2)) / log(analytic_plaq_SU2(β_1)))
end

function LCP_a_U2(a_1, β_1, β_2)
    return a_1 * sqrt(log(analytic_plaq_U2(β_2)) / log(analytic_plaq_U2(β_1)))
end

# test_image = plot(2:1:64, [test_F(2.0, beta_2, 32^2, 128^2, "U(2)") for beta_2 in 2:1:64])
# display(test_image)

# display(plot_LCP_beta(32,128,"U(1)"))

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\LCP_beta_U1.pdf")



# a_2 = a_factor ⋅ a_1
function LCP_beta_U1(a_factor, β_1)
    return find_zero(β_2 -> analytic_plaq_U1(β_2) - analytic_plaq_U1(β_1)^(a_factor^2), [β_1, 1.1*β_1/a_factor^2])
end

function LCP_beta_SU2(a_factor, β_1)
    return find_zero(β_2 -> analytic_plaq_SU2(β_2) - analytic_plaq_SU2(β_1)^(a_factor^2), [β_1, 1.1*β_1/a_factor^2])
end

function LCP_beta_U2(a_factor, β_1)
    return find_zero(β_2 -> analytic_plaq_U2(β_2) - analytic_plaq_U2(β_1)^(a_factor^2), [β_1, 1.1*β_1/a_factor^2])
end

function LCP_beta_old(a_factor, β_1)
    return β_1/a_factor^2
end
# LCP_beta_U1(0.5,2.0)

# f(β_1, β_2, a_factor) = analytic_plaq_U1(β_2) - analytic_plaq_U1(β_1)^(a_factor^2)
# display(plot(Vector(1:0.1:16), [f(1.0, x, 0.25) for x in 1:0.1:16]))



# function LCP_a_U1_largebeta(a_1, β_1, β_2)




function plot_LCP_a(β1,β2)
    sepa = 0.1
    a1 = 1.0
    image_LCP = plot(
        Vector(β1:sepa:β2), 
        [a1*sqrt(β1/beta) for beta in β1:sepa:β2],
        color = :grey25,
        linestyle = :dash,
        label = L"$a_2 = a_1 \cdot \sqrt{\beta_1/\beta_2} $" ,
        # title = "Ratio of Lattice Spacings in 2D for β₁ = $β1",
        xlabel = L"\beta_2",
        ylabel = L"a_2/a_1",
        xaxis = :log2,
        yaxis = :log,
        # xticks = β1:bla:β2
        xticks = [2^n for n = 1:2:9],
        # foreground_color_legend = nothing
    )
    image_LCP = plot!(
        Vector(β1:sepa:β2), 
        [LCP_a_U1(a1,β1,beta) for beta in β1:sepa:β2] .+ sqrt(β1/256) .- LCP_a_U1(a1,β1,256),
        label = L"$a_2$ from $\sigma$ in U(1)",
        color = palette(:default)[1]
    )
    image_LCP = plot!(
        Vector(β1:sepa:β2), 
        [LCP_a_SU2(a1,β1,beta) for beta in β1:sepa:β2] .+ sqrt(β1/256) .- LCP_a_SU2(a1,β1,256),
        label = L"$a_2$ from $\sigma$ in SU(2)",
        color = palette(:default)[2]
    )
    image_LCP = plot!(
        Vector(β1:sepa:β2), 
        [LCP_a_U2(a1,β1,beta) for beta in β1:sepa:β2] .+ sqrt(β1/256) .- LCP_a_U2(a1,β1,256),
        label = L"$a_2$ from $\sigma$ in U(2)",
        color = palette(:default)[3]
    )
    image_LCP = plot!(
        [256],
        seriestype = :vline,
        color = :red,
        linestyle = :dashdot,
        label = "LCPs laid on top of \n each other at β₂ = 256",
        foreground_color_legend = nothing,
        legend = :bottomleft
    )
    image_LCP = plot!(
        tickfontsize = 11,
        labelfontsize = 16,
        legendfontsize = 9,
        background_color_legend = nothing
    )
    return image_LCP
end

# plot_LCP_a(2.0, 512)

function susc_LCP_U1(β_1, lower_a_factor)
    x_vals = Vector(lower_a_factor:0.01:0.999)
    image_susc = plot(
        x_vals.^2, 
        [analytic_susc_U1(β) for β in [LCP_beta_U1(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\ln\langle P_{xt}\rangle (\beta_2) }{ \ln\langle P_{xt}\rangle (\beta_1)} } $"
    )
    image_susc = plot!(
        x_vals.^2, 
        [analytic_susc_U1(β) for β in [LCP_beta_old(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\beta_1}{\beta_2}} $",
        linestyle = :dash)
    image_susc = plot!(
        # title = "Cont. Limit of Topological Susc. \n in 2D U(1) with β₁ = $β_1",
        xlabel = L"$(a_2 / a_1)^2$",
        foreground_color_legend = nothing
    )
    image_susc = scatter!(
        [1.0, 0.0],
        [analytic_susc_U1(β_1), 0.0],
        markershape = :star,
        markersize = 7,
        color = :white,
        label = L"$\chi_{top}$ for $a_2/a_1 \in \{0,1\}$"
    )
    image_susc = plot!(
        tickfontsize = 11,
        labelfontsize = 16,
        legendfontsize = 9,
        background_color_legend = nothing
    )
    display(image_susc)
end

function susc_LCP_U2(β_1, lower_a_factor)
    x_vals = Vector(lower_a_factor:0.01:0.999)
    image_susc = plot(
        x_vals.^2, 
        [analytic_susc_U2(β) for β in [LCP_beta_U2(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\ln\langle P_{xt}\rangle (\beta_2) }{ \ln\langle P_{xt}\rangle (\beta_1)} } $"
    )
    image_susc = plot!(
        x_vals.^2, 
        [analytic_susc_U2(β) for β in [LCP_beta_old(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\beta_1}{\beta_2}} $",
        linestyle = :dash)
    image_susc = plot!(
        # title = "Cont. Limit of Topological Susc. \n in 2D U(2) with β₁ = $β_1",
        xlabel = L"$(a_2 / a_1)^2$",
        foreground_color_legend = nothing
    )
    display(image_susc)
    image_susc = scatter!(
        [1.0, 0.0],
        [analytic_susc_U2(β_1), 0.0],
        markershape = :star,
        markersize = 7,
        color = :white,
        label = L"$\chi_{top}$ for $a_2/a_1 \in \{0,1\}$"
    )
    image_susc = plot!(
        tickfontsize = 11,
        labelfontsize = 16,
        legendfontsize = 9,
        background_color_legend = nothing
    )
end

function plaq_LCP_U1(β_1, lower_a_factor)
    x_vals = Vector(lower_a_factor:0.01:0.999)
    image_plaq = plot(
        x_vals.^2, 
        [analytic_plaq_U1(β) for β in [LCP_beta_U1(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\ln\langle P_{xt}\rangle (\beta_2) }{ \ln\langle P_{xt}\rangle (\beta_1)} } $"
    )
    image_plaq = plot!(
        x_vals.^2, 
        [analytic_plaq_U1(β) for β in [LCP_beta_old(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\beta_1}{\beta_2}} $",
        linestyle = :dash)
    image_plaq = plot!(
        # title = "Cont. Limit of the Plaquette \n in 2D U(1) with β₁ = $β_1",
        xlabel = L"$(a_2 / a_1)^2$",
        legend = :bottomleft,
        foreground_color_legend = nothing
    )
    image_plaq = scatter!(
        [1.0, 0.0],
        [analytic_plaq_U1(β_1), 1.0],
        markershape = :star,
        markersize = 7,
        color = :white,
        label = L"$P_{xt}$ for $a_2/a_1 \in \{0,1\}$"
    )
    image_plaq = plot!(
        tickfontsize = 11,
        labelfontsize = 16,
        legendfontsize = 9,
        background_color_legend = nothing
    )
    display(image_plaq)
end

function plaq_LCP_SU2(β_1, lower_a_factor)
    x_vals = Vector(lower_a_factor:0.01:0.999)
    image_plaq = plot(
        x_vals.^2, 
        [analytic_plaq_SU2(β) for β in [LCP_beta_U1(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\ln\langle P_{xt}\rangle (\beta_2) }{ \ln\langle P_{xt}\rangle (\beta_1)} } $"
    )
    image_plaq = plot!(
        x_vals.^2, 
        [analytic_plaq_SU2(β) for β in [LCP_beta_old(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\beta_1}{\beta_2}} $",
        linestyle = :dash)
    image_plaq = plot!(
        # title = "Cont. Limit of the Plaquette \n in 2D U(2) with β₁ = $β_1",
        xlabel = L"$(a_2 / a_1)^2$",
        legend = :topright,
        foreground_color_legend = nothing
    )
    image_plaq = scatter!(
        [1.0, 0.0],
        [analytic_plaq_SU2(β_1), 1.0],
        markershape = :star,
        markersize = 7,
        color = :white,
        label = L"$P_{xt}$ for $a_2/a_1 \in \{0,1\}$"
    )
    image_plaq = plot!(
        tickfontsize = 11,
        labelfontsize = 16,
        legendfontsize = 9,
        background_color_legend = nothing
    )
    display(image_plaq)
end

function plaq_LCP_U2(β_1, lower_a_factor)
    x_vals = Vector(lower_a_factor:0.01:0.999)
    image_plaq = plot(
        x_vals.^2, 
        [analytic_plaq_U2(β) for β in [LCP_beta_U1(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\ln\langle P_{xt}\rangle (\beta_2) }{ \ln\langle P_{xt}\rangle (\beta_1)} } $"
    )
    image_plaq = plot!(
        x_vals.^2, 
        [analytic_plaq_U2(β) for β in [LCP_beta_old(a,β_1) for a in x_vals]],
        label = L"$\frac{a_2}{a_1} = \sqrt{\frac{\beta_1}{\beta_2}} $",
        linestyle = :dash)
    image_plaq = plot!(
        # title = "Cont. Limit of the Plaquette \n in 2D U(2) with β₁ = $β_1",
        xlabel = L"$(a_2 / a_1)^2$",
        legend = :topright,
        foreground_color_legend = nothing
    )
    image_plaq = scatter!(
        [1.0, 0.0],
        [analytic_plaq_U2(β_1), 1.0],
        markershape = :star,
        markersize = 7,
        color = :white,
        label = L"$P_{xt}$ for $a_2/a_1 \in \{0,1\}$"
    )
    image_plaq = plot!(
        tickfontsize = 11,
        labelfontsize = 16,
        legendfontsize = 9,
        background_color_legend = nothing
    )
    display(image_plaq)
end





β_1 = 2.0

plot_LCP_a(β_1, 2^9)
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\master_plots\\LCP_beta_1_$β_1.later.pdf")

plaq_LCP_U2(β_1, 0.0561)
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\master_plots\\U2_plaq_LCP_beta_1_$β_1.pdf")

susc_LCP_U2(β_1, 0.0561)
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\master_plots\\U2_susc_LCP_beta_1_$β_1.pdf")

plaq_LCP_U1(β_1, 0.0561)
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\master_plots\\U1_plaq_LCP_beta_1_$β_1.pdf")

susc_LCP_U1(β_1, 0.0561)
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\master_plots\\U1_susc_LCP_beta_1_$β_1.pdf")

plaq_LCP_SU2(β_1, 0.0561)
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik_Uni\\Master_Thesis\\master_plots\\SU2_plaq_LCP_beta_1_$β_1.pdf")
