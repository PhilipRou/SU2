include("D:\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("D:\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("D:\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("D:\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")

using Plots
using LsqFit
using LaTeXStrings
using LinearAlgebra
using Roots
using DelimitedFiles
using Statistics


# insta_update() doesn't really work, but cooling the proposal
# only a couple of times (N_insta_cool) seems to do the trick
function insta_cool_update_U2!(U, step, β, N_insta_cool, acc)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U2(NX,NT,rand([-2,2])) .* U
    for i = 1:N_insta_cool
        chess_cool!(U_prop,step,β,[0],"U2")
    end
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end

# insta_update() doesn't really work, but cooling the proposal
# only a couple of times (N_insta_cool) seems to do the trick
# (which does not mean that the results obtained are correct
# in any way).
function insta_ran_cool_update_U2!(U, step, β, N_insta_cool, cool_degree_in_percent, acc)
    NX = size(U,2)
    NT = size(U,3)
    U_prop = insta_U2_naive(NX,NT,rand([-2,2])) .* U
    for i = 1:N_insta_cool
        ran_cool!(U_prop,cool_degree_in_percent,step,β,[0],"U2")
    end
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT*cool_degree_in_percent/100
    end
    return nothing
end

function insta_flow_update_U2!(U, N_stout, acc, Q)
    NX = size(U,2)
    NT = size(U,3)
    ρ  = 1/N_stout
    U_prop = stout(U,ρ)
    for i = 1:N_stout-1
        U_prop = stout(U_prop,ρ)
    end
    U_prop = insta_U2_naive(NX,NT,rand([-Q,Q])) .* U_prop
    for i = 1:N_stout
        U_prop = stout(U_prop,-ρ)
    end
    if rand() < exp(action(U,β) - action(U_prop,β)) # Definitely old minus new!!
        U[:,:,:] = U_prop[:,:,:]
        acc[1] += 2*NX*NT
    end
    return nothing
end





λ1 = Complex.([0 1 0; 1 0 0; 0 0 0])
λ2 = Complex.([0 -im 0; im 0 0; 0 0 0])
λ3 = Complex.([1 0 0; 0 -1 0; 0 0 0])
λ4 = Complex.([0 0 1; 0 0 0; 1 0 0])
λ5 = Complex.([0 0 -im; 0 0 0; im 0 0])
λ6 = Complex.([0 0 0; 0 0 1; 0 1 0])
λ7 = Complex.([0 0 0; 0 0 -im; 0 im 0])
λ8 = Complex.([1 0 0; 0 1 0; 0 0 -2]) / sqrt(3)

Λ = [λ1, λ2, λ3, λ4, λ5, λ6, λ7, λ8]

λ0 = Complex.([1 0 0; 0 1 0; 0 0 1])
λ_zero = Complex.([0 0 0; 0 0 0; 0 0 0])

N_x = N_t = 32

function ran_U3(ϵ)
    return exp(ϵ * im * sum((rand(8) .- 0.5).*Λ))
end

function action_U3(U, β)
    NX = size(U,2)
    NT = size(U,3)
    S = 3*NX*NT   # later generalization: N_colour * NT * (NX)^d_s
    for t = 1:NT
        for x = 1:NX
            S -= real(tr(plaq(U,x,t)))
        end
    end
    return β*S/3    # later generalization: β*S/N_colour
end

function insta_U3_attempt(N_x, N_t, Q)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((im*Q*t*2*π)/(sqrt(3)*N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp(-(im*Q*x*2*π)/(sqrt(3)*N_x) * λ8) for x = 1:N_x]
    return U
    # if Q%3 == 0
    #     U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * λ0 for x = 1:N_x, t = 1:N_t]
    #     U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    #     U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * λ0 for x = 1:N_x]
    # elseif (Q-1)%3 == 0
    #     U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((im*t*2*π)/(sqrt(3)*N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    #     U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    #     U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp(-im*x*2*π/(sqrt(3)*N_x) * λ8) for x = 1:N_x]
    # elseif (Q-2)%3 == 0
    #     U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((im*2*t*2*π)/(sqrt(3)*N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    #     U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    #     U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp(-im*2*x*2*π/(sqrt(3)*N_x) * λ8) for x = 1:N_x]
    # end
    return U
end

function insta_U3_w(N_x, N_t, Q, z)
    w = sqrt(3) * 2 * pi * (z-Q/3)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((-im*t*w)/(N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp((im*x*w)/(N_x) * λ8) for x = 1:N_x]
    return U
end

function cool!(U, μ, x, t, step, β, acc, group)
    new_coeffs = U[μ,x,t]
    if group == "SU2"
        new_coeffs = ran_SU2(step) * new_coeffs
    elseif group == "U2"
        new_coeffs = ran_U2(step) * new_coeffs
    elseif group == "U3"
        new_coeffs = ran_U3(step) * new_coeffs
    end
    staple_d = staple_dag(U,μ,x,t)
    S_old = real(tr(U[μ,x,t] * staple_d)) # ❗ technically missing global factors like β,
    S_new = real(tr(new_coeffs * staple_d)) # but we don't need them here
    if S_old < S_new
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function insta_U2_w(N_x, N_t, Q, z)
    w = π*(2*z-Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*w/N_x/N_t)*coeffs_Id_U2() - sin(t*w/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*w/N_x)*coeffs_Id_U2() + sin(x*w/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    return U
end

function insta_U2_opt(N_x, N_t, Q)
    w = π*(2*mod(Q,2)*π -Q)
    U = Array{coeffs_U2}(undef, 2, N_x, N_t)
    U[1,:,:]       = [exp(-im*Q*t*π/N_x/N_t) * (cos(t*w/N_x/N_t)*coeffs_Id_U2() - sin(t*w/N_x/N_t)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [coeffs_Id_U2() for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*π/N_x) * (cos(x*w/N_x)*coeffs_Id_U2() + sin(x*w/N_x)*coeffs_U2(0.0*im, 0.0*im, 0.0*im, 1.0 + 0.0*im)) for x = 1:N_x]
    return U
end

function insta_action(β, N_c, N_x, N_t, Q, z)
    ReTrP = (N_c-1)*cos(2*π*z/N_x/N_t) + cos(2*π*(Q - (N_c-1)*z)/N_x/N_t ) 
    return β*N_x*N_t*(1 - ReTrP/N_c)
end

function insta_U3(N_x, N_t, Q)
    U = Array{Matrix}(undef, 2, N_x, N_t)
    w = 0.0
    if mod(Q%3,3) == 1
        w = -2*π/sqrt(3)
    elseif mod(Q%3,3) == 2
        w = 2*π/sqrt(3)
    end
    U[1,:,:]       = [exp(-(im*Q*t*2*π)/(3*N_x*N_t)) * exp((-im*t*w)/(N_x*N_t) * λ8) for x = 1:N_x, t = 1:N_t]
    U[2,:,1:N_t-1] = [λ0 for x = 1:N_x, t = 1:N_t-1]
    U[2,:,N_t]     = [exp(im*Q*x*2*π/(3*N_x)) * exp((im*x*w)/(N_x) * λ8) for x = 1:N_x]
    return U
end

function insta_action_min_U3(β, N_x, N_t, Q, q_plot)
    ReTrP = 0.0
    if mod(Q%3,3) == 0
        ReTrP = 3*cos(2*π/N_x/N_t * q_plot/3)
    elseif mod(Q%3,3) == 1
        ReTrP = 2*cos(2*π/N_x/N_t * (q_plot-1)/3) + cos(2*π/N_x/N_t * (q_plot+2)/3)
    elseif mod(Q%3,3) == 2
        ReTrP = 2*cos(2*π/N_x/N_t * (q_plot+1)/3) + cos(2*π/N_x/N_t * (q_plot-2)/3)
    end
    return β*N_x*N_t*(1 - ReTrP/3)
end

function insta_action_min(β, N_c, N_x, N_t, Q, q_plot)
    m = N_c*(round(Q/N_c, RoundNearestTiesAway) - Q/N_c)
    ReTrP = (N_c-1)*cos(2*π/N_x/N_t * (q_plot+m)/N_c) + cos(2*π/N_x/N_t *(q_plot - (N_c-1)*m)/N_c ) 
    return β*N_x*N_t*(1 - ReTrP/N_c)
end

function insta_action_min_U2(β, N_x, N_t, Q, q_plot)
    ReTrP = 2*cos(q_plot*π/N_x/N_t) 
    if isodd(Q)
        ReTrP *= cos(π/N_x/N_t)
    end
    return β*N_x*N_t*(1-ReTrP/2)
end

# plot(-6:0.01:6, [insta_action_min_U3(1,N_x,N_t,0,q) for q = -6:0.01:6])

# plot(-6:0.01:6, [insta_action_min(1,4,N_x,N_t,0,q) for q = -6:0.01:6])
# plot!(-6:0.01:6, [insta_action_min(1,4,N_x,N_t,4,q) for q = -6:0.01:6])





let
    Q_lower = -11
    Q_upper = 11
    z_lower = -4
    z_upper = 4
    image_insta = plot(
        title = L"$S/\beta$ for Local Minimum Solutions of 2D U(3)",
        xlabel = L"Top. Charge $Q$",
        legend = :top,
        xticks = Q_lower:Q_upper,
        # foreground_color_legend = :false
        # ylims = [-0.1,1.0]
    )
    for z = z_lower:z_upper
        image_insta = plot!(
            Q_lower:0.01:Q_upper, 
            [insta_action(1, 3, N_x, N_t, q, z) for q = Q_lower:0.01:Q_upper],
            label = :false,
            color = palette(:default)[1+z-z_lower],
        )
        if (z-z_lower)%3 == 0
            image_insta = scatter!(
                Q_lower:Q_upper, 
                [action_U3(insta_U3_w(N_x, N_t, q, z), 1) for q = Q_lower:Q_upper],
                label = "z = $z",
                markershape = :circle,
                alpha = 0.6,
                color = palette(:default)[1+z-z_lower],
            )
        elseif (z-z_lower)%3 == 1
            image_insta = scatter!(
                Q_lower:Q_upper, 
                [action_U3(insta_U3_w(N_x, N_t, q, z), 1) for q = Q_lower:Q_upper],
                label = "z = $z",
                markershape = :diamond,
                alpha = 0.6,
                color = palette(:default)[1+z-z_lower],
            )
        elseif (z-z_lower) %3 == 2
            image_insta = scatter!(
                Q_lower:Q_upper, 
                [action_U3(insta_U3_w(N_x, N_t, q, z), 1) for q = Q_lower:Q_upper],
                label = "z = $z",
                markershape = :star4,
                alpha = 0.6,
                color = palette(:default)[1+z-z_lower],
            )
        end
    end
    display(image_insta)
end
image_insta = plot!(xlims = [-10.5, 10.5], ylims = [0.01,0.5], legend = :bottomleft, background_color_legend = "rgba(100%,100%,100%,60%)", foreground_color_legend = :false)
image_insta = plot!(yaxis = :log)

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\local_minima_2D_U3_zoom.pdf")


let
    Q_lower = -10
    Q_upper = 10
    z_lower = -4
    z_upper = 4
    image_insta = plot(
        title = L"$S/\beta$ for Local Minimum Solutions of 2D U(2)",
        xlabel = L"Top. Charge $Q$",
        legend = :top,
        xticks = Q_lower:Q_upper,
        # ylims = 
    )
    for z = z_lower:z_upper
        image_insta = plot!(
            Q_lower:0.01:Q_upper, 
            [insta_action(1, 2, N_x, N_t, q, z) for q = Q_lower:0.01:Q_upper],
            label = :false,
            color = palette(:default)[1+z-z_lower],
        )
        if (z-z_lower)%2 == 0
            image_insta = scatter!(
                Q_lower:Q_upper, 
                [action(insta_U2_w(N_x, N_t, q, z), 1) for q = Q_lower:Q_upper],
                label = "z = $z",
                markershape = :rect,
                alpha = 0.6,
                color = palette(:default)[1+z-z_lower],
            )
        elseif (z-z_lower)%2 == 1
            image_insta = scatter!(
                Q_lower:Q_upper, 
                [action(insta_U2_w(N_x, N_t, q, z), 1) for q = Q_lower:Q_upper],
                label = "z = $z",
                markershape = :diamond,
                alpha = 0.6,
                color = palette(:default)[1+z-z_lower],
            )
        end
    end
    display(image_insta)
end
image_insta = plot!(xlims = [-10.5, 10.5], ylims = [-0.03,0.5], legend = :bottomleft, background_color_legend = "rgba(100%,100%,100%,60%)", foreground_color_legend = :false)
image_insta = plot!(yaxis = :log)

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\local_minima_2D_U2_zoom.pdf")



Q_bound = 10
image_insta = plot(
    title = "S/β for Local Minimum Solutions of 2D U(2) \n Q Taken as Continuous",
    xlabel = L"Top. Charge $Q$",
    legend = :top,
    xticks = -Q_bound:Q_bound
)
for z = -2:2
    image_insta = plot!(
        -Q_bound:0.01:Q_bound, 
        [action(insta_U2_w(N_x, N_t, q, z), 1) for q = -Q_bound:0.01:Q_bound],
        label = "z = $z",
    )
end
display(image_insta)

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\local_minima_2D_U2_Q_cont.pdf")



let
    q_low = -6
    q_up = 6
    versatz = 0.02
    image_einh = plot(
        title = "S/β for Q-sector-global Minimum Solutions in 2D U(3)",
        xticks = q_low:q_up,
        xlabel = "Top. Charge Q",
        legend = :top,
        # ylims = [0.0,0.5]
    )
    image_einh = plot!(
        q_low:0.01:q_up, 
        [insta_action_min_U3(1,N_x,N_t,0,q) for q = q_low:0.01:q_up], 
        color = palette(:default)[1],
        label = :false
    )
    image_einh = plot!(
        q_low:0.01:q_up, 
        [insta_action_min_U3(1,N_x,N_t,1,q) for q = q_low:0.01:q_up], 
        color = palette(:default)[2],
        linestyle = :dot,
        label = :false
    )
    image_einh = plot!(
        q_low+versatz:0.01:q_up, 
        [insta_action_min_U3(1,N_x,N_t,2,q) for q = q_low+versatz:0.01:q_up], 
        color = palette(:default)[3],
        linestyle = :dot,
        label = :false
    )
    image_einh = scatter!(
        q_low:3:q_up, 
        [action_U3(insta_U3(N_x, N_t, q),1) for q = q_low:3:q_up], 
        color = palette(:default)[1],
        markershape = :circle,
        label = L"$Q = 3\cdot z$"
    )
    image_einh = scatter!(
        q_low+1:3:q_up, 
        [action_U3(insta_U3(N_x, N_t, q),1) for q = q_low+1:3:q_up], 
        color = palette(:default)[2],
        markershape = :diamond,
        label = L"$Q = 3\cdot z + 1$"
    )
    image_einh = scatter!(
        q_low+2:3:q_up, 
        [action_U3(insta_U3(N_x, N_t, q),1) for q = q_low+2:3:q_up], 
        color = palette(:default)[3],
        markershape = :star4,
        label = L"$Q = 3\cdot z - 1$"
    )
    display(image_einh)
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\global_minima_2D_U3.pdf")



let
    q_low = -4
    q_up = 4
    versatz = 0.02
    image_einh = plot(
        title = "S/β for Q-sector-global Minimum Solutions in 2D U(2)",
        xticks = q_low:q_up,
        xlabel = "Top. Charge Q",
        legend = :top
    )
    image_einh = plot!(
        q_low:0.01:q_up, 
        [insta_action_min_U2(1,N_x,N_t,0,q) for q = q_low:0.01:q_up], 
        color = palette(:default)[1],
        label = :false
    )
    image_einh = plot!(
        q_low:0.01:q_up, 
        [insta_action_min_U2(1,N_x,N_t,1,q) for q = q_low:0.01:q_up], 
        color = palette(:default)[2],
        linestyle = :dot,
        label = :false
    )
    image_einh = scatter!(
        q_low:2:q_up, 
        [action(insta_U2(N_x, N_t, q),1) for q = q_low:2:q_up], 
        color = palette(:default)[1],
        markershape = :circle,
        label = L"$Q = 2\cdot z$"
    )
    image_einh = scatter!(
        q_low+1:2:q_up, 
        [action(insta_U2(N_x, N_t, q),1) for q = q_low+1:2:q_up], 
        color = palette(:default)[2],
        markershape = :diamond,
        label = L"$Q = 2\cdot z + 1$"
    )
    display(image_einh)
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\global_minima_2D_U2.pdf")



let
    q_low = -8
    q_up = 8
    versatz = 0.02
    image_einh = plot(
        title = "S/β for Q-sector-global Minimum Solutions in 2D U(4)",
        xticks = q_low:q_up,
        xlabel = "Top. Charge Q",
        legend = :top
    )
    image_einh = plot!(
        q_low:0.01:q_up, 
        [insta_action_min(1,4,N_x,N_t,0,q) for q = q_low:0.01:q_up], 
        color = palette(:default)[1],
        label = :false
    )
    image_einh = plot!(
        q_low:0.01:q_up, 
        [insta_action_min(1,4,N_x,N_t,1,q) for q = q_low:0.01:q_up], 
        color = palette(:default)[2],
        # linestyle = :dot,
        label = :false
    )
    image_einh = plot!(
        q_low+versatz:0.01:q_up, 
        [insta_action_min(1,4,N_x,N_t,3,q) for q = q_low+versatz:0.01:q_up], 
        color = palette(:default)[3],
        linestyle = :dash,
        label = :false
    )
    image_einh = plot!(
        q_low+versatz:0.01:q_up, 
        [insta_action_min(1,4,N_x,N_t,2,q) for q = q_low+versatz:0.01:q_up], 
        color = palette(:default)[4],
        # linestyle = :dot,
        label = :false
    )
    ############################################################################
    image_einh = scatter!(
        q_low:4:q_up, 
        [insta_action_min(1,4,N_x,N_t,0,q) for q = q_low:4:q_up], 
        color = palette(:default)[1],
        markershape = :circle,
        label = L"$Q = 4\cdot z$"
    )
    image_einh = scatter!(
        q_low+1:4:q_up, 
        [insta_action_min(1,4,N_x,N_t,1,q) for q = q_low+1:4:q_up], 
        color = palette(:default)[2],
        markershape = :diamond,
        label = L"$Q = 4\cdot z + 1$"
    )
    image_einh = scatter!(
        q_low+3:4:q_up, 
        [insta_action_min(1,4,N_x,N_t,3,q) for q = q_low+3:4:q_up], 
        color = palette(:default)[3],
        markershape = :star4,
        # markersize = 3,
        label = L"$Q = 4\cdot z - 1$"
    )
    image_einh = scatter!(
        q_low+2:4:q_up, 
        [insta_action_min(1,4,N_x,N_t,2,q) for q = q_low+2:4:q_up], 
        color = palette(:default)[4],
        markershape = :rect,
        markersize = 3,
        label = L"$Q = 4\cdot z + 2$"
    )
    display(image_einh)
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\global_minima_2D_U4.pdf")



let 
    N_c   = 2
    q_low = 0
    q_up  = 7
    z_low = -10
    z_up  = 10
    x_low = z_low - 0.2
    x_up  = z_up + 0.2
    y_frame = 0.02

    action_max = 0.0

    image_action_z = plot(
        title = "S/β for Local Minimum Solutions in 2D U($N_c)",
        x_ticks = z_low:z_up,
        xlabel = L"Integer $z$",
        legend = :top,
    )
    for q = q_low:q_up
        actions = [insta_action(1,N_c,N_x,N_t,q,z) for z = z_low:z_up]
        if mod(q,N_c) == 0
            image_action_z = scatter!(
                z_low:z_up,
                actions,
                label = latexstring("\$ Q = $q \$"),
                # label = :false,
                markershape = :circ,
                # markersize = 6,
                color = palette(:default)[1+q-q_low],
                alpha = 0.6,
            )
        elseif mod(q,N_c) == 1
            image_action_z = scatter!(
                z_low:z_up,
                actions,
                label = latexstring("\$ Q = $q \$"),
                # label = :false,
                markershape = :diamond,
                # markersize = 6,
                color = palette(:default)[1+q-q_low],
                alpha = 0.6,
            )
        elseif mod(q,N_c) == 2
            image_action_z = scatter!(
                z_low:z_up,
                actions,
                label = latexstring("\$ Q = $q \$"),
                # label = :false,
                markershape = :star4,
                # markersize = 6,
                color = palette(:default)[1+q-q_low],
                alpha = 0.6,
            )
        elseif mod(q,N_c) == 3
            image_action_z = scatter!(
                z_low:z_up,
                actions,
                label = latexstring("\$ Q = $q \$"),
                # label = :false,
                markershape = :rect,
                # markersize = 6,
                color = palette(:default)[1+q-q_low],
                alpha = 0.6,
            )
        end
    end
    for q = q_low:q_up
        actions = [insta_action(1,N_c,N_x,N_t,q,z) for z = x_low:0.01:x_up]
        if mod(q,N_c) == 0
            image_action_z = plot!(
                x_low:0.01:x_up,
                actions,
                color = palette(:default)[1+q-q_low],
                label = :false,
                # label = latexstring("\$ Q = $q \$"),
                alpha = 1
            )
        elseif mod(q,N_c) == 1
            image_action_z = plot!(
                x_low:0.01:x_up,
                actions,
                color = palette(:default)[1+q-q_low],
                label = :false,
                # label = latexstring("\$ Q = $q \$"),
                alpha = 1,
                # linestyle = :dash
            )
        elseif mod(q,N_c) == 2
            image_action_z = plot!(
                x_low:0.01:x_up,
                actions,
                color = palette(:default)[1+q-q_low],
                label = :false,
                # label = latexstring("\$ Q = $q \$"),
                alpha = 1,
                # linestyle = :dot,
                linewidth = 1.0
            )
        elseif mod(q,N_c) == 3
            image_action_z = plot!(
                x_low:0.01:x_up,
                actions,
                color = palette(:default)[1+q-q_low],
                label = :false,
                # label = latexstring("\$ Q = $q \$"),
                alpha = 1,
                # linestyle = :dashdot
            )
        end

        action_max = max(action_max, maximum(actions))
    end
    image_action_z = plot!(
        xlims = [x_low, x_up],
        ylims = [-y_frame, action_max+y_frame]
        # ylims = [-y_frame, 0.2]
    )
    display(image_action_z)
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\insta_actions_z_U2_zoomout.pdf")



let
    Q_lower = 0
    Q_upper = 7
    z_lower = -4
    z_upper = 4
    image_insta = plot(
        title = L"$S/\beta$ for Local Minimum Solutions of 2D U(3)",
        xlabel = L"Top. Charge $Q$",
        legend = :topleft,
        xticks = Q_lower:Q_upper,
        background_color_legend = "rgba(100%,100%,100%,60%)", 
        # foreground_color_legend = :false
        # foreground_color_legend = :false
        # ylims = [-0.1,1.0]
    )
    for z = z_lower:z_upper
        image_insta = plot!(
            Q_lower:0.01:Q_upper, 
            [insta_action(1, 3, N_x, N_t, q, z) for q = Q_lower:0.01:Q_upper],
            label = :false,
            color = palette(:default)[1+z-z_lower],
            alpha = 0.3
        )
        if (z-z_lower)%3 == 0
            image_insta = scatter!(
                Q_lower:Q_upper, 
                [insta_action(1, 3, N_x, N_t, q, z) for q = Q_lower:Q_upper],
                label = "z = $z",
                # markershape = :circle,
                markershape = :hline,
                markersize = 15,
                markerstrokewidth = 2,
                # alpha = 0.6,
                color = palette(:default)[1+z-z_lower],
            )
        elseif (z-z_lower)%3 == 1
            image_insta = scatter!(
                Q_lower:Q_upper, 
                [insta_action(1, 3, N_x, N_t, q, z) for q = Q_lower:Q_upper],
                label = "z = $z",
                # markershape = :diamond,
                markershape = :hline,
                markersize = 15,
                markerstrokewidth = 2,
                # alpha = 0.6,
                color = palette(:default)[1+z-z_lower],
            )
        elseif (z-z_lower) %3 == 2
            image_insta = scatter!(
                Q_lower:Q_upper, 
                [insta_action(1, 3, N_x, N_t, q, z) for q = Q_lower:Q_upper],
                label = "z = $z",
                # markershape = :star4,
                markershape = :hline,
                markersize = 15,
                markerstrokewidth = 2,
                # alpha = 0.6,
                color = palette(:default)[1+z-z_lower],
            )
        end
    end
    display(image_insta)
end

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\local_minima_2D_U3_spec.pdf")


let 
    N_c = 2
    q_low = 0
    q_up  = 7
    z_length = 4
    actions = zeros(2*z_length+1, length(q_low:q_up))
    for q = q_low:q_up
        z_min = floor(q/N_c)
        zs = Array(z_min-z_length:z_min+z_length)
        actions[:,1+q-q_low] = [insta_action(1,N_c,N_x,N_t,q,z) for z in zs]
        # println(q)
    end
    image_spec = plot(
        title = latexstring(" \$ S / \\beta \$ for Local Minimum Solutions of 2D U($N_c)\n The $(2*z_length+1) Lowest Actions per Sector"),
        xticks = q_low:q_up,
        xlabel = L"Top. Charge $Q$",
        legend = :false,
        background_color_legend = "rgba(100%,100%,100%,60%)"
    )
    for z = -z_length:z_length
        image_spec = scatter!(
            q_low:q_up,
            actions[1+z+z_length,:],
            markershape = :hline,
            markersize = 10,
            markerstrokewidth = 2,
            label = :false,
            # alpha = 0.5
        )
    end
    display(image_spec)
end

# for k = -5:5
#     for q = -6:6
#         bla = insta_U3_w(N_x, N_t, 3, 1);
#         S = action_U3(bla,1)
#         bla .*= ran_U3(N_x, N_t, 0.001);
#         # S = action_U3(bla,1)
#         for i = 1:100 chess_cool!(bla, 0.001, 1.0, [0], "U3") end
#         if action_U3(bla,1) < S
#             println("Problem at k = $k, q = $q: ΔS = ", action_U3(bla,1)-S)
#         end
#         # println("k = $k, q = $q, ΔS = ", action_U3(bla,1)-S)
#     end
# end

# function write_conf_U3(U, path)
#     NX = size(U,2)
#     NT = size(U,3)
#     bla = []
#     for μ = 1:2
#         for x = 1:NX
#             for t = 1:NT
#                 push!(bla, U[μ,x,t])
#             end
#         end
#     end
#     writedlm(path, bla, ',')
#     return nothing
# end

# function read_config_U3(path)
#     blub = readdlm(path, ',', ComplexF64)
#     L = Int(sqrt(size(blub,1)/2))
#     N = Int(sqrt(size(blub,2)))
#     return [reshape(blub[t+(x-1)*L+(μ-1)*L^2,:], N, N) for μ = 1:2, x = 1:L, t = 1:L] 
# end

# heatmap([real(tr(plaq(insta_U2(N_x, N_t, 2), x, t)))/2 for t = 1:N_t, x = 1:N_x])
# heatmap([real(tr(plaq(insta_U3(N_x, N_t, 0), x, t)))/3 for t = 1:N_t, x = 1:N_x])
# heatmap([real(tr(plaq(blub2, x, t)))/2 for t = 1:N_t, x = 1:N_x], 
#     xticks = vcat([1],Array(8:8:32)), 
#     yticks = vcat([1],Array(8:8:32)),
#     aspect_ratio = :equal,
#     size = (600,600), 
#     xlabel = L"$x$-coordinate",
#     ylabel = L"$t$-coordinate",
#     title = L"$\mathfrak{Re} Tr P_{xt}(n) / N$ of a cooled 2D U(2) config with $Q=2$",
#     rightmargin=8Plots.mm
# )
# heatmap([real(tr(plaq(bla1, x, t)))/3 for t = 1:N_t, x = 1:N_x], 
#     xticks = vcat([1],Array(8:8:32)), 
#     yticks = vcat([1],Array(8:8:32)),
#     aspect_ratio = :equal,
#     size = (600,600), 
#     xlabel = L"$x$-coordinate",
#     ylabel = L"$t$-coordinate",
#     title = L"$\mathfrak{Re} Tr P_{xt}(n) / N$ of a cooled 2D U(3) config with $Q=1$",
#     rightmargin=8Plots.mm
# )








