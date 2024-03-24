include("SU2_analyze_head.jl")
using QuadGK
# using ADerrors



# Note: test the upper bound of error E, quadgk(f,a,b) = (I,E)
function analytic_plaq_U2(β)
    numer(α) = besseli(0,β*cos(α)) + besseli(2,β*cos(α))
    denom(α) = 2*besseli(1,β*cos(α))/cos(α)
    return quadgk(numer,0,π)[1]/quadgk(denom,0,π)[1] - 1/β
end

function analytic_susc_U2(β)
    nasty(α)   = besseli(1,β*cos(α))/cos(α)
    nastier(α) = α^2 * besseli(1,β*cos(α))/cos(α)
    return quadgk(nastier,-π/2,π/2)[1] / quadgk(nasty,-π/2,π/2)[1] / π^2
end

# For 2D U(1) theory
# function test_chi(β)
#     fun(ϕ) = ϕ^2*exp(β*cos(ϕ))
#     return quadgk(fun,-π,π)[1] / (2π)^3 / besseli(0,β)
# end

# test_chi(10)*(10*80)

function creutz(means::Vector, a, b, c, d)
    return means[a]*means[d]/(means[b]*means[c])
end

function analytic_plaq(β)
    return (besseli(0,β) + besseli(2,β)) / (2*besseli(1,β)) - 1/β
end

function string_jack(means, b_size, a, b, c, d)
    N_blocks = Int(div(size(means,1), b_size, RoundDown))
    ratios = [creutz(means[i, :], a, b, c, d) for i = 1:N_blocks*b_size]
    string_mean = -log(mean(ratios))

    first_string = -log(mean(ratios[b_size+1:b_size*N_blocks]))
    last_string = -log(mean(ratios[1:b_size*(N_blocks-1)]))
    jack_strings = [first_string, last_string]
    for i = 2:N_blocks-1
        push!(jack_strings, -log(mean(vcat(ratios[1:(i-1)*b_size], ratios[i*b_size+1:b_size*N_blocks]))))
    end
    σ = sqrt((N_blocks-1) * mean((jack_strings .- string_mean).^2) )

    return [string_mean, σ]
end

const sp_fac = 0.49/sqrt(1.65)

function space_jack(means, b_size, a, b, c, d)
    N_blocks = Int(div(size(means,1), b_size, RoundDown))
    ratios = [creutz(means[i, :], a, b, c, d) for i = 1:N_blocks*b_size]
    space_mean = sp_fac * sqrt(-log(mean(ratios)))

    first_space = sp_fac * sqrt(-log(mean(ratios[b_size+1:b_size*N_blocks])))
    last_space = sp_fac * sqrt(-log(mean(ratios[1:b_size*(N_blocks-1)])))
    jack_spaces = [first_space, last_space]
    for i = 2:N_blocks-1
        push!(jack_spaces, sp_fac * sqrt(-log(mean(vcat(ratios[1:(i-1)*b_size], ratios[i*b_size+1:b_size*N_blocks])))))
    end
    σ = sqrt((N_blocks-1) * mean((jack_spaces .- space_mean).^2) )

    return [space_mean, σ]
end





L = 32
β = 2.0
N_t = L #+ i*16
N_x = L #+ i*16
hot = true
# ϵ   = 0.2 
n_stout = 0
ρ   = 0.12
sim_count = 1
loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
num_loops = length(loops)

# base_base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ" 
# base_path = string(base_base_path, "\\sim_count_$sim_count")
base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
mean_vals_path = string(base_path,"\\mean_vals.txt")

means = readdlm(mean_vals_path)
# means = real.(means)
num_means = size(means,1)
num_loops = size(means,2)

jack_means = []
jack_mean_errs = []
b_sizes = []
for i = 1:num_loops
    b_size = Int(round(2*auto_corr_time(means[:,i]) + 1, RoundUp))    
    bla = jackknife(means[:,i], b_size)#, 500)
    push!(jack_means, bla[1])
    push!(jack_mean_errs, bla[2])
    push!(b_sizes, b_size)
end

all_jack_means = jack_means
all_jack_errs = jack_mean_errs

x_lab = string.(loops)

# image = scatter(x_lab[int_start:int_end], jack_means[int_start:int_end], yerror = jack_mean_errs[int_start:int_end], label = "Conventional", markerstrokecolor = :auto)
image = scatter(
    x_lab[1:16], 
    all_jack_means[1:16], 
    yerror = all_jack_errs[1:16], 
    label = "Square Lat.",
    # colors = palette(:default)[1],
    markerstrokecolor = :auto, 
    # markershape = :square,
    legend = :top,
    foreground_color_legend = nothing,
    background_color_legend = nothing)
image = plot!(title = "Various Wilson Loops, 
Group: U(2), β = $β, N_t = N_x = $L", 
xlabel = "[R,T] in ⟨W(R,T)⟩ or name of loop")
# xticks = int_start-1:1:int_end+1)
display(image)

# bli = all_jack_means




L = 32
# for β in [12.0] # [2.0,4.0,6.0,8.0]
β = 12.0
N_t = L #+ i*16
N_x = L #+ i*16
hot = true
# ϵ   = 0.2 
n_stout = 0
ρ   = 0.12
sim_count = 16
loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
num_loops = length(loops)

# base_base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ" 
# base_path = string(base_base_path, "\\sim_count_$sim_count")
base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
Q_path = string(base_path, "\\top_charge.txt")
acc_path = string(base_path, "\\acceptances.txt")

charges = readdlm(Q_path)
# Qs = copy(charges)
susc = charges.^2 ./ L^2
# susc = Qs.^2 ./ L^2

b_size_Q = Int(round(2*auto_corr_time(charges)+1, RoundUp))
Q_mean, Q_err = round.(jackknife(charges, b_size_Q), digits = 4)
println("β: ", β, ", block size for Q:  ",  b_size_Q)

b_size = Int(round(2*auto_corr_time(susc)+1, RoundUp))
χ_mean, χ_err = round.(jackknife(susc, b_size), digits = 4)
χ_anal = round.(analytic_susc_U2(β), digits = 4)
println("        block size for Q²: ",  b_size)

skew = round(skewness(charges), digits = 3)
kurt = round(kurtosis(charges), digits = 3)

# bin_width = 2
image_Q = histogram(
    round.(Int,charges), 
    # round.(Int,Qs), 
    label = "⟨Q⟩ = $Q_mean ± $Q_err \n χ = $χ_mean ± $χ_err \n χ_anal = $χ_anal \n skew = $skew \n kurt = $kurt",
    # bins = 150,
    # bins = minimum(charges)-bin_width/2:bin_width:maximum(charges)-bin_width/2,
    normalize = :true,
    legend = :topleft,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    title = "Topological charge in SQUARE 2D U(2)
β = $β, N_t = N_x = $L, insta-update")
display(image_Q)
plot(charges)
# plot(Qs)

println(" ")

# Qs = vcat(Qs, charges)

acc = readdlm(acc_path)
last(acc[:,3])/32^2/2/6400


# Q_uw=uwreal(vec(charges),"daten")
# uwerr(Q_uw)
# Q_uw

# susc_uw=uwreal(vec(charges.^2)./L^2,"andere_daten")
# uwerr(susc_uw)
# susc_uw

# P_uw = uwreal(vec(means[:,1]),"daten")
# uwerr(P_uw)
# P_uw








L = 32
for β in [2.0,4.0,6.0,8.0]
N_t = 2*L #+ i*16
N_x = L #+ i*16
hot = true
# ϵ   = 0.2 
n_stout = 0
ρ   = 0.12
sim_count = 2
loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
num_loops = length(loops)

# base_base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ" 
# base_path = string(base_base_path, "\\sim_count_$sim_count")
base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\U2_data\\hex_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
Q_path = string(base_path, "\\top_charge.txt")

charges = readdlm(Q_path)
susc = charges.^2 ./ L^2

b_size_Q = Int(round(2*auto_corr_time(charges)+1, RoundUp))
Q_mean, Q_err = round.(jackknife(charges, b_size_Q), digits = 4)
println("β: ", β, ", block size for Q:  ",  b_size_Q)

b_size = Int(round(2*auto_corr_time(susc)+1, RoundUp))
χ_mean, χ_err = round.(jackknife(susc, b_size), digits = 4)
χ_anal = round.(analytic_susc_U2(β), digits = 4)
println("        block size for Q²: ",  b_size)

skew = round(skewness(charges), digits = 3)
kurt = round(kurtosis(charges), digits = 3)

bin_width = 2
hist_charges = round.(Int,charges)
image_Q = histogram(
    hist_charges, 
    label = "⟨Q⟩ = $Q_mean ± $Q_err \n χ = $χ_mean ± $χ_err \n χ_anal = $χ_anal \n skew = $skew \n kurt = $kurt",
    # bins = 150,
    # bins = minimum(charges)-bin_width/2:bin_width:maximum(charges)-bin_width/2,
    normalize = :true,
    legend = :topleft,
    foreground_color_legend = nothing,
    background_color_legend = nothing,
    title = "Topological charge in HEX. 2D U(2)
β = $β, N_t = 2⋅N_x = 2⋅$L, insta-update")
display(image_Q)

end

std(charges)



Q_uw=uwreal(vec(charges),"daten")
uwerr(Q_uw)
Q_uw
