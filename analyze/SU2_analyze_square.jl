include("SU2_analyze_head.jl")




#=
for L = 32:32:128
    # L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    β   = 6.0 #N_t*N_x/128
    hot = true
    ϵ   = 0.2 
    n_stout = 0
    ρ   = 0.12
    sim_count = 3
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_path = string(base_path,"\\mean_vals.txt")
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    # means = transpose(reshape(readdlm(mean_vals_path), (num_loops,:)))
    # means_mike = transpose(reshape(readdlm(mean_vals_mike_path), (num_loops,:)))
    # num_means = size(means,1)

    sim_count_2 = sim_count +1
    base_path_2 = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count_2"
    mean_vals_path_2 = string(base_path_2,"\\mean_vals.txt")
    mean_vals_mike_path_2 = string(base_path_2,"\\mean_vals_mike.txt")
    sim_count_3 = sim_count +2
    base_path_3 = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count_3"
    mean_vals_path_3 = string(base_path_3,"\\mean_vals.txt")
    mean_vals_mike_path_3 = string(base_path_3,"\\mean_vals_mike.txt")

    means_1 = transpose(reshape(readdlm(mean_vals_path), (num_loops,:)))
    means_mike_1 = transpose(reshape(readdlm(mean_vals_mike_path), (num_loops,:)))
    num_means_1 = size(means_1,1)

    means_2 = transpose(reshape(readdlm(mean_vals_path_2), (num_loops,:)))
    means_mike_2 = transpose(reshape(readdlm(mean_vals_mike_path_2), (num_loops,:)))
    means_3 = transpose(reshape(readdlm(mean_vals_path_3), (num_loops,:)))
    means_mike_3 = transpose(reshape(readdlm(mean_vals_mike_path_3), (num_loops,:)))

    means = Matrix{Float64}(undef, 3*size(means_1,1), size(means_1,2))
    means_mike = similar(means)
    for i = 1:num_loops
        means[:,i] = vcat(means_1[:,i], means_2[:,i], means_3[:,i])
        means_mike[:,i] = vcat(means_mike_1[:,i], means_mike_2[:,i], means_mike_3[:,i])
    end


    jack_means = []
    jack_mean_errs = []
    jack_means_mike = []
    jack_mean_errs_mike = []
    for i = 1:num_loops
        b_size = Int(round(2*auto_corr_time(means[:,i]) + 1, RoundUp))    
        bla = jackknife(means[:,i], b_size)#, 500)
        push!(jack_means, bla[1])
        push!(jack_mean_errs, bla[2])

        b_size = Int(round(2*auto_corr_time(means_mike[:,i]) + 1, RoundUp))    
        bla = jackknife(means_mike[:,i], b_size)#, 500)
        push!(jack_means_mike, bla[1])
        push!(jack_mean_errs_mike, bla[2])
    end

    int_start = 13
    int_end = 16

    x_lab = string.(loops)

    image = scatter(x_lab[int_start:int_end], jack_means[int_start:int_end], yerror = jack_mean_errs[int_start:int_end], label = "Conventional", markerstrokecolor = :auto)
    image = scatter!(x_lab[int_start:int_end], jack_means_mike[int_start:int_end], yerror = jack_mean_errs_mike[int_start:int_end], label = "C.Michael, NPB 259, 58, eq.(6)", markerstrokecolor = :auto)
    image = plot!(title = "Various ⟨W(R,T)⟩ with N_t = N_x = $N_t", 
    xlabel = "[R,T] in ⟨W(R,T)⟩")
    # xticks = int_start-1:1:int_end+1)
    display(image)

    # savefig("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\mikes_loops_zoomed_N_t_$N_t.pdf")

    image_err = plot(title = "Uncertainties of Various ⟨W(R,T)⟩ 
    with N_t = N_x = $N_t",
    xlabel = "[R,T] in ⟨W(R,T)⟩",
    legend = :right)
    # xticks = int_start-1:1:int_end+1)
    image_err = scatter!(x_lab[int_start:int_end], jack_mean_errs[int_start:int_end], label = "Conventional", markerstrokecolor = :auto)
    image_err = scatter!(x_lab[int_start:int_end], jack_mean_errs_mike[int_start:int_end], label = "C.Michael, NPB 259, 58, eq.(6)", markerstrokecolor = :auto)



    display(image_err)

    # savefig("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\mikes_loops_uncer_zoomed_N_t_$N_t.pdf")
end


for L = 32:32:128
# L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    β   = 6.0 #N_t*N_x/128
    hot = true
    ϵ   = 0.2 
    n_stout = 0
    ρ   = 0.12
    sim_count = 3
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_path = string(base_path,"\\mean_vals.txt")
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    # means = transpose(reshape(readdlm(mean_vals_path), (num_loops,:)))
    # means_mike = transpose(reshape(readdlm(mean_vals_mike_path), (num_loops,:)))
    # num_means = size(means,1)

    sim_count_2 = sim_count +1
    base_path_2 = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count_2"
    mean_vals_path_2 = string(base_path_2,"\\mean_vals.txt")
    mean_vals_mike_path_2 = string(base_path_2,"\\mean_vals_mike.txt")
    sim_count_3 = sim_count +2
    base_path_3 = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count_3"
    mean_vals_path_3 = string(base_path_3,"\\mean_vals.txt")
    mean_vals_mike_path_3 = string(base_path_3,"\\mean_vals_mike.txt")


    means_1 = transpose(reshape(readdlm(mean_vals_path), (num_loops,:)))
    means_mike_1 = transpose(reshape(readdlm(mean_vals_mike_path), (num_loops,:)))
    num_means_1 = size(means_1,1)

    means_2 = transpose(reshape(readdlm(mean_vals_path_2), (num_loops,:)))
    means_mike_2 = transpose(reshape(readdlm(mean_vals_mike_path_2), (num_loops,:)))
    means_3 = transpose(reshape(readdlm(mean_vals_path_3), (num_loops,:)))
    means_mike_3 = transpose(reshape(readdlm(mean_vals_mike_path_3), (num_loops,:)))

    means = Matrix{Float64}(undef, 3*size(means_1,1), size(means_1,2))
    means_mike = similar(means)
    for i = 1:num_loops
        means[:,i] = vcat(means_1[:,i], means_2[:,i], means_3[:,i])
        means_mike[:,i] = vcat(means_mike_1[:,i], means_mike_2[:,i], means_mike_3[:,i])
    end

    num_means = size(means,1)

    creutz_means = [creutz(means[i,:], j,j+1,j+2,j+3) for i = 1:num_means, j = 1:3:Int(num_loops-3)]
    creutz_means_mike = [creutz(means_mike[i,:], j,j+1,j+2,j+3) for i = 1:num_means, j = 1:3:Int(num_loops-3)]
    num_ratios = size(creutz_means,2)

    ratio_means = []
    ratio_mean_errs = []
    ratio_means_mike = []
    ratio_mean_errs_mike = []
    for i = 1:num_ratios
        b_size = Int(round(2*auto_corr_time(creutz_means[:,i]) + 1, RoundUp))    
        bla = jackknife(creutz_means[:,i], b_size)#, 500)
        push!(ratio_means, bla[1])
        push!(ratio_mean_errs, bla[2])

        b_size = Int(round(2*auto_corr_time(creutz_means_mike[:,i]) + 1, RoundUp))    
        bla = jackknife(creutz_means_mike[:,i], b_size)#, 500)
        push!(ratio_means_mike, bla[1])
        push!(ratio_mean_errs_mike, bla[2])
    end

    int_start = 1
    int_end = num_ratios
    x_lab = string.([[1,1], [2,2], [3,3], [4,4], [5,5]])

    # image = scatter(x_lab[int_start:int_end], ratio_means[int_start:int_endnum_ratios], yerror = ratio_mean_errs[int_start:int_end], label = "conventional", markerstrokecolor = :auto)
    image = scatter(x_lab[int_start:int_end], ratio_means_mike[int_start:int_end], yerror = ratio_mean_errs_mike[int_start:int_end], label = "⟨W(R,T)⟩ via Multihit [C.Michael, NPB 259, 58, eq.(6)]", markerstrokecolor = :auto)
    image = plot!(title = "Creutz Ratios with N_t = N_x = $N_t", 
    xlabel = "R and T in {⟨W(R,T)⟩ ⟨W(R+1,T+1)⟩} / {⟨W(R+1,T)⟩ ⟨W(R,T+1)⟩}")
    # xticks = 1:num_ratios)
    display(image)

    # savefig("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\ratios_N_t_$N_t.pdf")
end



β = 6.0
first_strings = []
first_string_errs_p = []
first_string_errs_m = []
for L = 32:32:128
    # L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    # β   = 6.0 #N_t*N_x/128
    hot = true
    ϵ   = 0.2 
    n_stout = 0
    ρ   = 0.12
    sim_count = 7
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_path = string(base_path,"\\mean_vals.txt")
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    # means = transpose(reshape(readdlm(mean_vals_path), (num_loops,:)))
    # means_mike = transpose(reshape(readdlm(mean_vals_mike_path), (num_loops,:)))
    # num_means = size(means,1)

    sim_count_2 = sim_count +1
    base_path_2 = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count_2"
    mean_vals_path_2 = string(base_path_2,"\\mean_vals.txt")
    mean_vals_mike_path_2 = string(base_path_2,"\\mean_vals_mike.txt")
    sim_count_3 = sim_count +2
    base_path_3 = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count_3"
    mean_vals_path_3 = string(base_path_3,"\\mean_vals.txt")
    mean_vals_mike_path_3 = string(base_path_3,"\\mean_vals_mike.txt")


    means_1 = transpose(reshape(readdlm(mean_vals_path), (num_loops,:)))
    means_mike_1 = transpose(reshape(readdlm(mean_vals_mike_path), (num_loops,:)))
    num_means_1 = size(means_1,1)

    means_2 = transpose(reshape(readdlm(mean_vals_path_2), (num_loops,:)))
    means_mike_2 = transpose(reshape(readdlm(mean_vals_mike_path_2), (num_loops,:)))
    means_3 = transpose(reshape(readdlm(mean_vals_path_3), (num_loops,:)))
    means_mike_3 = transpose(reshape(readdlm(mean_vals_mike_path_3), (num_loops,:)))

    means = Matrix{Float64}(undef, 3*size(means_1,1), size(means_1,2))
    means_mike = similar(means)
    for i = 1:num_loops
        means[:,i] = vcat(means_1[:,i], means_2[:,i], means_3[:,i])
        means_mike[:,i] = vcat(means_mike_1[:,i], means_mike_2[:,i], means_mike_3[:,i])
    end
    
    num_means = size(means,1)

    creutz_means = [creutz(means[i,:], j,j+1,j+2,j+3) for i = 1:num_means, j = 1:3:Int(num_loops-3)]
    creutz_means_mike = [creutz(means_mike[i,:], j,j+1,j+2,j+3) for i = 1:num_means, j = 1:3:Int(num_loops-3)]
    num_ratios = size(creutz_means,2)

    ratio_means = []
    ratio_mean_errs = []
    ratio_means_mike = []
    ratio_mean_errs_mike = []
    for i = 1:num_ratios
        b_size = Int(round(2*auto_corr_time(creutz_means[:,i]) + 1, RoundUp))    
        bla = jackknife(creutz_means[:,i], b_size)#, 500)
        push!(ratio_means, bla[1])
        push!(ratio_mean_errs, bla[2])

        b_size = Int(round(2*auto_corr_time(creutz_means_mike[:,i]) + 1, RoundUp))    
        bla = jackknife(creutz_means_mike[:,i], b_size)#, 500)
        push!(ratio_means_mike, bla[1])
        push!(ratio_mean_errs_mike, bla[2])
    end

    # strings = zeros(num_ratios)
    # string_errs_p = zeros(num_ratios)
    # string_errs_m = zeros(num_ratios)
    # strings_mike = zeros(num_ratios)
    # string_errs_p_mike = zeros(num_ratios)
    # string_errs_m_mike = zeros(num_ratios)
    indices = []
    indices_mike = []
    for i = 1:num_ratios
        if ratio_means[i] - ratio_mean_errs[i] > 0.0
            push!(indices, i)
            # strings[i] = -log(ratio_means[i])
            # string_errs_p[i] = -log(ratio_means[i] + ratio_mean_errs[i])
            # string_errs_m[i] = -log(ratio_means[i] - ratio_mean_errs[i])
        end
        if ratio_means_mike[i] - ratio_mean_errs_mike[i] > 0.0
            push!(indices_mike, i)
            # strings_mike[i] = -log(ratio_means_mike[i])
            # string_errs_p_mike[i] = -log(ratio_means_mike[i] + ratio_mean_errs_mike[i])
            # string_errs_m_mike[i] = -log(ratio_means_mike[i] - ratio_mean_errs_mike[i])
        end
    end


    strings = zeros(num_ratios)
    string_errs_p = zeros(num_ratios)
    string_errs_m = zeros(num_ratios)
    strings_mike = zeros(num_ratios)
    string_errs_p_mike = zeros(num_ratios)
    string_errs_m_mike = zeros(num_ratios)

    for i in indices
        strings[i] = -log(ratio_means[i])
        string_errs_p[i] = abs(-log(ratio_means[i] + ratio_mean_errs[i]) - strings[i])
        string_errs_m[i] = abs(-log(ratio_means[i] - ratio_mean_errs[i]) - strings[i])
    end
    for i in indices_mike
        strings_mike[i] = -log(ratio_means_mike[i])
        string_errs_p_mike[i] = abs(-log(ratio_means_mike[i] + ratio_mean_errs_mike[i]) - strings[i])
        string_errs_m_mike[i] = abs(-log(ratio_means_mike[i] - ratio_mean_errs_mike[i]) - strings[i])
    end

    # int_start = 1
    # int_end = num_ratios
    x_lab = string.([[1,1], [2,2], [3,3], [4,4], [5,5]])
    indices_mike = [1,2,3]

    # image = scatter(x_lab[indices], strings[indices], yerror = (string_errs_m[indices], string_errs_p[indices]), label = "conventional", markerstrokecolor = :auto)
    image = scatter(x_lab[indices_mike], strings_mike[indices_mike], yerror = (string_errs_m_mike[indices_mike], string_errs_p_mike[indices_mike]), label = "⟨W(R,T)⟩ via Multihit [C.Michael, NPB 259, 58, eq.(6)]", markerstrokecolor = :auto)
    image = plot!(title = "String Tensions from C.-ratios with N_t = N_x = $N_t", 
    xlabel = "R and T in -log{⟨W(R,T)⟩ ⟨W(R+1,T+1)⟩} / {⟨W(R+1,T)⟩ ⟨W(R,T+1)⟩}")
    # xticks = 1:num_ratios)
    display(image)

    # savefig("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\strings_zoomed_N_t_$N_t.pdf")
    push!(first_strings, strings_mike[1])
    push!(first_string_errs_m, string_errs_m_mike[1])
    push!(first_string_errs_p, string_errs_p_mike[1])

end


string_sqrts = sqrt.(first_strings)
string_sqrt_errs_p = abs.(sqrt.(first_strings .+ first_string_errs_p) .- sqrt.(first_strings))
string_sqrt_errs_m = abs.(sqrt.(first_strings .- first_string_errs_m) .- sqrt.(first_strings))
image_sqrts = scatter(
    32:32:128, 
    string_sqrts, 
    yerror = (string_sqrt_errs_m,string_sqrt_errs_p),
    xlabel = "N_t = N_x",
    title = "Square Root of String Tension 
    using ⟨W(1,1)⟩, ..., ⟨W(2,2)⟩, β = $β",
    label = :false,
    markerstrokecolor = :auto
)

=#


for L = 32:32:128
# L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    β   = 8.0 #N_t*N_x/128
    hot = true
    ϵ   = 0.2 
    n_stout = 0
    ρ   = 0.12
    sim_count = 1
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_path = string(base_path,"\\mean_vals.txt")
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    means = readdlm(mean_vals_path)
    means_mike = readdlm(mean_vals_mike_path)
    num_means = size(means,1)
    num_loops = size(means,2)

    jack_means = []
    jack_mean_errs = []
    jack_means_mike = []
    jack_mean_errs_mike = []
    for i = 1:num_loops
        b_size = Int(round(2*auto_corr_time(means[:,i]) + 1, RoundUp))    
        bla = jackknife(means[:,i], b_size)#, 500)
        push!(jack_means, bla[1])
        push!(jack_mean_errs, bla[2])

        b_size = Int(round(2*auto_corr_time(means_mike[:,i]) + 1, RoundUp))    
        bla = jackknife(means_mike[:,i], b_size)#, 500)
        push!(jack_means_mike, bla[1])
        push!(jack_mean_errs_mike, bla[2])
    end

    int_start = 13
    int_end = 16

    x_lab = string.(loops)

    image = scatter(x_lab[int_start:int_end], jack_means[int_start:int_end], yerror = jack_mean_errs[int_start:int_end], label = "Conventional", markerstrokecolor = :auto)
    image = scatter!(x_lab[int_start:int_end], jack_means_mike[int_start:int_end], yerror = jack_mean_errs_mike[int_start:int_end], label = "C.Michael, NPB 259, 58, eq.(6)", markerstrokecolor = :auto)
    image = plot!(title = "Various ⟨W(R,T)⟩ 
    with N_t = N_x = $N_t, β = $β", 
    xlabel = "[R,T] in ⟨W(R,T)⟩")
    # xticks = int_start-1:1:int_end+1)
    display(image)

    # savefig("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\beta_$β\\mikes_loops_zoomed_beta_$β._N_t_$N_t.pdf")

    image_err = plot(title = "Uncertainties of Various ⟨W(R,T)⟩ 
    with N_t = N_x = $N_t, β = $β",
    xlabel = "[R,T] in ⟨W(R,T)⟩",
    legend = :left)
    # xticks = int_start-1:1:int_end+1)
    image_err = scatter!(x_lab[int_start:int_end], jack_mean_errs[int_start:int_end], label = "Conventional", markerstrokecolor = :auto)
    image_err = scatter!(x_lab[int_start:int_end], jack_mean_errs_mike[int_start:int_end], label = "Via [C.Michael, NPB 259, 58, eq.(6)]", markerstrokecolor = :auto)



    display(image_err)

    # savefig("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\beta_$β\\mikes_loops_uncer_zoomed_beta_$β.N_t_$N_t.pdf")
end


for L = 32:32:128
# L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    β   = 8.0 #N_t*N_x/128
    hot = true
    ϵ   = 0.2 
    n_stout = 0
    ρ   = 0.12
    sim_count = 1
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_path = string(base_path,"\\mean_vals.txt")
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    means = readdlm(mean_vals_path)
    means_mike = readdlm(mean_vals_mike_path)
    num_means = size(means,1)
    num_loops = size(means,2)

    creutz_means = [creutz(means[i,:], j,j+1,j+2,j+3) for i = 1:num_means, j = 1:3:Int(num_loops-3)]
    creutz_means_mike = [creutz(means_mike[i,:], j,j+1,j+2,j+3) for i = 1:num_means, j = 1:3:Int(num_loops-3)]
    num_ratios = size(creutz_means,2)

    ratio_means = []
    ratio_mean_errs = []
    ratio_means_mike = []
    ratio_mean_errs_mike = []
    for i = 1:num_ratios
        b_size = Int(round(2*auto_corr_time(creutz_means[:,i]) + 1, RoundUp))    
        bla = jackknife(creutz_means[:,i], b_size)#, 500)
        push!(ratio_means, bla[1])
        push!(ratio_mean_errs, bla[2])

        b_size = Int(round(2*auto_corr_time(creutz_means_mike[:,i]) + 1, RoundUp))    
        bla = jackknife(creutz_means_mike[:,i], b_size)#, 500)
        push!(ratio_means_mike, bla[1])
        push!(ratio_mean_errs_mike, bla[2])
    end

    int_start = 1
    int_end = num_ratios
    x_lab = string.([[1,1], [2,2], [3,3], [4,4], [5,5]])

    # image = scatter(x_lab[int_start:int_end], ratio_means[int_start:int_endnum_ratios], yerror = ratio_mean_errs[int_start:int_end], label = "conventional", markerstrokecolor = :auto)
    image = scatter(x_lab[int_start:int_end], ratio_means_mike[int_start:int_end], yerror = ratio_mean_errs_mike[int_start:int_end], label = "⟨W(R,T)⟩ via Multihit [C.Michael, NPB 259, 58, eq.(6)]", markerstrokecolor = :auto)
    image = plot!(title = "Creutz Ratios 
    with N_t = N_x = $N_t, β = $β", 
    xlabel = "R and T in {⟨W(R,T)⟩ ⟨W(R+1,T+1)⟩} / {⟨W(R+1,T)⟩ ⟨W(R,T+1)⟩}")
    # xticks = 1:num_ratios)
    display(image)

    # savefig("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\beta_$β\\ratios_beta_$β._N_t_$N_t.pdf")
end




β = 8.0
first_strings = []
first_string_errs_p = []
first_string_errs_m = []
for L = 32:32:128
    # L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    # β   = 4.0 #N_t*N_x/128
    hot = true
    ϵ   = 0.2 
    n_stout = 0
    ρ   = 0.12
    sim_count = 1 # 7 for β = 6.0
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Physik Uni\\julia_projects\\SU2\\data\\N_t_$N_t.N_x_$N_x._beta_$β._eps_$ϵ\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_path = string(base_path,"\\mean_vals.txt")
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    means = readdlm(mean_vals_path)
    means_mike = readdlm(mean_vals_mike_path)
    num_means = size(means,1)
    num_loops = size(means,2)

    creutz_means = [creutz(means[i,:], j,j+1,j+2,j+3) for i = 1:num_means, j = 1:3:Int(num_loops-3)]
    creutz_means_mike = [creutz(means_mike[i,:], j,j+1,j+2,j+3) for i = 1:num_means, j = 1:3:Int(num_loops-3)]
    num_ratios = size(creutz_means,2)

    ratio_means = []
    ratio_mean_errs = []
    ratio_means_mike = []
    ratio_mean_errs_mike = []
    for i = 1:num_ratios
        b_size = Int(round(2*auto_corr_time(creutz_means[:,i]) + 1, RoundUp))    
        bla = jackknife(creutz_means[:,i], b_size)#, 500)
        push!(ratio_means, bla[1])
        push!(ratio_mean_errs, bla[2])

        b_size = Int(round(2*auto_corr_time(creutz_means_mike[:,i]) + 1, RoundUp))    
        bla = jackknife(creutz_means_mike[:,i], b_size)#, 500)
        push!(ratio_means_mike, bla[1])
        push!(ratio_mean_errs_mike, bla[2])
    end

    # strings = zeros(num_ratios)
    # string_errs_p = zeros(num_ratios)
    # string_errs_m = zeros(num_ratios)
    # strings_mike = zeros(num_ratios)
    # string_errs_p_mike = zeros(num_ratios)
    # string_errs_m_mike = zeros(num_ratios)
    indices = []
    indices_mike = []
    for i = 1:num_ratios
        if ratio_means[i] - ratio_mean_errs[i] > 0.0
            push!(indices, i)
            # strings[i] = -log(ratio_means[i])
            # string_errs_p[i] = -log(ratio_means[i] + ratio_mean_errs[i])
            # string_errs_m[i] = -log(ratio_means[i] - ratio_mean_errs[i])
        end
        if ratio_means_mike[i] - ratio_mean_errs_mike[i] > 0.0
            push!(indices_mike, i)
            # strings_mike[i] = -log(ratio_means_mike[i])
            # string_errs_p_mike[i] = -log(ratio_means_mike[i] + ratio_mean_errs_mike[i])
            # string_errs_m_mike[i] = -log(ratio_means_mike[i] - ratio_mean_errs_mike[i])
        end
    end


    strings = zeros(num_ratios)
    string_errs_p = zeros(num_ratios)
    string_errs_m = zeros(num_ratios)
    strings_mike = zeros(num_ratios)
    string_errs_p_mike = zeros(num_ratios)
    string_errs_m_mike = zeros(num_ratios)

    for i in indices
        strings[i] = -log(ratio_means[i])
        string_errs_p[i] = abs(-log(ratio_means[i] + ratio_mean_errs[i]) - strings[i])
        string_errs_m[i] = abs(-log(ratio_means[i] - ratio_mean_errs[i]) - strings[i])
    end
    for i in indices_mike
        strings_mike[i] = -log(ratio_means_mike[i])
        string_errs_p_mike[i] = abs(-log(ratio_means_mike[i] + ratio_mean_errs_mike[i]) - strings_mike[i])
        string_errs_m_mike[i] = abs(-log(ratio_means_mike[i] - ratio_mean_errs_mike[i]) - strings_mike[i])
    end

    # int_start = 1
    # int_end = num_ratios
    x_lab = string.([[1,1], [2,2], [3,3], [4,4], [5,5]])
    indices_mike = [1,2,3]

    # image = scatter(x_lab[indices], strings[indices], yerror = (string_errs_m[indices], string_errs_p[indices]), label = "conventional", markerstrokecolor = :auto)
    image = scatter(x_lab[indices_mike], strings_mike[indices_mike], yerror = (string_errs_m_mike[indices_mike], string_errs_p_mike[indices_mike]), label = "⟨W(R,T)⟩ via Multihit [C.Michael, NPB 259, 58, eq.(6)]", markerstrokecolor = :auto)
    image = plot!(title = "String Tensions from C.-ratios 
    with N_t = N_x = $N_t, β = $β", 
    xlabel = "R and T in -log{⟨W(R,T)⟩ ⟨W(R+1,T+1)⟩} / {⟨W(R+1,T)⟩ ⟨W(R,T+1)⟩}")
    # xticks = 1:num_ratios)
    display(image)

    # savefig("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\beta_$β\\strings_zoomed_β_$β._N_t_$N_t.pdf")
    push!(first_strings, strings_mike[1])
    push!(first_string_errs_m, string_errs_m_mike[1])
    push!(first_string_errs_p, string_errs_p_mike[1])

end


string_sqrts = sqrt.(first_strings)
string_sqrt_errs_p = abs.(sqrt.(first_strings .+ first_string_errs_p) .- sqrt.(first_strings))
string_sqrt_errs_m = abs.(sqrt.(first_strings .- first_string_errs_m) .- sqrt.(first_strings))
image_sqrts = scatter(
    32:32:128, 
    xticks = 32:32:128,
    string_sqrts, 
    yerror = (string_sqrt_errs_m,string_sqrt_errs_p),
    xlabel = "N_t = N_x",
    title = "Square Root of the String Tension
    using ⟨W(1,1)⟩, ..., ⟨W(2,2)⟩, β = $β",
    label = :false,
    markerstrokecolor = :auto
)



# writedlm("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\beta_$β\\string_sqrts.txt", string_sqrts)
# writedlm("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\beta_$β\\string_sqrt_errs_m.txt", string_sqrt_errs_m)
# writedlm("C:\\Physik Uni\\julia_projects\\SU2\\data\\creutz_ratios\\beta_$β\\string_sqrt_errs_p.txt", string_sqrt_errs_p) 











