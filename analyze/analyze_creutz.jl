include("SU2_analyze_head.jl")
include("SU2_jackknives.jl")

function creutz(means::Vector, a, b, c, d)
    return means[a]*means[d]/(means[b]*means[c])
end

function analytic_plaq(β)
    return (besseli(0,β) + besseli(2,β)) / (2*besseli(1,β)) - 1/β
end

function analytic_plaq_approx(β)
    b = big(β)
    return -1/b + (2 - 14/(8*b)) / (2 * (1-3/(8*b)))
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





for L = 32:32:128
    swil_sq_str = []
    swil_hex_str = []
    swil_anal_str = []
    plaqs_sq = []
    plaqs_sq_err = []
    plaqs_hex = []
    plaqs_hex_err = []
    plaqs_anal = []
    for β in [2.0, 4.0, 6.0, 8.0]
        # L = 32
        # β   = 2.0 #N_t*N_x/128
        N_t = L #+ i*16
        N_x = L #+ i*16
        hot = true
        # ϵ   = 0.2 
        n_stout = 0
        ρ   = 0.12
        sim_count = 4
        loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4]]#, [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
        num_loops = length(loops)

        # base_base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ" 
        # base_path = string(base_base_path, "\\sim_count_$sim_count")
        base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
        # mean_vals_path = string(base_path,"\\mean_vals.txt")
        mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")
        L_loop_means_path = string(base_path,"\\L_loop_means.txt")
        edge_loop_means_path = string(base_path,"\\edge_loop_means.txt")
        rhomb_means_path = string(base_path,"\\rhomb_means.txt")
        half_rhomb_means_path = string(base_path,"\\half_rhomb_means.txt")


        # means = readdlm(mean_vals_path)
        means_mike = readdlm(mean_vals_mike_path)[:,1:num_loops]
        num_means = size(means_mike,1)
        num_loops = size(means_mike,2)

        

        jack_means_mike = []
        jack_mean_errs_mike = []
        for i = 1:num_loops
            # b_size = Int(round(2*auto_corr_time(means[:,i]) + 1, RoundUp))    
            # bla = jackknife(means[:,i], b_size)#, 500)
            # push!(jack_means, bla[1])
            # push!(jack_mean_errs, bla[2])

            b_size = Int(round(2*auto_corr_time(means_mike[:,i]) + 1, RoundUp))    
            bla = jackknife(means_mike[:,i], b_size)#, 500)
            push!(jack_means_mike, bla[1])
            push!(jack_mean_errs_mike, bla[2])
        end

        edge_loop_means = readdlm(edge_loop_means_path)
        edge_b_size = Int(round(2*auto_corr_time(edge_loop_means) + 1, RoundUp))
        edge_jack = jackknife(edge_loop_means,edge_b_size)
        edge_mean = edge_jack[1]
        edge_err = edge_jack[2]

        L_loop_means = readdlm(L_loop_means_path)
        L_b_size = Int(round(2*auto_corr_time(L_loop_means) + 1, RoundUp))
        L_jack = jackknife(L_loop_means,L_b_size)
        L_mean = L_jack[1]
        L_err = L_jack[2]

        half_rhomb_means = readdlm(half_rhomb_means_path)
        half_rhomb_b_size = Int(round(2*auto_corr_time(half_rhomb_means) + 1, RoundUp))
        half_rhomb_jack = jackknife(half_rhomb_means,half_rhomb_b_size)
        half_rhomb_mean = half_rhomb_jack[1]
        half_rhomb_err = half_rhomb_jack[2]

        rhomb_means = readdlm(rhomb_means_path)
        rhomb_b_size = Int(round(2*auto_corr_time(rhomb_means) + 1, RoundUp))
        rhomb_jack = jackknife(rhomb_means,rhomb_b_size)
        rhomb_mean = rhomb_jack[1]
        rhomb_err = rhomb_jack[2]

        all_jack_means = vcat(jack_means_mike, edge_mean, L_mean, half_rhomb_mean, rhomb_mean)
        all_jack_errs = vcat(jack_mean_errs_mike, edge_err, L_err, half_rhomb_err, rhomb_err)

        x_lab = string.(loops)
        push!(x_lab, "edge", "L", "h-rh.", "rh.")

        # image = scatter(x_lab[int_start:int_end], jack_means[int_start:int_end], yerror = jack_mean_errs[int_start:int_end], label = "Conventional", markerstrokecolor = :auto)
        image = scatter(
            x_lab[1:12], 
            all_jack_means[1:12], 
            yerror = all_jack_errs[1:12], 
            label = "Square, one-link int.", # \n (no h-rh. or rh.)",
            color = cb_purple,# palette(:default)[1],
            markerstrokecolor = cb_purple, # :auto,  
            markershape = :square,
            markersize = 6,
            legend = :top,
            # foreground_color_legend = nothing,
            background_color_legend = nothing,
            xlabel = latexstring("\$[R,T]\$ or name of loop"),
            ylabel = latexstring("\$\\langle W(R,T) \\rangle\$"),
            tickfontsize = 9,
            labelfontsize = 14,
            legendfontsize = 10,
            )
        # image = plot!(title = "Various Wilson Loops, β = $β
        # Square: N_t = N_x = $L, Hex: N_t = 2⋅N_x = 2⋅$L", 
        # xticks = int_start-1:1:int_end+1)
        # display(image)



        # ⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔⬣⎔



        sim_count = 3
        N_t = 2*L
        base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\hex_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
        mean_vals_path = string(base_path,"\\mean_vals.txt")
        L_loop_means_path = string(base_path,"\\L_loop_means.txt")
        edge_loop_means_path = string(base_path,"\\edge_loop_means.txt")
        rhomb_half_means_path = string(base_path, "\\rhomb_half_loop_means.txt")
        rhomb_means_path = string(base_path, "\\rhomb_loop_means.txt")

        means = readdlm(mean_vals_path)
        jack_means = []
        jack_mean_errs = []
        for i = 1:num_loops
            b_size = Int(round(2*auto_corr_time(means[:,i]) + 1, RoundUp))    
            bla = jackknife(means[:,i], b_size)#, 500)
            push!(jack_means, bla[1])
            push!(jack_mean_errs, bla[2])
        end

        edge_loop_means = 2 .* readdlm(edge_loop_means_path)
        edge_b_size = Int(round(2*auto_corr_time(edge_loop_means) + 1, RoundUp))
        edge_jack = jackknife(edge_loop_means,edge_b_size)
        edge_mean_hex = edge_jack[1]
        edge_err_hex = edge_jack[2]

        L_loop_means = 2 .* readdlm(L_loop_means_path)
        L_b_size = Int(round(2*auto_corr_time(L_loop_means) + 1, RoundUp))
        L_jack = jackknife(L_loop_means,L_b_size)
        L_mean_hex = L_jack[1]
        L_err_hex = L_jack[2]

        rhomb_half_loop_means = 2 .* readdlm(rhomb_half_means_path)
        rhomb_half_b_size = Int(round(2*auto_corr_time(rhomb_half_loop_means) + 1, RoundUp))
        rhomb_half_jack = jackknife(rhomb_half_loop_means, rhomb_half_b_size)
        rhomb_half_mean = rhomb_half_jack[1]
        rhomb_half_err = rhomb_half_jack[2]

        rhomb_loop_means = 2 .* readdlm(rhomb_means_path)
        rhomb_b_size = Int(round(2*auto_corr_time(rhomb_loop_means) + 1, RoundUp))
        rhomb_jack = jackknife(rhomb_loop_means,rhomb_b_size)
        rhomb_mean = rhomb_jack[1]
        rhomb_err = rhomb_jack[2]

        # push!(x_lab, "h-rh.", "rh.")

        all_jack_means_hex = vcat(jack_means, edge_mean_hex, L_mean_hex, rhomb_half_mean, rhomb_mean)
        all_jack_errs_hex = vcat(jack_mean_errs, edge_err_hex, L_err_hex, rhomb_half_err, rhomb_err)

        image = scatter!(
            x_lab[1:12], 
            all_jack_means_hex[1:12], 
            yerror = all_jack_errs_hex[1:12], 
            label = "Hexagonal", 
            color = cb_pink; # palette(:default)[2], 
            markerstrokecolor = cb_pink, # :auto,
            markershape = :hexagon,
            markersize = 6,
            )





            # sim_count = 2
            # N_t = L
            # base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
            # mean_vals_path = string(base_path,"\\mean_vals.txt")
            # L_loop_means_path = string(base_path,"\\L_loop_means.txt")
            # edge_loop_means_path = string(base_path,"\\edge_loop_means.txt")
            # rhomb_half_means_path = string(base_path, "\\rhomb_half_loop_means.txt")
            # rhomb_means_path = string(base_path, "\\rhomb_loop_means.txt")
        
            # means = readdlm(mean_vals_path)
            # jack_means = []
            # jack_mean_errs = []
            # for i = 1:num_loops
            #     b_size = Int(round(2*auto_corr_time(means[:,i]) + 1, RoundUp))    
            #     bla = jackknife(means[:,i], b_size)#, 500)
            #     push!(jack_means, bla[1])
            #     push!(jack_mean_errs, bla[2])
            # end
        
            # edge_loop_means = readdlm(edge_loop_means_path)
            # edge_b_size = Int(round(2*auto_corr_time(edge_loop_means) + 1, RoundUp))
            # edge_jack = jackknife(edge_loop_means,edge_b_size)
            # edge_mean_hex = edge_jack[1]
            # edge_err_hex = edge_jack[2]
        
            # L_loop_means = readdlm(L_loop_means_path)
            # L_b_size = Int(round(2*auto_corr_time(L_loop_means) + 1, RoundUp))
            # L_jack = jackknife(L_loop_means,L_b_size)
            # L_mean_hex = L_jack[1]
            # L_err_hex = L_jack[2]
        
            # rhomb_half_loop_means = readdlm(rhomb_half_means_path)
            # rhomb_half_b_size = Int(round(2*auto_corr_time(rhomb_half_loop_means) + 1, RoundUp))
            # rhomb_half_jack = jackknife(rhomb_half_loop_means, rhomb_half_b_size)
            # rhomb_half_mean = rhomb_half_jack[1]
            # rhomb_half_err = rhomb_half_jack[2]
        
            # rhomb_loop_means = readdlm(rhomb_means_path)
            # rhomb_b_size = Int(round(2*auto_corr_time(rhomb_loop_means) + 1, RoundUp))
            # rhomb_jack = jackknife(rhomb_loop_means,rhomb_b_size)
            # rhomb_mean = rhomb_jack[1]
            # rhomb_err = rhomb_jack[2]
            
            # all_jack_means_hex_2 = vcat(jack_means, edge_mean_hex, L_mean_hex, rhomb_half_mean, rhomb_mean)
            # all_jack_errs_hex_2 = vcat(jack_mean_errs, edge_err_hex, L_err_hex, rhomb_half_err, rhomb_err)
        
            # image = scatter!(
            #     x_lab, 
            #     all_jack_means_hex_2, 
            #     yerror = all_jack_errs_hex_2, 
            #     label = "Hexagonal, N_x : N_t = 1:1", 
            #     color = palette(:default)[3], 
            #     markerstrokecolor = :auto,
            #     markershape = :hexagon)





        

        loop_sizes = [1, 2, 2, 4, 6, 6, 9, 12, 12, 16, 3, 4, 3, 4]
        jack_means_anal = [2*analytic_plaq(β)^l for l in loop_sizes]

        image = scatter!(
            x_lab[1:12], 
            jack_means_anal[1:12], 
            label = "Analytic", 
            color = :black,# :red, # palette(:default)[3], 
            # markerstrokecolor = :auto,
            markershape = :hline, 
            markersize = 10
            )


        display(image)
        # mkpath("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_loops\\beta_$β")
        # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_loops\\beta_$β\\compare_loops_many_beta_$β._L_$L.pdf")


        #=
        int_zoom_1 = collect(5:10)

        image_zoom_1 = scatter(
            x_lab[int_zoom_1],
            all_jack_means[int_zoom_1], 
            yerror = all_jack_errs[int_zoom_1], 
            label = "Square",
            markerstrokecolor = :auto,
            markershape = :diamond)
            # legend = :topright)
        image_zoom_1 = plot!(title = "Various Wilson Loops, β = $β
        Square: N_t = N_x = $L,   Hex: N_t = 2⋅N_x = 2⋅$L", 
        xlabel = "[R,T] in ⟨W(R,T)⟩ or name of loop",
        foreground_color_legend = nothing,
        background_color_legend = nothing)

        image_1 = scatter!(
            x_lab[int_zoom_1],
            all_jack_means_hex[int_zoom_1], 
            yerror = all_jack_errs_hex[int_zoom_1], 
            label = "Hexagonal",
            markerstrokecolor = :auto,
            markershape = :hexagon)

        image_1 = scatter!(
            x_lab[int_zoom_1],
            jack_means_anal[int_zoom_1], 
            # yerror = all_jack_errs_hex[int_zoom_1], 
            label = "Analytical",
            # markerstrokecolor = :auto,
            color = :red,
            markershape = :hline, 
            markersize = 10)

        display(image_zoom_1)

        # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_loops\\beta_$β\\compare_loops_zoom_1_beta_$β._L_$L.pdf")


        int_zoom_2 = [4,5,6,11,12]
        image_2 = scatter(
            x_lab[int_zoom_2],
            all_jack_means[int_zoom_2], 
            yerror = all_jack_errs[int_zoom_2], 
            label = "Square \n (no h-rh. or rh.)",
            markerstrokecolor = :auto,
            markershape = :diamond)
            # legend = :topright)
        image_2 = plot!(title = "Various Wilson Loops, β = $β
        Square: N_t = N_x = $L,   Hex: N_t = 2⋅N_x = 2⋅$L", 
        xlabel = "[R,T] in ⟨W(R,T)⟩ or name of loop",
        foreground_color_legend = nothing,
        background_color_legend = nothing)

        int_zoom_2 = [4,5,6,11,12,13,14]
        image_2 = scatter!(
            x_lab[int_zoom_2],
            all_jack_means_hex[int_zoom_2], 
            yerror = all_jack_errs_hex[int_zoom_2], 
            label = "Hexagonal",
            markerstrokecolor = :auto,
            markershape = :hexagon)

        image_2 = scatter!(
            x_lab[int_zoom_2],
            jack_means_anal[int_zoom_2], 
            # yerror = all_jack_errs_hex[int_zoom_1], 
            label = "Analytical",
            # markerstrokecolor = :auto,
            color = :red,
            markershape = :hline, 
            markersize = 10)

        display(image_2)
        # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_loops\\beta_$β\\compare_loops_zoom_2_beta_$β._L_$L.pdf")


        # all_jack_means
        # all_jack_errs
        # all_jack_means_hex
        # all_jack_errs_hex

        println("β = $β")
        println("[1,1]:      ", all_jack_means[1], " ± ", all_jack_errs[1])
        println("[1,1]_hex:  ", all_jack_means_hex[1], " ± ", all_jack_errs_hex[1])
        println("[2,2]:      ", all_jack_means[4], " ± ", all_jack_errs[4])
        println("[2,2]_hex:  ", all_jack_means_hex[4], " ± ", all_jack_errs_hex[4])
        println("[2,3]:      ", all_jack_means[5], " ± ", all_jack_errs[5])
        println("[2,3]_hex:  ", all_jack_means_hex[5], " ± ", all_jack_errs_hex[5])
        println("[3,3]:      ", all_jack_means[7], " ± ", all_jack_errs[7])
        println("[3,3]_hex:  ", all_jack_means_hex[7], " ± ", all_jack_errs_hex[7])
        println("h-rhomb:    ", all_jack_means[13], " ± ", all_jack_errs[13])
        println("h-rhomb_hex:", all_jack_means_hex[13], " ± ", all_jack_errs_hex[13])
        println("rhomb:      ", all_jack_means[14], " ± ", all_jack_errs[14])
        println("rhomb_hex:  ", all_jack_means_hex[14], " ± ", all_jack_errs_hex[14])
        println(" ")
        =#

        
        fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\plots\\WRT_square_hex_anal\\WRT_sha_beta_$(β)_L_$L.pdf"
        # savefig(fig_path)
        tab_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\tabellen\\WRT_square_hex_anal\\WRT_sha_beta_$(β)_L_$L.txt"
        err_digs  = 2
        anal_digs_pre  = 1
        anal_digs_post = 7
        bla = open(tab_path, "w")
        write(bla, "\\begin{table}[H]\n\t\\centering\n\t\\hline\n\t\\begin{tabular}{|c|c|c|c|}\n")
        write(bla, "\t\t loop &\t square &\t hex &\t analyt. \t \\\\\\hline\n")
        for i = 1:12
            sq_dat   = format_x_err(all_jack_means[i], all_jack_errs[i], err_digs)
            hex_dat  = format_x_err(all_jack_means_hex[i], all_jack_errs_hex[i], err_digs)
            # anal_dat = round(jack_means_anal[i], digits = anal_digs)
            anal_dat = @sprintf("%*.*f", anal_digs_pre, anal_digs_post, jack_means_anal[i])
            write(bla, "\t\t $(x_lab[i]) &\t $sq_dat &\t $hex_dat &\t $anal_dat\t \\\\\\hline\n")
            if i == 1
                sq_swil   = format_x_err((1-all_jack_means[i]/2)*β/4, all_jack_errs[i]/2*β/4, err_digs)
                hex_swil  = format_x_err((1-all_jack_means_hex[i]/2)*β/4, all_jack_errs_hex[i]/2*β/4, err_digs)
                anal_swil = @sprintf("%*.*f", anal_digs_pre, anal_digs_post, (1-jack_means_anal[i]/2)*β/4)
                push!(swil_sq_str, sq_swil)
                push!(swil_hex_str, hex_swil)
                push!(swil_anal_str, anal_swil)
                push!(plaqs_sq, all_jack_means[i])
                push!(plaqs_sq_err, all_jack_errs[i])
                push!(plaqs_hex, all_jack_means_hex[i])
                push!(plaqs_hex_err, all_jack_errs_hex[i])
                push!(plaqs_anal, jack_means_anal[i])
            end
        end
        write(bla,"\t\\end{tabular}\n\t\\caption{Caption}\n\t\\label{tab:my_label}\n\\end{table}")
        close(bla)
    end

    betas = Vector(2.0:2.0:8.0)
    beta_range_anal = 1.0:0.1:700
    slop = ((1-analytic_plaq(700))*700/4 - (1-analytic_plaq(699))*699/4) / (1/700-1/699)
    image_plaq_anal = plot(
        1 ./ beta_range_anal,
        # (1 .- analytic_plaq.(beta_range_anal)) .* beta_range_anal ./4 ,
        [(1-analytic_plaq(β))*β/4 for β in beta_range_anal],
        label = :false,
        color = :black
    )
    image_plaq_anal = plot!(
        0:0.05:0.3,
        [3/8 + slop*x for x in 0:0.05:0.3],
        label = :false, # "Slope at cont. value",
        color = cb_red
    )
    image_plaq_anal = scatter!(
        1 ./ betas,
        (1 .- plaqs_sq./2).*betas./4,
        yerror = plaqs_sq_err.*betas./8, 
        label = "Square", 
        color = cb_purple, 
        markerstrokecolor = cb_purple, # :auto, 
        markershape = :square,
        markersize = 7,
        legend = :topright,
        # foreground_color_legend = nothing,
        background_color_legend = nothing,
        xlabel = latexstring("\$1/\\beta = (ag)^2/4\$"),
        ylabel = latexstring("\$ s_\\mathrm{wil}/g^2 \$"),
        tickfontsize = 10,
        labelfontsize = 14,
        legendfontsize = 10,
        # msw = 10,
    )
    image_plaq_anal = scatter!(
        1 ./ betas,
        (1 .- plaqs_hex./2).*betas./4,
        yerror = plaqs_hex_err.*betas./8,            
        label = "Hexagonal", 
        color = cb_pink, 
        markerstrokecolor = cb_pink, # :auto, 
        markershape = :hexagon,
        markersize = 7,
    )
    image_plaq_anal = plot!(
        [0.5],
        [0.3],
        label = "Analytic",
        color = :black
    )
    image_plaq_anal = scatter!(
        [0.0],
        [3/8],
        markershape = :star5,
        markersize = 7,
        color = :white,
        label = latexstring("Cont. value \$3/8\$")
    )
    image_plaq_anal = plot!(
        [0.0],
        [3/8],
        label = "Slope at cont. value",
        color = cb_red,
    )
        
    display(image_plaq_anal)
    fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\plots\\plaq_square_hex_anal\\swil_sha_L_$L.pdf"
    # savefig(fig_path)

    tab_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\tabellen\\plaq_square_hex_anal\\swil_sha_L_$L.txt"
    bla = open(tab_path, "w")
    write(bla, "\\begin{table}[H]\n\t\\centering\n\t\\hline\n\t\\begin{tabular}{|c|c|c|c|}\n")
    write(bla, "\t\t \$\\beta\$ &\t square &\t hex &\t analyt. \t \\\\\\hline\n")
    for beta_ind = 1:4
        write(bla, "\t\t $(2.0*beta_ind) &\t $(swil_sq_str[beta_ind]) &\t $(swil_hex_str[beta_ind]) &\t $(swil_anal_str[beta_ind])\t \\\\\\hline\n")
    end
    write(bla,"\t\\end{tabular}\n\t\\caption{Caption}\n\t\\label{tab:my_label}\n\\end{table}")
    close(bla)
end

# \begin{table}[H]
#     \centering
#     \begin{tabular}{|c|c|c|}
#          &  \\
#          & 
#     \end{tabular}
#     \caption{Caption}
#     \label{tab:my_label}
# \end{table}




hex_best = []
hex_best_err = []
square_best = []
square_best_err = []

L = 32
# for L in [32] # 32:32:128
for β in [2.0, 4.0, 6.0, 8.0]
    # L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    n_stout = 0
    ρ   = 0.12
    sim_count = 4
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    means = readdlm(mean_vals_mike_path)
    num_means = size(means,1)
    num_loops = size(means,2)

    ratio_means = []
    ratio_mean_errs  = []
    for j = 1:3:4
        # @show j
        b_size = maximum([Int(round(2*auto_corr_time(means[:,k]) + 1, RoundUp)) for k in [j,j+1,j+2,j+3]])
        bla = jack_string_creutz(means,  j,j+1,j+2,j+3, b_size)
        push!(ratio_means, bla[1])
        push!(ratio_mean_errs, bla[2])
    end

    int_start = 1
    int_end = 2
    x_lab = string.([[1,1], [2,2], [3,3], [4,4], [5,5], [6,6]])

    image = scatter(
        x_lab[int_start:int_end], 
        ratio_means[int_start:int_end], 
        yerror = ratio_mean_errs[int_start:int_end], 
        label = "Square", 
        markerstrokecolor = :auto,
        markershape = :diamond
    )




    N_t = 2*L #+ i*16
    N_x = L #+ i*16
    hot = true
    # ϵ   = 0.2 
    n_stout = 0
    ρ   = 0.12
    sim_count = 3
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4]]
    num_loops = length(loops)

    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2_data\\hex_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_hex_path = string(base_path,"\\mean_vals.txt")

    means_hex = readdlm(mean_vals_hex_path)
    num_means = size(means_hex,1)
    num_loops = size(means_hex,2)

    ratio_means_hex = []
    ratio_mean_errs_hex  = []
    for j = 1:3:4
        b_size = maximum([Int(round(2*auto_corr_time(means_hex[:,k]) + 1, RoundUp)) for k in [j,j+1,j+2,j+3]])
        bla = jack_string_creutz(means_hex,  j,j+1,j+2,j+3, b_size)
        push!(ratio_means_hex, bla[1])
        push!(ratio_mean_errs_hex, bla[2])
    end

    int_end = 2

    image = scatter!(
        x_lab[int_start:int_end], 
        ratio_means_hex[int_start:int_end], 
        yerror = ratio_mean_errs_hex[int_start:int_end], 
        label = "Hexagonal", 
        markerstrokecolor = :auto,
        markershape = :hexagon)

    image = hline!([-log(analytic_plaq(β))], label = "Analytic", color = :red)
    
    image = plot!(title = "Creutz Ratios with β = $β,
    Sq.: N_t = N_x = $L, Hex.: N_t = 2⋅N_x = 2⋅$L", 
    xlabel = "R and T in {⟨W(R,T)⟩ ⟨W(R+1,T+1)⟩} / {⟨W(R+1,T)⟩ ⟨W(R,T+1)⟩}")
    display(image)

    # mkpath("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_creutz\\beta_$β")
    # savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_creutz\\beta_$β\\ratios_beta_$β._N_x_$L.pdf")

    push!(hex_best, ratio_means_hex[1])
    push!(hex_best_err, ratio_mean_errs_hex[1])
    push!(square_best, ratio_means[1])
    push!(square_best_err, ratio_mean_errs[1])
end
# end

# L = 32
let
    betas = [2.0, 4.0, 6.0, 8.0]
    betas_plot = Vector(1.0:0.1:700) 
    slop = - (log(analytic_plaq(700))*700/4 - log(analytic_plaq(699))*699/4 ) / (1/700-1/699)
    image_creutz_compare = plot(
        1 ./betas_plot,
        -log.(analytic_plaq.(betas_plot)) .* betas_plot ./4,
        label = :false,
        color = :black
    )
    image_creutz_compare = plot!(
        0:0.05:0.25,
        [3/8 + slop*x for x in 0:0.05:0.25],
        label = :false, # "Slope at cont. value",
        color = cb_red
    )
    image_creutz_compare = scatter!(
        1 ./betas,
        square_best.*betas./4,
        yerror = square_best_err.*betas./4,
        label = "Square",
        color = cb_purple,
        markerstrokecolor = cb_purple,
        markershape = :square,
        markersize = 7,
    )
    image_creutz_compare = scatter!(
        1 ./betas,
        hex_best.*betas./4,
        yerror = hex_best_err.*betas./4,
        label = "Hexagonal",
        color = cb_pink,
        markerstrokecolor = cb_pink,
        markershape = :hexagon,
        markersize = 7,
    )
    image_creutz_compare = plot!(
        [0.0],
        [3/8],
        label = "Analytic",
        color = :black
    )
    image_creutz_compare = plot!(
        xlabel = latexstring("\$1/\\beta = (ag)^2/4\$"),
        ylabel = latexstring("\$\\sigma/g^2\$"),
        tickfontsize = 9,
        labelfontsize = 14,
        legendfontsize = 10,
    )
    image_creutz_compare = scatter!(
        [0.0],
        [3/8],
        markershape = :star5,
        markersize = 7,
        color = :white,
        label = latexstring("Cont. value \$3/8\$")
    )
    image_creutz_compare = plot!(
        [0.0],
        [3/8],
        color = cb_red,
        label = "Slope at cont. value"
    )
    image_creutz_compare = plot!(
        legend = :bottom,
        background_color_legend = :false
    )
    display(image_creutz_compare)
    fig_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\plots\\string_square_hex_anal\\string_sha_L_$L.pdf"
    # savefig(fig_path)

    tab_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\Master_Thesis\\tabellen\\string_square_hex_anal\\string_sha_L_$L.txt"
    bla = open(tab_path, "w")
    write(bla, "\\begin{table}[H]\n\t\\centering\n\t\\hline\n\t\\begin{tabular}{|c|c|c|c|}\n")
    write(bla, "\t\t \$\\beta\$ &\t square &\t hex &\t analyt. \t \\\\\\hline\n")
    error_digs = 2
    anal_digs_pre  = 1 
    anal_digs_post = 7
    string_sq_str   = format_x_err.(square_best.*betas./4, square_best_err.*betas./4, error_digs)
    string_hex_str  = format_x_err.(hex_best.*betas./4, hex_best_err.*betas./4, error_digs)
    for beta_ind = 1:4
        string_anal_str = @sprintf("%*.*f", anal_digs_pre, anal_digs_post, -log(analytic_plaq(betas[beta_ind]))*betas[beta_ind]/4)
        write(bla, "\t\t $(2.0*beta_ind) &\t $(string_sq_str[beta_ind]) &\t $(string_hex_str[beta_ind]) &\t $string_anal_str\t \\\\\\hline\n")
    end
    write(bla,"\t\\end{tabular}\n\t\\caption{Caption}\n\t\\label{tab:my_label}\n\\end{table}")
    close(bla)
end
# end
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_creutz\\compare_fav_ratios_N_x_$L.pdf")











L = 128
betas = [2.0, 4.0, 6.0, 8.0]
strings = []
strings_hex = []
string_errs = []
string_hex_errs = []
for β in betas
    # L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    n_stout = 0
    ρ   = 0.12
    sim_count = 4
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    means = readdlm(mean_vals_mike_path)
    num_means = size(means,1)
    num_loops = size(means,2)

    # string_vals = [-real.(log(creutz(means[i,1:4], 1, 2, 3, 4) + 0im)) for i = 1:num_means]
    # string_vals = [-(log(creutz(means[i,1:4], 1, 2, 3, 4))) for i = 1:num_means]
    # b_size = Int(round(2*auto_corr_time(string_vals) + 1, RoundUp))
    # string_mean, string_err = jackknife(string_vals, b_size)
    b_size = 1 + 2 * maximum([Int(round(auto_corr_time(means[i,:]), RoundUp)) for i = 1:num_means ])
    string_mean, string_err = string_jack(means, b_size, 1, 2, 3, 4)



    N_t = 2*L #+ i*16
    N_x = L #+ i*16
    hot = true
    n_stout = 0
    ρ   = 0.12
    sim_count = 3
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4]]
    num_loops = length(loops)

    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_hex_path = string(base_path,"\\mean_vals.txt")

    means_hex = readdlm(mean_vals_hex_path)
    num_means = size(means_hex,1)
    num_loops = size(means_hex,2)

    # string_vals_hex = [-real.(log(creutz(means_hex[i,1:4], 1, 2, 3, 4) + 0im)) for i = 1:num_means]
    # string_vals_hex = [-(log(creutz(means_hex[i,1:4], 1, 2, 3, 4))) for i = 1:num_means]
    # b_size = Int(round(2*auto_corr_time(string_vals_hex) + 1, RoundUp))
    # string_hex_mean, string_hex_err = jackknife(string_vals_hex, b_size)
    b_size = 1 + 2 * maximum([Int(round(auto_corr_time(means_hex[i,:]), RoundUp)) for i = 1:num_means ])
    string_hex_mean, string_hex_err = string_jack(means_hex, b_size, 1, 2, 3, 4)



    push!(strings, string_mean)
    push!(strings_hex, string_hex_mean)
    push!(string_errs, string_err)
    push!(string_hex_errs, string_hex_err)    
end



# betas = [2.0, 4.0, 6.0, 8.0]
image_strings = scatter(
    betas,
    strings,
    yerror = string_errs,
    label = "Square",
    markerstrokecolor = :auto,
    markershape = :diamond
)
image_strings = scatter!(
    betas,
    strings_hex,
    yerror = string_hex_errs,
    label = "Hexagonal",
    markerstrokecolor = :auto,
    markershape = :hexagon
)
betas_plot = Array(first(betas)-0.2:0.01:last(betas)+0.2) 
image_strings = plot!(
    betas_plot,
    -log.(analytic_plaq.(betas_plot)),
    label = "Analytic",
    # markerstrokecolor = :auto,
    # markershape = :hline,
    markersize = 10,
    color = :red
)
image_strings = plot!(
    title = "String Tensions using ⟨W(1,1)⟩, ..., ⟨W(2,2)⟩
    Square: N_t = N_x = $L,  Hex.: N_t = 2⋅N_x = 2⋅$L",
    xlabel = "β"
)
display(image_strings)


anals = [-log(analytic_plaq(betas[i])) for i = 1:4]
# abs_distances_sq = [minimum([abs(strings[i]+string_errs[i]-anals[i]), abs(strings[i]-string_errs[i]-anals[i])]) for i = 1:4]
# abs_distances_hex = [minimum([abs(strings_hex[i]+string_hex_errs[i]-anals[i]), abs(strings_hex[i]-string_hex_errs[i]-anals[i])]) for i = 1:4]
sigmas_sq = [round(abs((anals[i]-strings[i])/string_errs[i]), digits = 3) for i = 1:4 ]
sigmas_hex = [round(abs((anals[i]-strings_hex[i])/string_hex_errs[i]), digits = 3) for i = 1:4 ]


for i = 1:4
    println(i, enu_endings[i], " Square String: ", strings[i], " ± ", string_errs[i])
    println(i, enu_endings[i], " Hexag. String: ", strings_hex[i], " ± ", string_hex_errs[i])
    println(i, enu_endings[i], " Anal. String:  ", anals[i])
    println(i, enu_endings[i], " Square distance:  ", sigmas_sq[i], "σ")
    println(i, enu_endings[i], " Hexag. distance:  ", sigmas_hex[i], "σ")
    println(" ")
    if i == 4
        println(" ")
    end
end

info = []
for i = 1:4
    push!(info,"β = $(betas[i])")
    push!(info,"Square String: $(strings[i]) ± $(string_errs[i])" )
    push!(info,"Hexag. String: $(strings_hex[i]) ± $(string_hex_errs[i])" )
    push!(info,"Anal. String:  $(anals[i])" )
    push!(info,"Square distance: $(sigmas_sq[i])σ" )
    push!(info,"Hexag. distance: $(sigmas_hex[i])σ" )
    push!(info," ")
    if i == 4
        push!(info, " ")
    end
end

# bla = open("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_creutz\\compare_strings_N_x_$L.txt", "w")
# writedlm(bla, info)
# close(bla)
# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_creutz\\compare_strings_N_x_$L.pdf")










L = 32
betas = [2.0, 4.0, 6.0, 8.0]
spaces = []
spaces_hex = []
space_errs = []
space_hex_errs = []
for β in betas
    # L = 32
    N_t = L #+ i*16
    N_x = L #+ i*16
    n_stout = 0
    ρ   = 0.12
    sim_count = 4
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4], [4,5], [5,4], [5,5], [5,6], [6,5], [6,6]]
    num_loops = length(loops)

    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_mike_path = string(base_path,"\\mean_vals_mike.txt")

    means = readdlm(mean_vals_mike_path)
    num_means = size(means,1)
    num_loops = size(means,2)

    # string_vals = [-real.(log(creutz(means[i,1:4], 1, 2, 3, 4) + 0im)) for i = 1:num_means]
    # string_vals = [-(log(creutz(means[i,1:4], 1, 2, 3, 4))) for i = 1:num_means]
    # b_size = Int(round(2*auto_corr_time(string_vals) + 1, RoundUp))
    # string_mean, string_err = jackknife(string_vals, b_size)
    b_size = 1 + 2 * maximum([Int(round(auto_corr_time(means[i,:]), RoundUp)) for i = 1:num_means ])
    space_mean, space_err = space_jack(means, b_size, 1, 2, 3, 4)



    N_t = 2*L #+ i*16
    N_x = L #+ i*16
    hot = true
    n_stout = 0
    ρ   = 0.12
    sim_count = 3
    loops   = [[1,1], [1,2], [2,1], [2,2], [2,3], [3,2], [3,3], [3,4], [4,3], [4,4]]
    num_loops = length(loops)

    base_path = "C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\beta_$β\\N_t_$N_t.N_x_$N_x\\n_stout_$n_stout._rho_$ρ\\sim_count_$sim_count"
    mean_vals_hex_path = string(base_path,"\\mean_vals.txt")

    means_hex = readdlm(mean_vals_hex_path)
    num_means = size(means_hex,1)
    num_loops = size(means_hex,2)

    # string_vals_hex = [-real.(log(creutz(means_hex[i,1:4], 1, 2, 3, 4) + 0im)) for i = 1:num_means]
    # string_vals_hex = [-(log(creutz(means_hex[i,1:4], 1, 2, 3, 4))) for i = 1:num_means]
    # b_size = Int(round(2*auto_corr_time(string_vals_hex) + 1, RoundUp))
    # string_hex_mean, string_hex_err = jackknife(string_vals_hex, b_size)
    b_size = 1 + 2 * maximum([Int(round(auto_corr_time(means_hex[i,:]), RoundUp)) for i = 1:num_means ])
    space_hex_mean, space_hex_err = space_jack(means_hex, b_size, 1, 2, 3, 4)



    push!(spaces, space_mean)
    push!(spaces_hex, space_hex_mean)
    push!(space_errs, space_err)
    push!(space_hex_errs, space_hex_err)    
end

# betas = [2.0, 4.0, 6.0, 8.0]
image_spaces = scatter(
    betas,
    spaces,
    yerror = space_errs,
    label = "Square",
    markerstrokecolor = :auto,
    markershape = :diamond
)
image_spaces = scatter!(
    betas,
    spaces_hex,
    yerror = space_hex_errs,
    label = "Hexagonal",
    markerstrokecolor = :auto,
    markershape = :hexagon
)
betas_plot = Array(first(betas)-0.2:0.01:last(betas)+0.2) 
image_spaces = plot!(
    betas_plot,
    sp_fac .* sqrt.( -log.(analytic_plaq.(betas_plot))),
    label = "Analytic",
    # markerstrokecolor = :auto,
    # markershape = :hline,
    markersize = 10,
    color = :red
)
image_spaces = plot!(
    title = "Spacings using ⟨W(1,1)⟩, ..., ⟨W(2,2)⟩
    Square: N_t = N_x = $L,  Hex.: N_t = 2⋅N_x = 2⋅$L",
    xlabel = "β"
)
display(image_spaces)

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_creutz\\compare_spacings_N_x_$L.pdf")



















# spacing = sp_fac .* sqrt.(strings)
# spacing_errs_p = sp_fac .* abs.(sqrt.(strings .+ string_errs) .- sqrt.(strings))
# spacing_errs_m = sp_fac .* abs.(sqrt.(strings .- string_errs) .- sqrt.(strings))

# hex_spacing = sp_fac .* sqrt.(strings_hex)
# hex_spacing_errs_p = sp_fac .* abs.(sqrt.(strings_hex .+ string_hex_errs) .- sqrt.(strings_hex))
# hex_spacing_errs_m = sp_fac .* abs.(sqrt.(strings_hex .- string_hex_errs) .- sqrt.(strings_hex))

# image_spacing = scatter(
#     betas,
#     spacing,
#     yerror = (spacing_errs_m, spacing_errs_p),
#     label = "Square",
#     markerstrokecolor = :auto,
#     markershape = :diamond
# )
# image_spacing = scatter!(
#     betas,
#     hex_spacing,
#     yerror = (hex_spacing_errs_m, hex_spacing_errs_p),
#     label = "Hexagonal",
#     markerstrokecolor = :auto,
#     markershape = :hexagon
# )
# betas_plot = Array(first(betas)-0.2:0.01:last(betas)+0.2) 
# image_spacing = plot!(
#     betas_plot,
#     sp_fac .* sqrt.(-log.(analytic_plaq.(betas_plot))),
#     label = "Analytic",
#     # markerstrokecolor = :auto,
#     # markershape = :hline,
#     markersize = 10,
#     color = :red
# )
# image_spacing = plot!(
#     title = "Sommer Spacings using ⟨W(1,1)⟩, ..., ⟨W(2,2)⟩
#     Sq: N_t = N_x = $L,  Hex: N_t = 2⋅N_x = 2⋅$L",
#     xlabel = "β"
# )
# display(image_spacing)

#=
println("1st Square String: ", strings[1], " ± ", string_errs[1])
println("1st Hexag. String: ", strings_hex[1], " ± ", string_hex_errs[1])
println("1st Anal. String:  ", -log(analytic_plaq(betas[1])))
println(" ")
println("2nd Square string: ", strings[2], " ± ", string_errs[2])
println("2nd Hexag. string: ", strings_hex[2], " ± ", string_hex_errs[2])
println("2nd Anal. string:  ", -log(analytic_plaq(betas[2])))
println(" ")
println(" ")
=#

# savefig("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\compare_creutz\\compare_spacings_N_x_$L.pdf")




function cosnx(n,x)
    bla = 0.0
    for i = 0:Int(round(Int,n/2,RoundDown))+1
        for j = 0:i
            bla += (-1)^(i-j) * binomial(n,2*i) * binomial(i,j) * cos(x)^(n-2*(i-j)) 
        end
    end
    return bla
end

# nrand = rand(1:100)
xrand = 2*pi*rand()
cos(50*xrand) - cosnx(50, xrand)