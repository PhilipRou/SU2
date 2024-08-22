include("SU2_analyze_head.jl")





for L in [16,32,64,96]
    # L = 16
        β = 6.0
        N_x = L
        N_t = L
    
        acceptances_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_chess_N_t_$N_t._N_x_$N_x.txt")
        acceptances_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_lex_N_t_$N_t._N_x_$N_x.txt")
        acceptances_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_ran_N_t_$N_t._N_x_$N_x.txt")
    
        actions_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\actions_chess_N_t_$N_t._N_x_$N_x.txt")
        actions_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\actions_lex_N_t_$N_t._N_x_$N_x.txt")
        actions_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\actions_ran_N_t_$N_t._N_x_$N_x.txt")
        
        n_meas = length(actions_chess)
    
        acc_chess = round(last(acceptances_chess) /2/N_t/N_x/n_meas, digits = 3)
        acc_lex = round(last(acceptances_lex) /2/N_t/N_x/n_meas, digits = 3)
        acc_ran = round(last(acceptances_ran) /2/N_t/N_x/n_meas, digits = 3)
    
        auto_time_chess = []
        auto_time_lex = []
        auto_time_ran = []
        for start = 50:50:n_meas-50
            push!(auto_time_chess, auto_corr_time(actions_chess[start:n_meas]))
            push!(auto_time_lex, auto_corr_time(actions_lex[start:n_meas]))
            push!(auto_time_ran, auto_corr_time(actions_ran[start:n_meas]))
        end
    
        image_auto_time = plot(50:50:n_meas-50, auto_time_chess, label = "chess, acc.rate = $acc_chess")
        image_auto_time = plot!(50:50:n_meas-50, auto_time_lex, label = "lexicog., acc.rate = $acc_lex")
        image_auto_time = plot!(50:50:n_meas-50, auto_time_ran, label = "ran., acc.rate = $acc_ran")
        image_auto_time = plot!(title = "Int. Auto. Time, Different Metr. Updates, Square
        N_t = N_x = $L,  β = $β,  ϵ = 0.2",
        xlabel = "MC time t at which integration started")
    
        display(image_auto_time)
        # savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\auto_time_N_t_$N_t._N_x_$N_x.pdf")
    end
    
    
    
    for L in [16,32,64,96]
    # L = 16
        β = 6.0
        N_x = L
        N_t = L
        
        acceptances_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_chess_N_t_$N_t._N_x_$N_x.txt")
        acceptances_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_lex_N_t_$N_t._N_x_$N_x.txt")
        acceptances_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_ran_N_t_$N_t._N_x_$N_x.txt")
    
        actions_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\actions_chess_N_t_$N_t._N_x_$N_x.txt")
        actions_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\actions_lex_N_t_$N_t._N_x_$N_x.txt")
        actions_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\actions_ran_N_t_$N_t._N_x_$N_x.txt")
        
        n_meas = length(actions_chess)
    
        acc_chess = round(last(acceptances_chess) /2/N_t/N_x/n_meas, digits = 3)
        acc_lex = round(last(acceptances_lex) /2/N_t/N_x/n_meas, digits = 3)
        acc_ran = round(last(acceptances_ran) /2/N_t/N_x/n_meas, digits = 3)
    
        auto_chess = []
        auto_lex = []
        auto_ran = []
        # start = 200
        for τ = 0:30
            push!(auto_chess, auto_corr_norm(actions_chess,τ))
            push!(auto_lex, auto_corr_norm(actions_lex,τ))
            push!(auto_ran, auto_corr_norm(actions_ran,τ))
        end
    
        image_time = plot(0:30, auto_chess, label = "chess, acc.rate = $acc_chess")
        image_time = plot!(0:30, auto_lex, label = "lexicog., acc.rate = $acc_lex")
        image_time = plot!(0:30, auto_ran, label = "ran, acc.rate = $acc_ran")
        image_auto = plot!(title = "Norm. Autocorr. Γ(t), Different Metr. Updates, Square
        N_t = 2⋅N_x = 2⋅$L,  β = $β,  ϵ = 0.2",
        xlabel = "MC time t")
    
        display(image_auto)
        # savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β._eps_0.2\\auto_corr_N_t_$N_t._N_x_$N_x.pdf")
    end



for L in [16,32,64,96]
    # L = 16
        β = 6.0
        N_x = L
        N_t = L
    
        eps_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\eps_chess_$N_t._N_x_$N_x.txt")
        eps_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\eps_lex_$N_t._N_x_$N_x.txt")
        eps_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\eps_ran_$N_t._N_x_$N_x.txt")
    
        actions_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\actions_chess_N_t_$N_t._N_x_$N_x.txt")
        actions_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\actions_lex_N_t_$N_t._N_x_$N_x.txt")
        actions_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\actions_ran_N_t_$N_t._N_x_$N_x.txt")
        
        n_meas = length(actions_chess)
    
        ϵ_chess = round(eps_chess[1], digits = 3)
        ϵ_lex = round(eps_lex[1], digits = 3)
        ϵ_ran = round(eps_ran[1], digits = 3)
    
        auto_time_chess = []
        auto_time_lex = []
        auto_time_ran = []
        for start = 50:50:n_meas-50
            push!(auto_time_chess, auto_corr_time(actions_chess[start:n_meas]))
            push!(auto_time_lex, auto_corr_time(actions_lex[start:n_meas]))
            push!(auto_time_ran, auto_corr_time(actions_ran[start:n_meas]))
        end
    
        image_auto_time = plot(50:50:n_meas-50, auto_time_chess, label = "chess, ϵ = $ϵ_chess")
        image_auto_time = plot!(50:50:n_meas-50, auto_time_lex, label = "lexicog., ϵ = $ϵ_lex")
        image_auto_time = plot!(50:50:n_meas-50, auto_time_ran, label = "ran., ϵ = $ϵ_ran")
        image_auto_time = plot!(title = "Int. Auto. Time, Different Metr. Updates, Square
        N_t = N_x = $L,  β = $β,  acc.rate ≈ 90%",
        xlabel = "MC time t at which integration started")
    
        display(image_auto_time)
        # savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\auto_time_N_t_$N_t._N_x_$N_x.pdf")
    end
    
    
    
    for L in [16,32,64,96]
    # L = 16
        β = 6.0
        N_x = L
        N_t = L
    
        eps_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\eps_chess_$N_t._N_x_$N_x.txt")
        eps_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\eps_lex_$N_t._N_x_$N_x.txt")
        eps_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\eps_ran_$N_t._N_x_$N_x.txt")
    
        actions_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\actions_chess_N_t_$N_t._N_x_$N_x.txt")
        actions_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\actions_lex_N_t_$N_t._N_x_$N_x.txt")
        actions_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\actions_ran_N_t_$N_t._N_x_$N_x.txt")
        
        n_meas = length(actions_chess)
    
        ϵ_chess = round(eps_chess[1], digits = 3)
        ϵ_lex = round(eps_lex[1], digits = 3)
        ϵ_ran = round(eps_ran[1], digits = 3)
    
        auto_chess = []
        auto_lex = []
        auto_ran = []
        # start = 200
        τ_max = 100
        for τ = 0:τ_max
            push!(auto_chess, auto_corr_norm(actions_chess,τ))
            push!(auto_lex, auto_corr_norm(actions_lex,τ))
            push!(auto_ran, auto_corr_norm(actions_ran,τ))
        end
    
        image_time = plot(0:τ_max, auto_chess, label = "chess, ϵ = $ϵ_chess")
        image_time = plot!(0:τ_max, auto_lex, label = "lexicog., ϵ = $ϵ_lex")
        image_time = plot!(0:τ_max, auto_ran, label = "ran, ϵ = $ϵ_ran")
        image_auto = plot!(title = "Norm. Autocorr. Γ(t), Different Metr. Updates, Square
        N_t = N_x = $L,  β = $β,  acc.rate ≈ 90%",
        xlabel = "MC time t")
    
        display(image_auto)
        # savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\auto_corr\\beta_$β\\auto_corr_N_t_$N_t._N_x_$N_x.pdf")
    end


# using SpecialFunctions

# function analytic_plaq(β)
#     return (besseli(0,β) + besseli(2,β)) / (2*besseli(1,β)) - 1/β
# end

# analytic_plaq(8.0)

# betas = [2.0, 4.0, 6.0, 8.0]
# plaq_means = []
# plaq_errs = []
# # plaq_exps = []
# for β in betas
#     L = 128 # [32, 64, 96, 128]
#     means = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\beta_$β\\N_t_$L.N_x_$L\\n_stout_0._rho_0.12\\sim_count_1\\mean_vals_mike.txt")
#     plaqs = 0.5 .* means[:,1]
#     b_size = Int(round(2*auto_corr_time(plaqs) + 1, RoundUp))    
#     push!(plaq_means,  mean(plaqs))
#     push!(plaq_errs, jackknife(plaqs, b_size)[2])
#     # push!(plaq_exps, 2*analytic_plaq(β))
# end

# x_vals = Array(1.0:0.01:9.0)
# plaq_exps = [analytic_plaq(β) for β in x_vals]
# image_plaqs = scatter(betas, plaq_means, yerror = plaq_errs, label = "Measurements: 32²-lattice")
# image_plaqs = plot!(x_vals, plaq_exps, label = "Analytic")
# image_plaqs = plot!(
#     title = "Comparison of Analytic Result and Measurements
#     of ⟨Pₓₜ⟩ in 2D Pure SU(2) Gauge Theory",
#     xlabel = "β",
#     xticks = [2.0, 4.0, 6.0, 8.0]
# )

# savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\square_data\\plaq_square_analytic.pdf")























for L in [16,32,64,96]
# L = 16
    β = 6.0
    N_x = L
    N_t = 2*L

    acceptances_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_chess_N_t_$N_t._N_x_$N_x.txt")
    acceptances_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_lex_N_t_$N_t._N_x_$N_x.txt")
    acceptances_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_ran_N_t_$N_t._N_x_$N_x.txt")

    actions_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\actions_chess_N_t_$N_t._N_x_$N_x.txt")
    actions_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\actions_lex_N_t_$N_t._N_x_$N_x.txt")
    actions_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\actions_ran_N_t_$N_t._N_x_$N_x.txt")
    
    n_meas = length(actions_chess)

    acc_chess = round(last(acceptances_chess) / 2 / 0.75 / N_t / N_x / n_meas, digits = 3)
    acc_lex = round(last(acceptances_lex) / 2 / 0.75 / N_t / N_x / n_meas, digits = 3)
    acc_ran = round(last(acceptances_ran) / 2 / 0.75 / N_t / N_x / n_meas, digits = 3)

    auto_time_chess = []
    auto_time_lex = []
    auto_time_ran = []
    for start = 50:50:n_meas-50
        push!(auto_time_chess, auto_corr_time(actions_chess[start:n_meas]))
        push!(auto_time_lex, auto_corr_time(actions_lex[start:n_meas]))
        push!(auto_time_ran, auto_corr_time(actions_ran[start:n_meas]))
    end

    image_auto_time = plot(50:50:n_meas-50, auto_time_chess, label = "chess, acc.rate = $acc_chess")
    image_auto_time = plot!(50:50:n_meas-50, auto_time_lex, label = "lexicog., acc.rate = $acc_lex")
    image_auto_time = plot!(50:50:n_meas-50, auto_time_ran, label = "ran., acc.rate = $acc_ran")
    image_auto_time = plot!(title = "Int. Auto. Time, Different Metrop. Updates, Hex.
    N_t = 2⋅N_x = 2⋅$L,  β = $β,  ϵ = 0.2",
    xlabel = "MC time t at which integration started")

    display(image_auto_time)
    # savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\auto_time_N_t_$N_t._N_x_$N_x.pdf")
end



for L in [16,32,64,96]
# L = 16
    β = 6.0
    N_x = L
    N_t = 2*L
    
    acceptances_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_chess_N_t_$N_t._N_x_$N_x.txt")
    acceptances_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_lex_N_t_$N_t._N_x_$N_x.txt")
    acceptances_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\acceptances_ran_N_t_$N_t._N_x_$N_x.txt")

    actions_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\actions_chess_N_t_$N_t._N_x_$N_x.txt")
    actions_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\actions_lex_N_t_$N_t._N_x_$N_x.txt")
    actions_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\actions_ran_N_t_$N_t._N_x_$N_x.txt")
    
    n_meas = length(actions_chess)

    acc_chess = round(last(acceptances_chess) / 2 / 0.75 / N_t / N_x / n_meas, digits = 3)
    acc_lex = round(last(acceptances_lex) / 2 / 0.75 / N_t / N_x / n_meas, digits = 3)
    acc_ran = round(last(acceptances_ran) / 2 / 0.75 / N_t / N_x / n_meas, digits = 3)

    auto_chess = []
    auto_lex = []
    auto_ran = []
    # start = 200
    for τ = 0:30
        push!(auto_chess, auto_corr_norm(actions_chess,τ))
        push!(auto_lex, auto_corr_norm(actions_lex,τ))
        push!(auto_ran, auto_corr_norm(actions_ran,τ))
    end

    image_time = plot(0:30, auto_chess, label = "chess, acc.rate = $acc_chess")
    image_time = plot!(0:30, auto_lex, label = "lexicog., acc.rate = $acc_lex")
    image_time = plot!(0:30, auto_ran, label = "ran, acc.rate = $acc_ran")
    image_auto = plot!(title = "Norm. Autocorr. Γ(t), Different Metrop. Updates, Hex.
    N_t = 2⋅N_x = 2⋅$L,  β = $β,  ϵ = 0.2",
    xlabel = "MC time t")

    display(image_auto)
    # savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β._eps_0.2\\auto_corr_N_t_$N_t._N_x_$N_x.pdf")
end










for L in [16,32,64,96]
# L = 16
    β = 6.0
    N_x = L
    N_t = 2*L

    eps_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\eps_chess_$N_t._N_x_$N_x.txt")
    eps_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\eps_lex_$N_t._N_x_$N_x.txt")
    eps_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\eps_ran_$N_t._N_x_$N_x.txt")

    actions_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\actions_chess_N_t_$N_t._N_x_$N_x.txt")
    actions_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\actions_lex_N_t_$N_t._N_x_$N_x.txt")
    actions_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\actions_ran_N_t_$N_t._N_x_$N_x.txt")
    
    n_meas = length(actions_chess)

    ϵ_chess = round(eps_chess[1], digits = 3)
    ϵ_lex = round(eps_lex[1], digits = 3)
    ϵ_ran = round(eps_ran[1], digits = 3)

    auto_time_chess = []
    auto_time_lex = []
    auto_time_ran = []
    for start = 50:50:n_meas-50
        push!(auto_time_chess, auto_corr_time(actions_chess[start:n_meas]))
        push!(auto_time_lex, auto_corr_time(actions_lex[start:n_meas]))
        push!(auto_time_ran, auto_corr_time(actions_ran[start:n_meas]))
    end

    image_auto_time = plot(50:50:n_meas-50, auto_time_chess, label = "chess, ϵ = $ϵ_chess")
    image_auto_time = plot!(50:50:n_meas-50, auto_time_lex, label = "lexicog., ϵ = $ϵ_lex")
    image_auto_time = plot!(50:50:n_meas-50, auto_time_ran, label = "ran., ϵ = $ϵ_ran")
    image_auto_time = plot!(title = "Int. Auto. Time, Different Metrop. Updates, Hex.
    N_t = 2⋅N_x = 2⋅$L,  β = $β,  acc.rate ≈ 90%",
    xlabel = "MC time t at which integration started")

    display(image_auto_time)
    # savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\auto_time_N_t_$N_t._N_x_$N_x.pdf")
end



for L in [16,32,64,96]
# L = 16
    β = 6.0
    N_x = L
    N_t = 2*L

    eps_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\eps_chess_$N_t._N_x_$N_x.txt")
    eps_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\eps_lex_$N_t._N_x_$N_x.txt")
    eps_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\eps_ran_$N_t._N_x_$N_x.txt")

    actions_chess = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\actions_chess_N_t_$N_t._N_x_$N_x.txt")
    actions_lex = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\actions_lex_N_t_$N_t._N_x_$N_x.txt")
    actions_ran = readdlm("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\actions_ran_N_t_$N_t._N_x_$N_x.txt")
    
    n_meas = length(actions_chess)

    ϵ_chess = round(eps_chess[1], digits = 3)
    ϵ_lex = round(eps_lex[1], digits = 3)
    ϵ_ran = round(eps_ran[1], digits = 3)

    auto_chess = []
    auto_lex = []
    auto_ran = []
    # start = 200
    τ_max = 100
    for τ = 0:τ_max
        push!(auto_chess, auto_corr_norm(actions_chess,τ))
        push!(auto_lex, auto_corr_norm(actions_lex,τ))
        push!(auto_ran, auto_corr_norm(actions_ran,τ))
    end

    image_time = plot(0:τ_max, auto_chess, label = "chess, ϵ = $ϵ_chess")
    image_time = plot!(0:τ_max, auto_lex, label = "lexicog., ϵ = $ϵ_lex")
    image_time = plot!(0:τ_max, auto_ran, label = "ran, ϵ = $ϵ_ran")
    image_auto = plot!(title = "Norm. Autocorr. Γ(t), Different Metrop. Updates, Hex.
    N_t = 2⋅N_x = 2⋅$L,  β = $β,  acc.rate ≈ 90%",
    xlabel = "MC time t")

    display(image_auto)
    # savefig("D:\\Physik Uni\\julia_projects\\SU2\\data\\hex_data\\auto_corr\\beta_$β\\auto_corr_N_t_$N_t._N_x_$N_x.pdf")
end


