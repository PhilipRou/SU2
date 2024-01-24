include("test_head.jl")


function staple_test()
    U = gaugefield_SU2(N_t,N_x,true)
    μ = rand(1:2)
    t = rand(1:N_t)
    x = rand(1:N_x)
    @assert isapprox(staple(U,μ,t,x), adjoint(staple_dag(U,μ,t,x)))
    return true
end

staple_test()

function traceless_test()
    coeffs = ran_SU2(rand())
    @assert tr(coeffs - adjoint(coeffs)) == 0.0
    return true
end

traceless_test()

function exp_traceless_test()
    coeffs = ran_SU2(rand())
    coeffs_tr_less = coeffs - adjoint(coeffs)
    @assert isapprox(exp_traceless(coeffs), grp2coeffs(exp(coeffs2grp(coeffs_tr_less))))
    return true
end

exp_traceless_test()

#=
function old_VS_new_stout_test()
    function old_stout(U, n_stout, α)
        V = deepcopy(U)
        for i = 1:n_stout
            for μ = 1:2
                for trip = 1:2  # For the chequer board pattern
                    for t = 1:N_t
                        for x = (1+mod(t+trip,2)):2:N_x
                            # stap_link = staple(V,μ,t,x) * adjoint(V[μ,t,x])
                            # stap_link *= α*0.5
                            # V[μ,t,x] = exp_traceless(stap_link) * V[μ,t,x]
                            stap_link = coeffs2grp(staple(V,μ,t,x) * adjoint(V[μ,t,x]))
                            V[μ,t,x] = grp2coeffs(exp(α*0.5*(stap_link - adjoint(stap_link)))) * V[μ,t,x]
                        end
                    end
                end
            end
        end
        return V
    end
    function new_stout(U, n_stout, α)
        V = deepcopy(U)
        for i = 1:n_stout
            for μ = 1:2
                for trip = 1:2  # For the chequer board pattern
                    for t = 1:N_t
                        for x = (1+mod(t+trip,2)):2:N_x
                            stap_link = staple(V,μ,t,x) * adjoint(V[μ,t,x])
                            stap_link = α*0.5*stap_link
                            V[μ,t,x] = exp_traceless(stap_link) * V[μ,t,x]
                            # stap_link = coeffs2grp(staple(V,μ,t,x) * adjoint(V[μ,t,x]))
                            # V[μ,t,x] = grp2coeffs(exp(α*0.5*(stap_link - adjoint(stap_link)))) * V[μ,t,x]
                        end
                    end
                end
            end
        end
        return V
    end
    test_field = gaugefield_SU2(N_t,N_x,true)
    old = old_stout(test_field, 1, 0.1)
    new = new_stout(test_field, 1, 0.1)
    for x = 1:N_x
        for t = 1:N_t
            for μ = 1:2
                @assert isapprox(old[μ,t,x], new[μ,t,x])
            end
        end
    end
    return true
end

old_VS_new_stout_test()
=#


# ⭕⭕⭕ Not really a test function
function stout_test()
    U = gaugefield_SU2(N_t,N_x,true)
    V = stout(U,1,0.12)
    return true 
end

stout_test()