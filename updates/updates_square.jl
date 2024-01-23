include("updates_head.jl")



function staple_dag(U, μ, x, t)
    NX = size(U,3)
    NT = size(U,2)
    a = coeffs_SU2(0.0,0.0,0.0,0.0)
    b = coeffs_SU2(0.0,0.0,0.0,0.0)
    x_p = mod1(x+1, NX) # x%NX +1                 
    t_p = mod1(t+1, NT) # t%NT +1                 
    x_m = mod1(x-1, NX) # (x + NX -2)%NX +1   
    t_m = mod1(t-1, NT) # (t + NT -2)%NT +1   

    # 🐌 More efficient: only use adj_SU2 once 🐌 (but less human-readable, no?)
    if μ == 1
        a = U[2,x_p,t] * adj_SU2(U[1,x,t_p]) * adj_SU2(U[2,x,t])
        b = adj_SU2(U[2,x_p,t_m]) * adj_SU2(U[1,x,t_m]) * U[2,x,t_m]
    else #if μ == 2
        a = U[1,x,t_p] * adj_SU2(U[2,x_p,t]) * adj_SU2(U[1,x,t])
        b = adj_SU2(U[1,x_m,t_p]) * adj_SU2(U[2,x_m,t]) * U[1,x_m,t]
    end
    return a + b 
end

# ❌ Not in use anymore! Only for testing purposes ❌
function delta_S_gauge(U, μ, x, t, old_coeffs::coeffs_SU2, new_coeffs::coeffs_SU2, β)
    return β*0.5*(tr((old_coeffs - new_coeffs) * staple_dag(U,μ,x,t)))
end

function metro!(U, μ, x, t, step, β, acc)
    new_coeffs = ran_SU2(step) * U[μ,x,t]
    staple_d = staple_dag(U,μ,x,t)
    S_old = β*0.5*tr(U[μ,x,t] * staple_d)
    S_new = β*0.5*tr(new_coeffs * staple_d)
    if rand() < exp(S_new-S_old)
        U[μ,x,t] = new_coeffs
        acc[1] += 1
    end
    return nothing
end

function lexico_metro!(U, step, β, acc)
    NX = size(U,3)
    NT = size(U,2)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                metro!(U,μ,x,t,step,β,acc)
            end
        end
    end
    return nothing
end

function chess_metro!(U, step, β, acc)
    NX = size(U,3)
    NT = size(U,2)
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    metro!(U,μ,x,t,step, β, acc)
                end
            end
        end
    end
    return nothing
end

#
function ran_metro!(U, step, β, acc)
    NX = size(U,2)
    NT = size(U,3)
    coords = [[rand(1:2), rand(1:NX), rand(1:NT)] for i = 1:2*NX*NT]
    for i = 1:2*NX*NT
        μ, x, t = coords[i]
        metro!(U,μ,x,t,step,β,acc)
    end    
    return nothing
end

#
function overrelax!(U, μ, x, t)
    v = proj_SU2(staple_dag(U,μ,x,t))
    U[μ,x,t] = adj_SU2(v *  U[μ,x,t] * v)
    return nothing
end

#
function lexico_overrelax!(U, acc)
    NT = size(U,2)
    NX = size(U,3)
    for t = 1:NT
        for x = 1:NX
            for μ = 1:2
                overrelax!(U, μ, x, t)
            end
        end
    end
    acc[1] += 2*NX*NT
    return nothing
end

#
function chess_overrelax!(U, acc)
    NT = size(U,2)
    NX = size(U,3)
    for μ = 1:2
        for trip = 1:2
            for t = 1:NT
                for x = (1+mod(t+trip,2)):2:NX
                    overrelax!(U,μ,x,t)
                end
            end
        end
    end
    acc[1] += 2*NX*NT
    return nothing
end