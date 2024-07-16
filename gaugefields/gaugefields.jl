include("SU2_gaugefields.jl")
include("U2_gaugefields.jl")



function gaugefield(N_x::Int64, N_t::Int64, hot::Bool, group::String, lattice::String)
    if group == "SU2"
        if lattice == "square"
            return gaugefield_SU2(N_x, N_t, hot)
        elseif lattice == "hexagon"
            return hexfield_SU2(N_x, N_t, hot)
        else
            error(
                "For the gauge group SU(2) the supported keywords for \"lattice\" are:
                \"square\"
                \"hexagon\" "
            )
        end
    elseif group == "U2"
        if lattice == "square"
            return gaugefield_U2(N_x, N_t, hot)
        elseif lattice == "hexagon"
            return hexfield_U2(N_x, N_t, hot)
        else
            error(
                "For the gauge group U(2) the supported keywords for \"lattice\" are:
                \"square\"
                \"hexagon\" "
            )
        end
    else
        error(
            "The supported keywords for \"group\" are:
            \"SU2\" for the gauge group SU(2)
            \"U2\" for the gauge group U(2)"
        )
    end
end

# test_field = gaugefield(32, 32, true, "U2", "pentagon")

# Apply temporal gauge onto a square config
function temp_gauge(U, group)
    NX = size(U,2)
    NT = size(U,3)
    if group == "SU2"
        V = gaugefield_SU2(NX, NT, false)
        Ω_slice = [coeffs_Id_SU2() for x = 1:NX] 
        for t = 1:NT
            V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice,-1))
            Ω_slice = Ω_slice .* U[2,:,t]
        end
        V[2,:,NT] = Ω_slice
        return V
    elseif group == "U2"
        V = gaugefield_U2(NX, NT, false)
        Ω_slice = [coeffs_Id_U2() for x = 1:NX] 
        for t = 1:NT
            V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice,-1))
            Ω_slice = Ω_slice .* U[2,:,t]
        end
        V[2,:,NT] = Ω_slice
        return V
    end
end

function max_gauge(U, group)
    NX = size(U,2)
    NT = size(U,3)
    if group == "U2"
        V  = gaugefield_U2(NX, NT, false)
        Ω_slice = [coeffs_Id_U2()]
        for x = 1:NX-1
            next_el = last(Ω_slice) * U[1,x,1]
            push!(Ω_slice, next_el)
        end
        Ω_slice_copy = Ω_slice
        V[1,NX,1] = last(Ω_slice) * U[1,NX,1]
        for t = 2:NT
            Ω_slice = Ω_slice .* U[2,:,t-1]
            V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice, -1))
        end
        V[2,:,NT] = Ω_slice .* U[2,:,NT] .* adjoint.(Ω_slice_copy)
        return V
    elseif group == "SU2"
        V  = gaugefield_SU2(NX, NT, false)
        Ω_slice = [coeffs_Id_SU2()]
        for x = 1:NX-1
            next_el = last(Ω_slice) * U[1,x,1]
            push!(Ω_slice, next_el)
        end
        Ω_slice_copy = Ω_slice
        V[1,NX,1] = last(Ω_slice) * U[1,NX,1]
        for t = 2:NT
            Ω_slice = Ω_slice .* U[2,:,t-1]
            V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice, -1))
        end
        V[2,:,NT] = Ω_slice .* U[2,:,NT] .* adjoint.(Ω_slice_copy)
        return V
    end
end

function max_gauge(U, group, g)
    NX = size(U,2)
    NT = size(U,3)
    if group == "U2"
        V  = gaugefield_U2(NX, NT, false)
        # Ω_slice = [coeffs_Id_U2()]
        Ω_slice = [g]
        for x = 1:NX-1
            next_el = last(Ω_slice) * U[1,x,1]
            push!(Ω_slice, next_el)
        end
        Ω_slice_copy = Ω_slice
        V[1,NX,1] = last(Ω_slice) * U[1,NX,1] * adjoint(g)
        for t = 2:NT
            Ω_slice = Ω_slice .* U[2,:,t-1]
            V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice, -1))
        end
        V[2,:,NT] = Ω_slice .* U[2,:,NT] .* adjoint.(Ω_slice_copy)
        return V
    elseif group == "SU2"
        V  = gaugefield_SU2(NX, NT, false)
        # Ω_slice = [coeffs_Id_SU2()]
        Ω_slice = [g]
        for x = 1:NX-1
            next_el = last(Ω_slice) * U[1,x,1]
            push!(Ω_slice, next_el)
        end
        Ω_slice_copy = Ω_slice
        V[1,NX,1] = last(Ω_slice) * U[1,NX,1] * adjoint(g)
        for t = 2:NT
            Ω_slice = Ω_slice .* U[2,:,t-1]
            V[1,:,t] = Ω_slice .* U[1,:,t] .* adjoint.(circshift(Ω_slice, -1))
        end
        V[2,:,NT] = Ω_slice .* U[2,:,NT] .* adjoint.(Ω_slice_copy)
        return V
    end
end