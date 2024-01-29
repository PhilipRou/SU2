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

# test_field = gaugefield(32, 32, true, "U2", "hexagon")