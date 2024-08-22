using SpecialFunctions

include("D:\\Physik Uni\\julia_projects\\SU2\\gaugefields\\gaugefields.jl")
include("D:\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")


function mywrite(path, obs)
    bla = open(path,"a")
    writedlm(bla, obs)
    # writedlm(bla, "\n")
    close(bla)
    return nothing
end

function mywrite(path, obs::Array)
    bla = open(path,"a")
    writedlm(bla, transpose(obs))
    # writedlm(bla, "\n")
    close(bla)
    return nothing
end