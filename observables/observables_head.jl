using SpecialFunctions

include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")


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