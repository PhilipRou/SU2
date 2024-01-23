################################################################################
# using Plots
using Statistics
using BenchmarkTools
using LinearAlgebra
using DelimitedFiles

include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\gaugefields.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_hex.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\observables_cube.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_square.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_hex.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\updates\\updates_cube.jl")
include("C:\\Users\\proue\\OneDrive\\Desktop\\Physik Uni\\julia_projects\\SU2\\observables\\smearing.jl")
# include("SU2_d_dim.jl")





### Some basic params needed for some functions
N_t = 64
N_x = 64
β   = 6.0
acc = [0]

### Basic Boolean functions for our way of representing SU(2) elements
# An SU(2) matrix has to be unitary and have determinant = 1.0
function is_SU2(mat::Matrix)
    return isapprox(det(mat), 1.0) * isapprox(mat*adjoint(mat), [1 0; 0 1])
end

function is_SU2(X::coeffs_SU2)
    return isapprox(X.a^2 + X.b^2 + X.c^2 + X.d^2, 1.0)
end
