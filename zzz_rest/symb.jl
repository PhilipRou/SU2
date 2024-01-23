################################################################################
#=
using Symbolics
# Still suffers children's diseases. Documentation bad, e.g. solve_system_eq(),
# and simply using solve_system_eq() results in some error involving occursin(),
# which cannot keep track of variables of type {Num}...

# So Symbolics.jl does not come with a scalar product. Here is my attempt of 
# one. Note that V = [x,y] is of type Vector{Num}, but [x y] without the comma
# is a (1×2)-matrix (just a row).
function Base.:*(X::Vector{Num}, Y::Vector{Num})
    @assert length(X) == length(Y) "*(::Vector{Num}, ::Vector{Num}) needs vectors of same length"
    return sum([X[i] * Y[i] for i = 1:length(X)])
end

@variables x y
z = [x, y]
z*z

# Same for Vector{Number} and Vector{Num}...
function Base.:*(X::Vector{T}, Y::Vector{Num}) where T<:Number
    @assert length(X) == length(Y) "*(::Vector{Num}, ::Vector{Num}) needs vectors of same length"
    return sum([X[i] * Y[i] for i = 1:length(X)])
end

function Base.:*(X::Vector{T}, Y::Vector{S}) where {T<:Number, S<: Number}
    @assert length(X) == length(Y) "*(::Vector{Num}, ::Vector{Num}) needs vectors of same length"
    return sum([X[i] * Y[i] for i = 1:length(X)])
end

@variables a b c d e f
@variables θ
e_1 = [1,0,0]
e_2 = [cos(θ), a, b]
e_3 = [cos(θ), c, d]
e_4 = [cos(θ), e, f]

e_1 * e_2

equations = [
    e_1 * e_1 ~ 1,
    e_2 * e_2 ~ 1,
    e_3 * e_3 ~ 1,
    e_4 * e_4 ~ 1,
    e_1 * e_2 ~ cos(θ),
    e_1 * e_3 ~ cos(θ),
    e_1 * e_4 ~ cos(θ),
    e_2 * e_3 ~ cos(θ),
    e_2 * e_4 ~ cos(θ),
    e_3 * e_4 ~ cos(θ)
]

solve_system_eq(equations, [a, b, c, d, e, f])
=#


################################################################################
# using PyCall
# using SymPy

# sp = pyimport("SymPy")

# x,y = Sym("x,y")
