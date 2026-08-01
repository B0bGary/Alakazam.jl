# include("operator.jl")

#we'll have a derivative as a funciton of multiple coords. this will be used to check if the derivative will annihilate tensors
#we'll then have a subtypes Covariant and NthDeriv to specify if the the coords are across different terms, or repeated derivatives

"""
Generic derivative object. Behaves covariant, i.e. non-commuting.
"""
mutable struct Derivative <: CovDerivSuperType
    data::OperatorData
    coord::Vector{Coordinate}
    Derivative(name::String, operand::Tensorial, coord::Coordinate, args...) =
        new(OperatorData(name, operand, args...), Coordinate[coord])
    Derivative(name::String, operand::Tensorial, coord::Vector{Coordinate}, args...) =
        new(OperatorData(name, operand, args...), (coord))
    Derivative(name::String, coord::Vector{Coordinate}, args...) =
        new(OperatorData(name, zero(), args...), (coord))
    Derivative(name::String, coord::Coordinate, args...) =
        new(OperatorData(name, zero(), args...), Coordinate[coord])
    Derivative(d::OperatorData, coord::Vector{Coordinate}) = new(d, coord)
end

"""
Generic (anti)commuting partial derivatives.
"""
mutable struct PartialDerivative <: PartialSuperType
    data::OperatorData
    coord::Vector{Coordinate}   # with respect to; shared between copies, never mutated in place

    PartialDerivative(name::String, operand::Tensorial, coord::Coordinate, args...) =
        new(OperatorData(name, TensorExpression(operand), args...), Coordinate[coord])

    PartialDerivative(name::String, operand::Tensorial, coord::Vector{Coordinate}, args...) =
        new(OperatorData(name, operand, args...), coord)

    PartialDerivative(name::String, coord::Coordinate, args...) =
        new(OperatorData(name, zero(), args...), Coordinate[coord])

    PartialDerivative(name::String, coord::Vector{Coordinate}, args...) =
        new(OperatorData(name, zero(), args...), coord)

    PartialDerivative(d::OperatorData, coord::Vector{Coordinate}) = new(d, coord)
end

#shorthand
const PD = PartialDerivative
const Partial = PartialDerivative

