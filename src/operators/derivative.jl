# include("operator.jl")

#we'll have a derivative as a funciton of multiple coords. this will be used to check if the derivative will annihilate tensors
#we'll then have a subtypes Covariant and NthDeriv to specify if the the coords are across different terms, or repeated derivatives

"""
Generic derivative object. Behaves covariant, i.e. non-commuting.
"""
mutable struct Derivative <: CovDerivSuperType
    data::OperatorData
    coord::Set{Coordinate}#w/r to
    function Derivative(name::String, operand::Tensorial, coord::Coordinate, args...)
        new(OperatorData(name, operand, args...), Set{Coordinate}((coord,)))
    end
    function Derivative(name::String, operand::Tensorial, coord::Set{Coordinate}, args...)
        new(OperatorData(name, operand, args...), coord)
    end
    function Derivative(name::String, coord::Set{Coordinate}, args...)
        new(OperatorData(name,  zero(), args...), coord)
    end
    function Derivative(name::String, coord::Coordinate, args...)
        new(OperatorData(name,  zero(), args...), Set((coord,)))
    end
    function Derivative(name::String, coord::Vector{Coordinate}, args...)
        new(OperatorData(name,  zero(), args...), Set(coord))
    end
    function Derivative(name::String, operand::Tensorial, coord::Vector{Coordinate}, args...)
        new(OperatorData(name, operand, args...), Set(coord))
    end
    function Derivative(d::OperatorData, coord::Set{Coordinate})
        new(d, coord)
    end
    function Derivative(d::OperatorData, coord::Vector{Coordinate})
        new(d, Set(coord))
    end
end

"""
Generic (anti)commuting partial derivatives.
"""
mutable struct PartialDerivative <: PartialSuperType
    data::OperatorData
    coord::Set{Coordinate}#w/r to
    function PartialDerivative(name::String, operand::Tensorial, coord::Coordinate, args...)
        new(OperatorData(name, TensorExpression(operand), args...), Set{Coordinate}((coord,)))
    end
    function PartialDerivative(name::String, operand::Tensorial, coord::Set{Coordinate}, args...)
        new(OperatorData(name, operand, args...), coord)
    end
    function PartialDerivative(name::String, coord::Set{Coordinate}, args...)
        new(OperatorData(name,  zero(), args...), coord)
    end
    function PartialDerivative(name::String, coord::Vector{Coordinate}, args...)
        new(OperatorData(name,  zero(), args...), Set(coord))
    end
    function PartialDerivative(name::String, coord::Coordinate, args...)
        new(OperatorData(name,  zero(), args...), Set((coord,)))
    end
    function PartialDerivative(name::String, operand::Tensorial, coord::Vector{Coordinate}, args...)
        new(OperatorData(name, operand, args...), Set(coord))
    end
    function PartialDerivative(d::OperatorData, coord::Set{Coordinate})
        new(d, coord)
    end
    function PartialDerivative(d::OperatorData, coord::Vector{Coordinate})
        new(d, Set(coord))
    end
end

#shorthand
const PD = PartialDerivative
const Partial = PartialDerivative

