include("operator.jl")

abstract type DerivativeSuperType <: LeibnizOperatorSuperType end

abstract type PartialSuperType <: DerivativeSuperType end

abstract type CovDerivSuperType <: DerivativeSuperType end

#we'll have a derivative as a funciton of multiple coords. this will be used to check if the derivative will annihilate tensors
#we'll then have a subtypes Covariant and NthDeriv to specify if the the coords are across different terms, or repeated derivatives

mutable struct Derivative <:DerivativeSuperType
    name::String
    operand::TensorExpression
    coord::Vector{T} where {T<:Coordinate}
    operator_indices::Vector{Pair{Index, IndexPosition}}
    operand_indices::Vector{Pair{Index, IndexPosition}}
    indices::Vector{Pair{Index, IndexPosition}}
    function_of::Vector{K} where {K<:Coordinate}


    function Derivative(name::String, operand::TensorExpression, coord::Vector{T}, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) where {T<:Coordinate}
        inds=operand.summands[1][2].indices
        new(name, operand, coord, operator_indices, inds, vcat(operator_indices, inds),coord)
    end
    function Derivative(name::String, operand::TensorSuperType, coord::Vector{T}, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) where {T<:Coordinate}
        new(name, TensorExpression(operand), coord, operator_indices, operand.indices, vcat(operator_indices, operand.indices))
    end
    function Derivative(name::String, operand::TensorTerm, coord::Vector{T}, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) where {T<:Coordinate} 
       new(name, TensorExpression(operand), coord, operator_indices, operand.indices, vcat(operator_indices, operand.indices))
    
    end
       function Derivative(name::String, operand::TensorExpression, coord::T, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) where {T<:Coordinate}
        inds=operand.summands[1][2].indices
        new(name, operand, [coord], operator_indices, inds, vcat(operator_indices, inds),[coord],operand.weights)
    end
    function Derivative(name::String, operand::TensorSuperType, coord::T, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) where {T<:Coordinate}
        new(name, TensorExpression(operand), [coord], operator_indices, operand.indices, vcat(operator_indices, operand.indices),[coord],operand.weights)
    end
    function Derivative(name::String, operand::TensorTerm, coord::T, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) where {T<:Coordinate} 
       new(name, TensorExpression(operand), [coord], operator_indices, operand.indices, vcat(operator_indices, operand.indices),[coord],operand.weights)
    end
 end


