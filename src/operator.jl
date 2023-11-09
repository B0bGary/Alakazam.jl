# include("tensor_term.jl")

abstract type OperatorSuperType <: TensorSuperType end

abstract type LeibnizOperatorSuperType <: OperatorSuperType end


mutable struct Operator <:OperatorSuperType
    name::String
    operand::TensorExpression
    operator_indices::Vector{Pair{IndexSuperType, IndexPosition}} 
    operand_indices::Vector{Pair{IndexSuperType, IndexPosition}} 
    indices::Vector{Pair{IndexSuperType, IndexPosition}} 
    function Operator(name::String, operand::TensorExpression, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}())where {I<:IndexSuperType}
        inds=operand.summands[1][2].indices
        new(name, operand, operator_indices, inds, vcat(operator_indices, inds))
    end
    function Operator(name::String, operand::TensorSuperType, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}())where {I<:IndexSuperType}
        
        new(name, TensorExpression(operand), operator_indices, operand.indices, vcat(operator_indices, operand.indices))
    end
    function Operator(name::String, operand::TensorTerm, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}()) where {I<:IndexSuperType}
       new(name, TensorExpression(operand), operator_indices, operand.indices, vcat(operator_indices, operand.indices))
    end
 end
