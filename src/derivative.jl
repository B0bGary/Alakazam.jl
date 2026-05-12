# include("operator.jl")

abstract type DerivativeSuperType <: LeibnizOperatorSuperType end

abstract type PartialSuperType <: DerivativeSuperType end

abstract type CovDerivSuperType <: DerivativeSuperType end

#we'll have a derivative as a funciton of multiple coords. this will be used to check if the derivative will annihilate tensors
#we'll then have a subtypes Covariant and NthDeriv to specify if the the coords are across different terms, or repeated derivatives
mutable struct Derivative <: DerivativeSuperType
    name::String
    operand::TensorExpression
    coord::Vector{T} where {T<:Coordinate}
    operator_indices::Vector{Pair{IndexSuperType,IndexPosition}}
    operand_indices::Vector{Pair{IndexSuperType,IndexPosition}}
    indices::Vector{Pair{IndexSuperType,IndexPosition}}
    function_of::Vector{K} where {K<:Coordinate}
    weights::Dict{String,Number}
    tableaux::Vector{YoungTableau}




    # Single constructor method
    function Derivative(name::String, operand, coord::Vector{Coordinate}, args::Vararg{Any})
        operand = TensorExpression(operand)  # Ensures operand is a TensorExpression

        # Default values
        operator_indices = Vector{Pair{IndexSuperType,IndexPosition}}()
        function_of = Vector{Coordinate}()
        weights = Dict{String,Number}()
        tableaux = Vector{YoungTableau}()

        # Process variable-length arguments
        for arg in args
            if arg isa Vector{Pair{I,IndexPosition}} where I<:IndexSuperType
                append!(operator_indices, arg)
            elseif arg isa Pair{I,IndexPosition} where I<:IndexSuperType
                push!(operator_indices, arg)
            elseif arg isa Coordinate
                function_of = [arg]
            elseif arg isa Vector{Coordinate}
                function_of = arg
            elseif arg isa Dict{String,Number}
                weights = merge(weights, arg)
            elseif arg isa YoungTableau
                tableaux = [arg]
            elseif arg isa Vector{YoungTableau}
                tableaux = arg
            end
        end

        append!(function_of, get_operand_function_of(operand))

        # Ensure operand_indices are taken from the operand
        operand_indices = operand.summands[1][2].indices
        indices = vcat(operator_indices, operand_indices)

        # Construct and return the Derivative object
        return new(name, operand, coord, operator_indices, operand_indices, indices, function_of, weights, tableaux)
    end
end

mutable struct PartialDerivative <: PartialSuperType
    name::String
    operand::TensorExpression
    coord::Vector{T} where {T<:Coordinate}
    operator_indices::Vector{Pair{IndexSuperType,IndexPosition}}
    operand_indices::Vector{Pair{IndexSuperType,IndexPosition}}
    indices::Vector{Pair{IndexSuperType,IndexPosition}}
    function_of::Vector{K} where {K<:Coordinate}
    weights::Dict{String,Number}
    tableaux::Vector{YoungTableau}

    # Single constructor method
    function PartialDerivative(name::String, operand, coord::Vector{Coordinate}, args::Vararg{Any})
        operand = TensorExpression(operand)  # Ensures operand is a TensorExpression

        # Default values
        operator_indices = Vector{Pair{IndexSuperType,IndexPosition}}()
        function_of = Vector{Coordinate}()
        weights = Dict{String,Number}()
        tableaux = Vector{YoungTableau}()

        # Process variable-length arguments
        for arg in args
            if arg isa Vector{Pair{I,IndexPosition}} where I<:IndexSuperType
                append!(operator_indices, arg)
            elseif arg isa Pair{I,IndexPosition} where I<:IndexSuperType
                push!(operator_indices, arg)
            elseif arg isa Coordinate
                function_of = [arg]
            elseif arg isa Vector{Coordinate}
                function_of = arg
            elseif arg isa Dict{String,Number}
                weights = merge(weights, arg)
            elseif arg isa YoungTableau
                tableaux = [arg]
            elseif arg isa Vector{YoungTableau}
                tableaux = arg
            end
        end

        append!(function_of, get_operand_function_of(operand))

        # Ensure operand_indices are taken from the operand
        operand_indices = operand.summands[1][2].indices
        indices = vcat(operator_indices, operand_indices)

        # Construct and return the PartialDerivative object
        return new(name, operand, coord, operator_indices, operand_indices, indices, function_of, weights, tableaux)
    end
end

function get_operand_function_of(expr::TensorExpression)
    coords = Coordinate[]
    for (_, term) in expr.summands
        for t in term.terms
            if hasproperty(t, :function_of)
                append!(coords, t.function_of)
            end
        end
    end
    return coords
end


# function Derivative(name::String, operand::TensorExpression, coord::Vector{T}, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}()) where {T<:Coordinate} where {I<:IndexSuperType}
#     inds=operand.summands[1][2].indices
#     new(name, operand, coord, operator_indices, inds, vcat(operator_indices, inds),coord,Vector{YoungTableau}())
# end
# function Derivative(name::String, operand::TensorSuperType, coord::Vector{T}, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}()) where {T<:Coordinate} where {I<:IndexSuperType}
#     new(name, TensorExpression(operand), coord, operator_indices, operand.indices, vcat(operator_indices, operand.indices),coord,Vector{YoungTableau}())
# end
# function Derivative(name::String, operand::TensorTerm, coord::Vector{T}, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}()) where {T<:Coordinate} where {I<:IndexSuperType}
#    new(name, TensorExpression(operand), coord, operator_indices, operand.indices, vcat(operator_indices, operand.indices),coord,Vector{YoungTableau}())

# end
#    function Derivative(name::String, operand::TensorExpression, coord::T, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}()) where {T<:Coordinate} where {I<:IndexSuperType}
#     inds=operand.summands[1][2].indices
#     new(name, operand, [coord], operator_indices, inds, vcat(operator_indices, inds),[coord],operand.weights,Vector{YoungTableau}())
# end
# function Derivative(name::String, operand::TensorSuperType, coord::T, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}()) where {T<:Coordinate} where {I<:IndexSuperType}
#     new(name, TensorExpression(operand), [coord], operator_indices, operand.indices, vcat(operator_indices, operand.indices),[coord],operand.weights,Vector{YoungTableau}())
# end
# function Derivative(name::String, operand::TensorTerm, coord::T, operator_indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}()) where {T<:Coordinate}  where {I<:IndexSuperType}
#    new(name, TensorExpression(operand), [coord], operator_indices, operand.indices, vcat(operator_indices, operand.indices),[coord],operand.weights,Vector{YoungTableau}())
# end


