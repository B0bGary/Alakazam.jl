include("tensor_term.jl")


mutable struct Operator <:TensorSuperType
    name::String
    operand::TensorExpression
    operator_indices::Vector{Pair{Index, IndexPosition}}
    operand_indices::Vector{Pair{Index, IndexPosition}}
    indices::Vector{Pair{Index, IndexPosition}}
    
    function Operator(name::String, operand::TensorExpression, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}())
        inds=operand.summands[1][2].indices
        new(name, operand, operator_indices, inds, vcat(operator_indices, inds))
    end
    function Operator(name::String, operand::TensorSuperType, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}())
        
        new(name, TensorExpression(operand), operator_indices, operand.indices, vcat(operator_indices, operand.indices))
    end
    function Operator(name::String, operand::TensorTerm, operator_indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}())
        new(name, TensorExpression(operand), operator_indices, operand.indices, vcat(operator_indices, operand.indices))
    end
 end



 function Base.show(io::IO, T::Operator)
    print(io, T.name)
    if !isempty(T.indices)
        lowers=""
        uppers=""
        lowersempt::Bool=1
        uppersempt::Bool=1
        for (ind, pos) in T.indices
            if pos.position==0
                uppers=uppers*ind.name
                uppersempt=0
                lowers=lowers * "*"
            else
                uppers=uppers * "*"
                lowersempt=0
                lowers=lowers * ind.name 
            end
        end
        if !uppersempt
            print(io, "^{"*uppers*"}")
        end
        if !lowersempt
            print(io, "_{"*lowers*"}")
        end
    end
    print(io,"[[")
    print(io, T.operand)
    print(io,"]]")
end