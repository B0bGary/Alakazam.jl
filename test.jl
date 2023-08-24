
mutable struct Tens<:TensorSuperType
    name::String
    indices::Vector{Pair{I, IndexPosition}} where {I<:IndexSuperType}

    function Tens()
        new("T")
    end

    function Tens(name::String)
        new(name, Vector{Pair{IndexSuperType, IndexPosition}}()) 
    end

    function Tens(indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}()) where {I<:IndexSuperType}
        new("T2", indices)
    end
end