
# include("tensor_term.jl")
include("metric.jl")

abstract type KroneckerSuperType <: TotallySymmetricSuperType end



mutable struct KroneckerDelta<:KroneckerSuperType
    name::String
    indices::Vector{Pair{I, IndexPosition}}where {I<:IndexSuperType}
    weights::Dict{String,Number}
    function_of::Vector{Coordinate}
    tableaux::Vector{YoungTableau}

    function KroneckerDelta(indices::Vector{Pair{I, IndexPosition}}) where {I<:IndexSuperType}#lets you have scalars with no indices
        if indices[1][2]==indices[2][2]
            throw(InvalidKroneckerIndices())
        end
        if indices[1][1].index_set!=indices[2][1].index_set
            throw(InvalidKroneckerIndices())
        end
        new("δ",indices,Dict{String,Number}(),[], [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end
 end
