
# include("tensor_term.jl")
include("metric.jl")

abstract type KroneckerSuperType <: TotallySymmetricSuperType end



mutable struct KroneckerDelta<:KroneckerSuperType
    name::String
    indices::Vector{Pair{Index, IndexPosition}}
    weights::Dict{String,Number}
    function_of::Vector{Coordinate}
    function KroneckerDelta(indices::Vector{Pair{Index, IndexPosition}}) #lets you have scalars with no indices
        if indices[1][2]==indices[2][2]
            throw(InvalidKroneckerIndices())
        end
        if indices[1][1].index_set!=indices[2][1].index_set
            throw(InvalidKroneckerIndices())
        end
        new("δ",indices,Dict{String,Number}(),[])
    end
 end
