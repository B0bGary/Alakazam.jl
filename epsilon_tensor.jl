
include("tensor_term.jl")
include("metric.jl")

abstract type EpsilonTensorSuperType <: AntisymMetricSuperType end


#levi civita tensor
mutable struct EpsilonTensor<:EpsilonTensorSuperType
    name::String
    indices::Vector{Pair{Index, IndexPosition}}
    weights::Dict{String,Number}
    function_of::Vector{Coordinate}

    function EpsilonTensor(indices::Vector{Pair{Index, IndexPosition}}) #lets you have scalars with no indices
        if indices[1][2]!=indices[2][2]
            throw(InvalidEpsilonIndices())
        end
        if indices[1][1].index_set!=indices[2][1].index_set
            println(indices[1][1])
            println(indices[2][1])
            throw(InvalidEpsilonIndices())
        end
        new("ϵ",indices,Dict{String,Number}(),[])
    end
 end
