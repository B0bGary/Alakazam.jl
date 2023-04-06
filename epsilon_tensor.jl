
include("tensor_term.jl")
include("metric.jl")

abstract type EpsilonTensorSuperType <: AntisymMetricSuperType end


#levi civita tensor
mutable struct EpsilonTensor<:EpsilonTensorSuperType
    name::String
    indices::Vector{Pair{Index, IndexPosition}}
    function EpsilonTensor(indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) #lets you have scalars with no indices
        new("ϵ",indices)
    end
 end
