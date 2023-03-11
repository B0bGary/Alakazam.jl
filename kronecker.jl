
include("tensor_term.jl")
include("metric.jl")

abstract type KroneckerSuperType <: SymMetricSuperType end



mutable struct KroneckerDelta<:KroneckerSuperType
    name::String
    indices::Vector{Pair{Index, IndexPosition}}
    function KroneckerDelta(indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) #lets you have scalars with no indices
        new("δ",indices)
    end
 end
