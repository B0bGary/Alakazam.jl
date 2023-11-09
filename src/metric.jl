include("tensor_term.jl")


abstract type SymMetricSuperType <: TotallySymmetricSuperType end
abstract type AntisymMetricSuperType <: TotallyAntisymmetricSuperType end


mutable struct Metric <:SymMetricSuperType
    name::String
    indices::Vector{Pair{IndexSuperType, IndexPosition}}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}

    function Metric(name::String="η", indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}())where {I<:IndexSuperType}
        new(name, indices,Vector{Coordinate}(),Dict{String, Number}())
    end
    function Metric(indices::Vector{Pair{I, IndexPosition}}=Vector{Pair{IndexSuperType, IndexPosition}}())where {I<:IndexSuperType}
        new("η", indices,Vector{Coordinate}(),Dict{String, Number}())
    end
 end

