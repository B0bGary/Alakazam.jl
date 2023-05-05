include("tensor_term.jl")


abstract type SymMetricSuperType <: TotallySymmetricSuperType end
abstract type AntisymMetricSuperType <: TotallyAntisymmetricSuperType end


mutable struct Metric <:SymMetricSuperType
    name::String
    indices::Vector{Pair{Index, IndexPosition}}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}

    function Metric(name::String="η", indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}())
        new(name, indices,Vector{Coordinate}(),Dict{String, Number}())
    end
    function Metric(indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}())
        new("η", indices,Vector{Coordinate}(),Dict{String, Number}())
    end
 end

