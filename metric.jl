include("tensor_term.jl")

abstract type MetricSuperType <: TensorSuperType end
abstract type AntiSymMetricSuperType <: MetricSuperType end
abstract type SymMetricSuperType <: MetricSuperType end

mutable struct Metric <:SymMetricSuperType
    name::String
    indices::Vector{Pair{Index, IndexPosition}}

    function Metric(name::String="η", indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}())
        new(name, indices)
    end
 end

