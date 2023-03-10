include("tensor_term.jl")
mutable struct Metric <:TensorSuperType
    name::String
    indices::Vector{Pair{Index, IndexPosition}}

    function Metric(name::String="η", indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}())
        new(name, indices)
    end
 end

