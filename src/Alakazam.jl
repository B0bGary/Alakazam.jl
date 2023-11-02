module Alakazam

using Pkg
Pkg.instantiate()
using Symbolics
using OrderedCollections
using PyCall
using StatsBase
import SymbolicUtils
import Base.:*
import Base.:+
import Base.:copy
import Base.:/
import Base.:(==)
import Base.:-
import Base.:replace
import DataStructures
include("indices.jl")
include("young.jl")
include("custom_exceptions.jl")
include("tensor_term.jl")
include("tensor_ops.jl")
include("metric.jl")
include("operator.jl")
include("operator_ops.jl")
include("derivative.jl")
include("derivative_ops.jl")
include("kronecker.jl")
include("epsilon_tensor.jl")
include("tensor_replacements.jl")
include("simplify_tableaux.jl")
include("lightcone.jl")

@pyimport pip
run(`$(PyCall.python) -m pip install networkx`)
run(`$(PyCall.python) -m pip install matplotlib`)
# pyimport("pip").pip("install networkx")
@pyimport networkx as nx
@pyimport matplotlib.pyplot as plt

sort_order = Dict()

function custom_sort_order(v::Vector{String})
    global sort_order
    sort_order = Dict(zip(v, 1:length(v)))
end

end