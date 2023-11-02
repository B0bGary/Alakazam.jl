module Alakazam

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


export pop_metric, contract_metrics, eliminate_kronecker,eliminate_epsilon,term_weight,singleterm
export get_weighted_terms,drop_weight_above,×,update_indices,get_unused_indices
export get_free_indices,get_free_indices_nopos,get_dummies_positional,get_index_sets
export get_dummies,rename_dummies,update_ext_indices,update_dummies,simplify,merge_sum,grading,substitute
export to_latex,replace,hasIndices,operator_weight,same_operator,partially_project_component
export add_derivative_action,index_to_pattern,same_tensor_mod_pattern,apply_derivative
export perform_derivative_action, expand_derivative,product_rule,to_lc_expression

export IndexSet,Index,DummyPatternIndex,Coordinate,DummyPatternTensor,TensorSuperType
export TotallyAntisymmetricSuperType,TotallySymmetricSuperType,Tensor,SymmetricTensor
export AntisymmetricTensor,TensorTerm,TensorExpression,EpsilonTensor,EpsilonTensorSuperType
export KroneckerDelta,KroneckerSuperType,Operator, OperatorSuperType
export LeibnizOperatorSuperType,Metric,SymMetricSuperType,AntisymMetricSuperType
export EpsilonTensorSuperType,EpsilonTensor,YoungTableau,DerivativeSuperType
export PartialSuperType,CovDerivSuperType,Derivative


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

pyimport(pip)
run(`$(PyCall.python) -m pip install networkx`)
run(`$(PyCall.python) -m pip install matplotlib`)
# pyimport("pip").pip("install networkx")
nx=pyimport(networkx)
plt=pyimport(matplotlib.pyplot)

sort_order = Dict()

function custom_sort_order(v::Vector{String})
    global sort_order
    sort_order = Dict(zip(v, 1:length(v)))
end

end