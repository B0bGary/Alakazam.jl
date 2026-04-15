module Alakazam

#deps
using InteractiveUtils
using Symbolics
using OrderedCollections
using PyCall
using StatsBase
using Combinatorics
using Serialization
using SymbolicUtils
import Base.:*
import Base.:+
import Base.:copy
import Base.:/
import Base.:(==)
import Base.:-
import Base.:replace
import DataStructures


function ensure_installed(package::String)
    try
        pyimport(package)
    catch
        run(`$(PyCall.python) -m pip install $package`)
    end
end

# Ensure required packages are installed
ensure_installed("networkx")
ensure_installed("matplotlib")
# Import the modules

# println("Checking if nx is correctly imported:", nx)
const nx= PyNULL()
const plt= PyNULL()

function __init__()
    copy!(nx,pyimport("networkx"))
    copy!(plt,pyimport("matplotlib.pyplot"))
end

#usable functions and types


export pop_metric, contract_metrics, eliminate_kronecker,eliminate_epsilon,term_weight,singleterm
export get_weighted_terms,drop_weight_above,×,update_indices,get_unused_indices
export ≅
export get_free_indices,get_free_indices_nopos,get_dummies_positional,get_index_sets
export get_dummies,rename_dummies,update_ext_indices,update_dummies,simplify,merge_sum,grading,substitute
export to_latex,replace,hasIndices,operator_weight,same_operator,partially_project_component
export add_derivative_action,index_to_pattern,same_tensor_mod_pattern,apply_derivative
export perform_derivative_action, expand_derivative,product_rule,to_lc_expression

export IndexSet,Index,DummyPatternIndex,Coordinate,DummyPatternTensor,TensorSuperType
export Tensor,@Tensor,RiemannTensor,SymmetricTensor,AntisymmetricTensor
#,TotallyAntisymmetricSuperType,TotallySymmetricSuperType,SymmetricTensor
export TensorTerm,TensorExpression,EpsilonTensor,EpsilonTensorSuperType
export KroneckerDelta,KroneckerSuperType,Operator, OperatorSuperType
export LeibnizOperatorSuperType,Metric,SymMetricSuperType,AntisymMetricSuperType
export EpsilonTensorSuperType,EpsilonTensor,YoungTableau,DerivativeSuperType
export PartialSuperType,CovDerivSuperType,Derivative
export lower,upper,sort_order,custom_sort_order,zero
export project_component
export same_tensorterm_mod_syms, same_tensor_mod_syms
export hard_simplify
export hard_equals

export ⋆,star


export read_from_file, write_to_file


include("indices.jl")
include("young.jl")
include("custom_exceptions.jl")
include("tensor_term.jl")
# include("metric.jl")
# include("kronecker.jl")
# include("epsilon_tensor.jl")
include("tensor_ops.jl")
# include("metric.jl")
include("operator.jl")
include("operator_ops.jl")
include("derivative.jl")
include("derivative_ops.jl")
# include("kronecker.jl")
# include("epsilon_tensor.jl")
include("tensor_replacements.jl")
include("simplify_tableaux.jl")
include("lightcone.jl")
include("superspace.jl")




#allows user to define a custom canonical sort order (other than alphabetic) for tensors
sort_order = Dict()

function custom_sort_order(v::Vector{String})
    global sort_order
    sort_order = Dict(zip(v, 1:length(v)))
end


types = Union{TensorSuperType, TensorTerm, TensorExpression}

function write_to_file(name::String, t::T) where T<:types
    return serialize(name*".akzm", t)
end


function read_from_file(name::String) 
    return deserialize(name*".akzm")
end


# # Get all subtypes of TensorSuperType
# function all_subtypes(T)
#     function collect_subtypes(T)
#         result = [T]  # Include T itself
#         for subtype in subtypes(T)  # Use subtypes directly
#             push!(result, collect_subtypes(subtype)...)
#         end
#         return result
#     end

#     return collect_subtypes(T)
# end

# # Function to serialize all subtypes of TensorSuperType
# function serialize_subtypes()
#     # println("Trying to serialize subtypes")

#     for ttyp in types
#         # println("Handling type: $ttyp")

#         for st in all_subtypes(ttyp)
#             # println("Subtype: $st")

#             try
#                 # Dynamically define the method using `@eval`
#                 @eval begin
#                     StructTypes.StructType(::Type{$st}) = StructTypes.OrderedStruct()
#                 end
#                 # println("Registered type: $st")
#             catch e
#                 println("Serialisation Error: Failed to register type $st: $e")
#             end
#         end
#     end
# end

# serialize_subtypes()

# # Serialize function: Convert instance to JSON, storing type information
# function serialize(obj)
#     try
#         return JSON3.write(Dict("type" => string(typeof(obj)), "data" => obj))
#     catch e
#         println("Serialization failed: $e")
#         return nothing
#     end
# end

# # Deserialize function: Convert JSON back to the correct type automatically
# function deserialize(json_str::String)
#     try
#         json_data = JSON3.read(json_str, TensorExpression)  # Read JSON as a dictionary
#         type_name = json_data["type"]  # Extract the type name
#         obj_type = get_type_from_string(type_name)  # Resolve the type

#         if obj_type === nothing
#             println("Unknown type: $type_name")
#             return nothing
#         end

#         return JSON3.read(json_data["data"], obj_type)  # Deserialize into the correct type
#     catch e
#         println("Deserialization failed: $e")
#         return nothing
#     end
# end

# # Helper function: Map type name from string to actual Julia type
# function get_type_from_string(type_name::String)
#     try
#         type_symbol = Symbol(type_name)  # Convert string to symbol
#         return getfield(Main, type_symbol)  # Get the type from the Main module
#     catch e
#         println("Type resolution failed: $e")
#         return nothing
#     end
# end

# # function to_file(file_name::String, tens::T) where {T<:types}
# #     open(file_name*".json", "w") do file
# #         JSON3.write_json(file, tens)
# #     end
# # end




# # function from_file(io)
# #     dict = JSON3.read_json(io, Dict)

# #     type_symbol = Symbol(type_name)  # Convert string to symbol
# #     return getfield(TensorTypes, type_symbol)  


# #     if haskey(dict, :value)
# #         return TensorTerm(dict[:value])
# #     elseif haskey(dict, :expression)
# #         return TensorExpression(dict[:expression])
# #     end
# # end


end