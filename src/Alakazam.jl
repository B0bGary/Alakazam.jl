module Alakazam
#deps
using Symbolics
using OrderedCollections
using DataStructures
using PrecompileTools
using InteractiveUtils
# using StatsBase
using Latexify
using Combinatorics
using SymbolicUtils
using JSON
import Base.:*
import Base: zero
import Base.:+
import Base.:copy
import Base.:/
import Base.:(==)
import Base.:-
import Base.:replace
import Base.//
import Symbolics: simplify
import Symbolics: terms

##########################
#basic objs
##########################
export IndexSuperType, Index, IndexSet, IndexPosition, IndexGrading, DummyPatternIndex, Coordinate, @indices
#index functions
export get_index, det, in_default_pos, op_pos, dimension, safe_pop!
export raise_all, lower_all, flip_all
export SymorNum
#global properties
export SET_REGISTRY, COORD_REGISTRY, sort_order, get_sort_order, custom_sort_order!
export add_global_property!, get_global_values
export clear_global!, clear_commutators!, clear_sort_order!, clear_covariant!, clear_coords!, clear_indexsets!
export lower, upper, bosonic, fermionic, ↑, ↓
export show_function_of, function_of_on


export construct_tensorterm_graph, find_weighted_square, same_tensorterm_mod_syms_grade

export derivative_actions

export to_latex

##########################
#base tensors
export TensorData, TensorDataSuperType, Tensorial, TensorTerm, Tensor, TensorExpression, TensorSuperType
export SymmetricTensor, AntisymmetricTensor, RiemannTensor
export Boson, Fermion
export @Tensor, @T_str

##########################
#all supertypes
##########################
export MetricSuperType, SymMetricSuperType, AntisymMetricSuperType, HermMetricSuperType
export EpsilonTensorSuperType, KroneckerSuperType
export MatrixSuperType, PauliSuperType, GammaSuperType, Gamma5SuperType
export OperatorSuperType, LeibnizOperatorSuperType, DerivativeSuperType, PartialSuperType
export CovDerivSuperType, VariationalSuperType, LieBracketSuperType
export HomomorphismSuperType, AntimorphismSuperType, InvolutionSuperType
export CyclicSuperType, TraceSuperType
export NumberSet, LieAlgebraSuperType
export GeneratorDerivSuperType

##########################
#basic tensor properties
##########################
export data
#expression
export terms, coefficients, coef, summands, is_const, first_summand, first_coef, singleterm, hasIndices, isZeroTens, isTensorTerm, isTensor
export name, function_of, indices, weights, tableaux, grading, labels
export set_name!, set_indices!, set_indices_deep!, update_indices!, set_index!
export free_indices, free_indices_nopos, get_index_sets, dummies, dummies_positional
export isAntisymmetric, isSymmetric
# zero() tensor expr
export coef_string

#serialisation
export to_JSON, from_JSON, save_header, save, load, load_header, @load

##########################
#operators
##########################
export Operator, OperatorData, op_data, op_indices, operand, op_function_of, op_weights, op_tableaux, op_labels, operand_indices, op_grading
export set_operand!, set_op_indices!


#exceptions
export FreeIndexMismatch, OverusedIndices, InvalidKroneckerIndices, InvalidEpsilonIndices, OperatorPrecedenceError, OutOfIndices, InvalidReplacementPattern
export MixedWeightOperandTerms, SameLevelDummyError, InvalidDerivativeBasis, InvalidDifferentialForm, InvalidCommutatorDefinition


##########################
#tableau
##########################
export YoungTableau, generate_symmetric_tableau, generate_antisymmetric_tableau, generate_riemann_tableau, young_projection
export young_projector_size, multiterm_reduce

##########################
#special tensors
##########################
export Metric, KroneckerDelta, EpsilonTensor, Delta
export is_riemann_metric
export Mat, Pauli, Gamma, Gamma5
export reduce_matrix, contract_sigma, @mrule, @matrix
export contract_adjacent#needed for macro

export SymplecticForm

##########################
#basic ops
##########################
export pop_metric, contract_metrics, eliminate_kronecker, eliminate_epsilon
export ≅, cong, get_weighted_terms, term_weight, drop_weight_above, drop_weight_below, ×, append, rename_dummies, rename_dummies!
export get_unused_indices
export simplify, hard_simplify, simp_ops_n_sym, canonicalise, canonicalize,collect_terms

##########################
#comrels
export Commutator, remove_commutator!, get_commutator, @commutator, @anticommutator
export add_commutator!#needed for macro

##########################
#replacements
export replace_pattern

###########################
#diffgeo etc
export isForm, ⋆, star, ⨼, ι, wedge, ∧, contraction


##########################
#derivatives
##########################
export RegisteredExpansion, Derivative, PartialDerivative
export @action, Partial, PD, derive, expand_derivative, derive, coord, apply_derivative, product_rule, symmetrise_partials
export @covariant
export add_derivative_action#needed for macro
export derivative_algebra, set_operand_deep!, commute, CovariantCommutator
export commute_right

##########################
#liebracket
export Lie, @Lie, first, second, swap, expand_commutators

##########################
#variational
export EL, VariationalDerivative




##########################
#superspace
export Superspace, generate_grassmann_product, GrassmannStructure, independent_structures, explain_structure, render_expression, generate_independant_structure
export generate_superfield, generate_superfield_naive, project_to_component

export measure, integrate

export render_compact


export to_dumless_same_tensor_mod_syms #for tests...



########

#### algebras
export generate_conformal_algebra, LieGenerator
export SLAlgebra, SOAlgebra, SUAlgebra, SpAlgebra
export epsilon, metric





export symmetrise_partials#fix this!!!!


export CoeffPair, fresh_indices
export @syms #reexport from symbolics...

# ###########still needs to be checked and tested
# export get_weighted_terms, drop_weight_above, ×, update_indices, get_unused_indices
# export ≅
# export free_indices, free_indices_nopos, dummies_positional, get_index_sets
# export dummies, rename_dummies, update_ext_indices, update_dummies, simplify, merge_sum, grading
# export to_latex, replace, hasIndices, operator_weight, same_operator, partially_project_component
# export add_derivative_action, index_to_pattern, same_tensor_mod_pattern, apply_derivative
# export perform_derivative_action, expand_derivative, product_rule, to_lc_expression

# export IndexSuperType, IndexSet, Index, DummyPatternIndex, Coordinate, TensorSuperType
# export Tensor, @Tensor, RiemannTensor, SymmetricTensor, AntisymmetricTensor
# #,TotallyAntisymmetricSuperType,TotallySymmetricSuperType,SymmetricTensor
# export TensorTerm, TensorExpression, EpsilonTensor, EpsilonTensorSuperType
# export KroneckerDelta, KroneckerSuperType, Operator, OperatorSuperType
# export LeibnizOperatorSuperType, Metric, SymMetricSuperType, AntisymMetricSuperType
# export EpsilonTensorSuperType, EpsilonTensor, YoungTableau, DerivativeSuperType
# export PartialSuperType, CovDerivSuperType, Derivative, PartialDerivative
# export lower, upper, sort_order, custom_sort_order!,  zero
# export project_component
# export same_tensorterm_mod_syms, same_tensor_mod_syms
# export hard_simplify
# export hard_equals

#core typelist
include("core/supertyping.jl")
#core components
#algebra
include("algebra/algebras.jl")
include("core/indices.jl")
include("tensors/young.jl")
include("core/custom_exceptions.jl")
#base components
include("tensors/tensorial.jl")
include("operators/operator.jl")
include("tensors/tensor_arithmetic.jl")
include("tensors/metric.jl")
include("commutation/matrices.jl")
include("commutation/commutators.jl")
include("tensors/tensor_props.jl")
include("tensors/tensor_ops.jl")
include("diffgeo/forms.jl")
include("operators/operator_ops.jl")
include("operators/derivative.jl")
include("operators/derivative_ops.jl")
include("core/tensor_replacements.jl")
include("core/simplify_tableaux.jl")
include("operators/variations.jl")
#lie bracketss
include("operators/lie.jl")
#superspace
include("superspace/superspace.jl")
include("superspace/superfield.jl")
# include("superspace/superfield2.jl")
include("superspace/grassmann_structures.jl")
#optional lightcone
include("superspace/lightcone.jl")
# using .Lightcone
include("algebra/generators.jl")



"""Should we print the function_of coordinates? Includes both @show and to_latex printing."""
function_of_on=false
show_function_of(b::Bool) = begin
    global function_of_on
    function_of_on=b
end

#allows user to define a custom canonical sort order (other than alphabetic) for tensors
const sort_order = Dict{String,Int}()
# const NCOPY = Ref(0); const NSET = Ref(0)
# export NCOPY,NSET

"""Global user workspace values. Collated for serialisation. User can add to this list."""
const global_values=Dict("sort_order"=>sort_order, "coordinates"=>COORD_REGISTRY, "index_sets"=>SET_REGISTRY, "covariant_expansions"=>DERIVATIVE_EXPANSIONS,
    "commutators"=>commutators, "derivative_actions"=>derivative_actions, "derivative_commutators"=>derivative_algebra)

"""
Clears all global registrations, like commutators, sort order, index sets...
"""
function clear_global!()
    global global_values
    for (k, v) in global_values
        empty!(v)
    end
end

function get_global_values()
    return global_values
end
function clear_sort_order!()
    empty!(sort_order)
end
function custom_sort_order!(v::Vector{String})
    empty!(sort_order)
    merge!(sort_order, Dict(zip(v, eachindex(v))))
    return sort_order
end
function custom_sort_order!(v::AbstractDict{String,Int})
    empty!(sort_order)
    merge!(sort_order, v)
    return sort_order
end
function get_sort_order()
    return sort_order
end


include("core/serialisation.jl")


"""
Add something to the list of global properties associated with the environment. 
These will be serialised when exporting to JSON. 
The user may have to define how to import if there are nontrivial structures.
"""
function add_global_property!(key::String, value)
    global_values[key]=value
end


@setup_workload begin
    @compile_workload begin
        include("core/precompile.jl")
    end
end



const BANNER_DIR = joinpath(@__DIR__, "..", "data")
function show_banner()

    width = try
        parse(Int, get(ENV, "COLUMNS", string(displaysize(stdout)[2])))
    catch
        displaysize(stdout)[2]
    end

    file = if width >= 126
        joinpath(BANNER_DIR, "alakazam_banner_right.ansi")
    else
        joinpath(BANNER_DIR, "alakazam_banner_bottom.ansi")
    end

    print(read(file, String))
    print("\e[0m")

end
function __init__()

    if isinteractive()
        show_banner()
    end

end
end
