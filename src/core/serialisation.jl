################################################################
# extra complicated functions to handle symbolics coefficients
##############################################################

# ---------- version-robust predicates ----------
const _iscall = isdefined(SymbolicUtils, :iscall) ? SymbolicUtils.iscall : SymbolicUtils.istree
const _issym = isdefined(SymbolicUtils, :issym) ? SymbolicUtils.issym : (x -> x isa SymbolicUtils.Sym)

function _resolve_module(name::Symbol)
    name === :Main && return Main
    name === :Base && return Base
    name === :Core && return Core
    for (_, mod) in Base.loaded_modules
        nameof(mod) === name && return mod
    end
    error("Module `$name` is not loaded.")
end

# Resolve a bare, undotted name (e.g. `Int64`, `sin`) by searching the
# obvious places first, then every loaded module's bindings.
function _resolve_simple_name(name::Symbol)
    for mod in (Core, Base, Main)
        isdefined(mod, name) && return getfield(mod, name)
    end
    for (_, mod) in Base.loaded_modules
        isdefined(mod, name) && return getfield(mod, name)
    end
    error("Could not resolve `$name` in any loaded module.")
end

function _resolve_qualified(path::AbstractString)
    parts = Symbol.(split(path, "."))
    if length(parts) == 1
        return _resolve_simple_name(parts[1])
    end
    obj = _resolve_module(parts[1])
    for p in parts[2:end]
        obj = getfield(obj, p)
    end
    obj
end

function _resolve_type_expr(ex)
    if ex isa Symbol
        _resolve_simple_name(ex)
    elseif ex isa Expr && ex.head === :(.)
        _resolve_qualified(string(ex))
    elseif ex isa Expr && ex.head === :curly
        base = _resolve_type_expr(ex.args[1])
        base{(_resolve_type_expr(a) for a in ex.args[2:end])...}
    else
        error("Unsupported type expression: $ex")
    end
end

_resolve_type(s::AbstractString) = _resolve_type_expr(Meta.parse(s))
_qualified_name(f) = string(nameof(parentmodule(f)), ".", nameof(f))

# structural (registry-free) tree serialize/deserialize, used as the "value" payload
# function _tree_to_dict(x)
#     if _issym(x)
#         Dict("kind"=>"sym", "name"=>String(nameof(x)),
#             "symtype"=>string(SymbolicUtils.symtype(x)))
#     elseif _iscall(x)
#         Dict("kind"=>"term", "op"=>_qualified_name(SymbolicUtils.operation(x)),
#             "args"=>[_tree_to_dict(a) for a in SymbolicUtils.arguments(x)])
#     elseif x isa Number
#         Dict("kind"=>"number", "type"=>string(typeof(x)), "value"=>x)
#     else
#         error("Don't know how to serialise $(typeof(x)): issym=$(_issym(x)) iscall=$(_iscall(x)) isnumber=$(x isa Number)")
#     end
# end
#latest symbolicsutils version
function _tree_to_dict(x)
    if _issym(x)
        Dict("kind"=>"sym", "name"=>String(nameof(x)),
            "symtype"=>string(SymbolicUtils.symtype(x)))
    elseif SymbolicUtils.isconst(x)
        v = SymbolicUtils.unwrap_const(x)
        Dict("kind"=>"number", "type"=>string(typeof(v)), "value"=>v)
    elseif _iscall(x)
        Dict("kind"=>"term", "op"=>_qualified_name(SymbolicUtils.operation(x)),
            "args"=>[_tree_to_dict(a) for a in SymbolicUtils.arguments(x)])
    elseif x isa Number
        Dict("kind"=>"number", "type"=>string(typeof(x)), "value"=>x)
    else
        error("Don't know how to serialise $(typeof(x))")
    end
end

function _tree_from_dict(d)
    kind = d["kind"]
    if kind == "sym"
        # SymbolicUtils.Sym{_resolve_type(d["symtype"])}(Symbol(d["name"]))
        SymbolicUtils.Sym{SymbolicUtils.SymReal}(Symbol(d["name"]);
            type=_resolve_type(d["symtype"]),
            shape=UnitRange{Int64}[])
    elseif kind == "term"
        f = _resolve_qualified(d["op"])
        f((_tree_from_dict(a) for a in d["args"])...)
    elseif kind == "number"
        T = _resolve_type(d["type"])
        v = d["value"]

        if T <: Complex
            return T(v["re"], v["im"])
        else
            return convert(T, v)
        end
    else
        error("Unknown kind $kind")
    end
end


# ░██████████ ░██    ░██ ░█████████    ░██████   ░█████████  ░██████████
# ░██          ░██  ░██  ░██     ░██  ░██   ░██  ░██     ░██     ░██
# ░██           ░██░██   ░██     ░██ ░██     ░██ ░██     ░██     ░██
# ░█████████     ░███    ░█████████  ░██     ░██ ░█████████      ░██
# ░██           ░██░██   ░██         ░██     ░██ ░██   ░██       ░██
# ░██          ░██  ░██  ░██          ░██   ░██  ░██    ░██      ░██
# ░██████████ ░██    ░██ ░██           ░██████   ░██     ░██     ░██


"""
Save an Alakazam object to "filename".json. Converted to canonical json form.
Overwrite the specified file."""
function save(obj, filename::String, save_headers::Bool=false)
    write(filename*".json", to_JSON(obj, save_headers))
end

to_JSON(x) = JSON.json(Dict("expression"=>to_dict(x)))

to_JSON(x, headers::Bool) = headers ? JSON.json(Dict("header"=>get_header(), "expression"=>to_dict(x))) : to_JSON(x)

"""
Get's global information to be serialised.
Includes coordinates, IndexSets, sort_order, and anything else defined by the user.
"""
get_header() = to_dict(global_values)
get_header_json() = JSON.json(Dict("header"=>get_header()))
save_header(filename::String) = write(filename*".json", get_header_json)


##########################################

#global vars
to_dict(x::Union{Number,AbstractString,Symbol,Bool,Char,Nothing}) =
    Dict("type"=>string(typeof(x)), "value"=>x)

to_dict(x::Union{AbstractArray,Tuple}) =
    Dict("type"=>string(typeof(x)), "value"=>[to_dict(v) for v in x])

function to_dict(x::Dict)
    Dict("type"=>string(typeof(x)), "value"=>Dict(k=>to_dict(v) for (k, v) in x))
end

# only reached for actual structs with real fields
function to_dict(x)
    Dict("type"=>string(typeof(x)),
        (string(f)=>to_dict(getfield(x, f)) for f in fieldnames(typeof(x)))...)
end


#basic elements
to_dict(x::IndexSuperType) = Dict(
    "type" => string(typeof(x)),
    "name" => x.name,
    "index_set" => x.index_set.name,#to_dict(x.index_set)
)

"""
Skips storing the array of contained indices, otherwise there would be a loop.
"""
to_dict(x::IndexSet) = Dict(
    "type" => string(typeof(x)),
    "indices" => [string(i) for i in x.indices],
    (string(f) => getfield(x, f) for f in fieldnames(typeof(x)) if f != :indices)...
)

to_dict(x::Pair{I,IndexPosition}) where I<:IndexSuperType = Dict("index"=>to_dict(x[1]), "position"=>x[2].position)

to_dict(x::Coordinate) = Dict(
    "type" => string(typeof(x)),
    "name" => x.name,
    "index_set" => [set.name for set in x.index_set]
)

to_dict(x::YoungTableau) = return Dict(
    "type" => string(typeof(x)),
    "row_lengths" => x.row_lengths,
    "index_assignment" => x.index_assignment
)


###################
# basic tensordata
#######################
treat_field(x, f) = f
treat_field(x, f::Set{YoungTableau}) = [to_dict(t) for t in tableaux(x)]
treat_field(x, f::Vector{Pair{I,IndexPosition}}) where I<:IndexSuperType = [to_dict(t) for t in indices(x)]
treat_field(x, f::TensorExpression) = to_dict(f)
treat_field(x, f::TensorDataSuperType) = x isa OperatorSuperType ? to_dict_op(f) : to_dict(f)
# treat_field(x, f::Set{Coordinate}) = [to_dict(c) for c in f]
treat_field(x, f::Vector{Coordinate}) = [to_dict(c) for c in f]

function to_dict(x::TensorDataSuperType)
    out = Dict(
        "type" => string(typeof(x)),
        (string(f) => treat_field(x, getfield(x, f)) for f in fieldnames(typeof(x)))...)
    return out
end

"""
Don't save overall indices etc for ops, only individual parts from which this is computed from
"""
function to_dict_op(x::TensorDataSuperType)
    out = Dict(
        "type" => string(typeof(x)),
        (string(f) => treat_field(x, getfield(x, f)) for f in fieldnames(typeof(x)) if !(string(f) in ["indices", "name", "function_of", "weights", "tableaux", "labels"]))...)
    return out
end

#### op data
function to_dict(x::OperatorData)
    out = Dict(
        "type" => string(typeof(x)),
        (string(f) => treat_field(x, getfield(x, f)) for f in fieldnames(typeof(x)))...)
    return out
end

function to_dict(x::OperatorSuperType)
    out = Dict(
        "type" => string(typeof(x)),
        (string(f) => treat_field(x, getfield(x, f)) for f in fieldnames(typeof(x)))...)
    return out
end



########################################
#basic tensors and TEXPs
function to_dict(x::TensorSuperType)
    out = Dict(
        "type" => string(typeof(x)), "data"=>to_dict(data(x)))
    return out
end

function to_dict(x::TensorTerm)
    return Dict("type"=>string(typeof(x)), "terms"=>[to_dict(term) for term in terms(x)])
end

function treat_coef(x::SymbolicUtils.BasicSymbolic)
    Dict("type" => "SymbolicUtils.BasicSymbolic", "value" => _tree_to_dict(x))
end
function treat_coef(x::Number)
    Dict("type"=>string(typeof(x)), "value"=>x)
end

function to_dict(x::TensorExpression)
    return Dict("type"=>string(typeof(x)), "summands"=>[Dict("coef"=>treat_coef(coef(t)), "summand"=>to_dict(term(t))) for t in summands(x)])
end

function to_dict(x::Dict{DerivativeSuperType,Dict{TensorSuperType,TensorExpression}})
    Dict(
        "type" => string(typeof(x)),
        "value" => [
            Dict(
                "derivative" => to_dict(deriv),
                "actions" => [
                    Dict(
                        "tensor" => to_dict(tensor),
                        "expression" => to_dict(expr)
                    )
                    for (tensor, expr) in actions
                ]
            )
            for (deriv, actions) in x
        ]
    )
end

# ░██████░███     ░███ ░█████████    ░██████   ░█████████  ░██████████
#   ░██  ░████   ░████ ░██     ░██  ░██   ░██  ░██     ░██     ░██
#   ░██  ░██░██ ░██░██ ░██     ░██ ░██     ░██ ░██     ░██     ░██
#   ░██  ░██ ░████ ░██ ░█████████  ░██     ░██ ░█████████      ░██
#   ░██  ░██  ░██  ░██ ░██         ░██     ░██ ░██   ░██       ░██
#   ░██  ░██       ░██ ░██          ░██   ░██  ░██    ░██      ░██
# ░██████░██       ░██ ░██           ░██████   ░██     ░██     ░██

# ─────────────────────────────────────────────────────────────
#  LOAD
#
#  `into` selects the module that reconstructed indices are bound
#  into (default `Main`, matching the previous behaviour).
#
#  It is a KEYWORD only on the public entry points (`load`,
#  `load_header`, `from_JSON`) for a clean API, and a trailing
#  POSITIONAL argument through the entire internal recursion.
#  Keyword args on the mutually-recursive reconstruct/from_dict
#  functions defeat inference and inflate compile time badly, so
#  the internals must stay positional.
# ─────────────────────────────────────────────────────────────

"""
Loads Alakazam objects from "filename".json.

`into` selects the module that reconstructed indices are bound into
(default `Main`, matching previous behaviour). From inside a module use
`into=@__MODULE__`, or the `@load` convenience macro.
"""
function load(filename::String, skip_header::Bool=false; into::Module=Main)
    json = read(filename*".json", String)
    dict = JSON.parse(json)
    if !skip_header && haskey(dict, "header")
        import_header(dict["header"], into)
        # println("Successfully imported header.")
    end
    return from_dict(dict["expression"], into)
end

"""Convenience macro: binds reconstructed indices into the calling module."""
macro load(args...)
    :(load($(map(esc, args)...); into=($__module__)))
end

function load_header(filename::String, skip_header::Bool=false; into::Module=Main)
    json = read(filename*".json", String)
    dict = JSON.parse(json)
    if !skip_header && haskey(dict, "header")
        import_header(dict["header"], into)
        # println("Successfully imported header.")
    end
end

function from_JSON(s; into::Module=Main)
    d = JSON.parse(s)
    return from_dict(d["expression"], into)
end

"""Extract data from the header. `into` is positional (internal)."""
function import_header(s, into::Module=Main)
    data = s["value"]
    if "index_sets" in keys(data) #always need to import index_sets first
        reconstruct_header(data["index_sets"], Val(Symbol("index_sets")), into)
    end
    for key in keys(data)
        key == "index_sets" ? continue : reconstruct_header(data[key], Val(Symbol(key)), into)
    end
end

#needed function for each header info
function reconstruct_header(s, ::Val{:coordinates}, into::Module=Main)
    data = s["value"]
    for key in keys(data)
        if !(key in keys(COORD_REGISTRY))
            from_dict(data[key], into)
        end
    end
end
function reconstruct_header(s, ::Val{:index_sets}, into::Module=Main)
    data = s["value"]
    for key in keys(data)
        if !(key in keys(SET_REGISTRY))
            from_dict(data[key], into)
        end
    end
end
function reconstruct_header(s, ::Val{:covariant_expansions}, into::Module=Main)
    data = s["value"]
    for key in keys(data)
        if !(key in keys(DERIVATIVE_EXPANSIONS))
            expx = data[key]
            cander = expx["canonical_deriv"]
            der = from_dict(cander, into)
            expan = expx["expansion"]
            expanded = from_dict(expan, into)
            register_expansion!(der, expanded)
        end
    end
end
function reconstruct_header(s, ::Val{:commutators}, into::Module=Main)
    data = s["value"]
    for key in keys(data)
        keypair = Meta.parse(key) |> eval
        if !(keypair in keys(commutators))
            comdat = data[key]
            reconcom = from_dict(comdat, into)
            add_commutator!(reconcom)
        end
    end
end
function reconstruct_header(s, ::Val{:sort_order}, into::Module=Main)
    data = s["value"]
    nd = Dict{String,Int}()
    for (k, v) in data
        nd[k] = v["value"]
    end
    custom_sort_order!(nd)
end
function reconstruct_header(s, ::Val{:derivative_actions}, into::Module=Main)
    for deriv_entry in s["value"]
        D = from_dict(deriv_entry["derivative"], into)

        for action in deriv_entry["actions"]
            exp = from_dict(action["tensor"], into)
            res = from_dict(action["expression"], into)

            add_derivative_action(D, exp, res)
        end
    end
end
function reconstruct_header(s, ::Val{:derivative_commutators}, into::Module=Main)
    data = s["value"]
    for key in keys(data)
        keypair = Meta.parse(key) |> eval
        if !(keypair in keys(commutators))
            comdat = data[key]
            reconcom = from_dict(comdat, into)
            add_commutator!(reconcom)
        end
    end
end



#sorting into dicts, vecs etc to do correct recursion
from_json_value(x, into::Module=Main) = x
# function from_json_value(x::Dict{String,<:Any}, into::Module=Main)
#     if haskey(x, "type")
#         return from_dict(x, into)
#     else
#         return Dict(k => from_json_value(v, into) for (k, v) in x)
#     end
# end
# from_json_value(x::Vector, into::Module=Main) = map(v -> from_json_value(v, into), x)

# function from_dict(d::Dict{String,<:Any}, into::Module=Main)
#     T = _resolve_type(d["type"])  # acceptable if trusted input
#     return reconstruct(T, d, into)
# end

function from_dict(d::AbstractDict, into::Module=Main)
    T = _resolve_type(d["type"])
    return reconstruct(T, d, into)
end

function from_json_value(x::AbstractDict, into::Module=Main)
    if haskey(x, "type")
        return from_dict(x, into)
    else
        return Dict(k => from_json_value(v, into) for (k, v) in x)
    end
end


from_field(::Type{T}, f, v, into::Module=Main) where T =
    from_json_value(v, into)

#special override for :indices, which is stored as a vector of dict pos/index pairs
# could be useful in other places or in user scope for custom types, hence the hackiness
# to use dispatch
function from_field(::Type{T}, ::Val{:indices}, v, into::Module=Main) where T
    return [from_json_value(x["index"], into)=>IndexPosition(x["position"]) for x in v]
end
function from_field(::Type{T}, ::Val{:metric_type}, v, into::Module=Main) where T
    return eval(Meta.parse(v))
end
function from_field(::Type{T}, ::Val{:default_position}, v, into::Module=Main) where T
    return (IndexPosition(v["position"]))
end
function from_field(::Type{T}, ::Val{:index_set}, v, into::Module=Main) where T
    return SET_REGISTRY[v]
end

function from_field(::Type{T}, ::Val{:function_of}, v, into::Module=Main) where T
    return [
        reconstruct(Coordinate, x, into)
        for x in v
    ]
end
# function from_field(::Type{T}, ::Val{:signature}, v, into::Module=Main) where T
#     return Int[x for x in v]
# end
function from_json_value(x::AbstractVector, into::Module=Main)
    out = map(v -> from_json_value(v, into), x)
    return isconcretetype(eltype(out)) ? out : identity.(out)   # narrow if possible
end
#generic reconstruction
function reconstruct(T, d, into::Module=Main)
    args = map(fieldnames(T)) do f
        from_field(T, Val(f), d[string(f)], into)
    end
    return T(args...)
end


function reconstruct(::Type{IndexSet}, d, into::Module=Main)
    args = [
        from_field(IndexSet, Val(f), d[string(f)], into)
        for f in fieldnames(IndexSet) if f != :indices
    ]
    iset = IndexSet(args...)
    for name in d["indices"]
        register_indices(name, iset, into)
    end

    return iset
end

# THE ONLY function that evals; targets the chosen module.
function register_indices(name::String, iset, into::Module=Main)
    sym = Symbol(name)
    idx = Index(name, iset)

    Core.eval(into, :($sym = $idx))

    return idx
end

#special cases
function reconstruct(::Type{YoungTableau}, d, into::Module=Main)
    return YoungTableau(
        from_json_value(d["row_lengths"], into),
        from_json_value(d["index_assignment"], into)
    )
end

"""Since indices are globally shared, special logic."""
function reconstruct(T::Type{I}, d, into::Module=Main) where I<:IndexSuperType
    sn = d["index_set"]
    if sn in keys(SET_REGISTRY)
        set = SET_REGISTRY[sn]
        iN = d["name"]
        idx = findfirst(x -> name(x) == iN, set.indices)
        if idx !== nothing
            return set.indices[idx]
        else
            args = [
                from_field(I, Val(f), d[string(f)], into)
                for f in fieldnames(I) if f != :grading
            ]
            return Index(args...)
        end
    else
        error("Missing index set "*sn*". Try importing from header.")
    end
end

"""Fix logic for already created inds"""
function reconstruct(::Type{c}, d, into::Module=Main) where c<:Coordinate
    name = d["name"]

    if haskey(COORD_REGISTRY, name)
        return COORD_REGISTRY[name]
    end
    sets = [
        SET_REGISTRY[s]
        for s in d["index_set"]
    ]

    return Coordinate(name, sets)
end
"""Reconstruct Operators"""
reconstruct(::Type{OperatorData}, d, into::Module=Main) =
    OperatorData(from_json_value(d["operator_data"], into), from_json_value(d["operand"], into))


function reconstruct(::Type{N}, d, into::Module=Main) where N<:Number
    v = d["value"]

    if N <: Complex
        return N(v["re"], v["im"])
    end

    return convert(N, v)
end

function reconstruct(::Type{N}, d, into::Module=Main) where N<:SymbolicUtils.BasicSymbolic
    return _tree_from_dict(d["value"])
end

reconstruct(::Type{TensorTerm}, d, into::Module=Main) =
    TensorTerm(TensorSuperType[from_json_value(t, into) for t in d["terms"]])

function reconstruct(::Type{TensorExpression}, d, into::Module=Main)
    terms = Pair{SymorNum,TensorTerm}[]

    for s in d["summands"]
        c = from_json_value(s["coef"], into)
        t = from_json_value(s["summand"], into)

        @assert c isa SymorNum "Bad coefficient $(typeof(c)), $(c)"

        push!(terms, c => t)
    end

    return TensorExpression(terms)
end