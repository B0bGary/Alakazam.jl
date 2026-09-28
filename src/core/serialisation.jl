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
    elseif ex isa Number || ex isa Bool || ex isa Char
        # A VALUE type parameter, not a type: `SpinAlgebra{4,1}` carries its signature this way,
        # as does `Array{Int,2}`. Parsing yields the literal itself, which is already the parameter.
        ex
    elseif ex isa QuoteNode
        ex.value                      # a Symbol parameter, as in `Val{:x}`
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
# `get_header_json()`, not `get_header_json`: without the call this wrote the FUNCTION OBJECT and
# raised `MethodError: no method matching write(::IOStream, ::typeof(get_header_json))`, so this
# exported entry point could never have worked.
save_header(filename::String) = write(filename*".json", get_header_json())


##########################################

#global vars
to_dict(x::Union{Number,AbstractString,Symbol,Bool,Char,Nothing}) =
    Dict("type"=>string(typeof(x)), "value"=>x)

to_dict(x::Union{AbstractArray,Tuple}) =
    Dict("type"=>string(typeof(x)), "value"=>[to_dict(v) for v in x])

function to_dict(x::Dict)
    Dict("type"=>string(typeof(x)), "value"=>Dict(k=>to_dict(v) for (k, v) in x))
end

"""
Serialise one field of a struct.

Identical to `to_dict` except for REFERENTIAL values. An `IndexSet` or `Coordinate` is a shared,
registry-owned object: nesting a copy of it would both bloat the file and, worse, revive as a
distinct object, so any comparison against the registered one fails. Storing a reference keeps
identity intact.

This is applied by every struct path rather than by hand-written cases per type, so a user who
defines a type with a new field -- including one holding an index set -- is serialised correctly
without teaching the serialiser anything.
"""
_ser_field(v::IndexSet)   = Dict("__ref__" => "IndexSet",   "name" => v.name)
_ser_field(v::Coordinate) = Dict("__ref__" => "Coordinate", "name" => v.name)
# Anything JSON represents natively is stored as itself. Wrapping a String as
# {"type":"String","value":...} sends revival through `reconstruct(String, ...)`, which has no
# method -- the wrapper is only useful for values whose type cannot be inferred from the JSON.
_ser_field(v::Union{Number,AbstractString,Symbol,Bool,Char,Nothing}) = v
_ser_field(v::Union{AbstractArray,Tuple}) = [_ser_field(e) for e in v]
# A Type is stored as its name. Descending into one instead reaches Julia's own type internals
# (`Core.TypeName` and friends), which are self-referential -- that recursion is unbounded and
# overflows the stack. `from_field(::Val{:metric_type})` reads this back with `eval(Meta.parse(·))`.
_ser_field(v::Type) = string(v)
_ser_field(v) = to_dict(v)

# only reached for actual structs with real fields
function to_dict(x)
    Dict("type"=>string(typeof(x)),
        (string(f)=>_ser_field(getfield(x, f)) for f in fieldnames(typeof(x)))...)
end


#basic elements
to_dict(x::IndexSuperType) = Dict(
    "type" => string(typeof(x)),
    "name" => x.name,
    "index_set" => x.index_set.name,#to_dict(x.index_set)
)

"""
Skips storing the array of contained indices, otherwise there would be a loop.

Two conjugate index sets likewise point at each other, so nesting the partner would recurse without
end. The conjugate is therefore stored by NAME and relinked once every set exists -- the same
treatment `to_dict(::Coordinate)` already gives the sets a coordinate belongs to.
"""
to_dict(x::IndexSet) = Dict(
    "type" => string(typeof(x)),
    "indices" => [string(i) for i in x.indices],
    "conjugate" => (x.conjugate === nothing ? nothing : x.conjugate.name),
    # via `to_dict`, not a raw `getfield`: an algebra is a struct, and storing it raw flattens it
    # to its bare fields with the TYPE discarded, so it came back as an anonymous Dict that the
    # IndexSet constructor rejected outright.
    "algebra" => (x.algebra === nothing ? nothing : to_dict(x.algebra)),
    # by NAME, like the conjugate: a coordinate holds its sets, so nesting it would recurse. Relinked
    # by `import_header` once the coordinates exist.
    "metric_function" => [c.name for c in x.metric_function],
    # the remaining fields stay RAW: `IndexSet` has `from_field` handlers matched to those exact
    # forms (metric_type as a type name, default_position, grading), and wrapping them here would
    # hand those handlers a shape they do not expect.
    (string(f) => getfield(x, f)
     for f in fieldnames(typeof(x)) if f ∉ (:indices, :conjugate, :algebra, :metric_function))...
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
treat_field(x, f::Vector{Pair{I,IndexPosition}}) where I<:IndexSuperType = [to_dict(t) for t in f]
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
A tensor's data and its indices are separate objects now
"""
function _flat_tensor_dict(x)
    return Dict(
        "type" => "Alakazam.TensorData",
        "name" => name(x),
        "indices" => [to_dict(i) for i in indices(x)],
        "function_of" => [to_dict(c) for c in function_of(x)],
        "weights" => weights(x),
        "tableaux" => [to_dict(t) for t in tableaux(x)],
        "labels" => labels(x),
        # by type name, not by instance: NumberSets are singletons, so the name is the whole value
        "numberfield" => string(nameof(typeof(numberfield(x)))),
    )
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
    d = Dict("type" => string(typeof(x)), "data" => _flat_tensor_dict(x))
    extra = _extra_fields(x)
    isempty(extra) || (d["extra"] = extra)
    return d
end

"""
Fields a tensor type declares BEYOND the `data`/`indices` pair.

The flat layout above carries only those two, because almost every tensor has exactly them and
flattening keeps the common case cheap. A type with more -- `IdMatrix` records which space it is
the identity on -- would otherwise have that field silently dropped, and since `IdMatrix` defines
`==` and `hash` on it, two identities over different spaces would compare equal after a round trip.

Collected generically rather than by writing a bespoke `to_dict` per type, so a new type with an
extra field is handled without anyone remembering to teach the serialiser about it. The revival
side reads them back through `from_field`, so a field can still opt into special treatment by
defining a method there.
"""
function _extra_fields(x)
    out = Dict{String,Any}()
    for f in fieldnames(typeof(x))
        (f === :data || f === :indices) && continue
        out[string(f)] = _extra_value(getfield(x, f))
    end
    return out
end

# Referential values are stored by registry key, never nested: an IndexSet reached this way must
# come back as THE registered set, not a copy, or identity comparisons against it fail.
_extra_value(v::Nothing) = nothing
_extra_value(v) = _ser_field(v)

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
    # last: a set's metric_function names coordinates, which exist only once the loop above is done
    "index_sets" in keys(data) && _relink_metric_functions(data["index_sets"])
end

"""Resolve the coordinate names each stored index set's metric depends on. Files written before
metric_function existed have no entry and stay flat."""
function _relink_metric_functions(s)
    data = s["value"]
    for key in keys(data)
        names = get(data[key], "metric_function", nothing)
        (names === nothing || isempty(names) || !haskey(SET_REGISTRY, key)) && continue
        all(n -> haskey(COORD_REGISTRY, n), names) ||
            error("index set \"$key\" has a metric depending on $(join(names, ", ")), but not all of these coordinates were found.")
        metric_function!(SET_REGISTRY[key], [COORD_REGISTRY[n] for n in names])
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
    # second pass: every set now exists, so the conjugate names stored above can be resolved into
    # the mutual references they stand for. Done here rather than in the constructor because the
    # partner is generally not built yet at that point.
    for key in keys(data)
        entry = data[key]
        haskey(entry, "conjugate") || continue
        partner = from_json_value(entry["conjugate"], into)
        partner === nothing && continue
        (haskey(SET_REGISTRY, key) && haskey(SET_REGISTRY, partner)) || continue
        a, b = SET_REGISTRY[key], SET_REGISTRY[partner]
        a.conjugate === b && continue          # already linked from the partner's own entry
        conjugates!(a, b)
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
        # `derivative_algebra`, not `commutators`: a CovariantCommutator is stored in the former,
        # so guarding on the latter tested a dict this entry could never be in and the guard never
        # fired. Harmless while re-adding overwrites with an identical value, but it meant the
        # "already present" check did nothing at all.
        if !(keypair in keys(derivative_algebra))
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
    # a reference resolves to the registered object, never to a rebuilt copy
    if haskey(x, "__ref__")
        r = x["__ref__"]
        r == "IndexSet"   && return SET_REGISTRY[x["name"]]
        r == "Coordinate" && return COORD_REGISTRY[x["name"]]
        error("Unknown reference kind \"$r\" in serialised data")
    end
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
function from_field(::Type{T}, ::Val{:operator_indices}, v, into::Module=Main) where T
    return Pair{IndexTypes,IndexPosition}[from_json_value(x["index"], into)=>IndexPosition(x["position"]) for x in v]
end

function from_field(::Type{T}, ::Val{:metric_type}, v, into::Module=Main) where T
    return eval(Meta.parse(v))
end
"""
Revive a grading as an `IndexGrading` rather than as the bare `Int` the field stores.

`reconstruct(::Type{IndexSet}, ...)` feeds the field values positionally into `IndexSet(name,
args...)`, but that constructor dispatches on the TYPE of each argument, not its position. A bare
`Int` there matches the `range` branch, so a stored grading was being read as a dimension and the
grading itself silently stayed at its default -- a fermionic set came back bosonic, which changes
every Koszul sign computed from it.
"""
function from_field(::Type{IndexSet}, ::Val{:grading}, v, into::Module=Main)
    g = from_json_value(v, into)
    g isa IndexGrading && return g
    return IndexGrading(g)
end

"""
Deferred: the partner may not exist yet when this set is built, so the link is made in a second
pass by `reconstruct_header(::Val{:index_sets})` once all of them are present.
"""
from_field(::Type{IndexSet}, ::Val{:conjugate}, v, into::Module=Main) = nothing
"""Deferred for the same reason: coordinates are built after the sets, see `_relink_metric_functions`."""
from_field(::Type{IndexSet}, ::Val{:metric_function}, v, into::Module=Main) = nothing

"""Rebuild the algebra from its recorded type rather than the anonymous field dump."""
function from_field(::Type{IndexSet}, ::Val{:algebra}, v, into::Module=Main)
    v === nothing && return nothing
    return from_dict(v, into)
end

function from_field(::Type{T}, ::Val{:default_position}, v, into::Module=Main) where T
    return (IndexPosition(v["position"]))
end
function from_field(::Type{T}, ::Val{:index_set}, v, into::Module=Main) where T
    v === nothing && return nothing            # an optional space, e.g. an unattached IdMatrix
    v isa AbstractString && return SET_REGISTRY[v]   # how an Index stores its set
    return from_json_value(v, into)                  # a reference dict from `_ser_field`
end

"""A `NumberSet` round-trips as its type name; rebuild the singleton from it."""
function from_field(::Type{T}, ::Val{:numberfield}, v, into::Module=Main) where T
    v isa NumberSet && return v
    nm = v isa AbstractString ? String(v) : string(v)
    nm = last(split(nm, '.'))                    
    return getfield(Alakazam, Symbol(nm))()
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
"""
Labels round-trip awkwardly: JSON renders a `Vector{Pair{String,IndexPosition}}` as a map, so an
empty label list comes back as `{}` rather than `[]`. Accept either shape.
"""
_as_labels(v::AbstractDict) = LabelPair[String(k) => _as_position(p) for (k, p) in v]
_as_labels(v::AbstractVector) = LabelPair[_as_label_pair(x) for x in v]
_as_label_pair(x::Pair) = LabelPair(String(x[1]), _as_position(x[2]))
_as_label_pair(x::AbstractDict) = LabelPair(String(x["first"]), _as_position(x["second"]))
_as_position(p::IndexPosition) = p
# An `IndexPosition` written raw serialises as an OBJECT, `{"position": true}`, not as the bare
# boolean -- so a non-empty label list came back as a Dict here and `Bool(::Dict)` failed. Empty
# lists never hit this, which is why it went unnoticed.
_as_position(p::AbstractDict) = IndexPosition(Bool(p["position"]))
_as_position(p) = IndexPosition(Bool(p))

"""Rebuild a `TensorData` alone - used for an operator's own data block, which has no indices."""
reconstruct(::Type{TensorData}, d, into::Module=Main) = _tensor_parts(d, into)[1]

"""Rebuild a `TensorData` and its index vector from the flat serialised layout."""
function _tensor_parts(d, into::Module=Main)
    nm   = from_field(TensorData, Val(:name), d["name"], into)
    # an operator's TensorData block has no "indices" key: its indices live on the OperatorData
    inds = haskey(d, "indices") ? from_field(TensorData, Val(:indices), d["indices"], into) : IndexPair[]
    fo   = from_field(TensorData, Val(:function_of), d["function_of"], into)
    w    = from_field(TensorData, Val(:weights), d["weights"], into)
    tabs = from_field(TensorData, Val(:tableaux), d["tableaux"], into)
    labs = from_field(TensorData, Val(:labels), d["labels"], into)
    # absent in files written before the field existed, which are all real by construction
    nf   = haskey(d, "numberfield") ?
           from_field(TensorData, Val(:numberfield), d["numberfield"], into) : realfield()
    # route through the argument parser so each piece lands in its declared field type
    return build_tensor(nm,
        Vector{IndexPair}(inds),
        Vector{Coordinate}(fo),
        Dict{String,Number}(w),
        Vector{YoungTableau}(tabs),
        _as_labels(labs),
        nf)
end

"""Tensors serialise as one "data" block; split it back into data and indices."""
function reconstruct(T::Type{<:TensorSuperType}, d, into::Module=Main)
    parts = _tensor_parts(d["data"], into)
    haskey(d, "extra") || return T(parts...)
    # extras follow data/indices in declaration order, which is the order the constructors take
    extras = Any[from_field(T, Val(f), d["extra"][string(f)], into)
                 for f in fieldnames(T) if f !== :data && f !== :indices]
    return T(parts..., extras...)
end

"""
Operators are also `<:TensorSuperType` but serialise field-by-field (see `to_dict_op`), so they
keep the generic reconstruction rather than the tensor one above.
"""
function reconstruct(T::Type{<:OperatorSuperType}, d, into::Module=Main)
    args = map(fieldnames(T)) do f
        from_field(T, Val(f), d[string(f)], into)
    end
    return T(args...)
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
        #a field missing from an older file is left at its default: the constructor reads by type, not position
        for f in fieldnames(IndexSet) if f != :indices && haskey(d, string(f))
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
            # Build via the normal constructor
            return I(iN, set)
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
"""Reconstruct Operators.

The operator's own indices used to live inside its `TensorData`; they are now a separate field, so
they are read back separately and handed to the constructor alongside it.
"""
function reconstruct(::Type{OperatorData}, d, into::Module=Main)
    op_data = from_json_value(d["operator_data"], into)
    op_inds = from_field(OperatorData, Val(:operator_indices), d["operator_indices"], into)
    return OperatorData(op_data, op_inds, from_json_value(d["operand"], into))
end


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