const LabelPair = Pair{String,IndexPosition}
const _T_INDEX_VEC = Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
const _T_INDEX_PAIR = Pair{I,IndexPosition} where {I<:IndexSuperType}
const _EMPTY_FUNCTION_OF =  Coordinate[]
const _EMPTY_TABLEAUX = Set{YoungTableau}()
const _EMPTY_WEIGHTS = Dict{String,Number}()
const _EMPTY_LABELS = LabelPair[]
const _EMPTY_INDICES = IndexPair[]
@inline coordfree() = _EMPTY_FUNCTION_OF
@inline notableaux()::Set{YoungTableau} = _EMPTY_TABLEAUX
@inline noweights()::Dict{String,Number} = _EMPTY_WEIGHTS
@inline nolabels()::Vector{LabelPair} = _EMPTY_LABELS
# Shared default, so every real tensor points at one instance rather than allocating its own.
# Real rather than complex because conjugation then stops at the leaf for the common case; a tensor
# whose conjugate is a distinct object has to say so.
const REAL_FIELD = ℝ()
@inline realfield()::NumberSet = REAL_FIELD
"""
Content hash of a `TensorData`, covering every field.

Immutability makes `objectid` expensive here, it recomputes over the fields on every call, ~75 ns against ~1 ns for comparing this.

`tableaux` and `weights` are unordered containers, so they fold in through XOR; `function_of` is
kept in canonical order by `_insert_coord!` and `labels` is compared as a vector by `==`, so both
are hashed in sequence.
"""
function _tensordata_uid(name::String, function_of, weights, tableaux, labels,
                         numberfield)::UInt64
    h = hash(nameof(typeof(numberfield)), hash(name))
    for c in function_of
        h = hash(c, h)
    end
    for l in labels
        h = hash(l[2].position, hash(l[1], h))
    end
    x = UInt(0)
    for t in tableaux
        x ⊻= hash(t)
    end
    for (k, v) in weights
        x ⊻= hash(v, hash(k))
    end
    return hash(x, h)
end

"""
Mutual tensor data that all tensors should have in common. This makes it easier to add new features to tensors.
Includes everything about a tensor except its indices.
"""
struct TensorData <: TensorDataSuperType
    name::String
    function_of::Vector{Coordinate} #set if ok
    weights::Dict{String,Number}
    tableaux::Set{YoungTableau}
    labels::Vector{LabelPair}
    numberfield::NumberSet
    uid::UInt64
    function TensorData(name::String, function_of, weights, tableaux, labels,
                        numberfield::NumberSet = REAL_FIELD)
        new(name, function_of, weights, tableaux, labels, numberfield,
            _tensordata_uid(name, function_of, weights, tableaux, labels, numberfield))
    end
end

"""
    build_tensor(name, args...) -> (TensorData, Vector{IndexPair})

Parse the tensor constructor arguments into the shared data and the index vector.
"""
function build_tensor(name::String, args...)
        indices = Vector{IndexPair}()
        function_of = coordfree();
        has_coords = false
        tableaux = notableaux();
        has_tabs = false
        labels = nolabels();
        has_labels = false
        weights = noweights()
        numberfield = realfield()

        for arg in args
            if arg isa _T_INDEX_VEC
                append!(indices, arg)
            elseif arg isa _T_INDEX_PAIR
                push!(indices, arg)
            elseif arg isa Coordinate
                has_coords || (function_of=Coordinate[]; has_coords=true)
                _insert_coord!(function_of, arg)
            elseif arg isa AbstractVector{<:Coordinate}
                isempty(arg) && continue
                has_coords || (function_of=Coordinate[]; has_coords=true)
                for c in arg
                    _insert_coord!(function_of, c)
                end
            elseif arg isa Dict{String,Number}
                weights = arg
            elseif arg isa YoungTableau
                has_tabs || (tableaux=Set{YoungTableau}(); has_tabs=true)
                push!(tableaux, arg)
            elseif arg isa Vector{YoungTableau}
                isempty(arg) && continue
                has_tabs || (tableaux=Set{YoungTableau}(); has_tabs=true)
                union!(tableaux, arg)
            elseif arg isa String
                has_labels || (labels=LabelPair[]; has_labels=true)
                push!(labels, arg => lower)
            elseif arg isa Pair{String,IndexPosition}
                has_labels || (labels=LabelPair[]; has_labels=true)
                push!(labels, arg)
            elseif arg isa Vector{Pair{String,IndexPosition}}
                isempty(arg) && continue
                has_labels || (labels=LabelPair[]; has_labels=true)
                append!(labels, arg)
            elseif arg isa NumberSet
                numberfield = arg
            elseif arg isa Type{<:NumberSet}
                # accept the bare alias too, so `Tensor("psi", ℂ, ...)` works alongside `ℂ()`
                numberfield = arg()
            end
        end

        tabs = has_tabs ? reduce_tableaux(tableaux) : tableaux

        if !has_coords
            return (TensorData(name, function_of, weights, tabs, labels, numberfield), indices)
        end
        w = Dict{String,Number}()
        for c in function_of
            ;
            w[c.name] = 1;
        end
        return (TensorData(name, function_of, merge(w, weights), tabs, labels, numberfield),
                indices)
end

"""Rebuild a `TensorData` with one field replaced. Used by the rare mutating operations."""
_with_name(p::TensorData, nm::String) =
    TensorData(nm, p.function_of, p.weights, p.tableaux, p.labels, p.numberfield)
_with_tableaux(p::TensorData, tabs::Set{YoungTableau}) =
    TensorData(p.name, p.function_of, p.weights, tabs, p.labels, p.numberfield)
"""Rebuild with a different number field. The name is left alone: bar-toggling belongs to
conjugation, which decides it from the field rather than the other way round."""
_with_numberfield(p::TensorData, f::NumberSet) =
    TensorData(p.name, p.function_of, p.weights, p.tableaux, p.labels, f)

# Defined here rather than in young.jl because that file loads before TensorData exists.
function _add_tableau(p::TensorData, tab::YoungTableau)
    tabs = tableaux(p)
    if tabs === notableaux()
        tabs = Set{YoungTableau}()
        tabs.dict[tab] = nothing
        return _with_tableaux(p, tabs)
    end
    push!(tabs, tab)
    return _with_tableaux(p, tabs)
end
add_tableau!(T::TensorSuperType, tab::YoungTableau) = (T.data = _add_tableau(T.data, tab); T)


function _insert_coord!(v::Vector{Coordinate}, c::Coordinate)
    for (i, x) in enumerate(v)
        x == c && return v                    # already present
        if isless(c.name, x.name)
            insert!(v, i, c)
            return v
        end
    end
    push!(v, c)
    return v
end


#####################
####### primitive tensor types
#####################

#simplest tensor constructu
"""Basic Tensor object."""
mutable struct Tensor <: TensorSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Tensor(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Tensor(name::String, args...)
        return new(build_tensor(name, args...)...)
    end
end

"""Basic Fermion. Fermionic no matter the assigned indices. Identical to a tensor otherwise."""
mutable struct Fermion <: TensorSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Fermion(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Fermion(name::String, args...)
        return new(build_tensor(name, args...)...)
    end
end
"""Basic Boson. Bosonic no matter the assigned indices. Identical to a tensor otherwise."""
mutable struct Boson <: TensorSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Boson(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Boson(name::String, args...)
        return new(build_tensor(name, args...)...)
    end
end

function SymmetricTensor(name::String, args...)
    t = Tensor(name, args...)
    if length(indices(t))>1
        add_tableau!(t, generate_symmetric_tableau(indices(t)))
        # data(t).tableaux = reduce_tableaux(tableaux(t))
    end
    return t
end

function AntisymmetricTensor(name::String, args...)
    t = Tensor(name, args...)
    if length(indices(t))>1
        add_tableau!(t, generate_antisymmetric_tableau(indices(t)))
    end
    return t
end

function RiemannTensor(args...)
    return RiemannTensor("R", args...)
end

function RiemannTensor(name::String, args...)
    t = Tensor(name, args...)
    if length(indices(t)) != 4 #perhaps at some point want to fix this incase people have some sort of multi-metric gravity etc
        throw(ArgumentError("Riemann tensor must have exactly 4 indices"))
    end
    add_tableau!(t, generate_riemann_tableau())
    return t
end


#######
# common tableaux generation functions
######
# These tableaux are generated many times (eg in the metric), so we cache and share them
const _SYM_TABLEAU_CACHE = Dict{Int,YoungTableau}()
const _ASYM_TABLEAU_CACHE = Dict{Int,YoungTableau}()

generate_symmetric_tableau(indices::Vector) = generate_symmetric_tableau(length(indices))
generate_antisymmetric_tableau(indices::Vector) = generate_antisymmetric_tableau(length(indices))
function generate_symmetric_tableau(n::Integer)
    k = Int(n)
    return get!(() -> YoungTableau([k], collect(1:k)), _SYM_TABLEAU_CACHE, k)
end
function generate_antisymmetric_tableau(n::Integer)
    k = Int(n)
    return get!(() -> YoungTableau(fill(1, k), collect(1:k)), _ASYM_TABLEAU_CACHE, k)
end
function generate_riemann_tableau()
    return YoungTableau([2, 2], [1, 3, 2, 4])
end

######
# compound tensor terms and expressions
# ADD proper validity checking...
######
mutable struct TensorTerm <: Tensorial
    terms::Vector{TensorSuperType}# where {T<:TensorSuperType}
    indices::Vector{IndexPair}
    # function TensorTerm(terms::Vector{T} where {T<:TensorSuperType})
    #     inds = IndexPair[]
    #     outs=Vector{TensorSuperType}()
    #     for t in terms
    #         if name(t)==""
    #             continue
    #         end
    #         for i in indices(t)
    #             push!(inds, i)

    #         end
    #         push!(outs, t)
    #     end
    #     return new(outs, inds)
    # end
    function TensorTerm(terms::Vector{T} where {T<:TensorSuperType})
        nt = 0;
        ni = 0
        for t in terms
            isempty(name(t)) && continue
            nt += 1
            ni += length(indices(t))
        end
        outs = Vector{TensorSuperType}(undef, nt)
        inds = Vector{IndexPair}(undef, ni)
        j = 1;
        k = 1
        for t in terms
            isempty(name(t)) && continue
            @inbounds outs[j] = t;
            j += 1
            k=_copy_inds!(inds, k, t)
        end
        return new(outs, inds)
    end
    function TensorTerm(term::TensorSuperType)
        return new([term], indices(term))
    end
    function TensorTerm(A::TensorSuperType, B::TensorSuperType)
        ia = indices(A);
        ib = indices(B)
        n = length(ia) + length(ib)
        inds = Vector{IndexPair}(undef, n)
        @inbounds for k in eachindex(ia)
            ;
            inds[k] = ia[k];
        end
        @inbounds for k in eachindex(ib)
            ;
            inds[length(ia)+k] = ib[k];
        end
        new(TensorSuperType[A, B], inds)

    end
    function TensorTerm(term::TensorTerm)
        return term
    end
end
"""Helper for the constructor to use specific typing"""
function _copy_inds!(inds, k, t)
    ti = indices(t)
    @inbounds for m in eachindex(ti)
        inds[k] = ti[m];
        k += 1
    end
    return k
end



const CoeffPair = Pair{SymorNum,TensorTerm}

const _RAW = Val(:raw)
mutable struct TensorExpression <: Tensorial
    summands::Vector{CoeffPair}
    TensorExpression(v::Vector{CoeffPair}) = merge_sum(new(v))
    TensorExpression(::Val{:raw}, v::Vector{CoeffPair}) = new(v)   # caller guarantees already merged, no zero coeffs
    TensorExpression() = new(CoeffPair[CoeffPair(0, TensorTerm(TensorSuperType[]))])
end

# ---------------------------------------------------------------------------
# single argument
TensorExpression(tex::TensorExpression) = tex
# coefficient is 1 no merge
TensorExpression(T::TensorTerm) = TensorExpression(_RAW, CoeffPair[CoeffPair(1, T)])
TensorExpression(T::TensorSuperType) = TensorExpression(_RAW, CoeffPair[CoeffPair(1, TensorTerm(T))])
# constant
TensorExpression(n::SymorNum) = isZero(n) ? zero() : TensorExpression(_RAW, CoeffPair[CoeffPair(n, TensorTerm(TensorSuperType[]))])
# ---------------------------------------------------------------------------
# two arguments so only term equality matters
TensorExpression(A::TensorTerm, B::TensorTerm) = _merge2(CoeffPair[CoeffPair(1, A), CoeffPair(1, B)])
function TensorExpression(A::TensorSuperType, B::TensorSuperType)
    if A == B
        return TensorExpression(_RAW, CoeffPair[CoeffPair(2, TensorTerm(A))])
    end
    return TensorExpression(_RAW, CoeffPair[CoeffPair(1, TensorTerm(A)), CoeffPair(1, TensorTerm(B))])
end
function TensorExpression(A::TensorTerm, B::TensorSuperType)
    if length(terms(A)) == 1 && terms(A)[1] == B
        return TensorExpression(_RAW, CoeffPair[CoeffPair(2, A)])
    end
    return TensorExpression(_RAW, CoeffPair[CoeffPair(1, A), CoeffPair(1, TensorTerm(B))])
end
# ---------------------------------------------------------------------------
#only the new term can merge
# TensorExpression(A::TensorTerm, E::TensorExpression) = _insert_merged(copy(E.summands), 1, A)
# TensorExpression(A::TensorSuperType, E::TensorExpression) = _insert_merged(copy(E.summands), 1, TensorTerm(A))
TensorExpression(E::TensorExpression, A::TensorSuperType) = _insert_merged(copy(E.summands), 1, TensorTerm(A))
TensorExpression(E::TensorExpression, A::TensorTerm) = _insert_merged(copy(E.summands), 1, A)
# both sides may contain matching terms so general merge
# TensorExpression(E1::TensorExpression, E2::TensorExpression) = TensorExpression(CoeffPair[E1.summands; E2.summands])
#simplify assuming each individual expression is already merged
function TensorExpression(E1::TensorExpression, E2::TensorExpression)
    v = copy(E1.summands)
    n1 = length(v)
    @inbounds for (c, t) in E2.summands
        isZero(c) && continue
        found = false
        for k in 1:length(v)
            if v[k][2] == t
                nc = v[k][1] + c
                if isZero(nc)
                    deleteat!(v, k)
                else
                    v[k] = CoeffPair(nc, t)
                end
                found = true;
                break
            end
        end
        found || push!(v, CoeffPair(c, t))
    end
    isempty(v) && return zero()
    return TensorExpression(_RAW, v)
end
# ---------------------------------------------------------------------------
TensorExpression(n::SymorNum, T::TensorSuperType) = isZero(n) ? zero() : TensorExpression(_RAW, CoeffPair[CoeffPair(n, TensorTerm(T))])
TensorExpression(p::Pair{<:SymorNum,TensorTerm}) = isZero(p[1]) ? zero() : TensorExpression(_RAW, CoeffPair[CoeffPair(p[1], p[2])])
TensorExpression(p::Pair{<:SymorNum,<:TensorSuperType}) = isZero(p[1]) ? zero() : TensorExpression(_RAW, CoeffPair[CoeffPair(p[1], TensorTerm(p[2]))])
TensorExpression(v::Vector{<:Pair{<:SymorNum,TensorTerm}}) = TensorExpression(CoeffPair[CoeffPair(p[1], p[2]) for p in v])
TensorExpression(v::Vector{<:Pair{<:SymorNum,<:TensorSuperType}}) = TensorExpression(CoeffPair[CoeffPair(p[1], TensorTerm(p[2])) for p in v])
# ---------------------------------------------------------------------------
# vectors of terms
TensorExpression(v::Vector{<:TensorTerm}) = TensorExpression(CoeffPair[CoeffPair(1, t) for t in v])
TensorExpression(v::Vector{<:TensorSuperType}) = TensorExpression(CoeffPair[CoeffPair(1, TensorTerm(t)) for t in v])

function _prepend_merged(c, t::TensorTerm, v::Vector{CoeffPair})
    @inbounds for i in eachindex(v)
        if v[i][2] == t
            nc = v[i][1] + c
            isZero(nc) && (w=copy(v); deleteat!(w, i); return isempty(w) ? zero() : TensorExpression(_RAW, w))
            w = copy(v);
            w[i] = CoeffPair(nc, t)
            return TensorExpression(_RAW, w)
        end
    end
    return TensorExpression(_RAW, CoeffPair[CoeffPair(c, t); v])
end

TensorExpression(A::TensorTerm, E::TensorExpression) = _prepend_merged(1, A, E.summands)
TensorExpression(A::TensorSuperType, E::TensorExpression) = _prepend_merged(1, TensorTerm(A), E.summands)

# https://docs.julialang.org/en/v1/manual/metaprogramming/
Base.@doc "Usage: `@Tensor A^μν_ρ^τ{x,y,z}` where `A` is a symbol, `μ,ν,ρ,τ` are indices, and x,y,z are coordinates."
macro Tensor(expr::Expr)
    functionof = Symbol[]

    symbol = string(first(string(expr)))  # Unicode-safe

    indices = []

    # Handle indices attached to the tensor name, e.g. A_bc
    if length(string(expr.args[2])) > 1 &&
       string(expr.args[2])[nextind(string(expr.args[2]), 1)] == '_'

        push!(
            indices,
            [v => lower for v in
                            split(string(expr.args[2])[nextind(string(expr.args[2]), 1, 2):end], "")]...
        )
    else
        @assert length(string(expr.args[2])) == 1 "Tensor symbol must be a single character, not $(expr.args[2])"
    end

    head = expr.args[3]

    while true

        # Peel off trailing {x,y,...}
        if head isa Expr && head.head == :curly
            append!(functionof, Symbol.(head.args[2:end]))
            head = head.args[1]
        end

        if head isa Expr && head.args[1] == :^
            indices = parse_indices(head, indices, false)
            head = head.args[3]
        else
            break
        end
    end

    # Last index group may also have {x,y}
    if head isa Expr && head.head == :curly
        append!(functionof, Symbol.(head.args[2:end]))
        head = head.args[1]
    end

    indices = parse_indices(head, indices, true)

    functions_expr = Expr(:vect, map(esc, functionof)...)

    return quote
        $(esc(Symbol(symbol))) = Tensor(
            $symbol,
            [$([:($(esc(Symbol(idx.first))) => $(esc(Symbol(idx.second))))
                for idx in indices]...)],
            [],
            $functions_expr
        )
    end
end

function parse_indices(head, indices, leaf)
    str_head = [v for (i, v) in enumerate(string(leaf ? head : head.args[2]))]
    # @show str_head
    if length(str_head) == 1
        indices = [indices..., str_head[1] => upper]
    else
        underscore = findfirst('_' .== str_head)
        # underscore = underscore === nothing ? length(string(head)) : underscore
        if underscore != nothing
            underscoreless_str = [str_head[1:(underscore-1)]..., str_head[(underscore+1):end]...]
            indices = [indices..., [v => i ≥ underscore ? lower : upper for (i, v) in enumerate(underscoreless_str)]...]
        else
            indices = [indices..., [v => upper for (i, v) in enumerate(str_head)]...]
        end
    end
    return indices
end

Base.@doc"Usage: `@Tensor A_bc` where `A` is a symbol and `b`/`c` are indices."
macro Tensor(expr::Symbol)
    length(string(expr)) == 1 && return Tensor(string(expr))
    index_count = length(string(expr)) - 2
    # split the symbol A_bc at the underscore
    symbol, indices = split(string(expr), "_")
    return quote
        $(esc(Symbol(symbol))) = Tensor(
            $(string(symbol)),
            [$([:($(esc(Symbol(index))) => lower) for index in indices]...)]
        )
    end
end

macro T_str(expr::String)
    esc(:(@Tensor $expr))
end
"""
Validate one entry of a `{...}` group in a tensor string.

`{...}` lists the coordinates a tensor is a function of. Indices are written bare after `^`/`_`,
with no braces unlike `replace_pattern`'s own string syntax, where braces group indices.
"""
function _coord_arg(x, tok::AbstractString)
    x isa Coordinate && return x
    if x === nothing
        error("In a tensor string, {...} lists the coordinates the tensor depends on, but no " *
              "variable named `$tok` is defined. If `$tok` was meant to be an index, write it " *
              "bare after ^ or _ (e.g. \"A_$tok\" or \"A^$tok\"): index groups are not braced.")
    end
    if x isa IndexSuperType
        error("In a tensor string, {...} lists the coordinates the tensor depends on, but `$tok` " *
              "is an index, not a Coordinate. Did you mean to write it as an index, i.e. " *
              "\"A_$tok\" (lower) or \"A^$tok\" (upper)? Index groups are not braced in " *
              "T\"...\"/DT\"...\" strings.")
    end
    error("In a tensor string, {...} lists the coordinates the tensor depends on, but `$tok` is a " *
          "$(typeof(x)), not a Coordinate.")
end

const index_suffixes = ["dot", "hat", "bar"]
function next_index(expr, i, comma_mode)

    if comma_mode
        # read until comma or structural character
        j = i
        while j <= lastindex(expr)
            c = expr[j]
            if c == ','
                break
            elseif c in ['^', '_', '[', ']', '(', ')', '{', '}', '|']
                break
            end
            j = nextind(expr, j)
        end

        idx = Symbol(expr[i:prevind(expr, j)])

        # skip comma
        if j <= lastindex(expr) && expr[j] == ','
            j = nextind(expr, j)
        end

        return idx, prevind(expr, j)
    end

    name = string(expr[i])
    pos = i

    for suffix in index_suffixes
        j = pos
        ok = true

        for _ in 1:length(suffix)
            if j == lastindex(expr)
                ok = false
                break
            end
            j = nextind(expr, j)
        end

        if ok && expr[nextind(expr, pos):j] == suffix
            name *= suffix
            pos = j
            break
        end
    end

    return Symbol(name), pos
end
"""
Build one index slot from a parsed token.

In order of precedence: a leading `?` is a wildcard, becoming a `DummyPatternIndex` resolved at
match time; a `'quoted'` name is a label even when a variable of that name is in scope; a defined
identifier is that index; anything else (digits, undefined names) is a label.
"""
function slot_expr(idx::Symbol, pos)
    name = string(idx)
    startswith(name, "?") &&
        return :($(Pair{IndexTypes,IndexPosition})($(DummyPatternIndex)($(name[nextind(name, 1):end])), $pos))
    startswith(name, "'") &&
        return :($(Pair{String,IndexPosition})($(name[nextind(name, 1):end]), $pos))
    Base.isidentifier(name) ||
        return :($(Pair{String,IndexPosition})($name, $pos))
    return :($(esc(:(@isdefined $idx))) ?
             ($(esc(idx)) => $pos) :
             $(Pair{String,IndexPosition})($name, $pos))
end

"""
Shared body for the `@Tensor` and `T"..."` parsers.

Patterns are built by exactly the same parser as ordinary tensors, and so inherit labels, tableaux,
symmetry syntax and coordinate dependence. A wildcard needs no separate type: `?x` in a slot is a
`DummyPatternIndex`, and a wildcard name is a tensor whose name begins with `?`, which is what the
matcher tests.
"""
function _parse_tensor_string(expr::String)
    length(expr) == 1 && return Tensor(expr)

    nameend = firstindex(expr)
    # "?A" is a wildcard name; "??A" additionally adopts the matched object's data
    if startswith(expr, "??")
        nameend = nextind(expr, nextind(expr, nameend))
    elseif startswith(expr, "?")
        nameend = nextind(expr, nameend)
    end
    symbol = expr[firstindex(expr):nameend]

    # A tensor name is one character. Several objects of the same kind are told apart by labels `T"k_1^μ"` and `T"k_2^μ"`, not by longer names.
    @assert length(expr) > 1 && contains("_^([{?", expr[nextind(expr, nameend)]) "Tensor name should be one character."

    positions = []
    indices = Vector{Symbol}()
    currpos = upper
    pause = false
    asym = false
    sym = false
    slotnum = 0
    functionoflist = []
    funof = false
    tabstart = 0
    tabaccum = Vector{Int64}()
    tableaulist = Set{YoungTableau}()
    coordstring = ""
    comma_mode = false
    i = firstindex(expr)

    while i <= lastindex(expr)
        if i <= nameend
            i = nextind(expr, i)
            continue
        end
        char = string(expr[i])
        if char == "^"
            currpos = upper
            comma_mode = has_comma_before_structure(expr, nextind(expr, i))
        elseif char == "_"
            currpos = lower
            comma_mode = has_comma_before_structure(expr, nextind(expr, i))
        elseif char == "|"
            pause = !pause
        elseif char == "["
            tabaccum = Vector{Int64}()
            asym = true
            tabstart = slotnum + 1
        elseif char == "]"
            asym = false
            push!(
                tableaulist,
                YoungTableau(fill(1, slotnum - tabstart + 1), tabaccum)
            )
        elseif char == "("
            tabaccum = Vector{Int64}()
            sym = true
            tabstart = slotnum + 1
        elseif char == ")"
            sym = false
            push!(tableaulist, YoungTableau([length(tabaccum)], tabaccum))
        elseif char == "{"
            funof = true
            # elseif char == ","
            #     i = nextind(expr, i)
            #     continue
        else
            if funof
                if char == "," || char == "}"
                    push!(functionoflist, :($(_coord_arg)(
                        $(esc(:(@isdefined $(Symbol(coordstring))))) ? $(esc(Symbol(coordstring))) : nothing,
                        $coordstring)))
                    coordstring = ""

                    if char == "}"
                        funof = false
                    end
                else
                    coordstring *= char
                end

            else
                slotnum += 1

                if expr[i] == '''
                    # 'quoted' -> always a label, whatever else that name might mean
                    j = nextind(expr, i)
                    k = j
                    while k <= lastindex(expr) && expr[k] != '''
                        k = nextind(expr, k)
                    end
                    push!(indices, Symbol("'" * expr[j:prevind(expr, k)]))
                    push!(positions, currpos)
                    i = k
                elseif expr[i] == '?'
                    # ?x -> wildcard index; consume the marker together with its name. A suffix is
                    # folded into its combining mark here, since there is no scope variable to do it.
                    idx, newi = next_index(expr, nextind(expr, i), comma_mode)
                    push!(indices, Symbol("?" * apply_index_suffix(string(special_index(idx)))))
                    push!(positions, currpos)
                    i = newi
                else
                    idx, newi = next_index(expr, i, comma_mode)
                    push!(indices, special_index(idx))
                    push!(positions, currpos)
                    i = newi
                end

                if !pause && asym
                    push!(tabaccum, slotnum)
                elseif !pause && sym
                    push!(tabaccum, slotnum)
                end
            end
        end

        i = nextind(expr, i)
    end

    functions_expr = Expr(:vect, functionoflist...)

    # indplacement = quote
    #     [$([:($(esc(Symbol(idx))) => $pos) for (idx, pos) in zip(indices, positions)]...)]
    # end
    slots = [slot_expr(idx, pos) for (idx, pos) in zip(indices, positions)]

    return quote
         $(Tensor)(
            $symbol,
            $(slots...),
            [$tableaulist...],
            $functions_expr
        )
    end
end

macro Tensor(expr::String)
    return _parse_tensor_string(expr)
end


# function has_comma_outside_curly(expr::String)
#     depth = 0

#     for c in expr
#         if c == '{'
#             depth += 1
#         elseif c == '}'
#             depth -= 1
#         elseif c == ',' && depth == 0
#             return true
#         end
#     end

#     return false
# end
function has_comma_before_structure(expr, i)
    j = i
    while j <= lastindex(expr)
        c = expr[j]
        if c == '{' || c == '^' || c == '_'
            return false
        elseif c == ','
            return true
        end
        j = nextind(expr, j)
    end
    return false
end
macro Tensor_str(expr::String)
    return esc(:(@Tensor $expr))
end

"""Combining mark written for each recognised index suffix, matching what `uni_to_TeX` decodes."""
const _INDEX_SUFFIX_MARKS = ("dot" => '̇', "bar" => '̄', "hat" => '̂')

"""
    apply_index_suffix(name) -> String

Rewrite a suffixed index name into the combining-mark form used throughout, so `xdot` becomes `̇x`.

A named index picks this up for free: `A_adot` resolves `adot` in the caller's scope, where the user
wrote the mark directly. A wildcard has no such variable, so the parser applies the convention
itself and `?xdot` denotes the dotted index rather than one literally called "xdot".
"""
function apply_index_suffix(nm::AbstractString)
    for (suf, mark) in _INDEX_SUFFIX_MARKS
        if length(nm) > length(suf) && endswith(nm, suf)
            return string(mark) * chop(nm, tail=length(suf))
        end
    end
    return String(nm)
end

function special_index(idx::Symbol)
    if idx == :+
        return Symbol("ℓ⁺")
    elseif idx == :-
        return Symbol("ℓ⁻")
    elseif idx == :±
        return Symbol("ℓ⁺⁻")
    elseif idx == :∓
        return Symbol("ℓ⁻⁺")
    else
        return Symbol(idx)
    end
end

function special_index(idx::String)
    if idx == "+"
        return Symbol("ℓ⁺")
    elseif idx == "-"
        return Symbol("ℓ⁻")
    elseif idx == "±"
        return Symbol("ℓ⁺⁻")
    elseif idx == "∓"
        return Symbol("ℓ⁻⁺")
    else
        return Symbol(idx)
    end
end