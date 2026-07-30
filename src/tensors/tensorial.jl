const LabelPair = Pair{String,IndexPosition}
const _EMPTY_FUNCTION_OF = Set{Coordinate}()
const _EMPTY_TABLEAUX = Set{YoungTableau}()
const _EMPTY_WEIGHTS = Dict{String,Number}()
const _EMPTY_LABELS = LabelPair[]
@inline coordfree() = _EMPTY_FUNCTION_OF
@inline notableaux()::Set{YoungTableau} = _EMPTY_TABLEAUX
@inline noweights()::Dict{String,Number} = _EMPTY_WEIGHTS
@inline nolabels()::Vector{LabelPair} = _EMPTY_LABELS
"""
Mutual tensor data that all tensors should have in common. This makes it easier to add new features to tensors.
"""
mutable struct TensorData <: TensorDataSuperType
    name::String
    indices::Vector{IndexPair}
    function_of::Set{Coordinate} #set if ok
    weights::Dict{String,Number}
    tableaux::Set{YoungTableau}
    labels::Vector{LabelPair}
    function TensorData(name::String, args...)
        indices = Vector{IndexPair}()
        function_of = coordfree();
        has_coords = false
        tableaux = notableaux();
        has_tabs = false
        labels = nolabels();
        has_labels = false
        weights = noweights()

        for arg in args
            if arg isa Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
                append!(indices, arg)
            elseif arg isa Pair{I,IndexPosition} where {I<:IndexSuperType}
                push!(indices, arg)
            elseif arg isa Coordinate
                has_coords || (function_of=Set{Coordinate}(); has_coords=true)
                push!(function_of, arg)
            elseif arg isa AbstractVector{<:Coordinate}
                isempty(arg) && continue
                has_coords || (function_of=Set{Coordinate}(); has_coords=true)
                union!(function_of, arg)
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
            end
        end

        tabs = has_tabs ? reduce_tableaux(tableaux) : tableaux

        if !has_coords
            return new(name, indices, function_of, weights, tabs, labels)
        end
        w = Dict{String,Number}()
        for c in function_of
            ;
            w[c.name] = 1;
        end
        return new(name, indices, function_of, merge(w, weights), tabs, labels)
    end
    # raw constructor — stores fields directly, no processing. For copying.
    function TensorData(name::String, indices::Vector{IndexPair},
        function_of::Set{Coordinate},
        weights::Dict{String,Number},
        tableaux::Set{YoungTableau},
        labels::Vector{LabelPair})
        return new(name, indices, function_of, weights, tableaux, labels)
    end
end

#pattern matching struct... don't touch
mutable struct DummyPatternTensor <: TensorSuperType
    data::TensorData
    function DummyPatternTensor(name::String, args...)
        return new(TensorData(name, args...))
    end
    function DummyPatternTensor(d::TensorData)
        return new(d)
    end
end



#####################
####### primitive tensor types
#####################

#simplest tensor constructu
"""Basic Tensor object."""
mutable struct Tensor <: TensorSuperType
    data::TensorData
    function Tensor(name::String, args...)
        return new(TensorData(name, args...))
    end
    function Tensor(d::TensorData)
        return new(d)
    end
end

"""Basic Fermion. Fermionic no matter the assigned indices. Identical to a tensor otherwise."""
mutable struct Fermion <: TensorSuperType
    data::TensorData
    function Fermion(name::String, args...)
        return new(TensorData(name, args...))
    end
    function Fermion(d::TensorData)
        return new(d)
    end
end
"""Basic Boson. Bosonic no matter the assigned indices. Identical to a tensor otherwise."""
mutable struct Boson <: TensorSuperType
    data::TensorData
    function Boson(name::String, args...)
        return new(TensorData(name, args...))
    end
    function Boson(d::TensorData)
        return new(d)
    end
end

function SymmetricTensor(name::String, args...)
    t = Tensor(name, args...)
    if length(indices(t))>1
        add_tableau!(data(t), generate_symmetric_tableau(indices(t)))
        # data(t).tableaux = reduce_tableaux(tableaux(t))
    end
    return t
end

function AntisymmetricTensor(name::String, args...)
    t = Tensor(name, args...)
    if length(indices(t))>1
        add_tableau!(data(t), generate_antisymmetric_tableau(indices(t)))
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
    add_tableau!(data(t), generate_riemann_tableau())
    return t
end


#######
# common tableaux generation functions
######
function generate_symmetric_tableau(indices::Vector)
    n = length(indices)
    return YoungTableau([n], collect(1:n))
end
function generate_antisymmetric_tableau(indices::Vector)
    n = length(indices)
    return YoungTableau(fill(1, n), collect(1:n))
end
function generate_symmetric_tableau(n::Integer)
    return YoungTableau([n], collect(1:n))
end
function generate_antisymmetric_tableau(n::Integer)
    return YoungTableau(fill(1, n), collect(1:n))
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

macro Tensor(expr::String)
    length(expr) == 1 && return Tensor(expr)

    symbol = string(first(expr))

    @assert length(expr) > 1 && contains("_^([{", expr[nextind(expr, firstindex(expr))]) "Tensor name should be one character."

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
        if i == firstindex(expr)
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
                    push!(functionoflist, :($(esc(Symbol(coordstring)))))
                    coordstring = ""

                    if char == "}"
                        funof = false
                    end
                else
                    coordstring *= char
                end

            else
                slotnum += 1

                idx, newi = next_index(expr, i, comma_mode)

                push!(indices, special_index(idx))
                push!(positions, currpos)

                i = newi

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

    indplacement = quote
        [$([:($(esc(Symbol(idx))) => $pos) for (idx, pos) in zip(indices, positions)]...)]
    end

    varname = Symbol(first(expr))
    # return :(Tensor($symbol, [$indplacement...], [$tableaulist...], $functions_expr))
    return quote
        $(esc(varname)) = Tensor(
            $symbol,
            [$indplacement...],
            [$tableaulist...],
            $functions_expr
        )
    end
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