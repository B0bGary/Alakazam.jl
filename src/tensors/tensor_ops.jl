set_default(convert_unicode=true)


##############################################################################################################################
##############################################################################################################################
# generic application of a function to tensorial objs. this sort of pattern occurs a lot
# i will progressively move to this syntax in future
transform(f, exp::TensorExpression) = sum(c * transform(f, t) for (c, t) in summands(exp))
transform(f, tt::TensorTerm) = prod(transform(f, t) for t in terms(tt))
transform(f, t::Tensorial) = f(t)




##############################################################################################################################
#equality checking
##############################################################################################################################
# == is a superficial check, ≅ is a deep check that checks for equivalence mod symmetries and dummy index relabeling. 
#hard_equals is a version of == that checks for equivalence mod symmetries and dummy index relabeling

#change to use our tabelau code!
"""
Two tensors are equal when they are the same type carrying the same name, indices,
coordinate dependence, tableaux and labels.
"""
function ==(a::TensorSuperType, b::TensorSuperType)#name first since easiest discriminator
    if (name(a) == name(b) && typeof(a) == typeof(b) && indices(a) == indices(b) && _same_coords(function_of(a), function_of(b)) && tableaux(a)==tableaux(b) && labels(a)==labels(b))
        # println("apparnetly equl inds")
        return true
    else
        # println("false at tens level")
        return false
    end
end

function Base.hash(d::TensorSuperType, h::UInt)::UInt
    h = hash(typeof(d), h)
    h = hash(name(d), h)
    @inbounds for p in indices(d)
        h = hash(p[1].uid, hash(p[2].position, h))
    end
    # labels are ordered and are compared as vectors by `==`, so they are hashed in order too.
    @inbounds for l in labels(d)
        h = hash(l[1], hash(l[2].position, h))
    end
    x = UInt(0)
    for c in function_of(d)
        ;
        x ⊻= hash(c);
    end
    for t in tableaux(d)
        ;
        x ⊻= hash(t);
    end
    return hash(x, h)
end
# add single term check to make things less expensive
# isequal, not ==
==(a::TensorExpression, b::TensorTerm) = singleterm(a) && isequal(first_coef(a), 1) && (first_summand(a) == b)
==(a::TensorTerm, b::TensorExpression) = b == a
==(a::TensorSuperType, b::TensorExpression) = b == a
function ==(a::TensorExpression, b::TensorSuperType)
    tens = isTensor(a)
    if (tens !== nothing)
        return tens == b
    else
        return false
    end
end
function ==(a::TensorSuperType, b::TensorTerm)
    tens = isTensor(b)
    if (tens !== nothing)
        a == tens
    else
        return false
    end
end
==(a::TensorTerm, b::TensorSuperType) = b == a
function ==(a::TensorTerm, b::TensorTerm)
    res = terms(a) == terms(b)
    return res
end
#hashing for faster merge_sum
function Base.hash(t::TensorTerm, h::UInt)::UInt
    for term in terms(t)
        h = hash(term, h)
    end
    return h
end

function ==(a::TensorExpression, b::TensorExpression)
    if length(a.summands) != length(b.summands)
        (length(a.summands)==0 && length(b.summands)==1 && coeff_equal(first_coef(b), 0)) || (length(b.summands)==0 && length(a.summands)==1 && coeff_equal(first_coef(a), 0)) ? (return true) : (return false)
    end
    aa = a.summands
    bb = b.summands
    if same_summands(aa, bb)
        return true
    end
    return false

end

==(N::SymorNum, T::Tensorial) = (T == N)
==(T::Tensorial, N::SymorNum) = TensorExpression(T) == TensorExpression(N)

coeff_equal(a::Number, b::Number) =
    a isa Num || b isa Num ? coeff_equal(Symbolics.unwrap(a), Symbolics.unwrap(b)) : a == b
coeff_equal(a, b) = isequal(a, b) || _sym_equal(a, b)   # symbolic fallback
# @noinline _sym_equal(a, b) = Symbolics.simplify(Symbolics.expand(a - b)) === 0 #1.11 version
@noinline function _sym_equal(a, b)#1.12 version
    d = Symbolics.simplify(Symbolics.expand(a - b))
    SymbolicUtils.isconst(d) && return iszero(SymbolicUtils.unwrap_const(d))
    return d isa Number && iszero(d)
end

function same_summands(a, b)
    # @show a
    # @show b
    if isequal(a, b)#naive check
        return true
    end
    if length(a) == 1 && isempty(terms(a[1][2])) #already, the leading coeffs are different to get here
        # @show "false here"
        return false
    end
    #currently doesnt do operator checks...
    # if singleterm(a.second) isa OperatorSuperType && singleterm(b.second) isa OperatorSuperType
    #     return singleterm(a.second) == singleterm(b.second)
    # end
    # eq(x, y) = ((is_equal(x.first, y.first) && (singleterm(x.second) && singleterm(y.second) && !(singleterm(x.second) isa OperatorSuperType) && !(singleterm(y.second) isa OperatorSuperType) ? same_tensor_mod_syms(x.second.terms[1], y.second.terms[1])[1] : same_tensorterm_mod_syms(TensorTerm(x.second), TensorTerm(y.second))))

    eq(x, y) = ((coeff_equal(x.first, y.first) && (singleterm(x.second) && singleterm(y.second) ? same_tensor_mod_syms(x.second.terms[1], y.second.terms[1])[2] == 1 : same_tensorterm_mod_syms(TensorTerm(x.second), TensorTerm(y.second))))
                ||
                coeff_equal(x.first, -y.first) && same_tensorterm_mod_syms_grade(TensorTerm(x.second), TensorTerm(y.second)) == (true, -1)
    )
    length(a) == length(b) || return false

    used = falses(length(b))
    for x in a
        found = false
        for j in eachindex(b)
            if !used[j] && eq(x, b[j])
                used[j] = true
                found = true
                break
            end
        end
        found || return false
    end

    return true
end


function hard_equals(a::TensorTerm, b::TensorTerm)
    return same_tensorterm_mod_syms(a, b)
end

function hard_equals_graded(a::TensorTerm, b::TensorTerm; g1=nothing, g2=nothing,
    ws1=nothing, ws2=nothing)
    return same_tensorterm_mod_syms_grade(a, b; g1, g2, ws1, ws2)
end


cong(T1, T2) = T1≅T2

function ≅(a::TensorTerm, b::TensorTerm)
    return hard_equals(a, b)
end

≅(a::TensorSuperType, b::TensorExpression) = (b ≅ a)
≅(a::TensorExpression, b::TensorSuperType) = a ≅ TensorExpression(b)


_eq_one(x, y) = x == y      # specialises on concrete typeof(x)
function ==(a::Vector{T} where {T<:TensorSuperType}, b::Vector{T} where {T<:TensorSuperType})
    length(a) == length(b) || return false
    @inbounds for i in eachindex(a)
        _eq_one(a[i], b[i]) || return false
    end
    return true
end

≅(b::TensorExpression, a::TensorTerm) = (a ≅ b)
≅(a::TensorTerm, b::TensorExpression) = (TensorExpression(a) ≅ b)
≅(a::OperatorSuperType, b::OperatorSuperType) = ≅(TensorExpression(a), TensorExpression(b))
function ≅(a::TensorSuperType, b::TensorSuperType)#change to use our tabelau code!
    function_of(a)==function_of(b) || return false
    (same, grweight) = same_tensor_mod_syms(a, b)
    return (same && grweight == 1)
end
"""
`≅` is the deep equality test, so it opts into the multi-term (Young projector) comparison that
`hard_simplify` leaves off by default. Two expressions differing by a Bianchi-type identity are
therefore congruent even though neither `==` nor `hard_simplify` alone would merge them.
"""
function ≅(a::TensorExpression, b::TensorExpression)
    a = hard_merge(a)
    b = hard_merge(b)
    C=merge_sum_syms(a-b)
    m = hard_merge(symmetrise_partials(C); multiterm=true)
    return m == zero()
end

function ≅(a::OperatorSuperType, b::TensorExpression)
    ≅(TensorExpression(a), b)
end

##############################################################################################################################
##############################################################################################################################
### Weightings
##############################################################################################################################
##############################################################################################################################

"""
Extract the overall weight value from a TensorTerm for a given weight label
"""
function term_weight(T::TensorTerm, w::String)
    weight = 0
    for term in terms(T)
        if typeof(term) <: OperatorSuperType
            weight += operator_weight(term, w)
        elseif haskey(weights(term), w)
            weight += weights(term)[w]
        end
    end
    return weight
end


"""
Get all summands in a TensorExpression of a particular weight. For example, this may be terms that carry a particular polynomial weighting.
"""
function get_weighted_terms(T::TensorExpression, w::String, value::Number)
    new_exp = zero()
    for summand in summands(T)
        if term_weight(summand[2], w) == value
            new_exp += summand[1] * summand[2]
        end
    end
    return new_exp
end

"""
Get all summands in a TensorExpression of below particular weight. For example, this may get leading order terms in a Taylor expansion.
Returns only weights<=value. Assumes weights are non-negative.
"""
function drop_weight_above(T::TensorExpression, w::String, value::Number)
    new_exp = zero()
    i = 0
    for summand in summands(T)
        if term_weight(summand[2], w) <= value
            new_exp += summand[1] * summand[2]
        end
    end
    return new_exp
end

"""
Get all summands in a TensorExpression of above particular weight. For example, this may get terms with dominant powers.
Returns only weights>=value. Assumes weights are non-negative.
This should be changed to be better
"""
function drop_weight_below(T::TensorExpression, w::String, value::Number)
    new_exp = zero()
    for summand in summands(T)
        if term_weight(summand[2], w) >= value
            new_exp += summand[1] * summand[2]
        end
    end
    return new_exp
end

##############################################################################################################################
##############################################################################################################################
### Advanced Mathematical Operations
##############################################################################################################################
##############################################################################################################################

# juxtaposition/appension/smart times, basically appends a tensor to a term so that all of its indices are free. 
# useful if you want to mutliply something without doing contractions, or want T^3 without manually having to create objects with new indices

# function ×(T1::TensorTerm, T2::TensorSuperType)

append(T1, T2) = T1×T2
# end
×(T2::SymorNum, N::Union{TensorExpression,TensorTerm,TensorSuperType}) = T2*N

"""
I previously had here an operator precendence error. But, the user should just bracket correctly.
"""
×(T2::Union{TensorExpression,TensorTerm,TensorSuperType}, N::SymorNum) = T2*N


×(N1::SymorNum, N2::SymorNum) = TensorExpression(N1*N2)


function ×(T2::TensorExpression, T1::TensorExpression)
    out = CoeffPair[]
    sizehint!(out, length(T2.summands) * length(T1.summands))
    for summand2 in T2.summands
        for summand1 in T1.summands
            append!(out, summands((summand2[1] * summand1[1]) * (summand2[2] × summand1[2])))
        end
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end

function ×(T1::TensorExpression, T2::Union{TensorSuperType,TensorTerm})
    T1 = rename_dummies(T1)
    out = CoeffPair[]
    sizehint!(out, length(T1.summands))
    for summand in T1.summands
        append!(out, summands(summand[1] * (summand[2] × T2)))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end

function ×(T2::Union{TensorSuperType,TensorTerm}, T1::TensorExpression)
    T1 = rename_dummies(T1)
    out = CoeffPair[]
    sizehint!(out, length(T1.summands))
    for summand in T1.summands
        append!(out, summands(summand[1] * (T2 × summand[2])))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end


function ×(T1::Union{TensorTerm,TensorSuperType}, T2::Union{TensorTerm,TensorSuperType})
    taken = OrderedSet{IndexSuperType}()
    collect_indices!(taken, T1)
    return T1 * rename_dummies_excluded(T2, taken, ())      #free indices too: × contracts nothing
end

"""
    rename_dummies_excluded!(T, taken, keep=free indices of T) -> T

Rename, in place, the indices `T` uses, inside operands too, that are in `taken` and not in `keep`, to
ones neither uses, each occurrence alike. By default only dummies are renamed, so the value is
unchanged; pass `keep=()` to rename clashing free indices as well. `taken` may hold indices or index
pairs. Lightcone indices, which repeat freely, and pattern indices are left alone.

`T` is edited, and with it any tensor it shares with other expressions (`TensorExpression(A)` holds
`A` itself), so pass a copy of anything you do not own: `rename_dummies_excluded!(copy(E), taken)`.
"""
function rename_dummies_excluded!(T::Tensorial, taken, keep=nothing)
    any(q -> q isa Pair, taken) && (taken = [q isa Pair ? q[1] : q for q in taken])
    keep === nothing && (keep = Set(p[1] for p in free_indices_ordered(T)))
    used = OrderedSet{IndexSuperType}()
    collect_indices!(used, T)
    claimed = nothing                         #built at the first clash, which there usually is not
    rep = Dict{IndexPair,IndexPair}()
    for q in used
        (isLightcone(q) || q isa DummyPatternIndex || !(q in taken) || q in keep) && continue
        claimed === nothing && (claimed = IndexTypes[r for r in Iterators.flatten((taken, used))])
        new = _fresh_index(q.index_set, claimed)
        push!(claimed, new)
        rep[q=>upper] = new => upper
        rep[q=>lower] = new => lower
    end
    isempty(rep) || set_indices_deep!(T, rep)
    return T
end

"""`rename_dummies_excluded!` on a copy, leaving `T` and everything it shares untouched."""
rename_dummies_excluded(T::Tensorial, taken, keep=nothing) = rename_dummies_excluded!(copy(T), taken, keep)

##############################################################################################################################
##############################################################################################################################
#indices
##############################################################################################################################
#update indices

update_indices!(T::TensorSuperType) = T

"""Below this many indices the pairwise duplicate scan beats building a count table."""
const _EXACT_SCAN_CUTOFF = 16

"""
Updates external indices of a tensor term, expression, or operator, by looking at the indices of all the terms it contains.
"""
function update_indices!(T::TensorTerm, n::Int=-1)
    n == -1 && (n = length(indices(T)))
    inds = Vector{Pair{IndexTypes,IndexPosition}}(undef, n)
    k = 1
    for t in terms(T)
        typeof(t) <: OperatorSuperType && update_indices!(t)
        k = _copy_inds!(inds, k, t)
    end

    # cheap prefilter: does any non-lightcone index appear more than once?
    mask = UInt64(0)
    nnl = 0
    @inbounds for p in inds
        isLightcone(p[1]) && continue
        nnl += 1
        mask |= _bit(p[1])
    end
    if count_ones(mask) != nnl
        if nnl <= _EXACT_SCAN_CUTOFF
            @inbounds for i in eachindex(inds)
                idx = inds[i][1]
                isLightcone(idx) && continue
                c = 1
                for j in (i+1):length(inds)
                    inds[j][1] == idx && (c += 1)
                end
                c > 2 && throw(OverusedIndices(T))
            end
        else
            counts = Dict{UInt64,Int}()
            sizehint!(counts, nnl)
            @inbounds for p in inds
                isLightcone(p[1]) && continue
                c = get(counts, p[1].uid, 0) + 1
                c > 2 && throw(OverusedIndices(T))
                counts[p[1].uid] = c
            end
        end
    end

    set_indices!(T, inds)
    return T
end

function update_indices!(T::SymorNum)
    return T
end
function update_indices!(T::TensorExpression)
    for summand in summands(T)
        update_indices!(summand[2])
    end
    return T
end

function get_unused_indices(T::Union{TensorSuperType,TensorTerm,TensorExpression}, IS::IndexSet)
    used = IndexTypes[]
    collect_indices!(used, T)                  # needs a Vector-accepting collect_indices!
    out = IndexTypes[]
    for i in IS.indices
        i isa DummyPatternIndex && continue    # inline filter_wildcard_dummies
        i in used && continue                   # inline setdiff
        push!(out, i)
    end
    return out
end
function filter_wildcard_dummies(I::Vector{IndexTypes})
    return filter(i -> !(i isa DummyPatternIndex), I)
end


function get_index_index(T::TensorSuperType)
    return OrderedSet([indices(T)[i][1] => i for i in 1:size(indices(T))[1]])
end


function index_index_to_dict(s::OrderedSet)
    return Dict(p[2] => p[1] for p in s)
end

function rename_dummies(T::TensorSuperType)
    T = copy(T)
    dumms = dummies(T)
    isempty(dumms) && return T

    index_sets = IndexSet[]
    for d in dumms
        d.index_set in index_sets || push!(index_sets, d.index_set)
    end

    all_indices = IndexTypes[]
    for i in indices(T)
        i[1] in all_indices || push!(all_indices, i[1])
    end

    replacements = Pair{IndexSuperType,IndexSuperType}[]
    inset = IndexTypes[]
    remain = IndexTypes[]
    ordered = IndexTypes[]

    for set in index_sets
        # indices of this set already on the tensor
        empty!(inset)
        for a in all_indices
            a.index_set === set && push!(inset, a)
        end

        # this set's non-wildcard indices, in declaration order
        empty!(ordered)
        for x in set.indices
            x isa DummyPatternIndex || push!(ordered, x)
        end

        # candidates: not already in use
        empty!(remain)
        for x in ordered
            x in inset || push!(remain, x)
        end

        for d in dumms
            d.index_set === set || continue
            for (ri, r) in enumerate(remain)
                if first_element(ordered, d, r) != d
                    push!(replacements, d => r)
                    deleteat!(remain, ri)
                    break
                end
            end
        end
    end

    isempty(replacements) && return T
    inds = indices(T)
    out = Vector{IndexPair}(undef, length(inds))
    @inbounds for i in eachindex(inds)
        x, y = inds[i]
        for p in replacements
            if x == p[1]
                x = p[2]
                break
            end
        end
        out[i] = x => y
    end
    set_indices!(T, out)

    return T
end
"""
Map the position of every free index through f, leaving dummy pairs untouched.
Dummies are skipped deliberately.
Which indices are dummies is a property of the whole term, not of any one tensor, so the set is
computed here and passed down.
"""
function _map_free_positions(tt::TensorTerm, f)
    dums = Set{UInt64}(d.uid for d in dummies(tt))
    return prod(set_indices_deep!(copy(t),
        IndexPair[(p[1].uid in dums ? p : (p[1] => f(p[2]))) for p in indices(t)])
                for t in terms(tt))
end
function _map_free_positions(ex::TensorExpression, f)
    isZeroTens(ex) && return ex
    return sum(c * _map_free_positions(tt, f) for (c, tt) in summands(ex))
end
_map_free_positions(t::TensorSuperType, f) = _map_free_positions(TensorTerm(t), f)

"""
    raise_all(T)

Return a copy of `T` with every free index moved to the upper position. Dummy pairs are left as
they are. See also [`lower_all`](@ref) and [`flip_all`](@ref).
"""
raise_all(T::Tensorial) = _map_free_positions(T, _ -> upper)

"""
    lower_all(T)

Return a copy of `T` with every free index moved to the lower position; the counterpart of
[`raise_all`](@ref).
"""
lower_all(T::Tensorial) = _map_free_positions(T, _ -> lower)

"""
    flip_all(T)

Return a copy of `T` with every free index position reversed, upper becoming lower and lower upper.
Dummy pairs are left untouched, as for [`raise_all`](@ref) and [`lower_all`](@ref).
"""
flip_all(T::Tensorial) = _map_free_positions(T, op_pos)

"""Reverse the positions of the free indices of an expression; an alias for [`flip_all`](@ref)."""
op_pos(tt::Tensorial) = flip_all(tt)




##############################################################################################################################
"""Update the external indices of TensorTerms and Operators from the true inner indices."""
function update_ext_indices!(T::TensorExpression)
    new_exp = zero()
    for tens in summands(T)
        if isZeroTens(tens[1])
            continue
        end
        update_ext_indices!(tens[2])
        new_exp += tens[1] * tens[2]
    end
    T = new_exp
end

function update_ext_indices!(T::TensorTerm)
    inds = Vector{Pair{IndexTypes,IndexPosition}}()
    for tens in terms(T)
        if typeof(tens) <: OperatorSuperType
            for summand in summands(operand(tens))
                for term in terms(summand[2])
                    if term isa OperatorSuperType
                        update_ext_indices!(term)
                    end
                end
                # set_operand_indices!(tens, indices(summands(operand(tens))[1][2]))
                set_indices!(tens, vcat(op_indices(tens), operand_indices(tens)))

            end
        end
        for i in indices(tens)
            push!(inds, i)

        end
    end
    set_indices!(T, inds)

end

#rename dummies function
"""Rename dummy Einstein summation indices to be more lexigraphically uniform across expressions."""
function rename_dummies(T::TensorTerm)
    T = copy(T)
    dumms = dummies(T)
    isempty(dumms) && return T

    index_sets = IndexSet[]
    for d in dumms
        d.index_set in index_sets || push!(index_sets, d.index_set)
    end

    all_indices = IndexTypes[]
    for t in terms(T), i in indices(t)
        i[1] in all_indices || push!(all_indices, i[1])
    end

    dumm_uids = Set{UInt64}()
    for d in dumms
        push!(dumm_uids, d.uid)
    end
    inset_uids = Set{UInt64}()
    remain = IndexTypes[]

    for set in index_sets
        # indices of this set already in use, excluding the dummies being replaced
        empty!(inset_uids)
        for a in all_indices
            a.index_set === set || continue
            a.uid in dumm_uids && continue
            push!(inset_uids, a.uid)
        end
        # candidates: this set's non-wildcard indices not already in use
        empty!(remain)
        for x in set.indices
            x isa DummyPatternIndex && continue
            x.uid in inset_uids && continue
            push!(remain, x)
        end
        # keyed by uid: an abstract key type makes every hash a dynamic dispatch
        replacements = Dict{UInt64,IndexSuperType}()
        cursor = 1                     # walk `remain` rather than popfirst!, which shifts it
        for d in dumms
            d.index_set === set || continue
            if cursor > length(remain)
                nid = generate_indices(set)
                nid === nothing && throw(OutOfIndices(set))
                push!(remain, nid)
            end
            replacements[d.uid] = remain[cursor]
            cursor += 1
        end
        isempty(replacements) || (T = rename_recursive!(T, replacements))
    end
    return T
end

function rename_recursive!(T::TensorExpression, replacements)
    s = summands(T)
    isempty(s) && return T
    out = CoeffPair[]
    sizehint!(out, length(s))
    for (coef, tt) in s
        renamed = rename_recursive!(tt, replacements)
        append!(out, summands(coef * renamed))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end

function rename_recursive!(T::TensorTerm, replacements)
    for term in terms(T)
        rename_recursive!(term, replacements)
    end
    update_ext_indices!(T)
    return T
end

function rename_recursive!(t::TensorSuperType, replacements)
    if t isa OperatorSuperType
        oi = op_indices(t)
        if any(p -> haskey(replacements, p[1].uid), oi)
            set_op_indices!(t, IndexPair[Pair(get(replacements, x.uid, x), y) for (x, y) in oi])
        end
        rename_recursive!(operand(t), replacements)
    else
        ti = indices(t)
        if any(p -> haskey(replacements, p[1].uid), ti)
            set_indices!(t, IndexPair[Pair(get(replacements, x.uid, x), y) for (x, y) in ti])
        end
    end
    return t
end


function rename_dummies(T::TensorExpression)
    s = summands(T)
    isempty(s) && return T
    out = CoeffPair[]
    sizehint!(out, length(s))
    for (coef, tt) in s
        renamed = rename_dummies(tt)
        append!(out, summands(coef * renamed))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end
function rename_dummies!(T::TensorExpression)
    s = summands(T)
    for i in eachindex(s)
        s[i] = CoeffPair(s[i][1], rename_dummies(s[i][2]))
    end
    return T
end


function update_dummies(T::TensorSuperType, reps::Pair{I,I}) where {I<:IndexSuperType}
    if typeof(T) <: OperatorSuperType
        update_dummies(T.operand, reps)
        set_op_indices!(T) = [Pair((x == reps[1] ? reps[2] : x), y) for (x, y) in op_indices(T)]
        set_operand_indices!(T) = [Pair((x == reps[1] ? reps[2] : x), y) for (x, y) in op_indices(T)]
    end

    if !isempty(indices(T))
        indices(T) = [Pair((x == reps[1] ? reps[2] : x), y) for (x, y) in indices(T)]
    end
end

function update_dummies(T::TensorTerm, reps::Pair{I,I}) where {I<:IndexSuperType}
    for t in terms(T)
        update_dummies(t, reps)
    end
end
function update_dummies(TEX::TensorExpression, reps::Pair{I,I}) where {I<:IndexSuperType}
    for summand in TEX.summands
        for t in summand[2].terms
            update_dummies(t, reps)
        end
    end
end

"""
Relabel indices of term so they don't use anything from extra_inds.
Uses the appension operation to relabel.
"""
function reserve_extra_indices(extra_inds::Vector{Pair{I,IndexPosition}}, trm) where I<:IndexSuperType
    isempty(extra_inds) && return trm
    tt = isTensorTerm(Tensor("_", extra_inds) × trm)
    tt === nothing && (tt = isTensorTerm(-Tensor("_", extra_inds) × trm))
    return TensorTerm(terms(tt)[2:end])
end


##############################################################################################################################
##############################################################################################################################
### index sorting and canonicalisation
##############################################################################################################################
##############################################################################################################################

"""Computes the distance between two configurations of indices, as the number of inversions needed to get from one to the other. 
Used as a cost function for finding the canonical representation of a tensor.
Penalised out of alphabetical order, and also by the number of swaps needed to get to the canonical order."""
function rep_distance(canon::Vector, config::Vector)
    inv = 0
    n = length(config)
    for i in 1:n
        for j in (i+1):n
            if config[i] > config[j]
                inv += 1
            end
        end
    end

    return inv
end

@inline _signbit(s::Int) = s > 0 ? 0x01 : 0x02

"""
Factors whose levels the metric of their sets does not move: a metric, etc.
Used only for signature checking, not for blocking metric raising/lowering
"""
@inline _keeps_levels(t::TensorSuperType) =
    t isa MetricSuperType || t isa KroneckerSuperType

"""
Does a dummy index contraction vanish due to symmetry against the symmetry of the metric?
One pass over the slots. Up to 64 of them, a mask of uid bits flags a slot that may repeat an earlier
index, and only those look back for it, so a factor with no trace, nearly all of them, allocates nothing
and compares nothing. Beyond 64 the mask fills, and a table of first slots takes its place.
"""
function _trace_vanishes(T::TensorSuperType, inds)
    tabs = tableaux(T)
    isempty(tabs) && return false
    (T isa OperatorSuperType || _keeps_levels(T)) && return false
    n = length(inds)
    if n <= 64
        mask = UInt64(0)
        @inbounds for q in 1:n
            idx = inds[q][1]
            b = _bit(idx)
            if mask & b != 0
                for p in 1:(q-1)
                    inds[p][1] === idx || continue
                    _traced_pair_vanishes(tabs, inds, p, q) && return true
                    break
                end
            end
            mask |= b
        end
    else
        firstat = Dict{UInt64,Int}()
        sizehint!(firstat, n)
        @inbounds for q in 1:n
            p = get!(firstat, inds[q][1].uid, q)
            p != q && _traced_pair_vanishes(tabs, inds, p, q) && return true
        end
    end
    return false
end

#the trace over slots p < q, which hold one index, by the sign its symmetry and its metric give
@inline function _traced_pair_vanishes(tabs, inds, p::Int, q::Int)
    @inbounds idx, pos = inds[p]
    @inbounds inds[q][2] === pos && return false
    #a lightcone component repeated is not summed over
    (isLightcone(idx) || idx isa DummyPatternIndex) && return false
    mt = idx.index_set.metric_type
    mt === nothing && return false
    m = mt <: SymMetricSuperType ? 1 : mt <: AntisymMetricSuperType ? -1 : 0
    for tab in tabs
        #adj_matrix holds the sign of each slot exchange the tableau makes a symmetry, 0 for none
        q <= size(tab.adj_matrix, 1) || continue
        @inbounds tab.adj_matrix[p, q] * m == -1 && return true
    end
    return false
end
"""
    find_canonical_rep(T::TensorSuperType) -> (sign, indices)

Find the canonical index ordering for `T` under the symmetries described by its Young tableaux, together with the sign picked up in reaching it.

Canonical form is the lexicographic ordering. Shortcut for fully symmetric or fully antisymmetric tensor using perm.
Otherwise the tableaux generate a group of allowed index permutations, and the function searches it breadth-first for a representative matching the lexicographic target.
"""
function find_canonical_rep(T::TensorSuperType)::Tuple{Int,Vector{IndexPair}}
    inds = indices(T)
    # Nothing to reorder below two indices
    length(inds) <= 1 && return (1, inds)
    perm = sortperm(inds)
    lexo = inds[perm]
    _trace_vanishes(T, inds) && return (0, lexo)

    isSymmetric(T) && return (1, lexo)
    isAntisymmetric(T) && return (permutation_sign(perm), lexo)

    start_inds = inds
    best_inds = start_inds
    best_cost = rep_distance(lexo, start_inds)
    best_sign = 1
    # generators depend only on T's tableaux, so they are cached on each YoungTableau on first use
    tabs = tableaux(T)
    if length(tabs) == 1
        tab = ensure_cycles!(first(tabs))
        row_cycles = tab.row_cycles
        col_cycles = tab.col_cycles
    else
        row_cycles = Vector{PermCycles}()
        col_cycles = Vector{PermCycles}()
        for tab in tabs
            ensure_cycles!(tab)
            append!(row_cycles, tab.row_cycles)
            append!(col_cycles, tab.col_cycles)
        end
    end
    queue = Tuple{Vector{IndexPair},Int}[(start_inds, 1)]
    head = 1
    # index pattern => bitmask of signs already seen (0x01 = +1, 0x02 = -1)
    seen = Dict{Vector{IndexPair},UInt8}()
    seen[copy(start_inds)] = 0x01
    while head <= length(queue)
        curr, sign = queue[head]
        head += 1
        cost = rep_distance(lexo, curr)
        cost == 0 && return (sign, curr)
        if cost < best_cost || (cost == best_cost && curr < best_inds)
            best_cost = cost
            best_inds = curr
            best_sign = sign
        end
        for cycle in row_cycles
            new_inds = apply_permutation_cycles(curr, cycle)
            b = _signbit(sign)
            m = get(seen, new_inds, 0x00)
            if m & b == 0x00
                seen[new_inds] = m | b
                push!(queue, (new_inds, sign))
            end
        end
        for cycle in col_cycles
            new_inds = apply_permutation_cycles(curr, cycle)
            new_sign = -sign
            b = _signbit(new_sign)
            m = get(seen, new_inds, 0x00)
            if m & b == 0x00
                seen[new_inds] = m | b
                push!(queue, (new_inds, new_sign))
            end
        end
    end

    return (best_sign, best_inds)
end

function permutation_sign(p)

    inv = 0
    n = length(p)

    for i in 1:n
        for j in (i+1):n
            inv += p[i] > p[j]
        end
    end

    return isodd(inv) ? -1 : 1
end
"""
Applies a permutation, given as disjoint cycles of slot positions, to `v`. Returns a new vector.
"""
function apply_permutation_cycles(v::Vector{T}, cycles::PermCycles) where {T}
    result = copy(v)
    @inbounds for c in cycles
        k = length(c)
        k < 2 && continue
        temp = result[c[k]]
        for i in k:-1:2
            result[c[i]] = result[c[i-1]]
        end
        result[c[1]] = temp
    end
    return result
end

"""Given a permutation cycle like "(1 5 2)", applies the corresponding permutation to the vector."""
function apply_permutation_cycle(v::Vector, cycle::String)
    result = copy(v)
    if cycle == "()"
        return result
    end
    # extract cycles like "(1 5 2)"
    cycles = collect(eachmatch(r"\(([^()]*)\)", cycle))

    for c in cycles
        # parse numbers inside this cycle
        nums = parse.(Int, split(c.captures[1]))
        k = length(nums)

        # rotate values according to cycle
        temp = result[nums[end]]
        for i in k:-1:2
            result[nums[i]] = result[nums[i-1]]
        end
        result[nums[1]] = temp
    end

    return result
end
"""
Sorts indices of a tensor into a canonical order, and applies the appropriate symmetries if the tensor is symmetric or antisymmetric. 
Swaps are based on permutations allowed by the Young Tableau
"""
function sort_symmetric(T::TensorSuperType)
    (sign, new_inds) = find_canonical_rep(T)
    sign == 0 && return zero()
    T_new = copy(T)
    set_indices!(T_new, new_inds)
    return sign * T_new
end

function sort_symmetric!(T::TensorSuperType)
    (sign, new_inds) = find_canonical_rep(T)
    set_indices!(T, new_inds)
    return (sign, T)
end


function simplify(T::TensorSuperType)#incase of operators.. (ops have their own)
    sort(T)
end


function simplify(TEx::TensorTerm)
    simplify(TensorExpression(TEx))
end


function simplify(TEx::TensorExpression)
    # Collect first
    exp=symmetrise_partials(collect_terms(TEx))
    exp = pop_metric(exp)
    exp = sort(exp)
    exp = eliminate_epsilon(exp)
    new_exp = simp_ops_n_sym(exp)
    new_exp = contract_metrics(new_exp, Metric)
    new_exp = eliminate_kronecker(new_exp)
    new_exp = collect_terms(new_exp)
    new_exp = rename_dummies!(new_exp)
    new_exp = sort(new_exp)
    return new_exp
end

canonicalise(T::Tensorial) = canonicalise(TensorExpression(T), true)
function canonicalise(TEx::TensorExpression, cp=false)
    cp==false && (TEx=copy(TEx))
    TEx=merge_sum_syms(TEx)
    return sort(TEx)
end
canonicalize(T) = canonicalise(T)


function simp_ops_n_sym(exp::TensorExpression)
    out = CoeffPair[]
    sizehint!(out, length(exp.summands))
    for summand in exp.summands
        sorted = Any[]
        vanishes = false
        for tens in summand[2].terms
            if typeof(tens) <: OperatorSuperType
                op = simp_ops_n_sym(operand(tens))
                if op == zero() || op === nothing
                    vanishes = true
                    break
                end
                tens = copy(tens)       #copy on write
                set_operand!(tens, op) #operator on nth term
            end
            st = sort_symmetric(tens)
            if isZeroTens(st)
                vanishes = true
                break
            end
            push!(sorted, st)
        end
        vanishes && continue
        append!(out, summands(TensorExpression(_factor_product(sorted, summand[1]))))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end


hard_simplify(T::TensorTerm; kw...) = hard_simplify(TensorExpression(T); kw...)
hard_simplify(T::TensorSuperType; kw...) = hard_simplify(TensorExpression(T); kw...)

"""
    hard_simplify(TEx; multiterm=false, multiterm_budget=DEFAULT_MULTITERM_BUDGET)

Merge terms up to symmetry and dummy relabelling, more aggressively than `simplify`.

By default only mono-term symmetries are used decided by graph isomorphism. Setting `multiterm=true` additionally
tries the Young projector on pairs that survive, which reaches identities spanning several terms
such as the first Bianchi identity. That test is exponential in the number of symmetric factors,
so it is opt-in and bounded by `multiterm_budget`; see [`multiterm_equal`](@ref).
"""
function hard_simplify(TEx::TensorExpression; multiterm::Bool=false,
    multiterm_budget::Int=DEFAULT_MULTITERM_BUDGET)
    return hard_merge(sort(collect_terms(TEx)); multiterm, multiterm_budget)
end

"""
Cheap test to see of terms are definitely different to avoid hard merge
"""
function term_fingerprint(T::TensorTerm)
    h = UInt(0)
    for t in terms(T)
        h ⊻= hash(length(indices(t)), hash(name(t)))
    end
    x = UInt(0)
    inds = indices(T)
    @inbounds for i in eachindex(inds)
        _is_dummy_at(inds, i) && continue
        x ⊻= hash(inds[i][1].uid, hash(inds[i][2].position))
    end
    return hash(x, hash(length(terms(T)), h))
end

function hard_merge(TEx::TensorExpression; multiterm::Bool=false,
    multiterm_budget::Int=DEFAULT_MULTITERM_BUDGET)
    src = TEx.summands
    n = length(src)
    n == 0 && return zero()

    out = Vector{CoeffPair}(undef, n)
    fps = Vector{UInt}(undef, n)          # fingerprint of each kept summand
    # the graph of each kept summand, built once here rather than twice inside every comparison
    gs = Vector{LabeledMultiGraph}(undef, n)
    ivs = Vector{UInt}(undef, n)          # isomorphism-invariant hash of that graph
    m = 0

    @inbounds for i in 1:n
        s = src[i]
        isZero(s[1]) && continue
        t = s[2]
        # Identically zero: some automorphism of the term sends it to minus itself. The single
        # transposition test runs first because it is cheap and catches most cases outright; the
        # general search is only reached when it finds nothing, and is budgeted so a term whose
        # automorphism group is too large to explore is reported as non-zero rather than guessed.
        _exchange_vanishes(t) && continue
        g = construct_tensorterm_graph(t)
        find_weighted_square(g) && continue

        fp = term_fingerprint(t)
        iv = graph_invariant(g)
        found = false
        for k in 1:m
            fps[k] == fp || continue      # cheap reject before the isomorphism search
            ivs[k] == iv || continue      # stronger one, far cheaper than vf2
            # every stored term already failed find_weighted_square above, so both flags are known
            eq, gr = hard_equals_graded(t, out[k][2]; g1=g, g2=gs[k], ws1=false, ws2=false)
            eq || continue
            c = gr * out[k][1] + s[1]
            if isZero(c)
                for q in k:(m-1)
                    out[q] = out[q+1]
                    fps[q] = fps[q+1]
                    gs[q] = gs[q+1]
                    ivs[q] = ivs[q+1]
                end
                m -= 1
            else
                out[k] = CoeffPair(c, t)   # same fingerprint as before, fps[k] still valid
                gs[k] = g
                ivs[k] = iv
            end
            found = true;
            break
        end

        found || (m += 1; out[m]=s; fps[m]=fp; gs[m]=g; ivs[m]=iv)
    end

    m == 0 && return zero()
    m < n && resize!(out, m)
    TEx.summands = out
    # multi-term symmetries relate whole sums rather than pairs of summands, so the projector runs once over the merged result rather than inside the loop above
    multiterm && return multiterm_reduce(TEx; budget=multiterm_budget)
    return TEx
end
##############################################################################################
# Canonical collection: merging summands that differ only by dummy names and factor order
##############################################################################################

"""
Number of distinct values in `col`, using `buf` as scratch so the count allocates nothing.
"""
function _ndistinct!(buf::Vector{UInt}, col::Vector{UInt})
    copyto!(buf, col)
    sort!(buf)
    c = 0
    @inbounds for i in eachindex(buf)
        (i == 1 || buf[i] != buf[i-1]) && (c += 1)
    end
    return c
end

"""
The slots of a monomial's factors, flattened: factor `i` holds slots `off[i]+1:off[i+1]`, and slot `k`
carries index `uid[k]` on factor `fac[k]`.
"""
function _flat_slots(ind::Vector{Vector{IndexPair}})
    n = length(ind)
    off = Vector{Int}(undef, n + 1)
    off[1] = 0
    @inbounds for i in 1:n
        off[i+1] = off[i] + length(ind[i])
    end
    uid = Vector{UInt64}(undef, off[n+1])
    fac = Vector{Int}(undef, off[n+1])
    @inbounds for i in 1:n, s in eachindex(ind[i])
        uid[off[i]+s] = ind[i][s][1].uid
        fac[off[i]+s] = i
    end
    return off, uid, fac
end

"""
Colour each factor of a monomial by its contraction pattern, so that two terms differing only by
dummy names and factor order get the same colours.

The slot table is flattened into parallel arrays and searched by linear scan. A monomial carries a
few dozen slots at most, and at that size a scan beats a hash map.
"""
function _contraction_colours(T::TensorTerm, cls=_slot_classes(T), orient=nothing,
                              ind=Vector{IndexPair}[indices(t) for t in terms(T)], flat=_flat_slots(ind))
    ts = terms(T)
    n = length(ts)
    off, uid, fac = flat
    S = off[n+1]

    scl = Vector{Int}(undef, S)                 # symmetry class of the slot
    psn = Vector{Bool}(undef, S)
    live = Vector{Bool}(undef, S)               # lightcone indices take no part
    @inbounds for i in 1:n, s in eachindex(ind[i])
        k = off[i] + s
        idxo, poso = ind[i][s]
        scl[k] = _slot_class(cls[i], s)
        #a pair the signature reorients is coloured without its levels, which it may exchange
        psn[k] = orient !== nothing && idxo.uid in orient ? false : poso.position
        live[k] = !isLightcone(idxo)
    end

    # occurrences of each index, and the other end when there are exactly two
    cnt = zeros(Int, S)
    partner = zeros(Int, S)
    @inbounds for k in 1:S
        live[k] || continue
        c = 0
        other = 0
        for j in 1:S
            (live[j] && uid[j] == uid[k]) || continue
            c += 1
            (j != k && other == 0) && (other = j)
        end
        cnt[k] = c
        partner[k] = other
    end
    #a dummy slot; the factor at its other end, 0 for none; and the label of that edge, the classes of
    #both slots, which refinement mixes with the colour it finds there
    dum = Vector{Bool}(undef, S)
    nbr = zeros(Int, S)
    edge = Vector{UInt}(undef, S)
    @inbounds for k in 1:S
        dum[k] = live[k] && cnt[k] == 2
        dum[k] || continue
        nbr[k] = fac[partner[k]]
        edge[k] = hash(scl[k], hash(scl[partner[k]], hash(psn[k], UInt(0))))
    end

    #slots are told apart by symmetry class, not position, and summed, so that a colour does not depend
    #on the order of slots a symmetry exchanges. Each slot's part is mixed by an outer hash before it is summed
    col = Vector{UInt}(undef, n)
    @inbounds for i in 1:n
        h = hash(name(ts[i]), hash(length(ind[i])))::UInt     # pin type for performance
        for s in eachindex(ind[i])
            k = off[i] + s
            if !dum[k]
                h += hash(hash(uid[k], hash(psn[k], hash(scl[k]))))     # free: the name matters
            else
                h += hash(hash(psn[k], hash(scl[k])))                   # dummy: the name does not
            end
        end
        col[i] = h
    end

    newcol = Vector{UInt}(undef, n)
    buf = Vector{UInt}(undef, n)
    _refine_colours!(col, newcol, buf, off, nbr, edge)
    #tie breaking
    joined = falses(n)
    @inbounds for k in 1:S
        nbr[k] != 0 && (joined[fac[k]] = true)
    end
    @inbounds while true
        v, cv = 0, typemax(UInt)
        for i in 1:n
            joined[i] || continue
            for j in 1:n
                (j != i && joined[j] && col[j] == col[i]) || continue
                (col[i] < cv || (col[i] == cv && i < v)) && (v = i; cv = col[i])
                break
            end
        end
        v == 0 && break
        if (prs = _closed_pairs(cv, col, joined, off, nbr)) !== nothing
            numbered = length(prs) > 1 && prs[1][2] != 0
            for k in eachindex(prs)
                p, q = prs[k]
                numbered && (col[p] = col[q] = hash(cv, UInt(k)))
                joined[p] = false
                q != 0 && (joined[q] = false)
            end
            continue
        end
        col[v] = hash(col[v], 0x9e3779b97f4a7c15)
        _refine_colours!(col, newcol, buf, off, nbr, edge)
    end
    return col
end

function _closed_pairs(cv, col, joined, off, nbr)
    n = length(col)
    partner = zeros(Int, n)
    lone = paired = false
    @inbounds for i in 1:n
        (joined[i] && col[i] == cv) || continue
        j = 0
        for k in off[i]+1:off[i+1]
            x = nbr[k]
            (x == 0 || x == i) && continue
            j == 0 ? (j = x) : (x == j || return nothing)
        end
        if j == 0
            lone = true
        else
            (joined[j] && col[j] == cv) || return nothing
            paired = true
        end
        partner[i] = j == 0 ? i : j
        lone && paired && return nothing
    end
    prs = Tuple{Int,Int}[]
    @inbounds for i in 1:n
        j = partner[i]
        j == 0 && continue
        if j == i
            push!(prs, (i, 0))
        else
            partner[j] == i || return nothing
            i < j && push!(prs, (i, j))
        end
    end
    return prs
end

"""
Rounds of Weisfeiler-Leman refinement, shared by the canonical key's factor colours and the isomorphism
graph's `wl_colours`: each vertex's colour absorbs those of its neighbours. The neighbours of vertex `i`
are `nbr[off[i]+1:off[i+1]]`, 0 standing for none, each reached along an edge labelled `lab`. Stops when
the partition is stable or discrete.
"""
function _refine_colours!(col, newcol, buf, off, nbr, lab)
    n = length(col)
    ndistinct = _ndistinct!(buf, col)
    for _ in 1:n
        ndistinct == n && break                 # every vertex apart already, nothing left to split
        @inbounds for i in 1:n
            acc = UInt(0)
            for k in off[i]+1:off[i+1]
                nbr[k] == 0 && continue
               
                acc += hash(hash(col[nbr[k]], lab[k]))
            end
            newcol[i] = hash(acc, col[i])
        end
        copyto!(col, newcol)
        nd = _ndistinct!(buf, col)
        nd == ndistinct && break                # partition stable, further rounds are noise
        ndistinct = nd
    end
    return col
end

"""
Sign from permuting graded factors: `(-1)^(odd-graded inversions)`.

Only odd-graded factors can contribute, so those are gathered first and the inversion count runs
over just them. The naive form is O(n²) over *all* factor pairs with a `grading` call inside the inner loop, and for an all-bosonic term every one of
those iterations is wasted.
"""
function _graded_perm_sign(ts, perm)
    odd = Int[]
    @inbounds for i in eachindex(ts)
        grading(ts[i]) == 1 && push!(odd, i)
    end
    length(odd) < 2 && return 1

    newpos = Vector{Int}(undef, length(perm))       # factor index => position after permuting
    @inbounds for k in eachindex(perm)
        newpos[perm[k]] = k
    end
    s = 0
    @inbounds for a in eachindex(odd), b in (a+1):length(odd)
        newpos[odd[a]] > newpos[odd[b]] && (s += 1)
    end
    return iseven(s) ? 1 : -1
end

"""
Can the factors of this term be freely permuted?

False as soon as any pair is subject to a commutation relation - matrices with a registered
commutator do not simply reorder, they generate extra terms, which is `bubble_sort_with_grading!`'s
job rather than ours.
"""
function _freely_reorderable(ts)
    nmat = 0
    @inbounds for t in ts
        t isa MatrixSuperType && (nmat += 1)
    end
    nmat <= 1 && return true                 # a commutator needs two matrices
    @inbounds for i in eachindex(ts), j in (i+1):length(ts)
        can_swap(ts[i], ts[j]) || return false
    end
    return true
end

"""
    canonical_form(T::TensorTerm) -> (sign, term)

Normalise `T` by reordering commuting factors and renaming dummy indices, such that
`sign * term == T`.
"""
function canonical_form(T::TensorTerm)
    ts = terms(T)
    n = length(ts)
    n <= 1 && return (1, rename_dummies(T))
    _freely_reorderable(ts) || return (1, rename_dummies(T))
    col = _contraction_colours(T)
    perm = collect(1:n)
    sort!(perm, by=i -> (col[i], i))       # ties keep source order, so this is deterministic
    if perm == 1:n
        return (1, rename_dummies(T))
    end
    sgn = _graded_perm_sign(ts, perm)
    reordered = TensorTerm(TensorSuperType[ts[p] for p in perm])
    return (sgn, rename_dummies(reordered))
end

"""
The exchange sign of any two slots of `t`: `-1` totally antisymmetric, `+1` totally symmetric, and
`0` for anything else. Mixed (Young) symmetry has no single exchange sign, so those factors carry
no symmetry step.
"""
function _exchange_sign(t::TensorSuperType)
    length(indices(t)) < 2 && return 0
    isAntisymmetric(t) && return -1
    isSymmetric(t) && return 1
    return 0
end

"""
Whether two factors are copies of the same tensor.

The `TensorData` uid screens pairs cheaply, but a match still falls through to the full structural
comparison. Non-matching pairs (almost all of them) exit on the uid, so the exact comparison runs
only for pairs that were going to be examined anyway.

Operators carry their own data block, so they skip the uid and go straight to the fields.
"""
function _same_tensor_kind(a::TensorSuperType, b::TensorSuperType)
    typeof(a) === typeof(b) || return false
    if !(a isa OperatorSuperType) && !(b isa OperatorSuperType)
        data(a).uid == data(b).uid || return false
    end
    return name(a) == name(b) && labels(a) == labels(b) &&
           tableaux(a) == tableaux(b) && function_of(a) == function_of(b) &&
           _extra_identity(a) == _extra_identity(b)
end

"Sentinel for an index the endpoint table has not seen; count 0 is never a dummy."
const _NO_ENDS = (0, 0, 0, 0, 0)

"""
    _exchange_vanishes(tt) -> Bool

True when exchanging two copies of one tensor maps the term to minus itself, forcing it to zero.

Exchanging two rank-`r` factors is `r` transpositions at once, one per slot, so it is a symmetry of
the term only if *every* slot closes a cycle: the dummies in slot `k` of the two factors must land
on the same partner factor, in two slots that partner relates by a symmetry. Each closed slot
contributes that symmetry's sign, and the exchange costs `(-1)^grading` on top. The term vanishes
when the product is `-1`.
"""
function _exchange_vanishes(tt::TensorTerm)
    ts = terms(tt)
    n = length(ts)
    n < 2 && return false

    # Endpoints of every index, packed as (count, factor1, slot1, factor2, slot2). A dummy has
    # exactly two and a free index one, so the whole record fits in a tuple, which the Dict stores
    # inline. The count saturates at 3, meaning "more than two", which is not a dummy either way.
    nslots = 0
    @inbounds for t in ts
        nslots += length(indices(t))
    end
    ends = Dict{UInt64,NTuple{5,Int}}()
    sizehint!(ends, nslots)
    @inbounds for i in 1:n, (k, (idx, _)) in enumerate(indices(ts[i]))
        e = get(ends, idx.uid, _NO_ENDS)
        ends[idx.uid] = e[1] == 0 ? (1, i, k, 0, 0) :
                        e[1] == 1 ? (2, e[2], e[3], i, k) :
                        (3, e[2], e[3], e[4], e[5])
    end

    @inbounds for i in 1:(n-1), j in (i+1):n
        ri = length(indices(ts[i]))
        ri == length(indices(ts[j])) || continue
        ri == 0 && continue
        _same_tensor_kind(ts[i], ts[j]) || continue
        grading(ts[i]) == grading(ts[j]) || continue
        can_swap(ts[i], ts[j]) || continue
        # every factor between them must also be crossable, twice over
        blocked = false
        for m in (i+1):(j-1)
            if !can_swap(ts[i], ts[m]) || !can_swap(ts[j], ts[m])
                blocked = true;
                break
            end
        end
        blocked && continue

        total = iseven(grading(ts[i])) ? 1 : -1
        ok = true
        for k in 1:ri
            # the two factors must agree on index level slot by slot: swapping an upper index for
            # a lower one is not a relabelling, so such a pair is not exchangeable at all
            indices(ts[i])[k][2] == indices(ts[j])[k][2] || (ok=false; break)
            u = indices(ts[i])[k][1].uid
            v = indices(ts[j])[k][1].uid
            eu = get(ends, u, _NO_ENDS)
            ev = get(ends, v, _NO_ENDS)
            # both must be dummies: a free index cannot be relabelled. A count other than 2 covers
            # the absent case too, so no separate nothing check is needed
            (eu[1] != 2 || ev[1] != 2) && (ok=false; break)
            pu = (eu[2] == i && eu[3] == k) ? (eu[4], eu[5]) : (eu[2], eu[3])
            pv = (ev[2] == j && ev[3] == k) ? (ev[4], ev[5]) : (ev[2], ev[3])
            pu[1] == pv[1] || (ok=false; break)      # must land on the same partner factor
            (pu[1] == i || pu[1] == j) && (ok=false; break)   # partner is one of the pair
            part = ts[pu[1]]
            sgn = isAntisymmetric(part) ? -1 : isSymmetric(part) ? 1 : 0
            sgn == 0 && (ok=false; break)            # partner relates those slots by no symmetry
            total *= sgn
        end
        ok && total == -1 && return true
    end
    return false
end

"""
Whatever a factor is told apart by that its shared `TensorData` does not hold.

`TensorData` is interned, so two factors of the same name and type share one object. Some matricies have other identifying properties
"""
_extra_identity(::TensorSuperType) = nothing

"""
Identity of a factor for signature purposes: everything about it except its indices.

Keyed on the shared `TensorData`'s content `uid`, not on `objectid`: ~75x cheaper to compare, since `objectid` on an
immutable struct recomputes over the fields every call.
"""
_factor_identity(t::TensorSuperType) =
    t isa OperatorSuperType ? hash(name(t)) :
    hash(_extra_identity(t), hash(typeof(t), data(t).uid))

"""
Emit the slot codes for one index vector.
"""
function _sig_slots!(key::Vector{UInt64}, ti::Vector{IndexPair},
    dummyno::Dict{UInt64,UInt64}, dums::Set{UInt64}, orient=nothing, flipped=nothing, fsgn=nothing, order=nothing)
    @inbounds for j in eachindex(ti)
        idx, pos = ti[order === nothing ? j : order[j]]
        u = idx.uid
        level = UInt64(pos.position)
        if u in dums
            reorients = orient !== nothing && u in orient
            first_seen = reorients && !haskey(dummyno, u)
            #numbered on first appearance, and the set goes in with the number, since a trace over
            #one set is not a trace over another; tagged so it cannot alias a uid
            code = get(dummyno, u, UInt64(0))        #a code is never 0, carrying the tag
            if code == 0
                code = hash(idx.index_set.name, UInt64(length(dummyno) + 1)) | 0x8000_0000_0000_0000
                dummyno[u] = code
            end
            push!(key, code)
            if reorients
                #the pair is read raised end first: met lowered, both its levels are exchanged, at the
                #metric's exchange sign
                if first_seen && pos == lower
                    push!(flipped, u)
                    idx.index_set.metric_type <: AntisymMetricSuperType && (fsgn[] = -fsgn[])
                end
                u in flipped && (level = level == 0 ? UInt64(1) : UInt64(0))
            end
        else
            push!(key, u)
        end
        push!(key, level)
    end
    return key
end

"""
The dummies of a term whose pair the canonical signature may read raised end first, exchanging the
levels of both ends at the cost of the metric's exchange sign.
"""
function _orientable_dummies(ts::Vector{TensorSuperType}, dums::Set{UInt64})
    isempty(dums) && return nothing
    seen = sizehint!(UInt64[], 2 * length(dums))      #eligible occurrences; a term has few dummies
    sets = sizehint!(IndexSet[], 2 * length(dums))
    symonly = sizehint!(Bool[], 2 * length(dums))     #a derivative's own index needs a symmetric metric
    for t in ts
        _eligible_slots!(seen, sets, symonly, t, dums)     #one dispatch per factor, its type known inside
    end
    out = nothing
    lastS, sym, anti = nothing, false, false
    @inbounds for j in eachindex(seen)
        u = seen[j]
        count(==(u), seen) == 2 || continue
        k = findnext(==(u), seen, j + 1)
        k === nothing && continue                     #the pair is taken at its first end
        #both ends share a set. Subtype tests on a stored type are slow and a term has few sets, so each
        #run of one set is classified once
        S = sets[j]
        if S !== lastS
            mt = S.metric_type
            lastS, sym = S, mt !== nothing && mt <: SymMetricSuperType
            anti = !sym && mt !== nothing && mt <: AntisymMetricSuperType
        end
        (sym || (anti && !symonly[j] && !symonly[k])) || continue
        out === nothing && (out = UInt64[])
        push!(out, u)
    end
    return out
end

#the occurrences on factor `t` that `_orientable_dummies` may reorient, before their metric is checked.
#`encl` holds the derivatives whose operand `t` sits in: an occurrence there counts only if its set's
#metric passes every one of them, so that it can be carried out to where the other end is
function _eligible_slots!(seen::Vector{UInt64}, sets::Vector{IndexSet}, symonly::Vector{Bool},
                          t::TensorSuperType, dums::Set{UInt64}, encl=())
    _keeps_levels(t) && return seen
    isop = t isa OperatorSuperType
    isop && !(t isa PartialSuperType || t isa CovDerivSuperType) && return seen
    inds = indices(t)
    own = isop ? op_indices(t) : inds
    @inbounds for s in eachindex(own)
        idx = own[s][1]
        #an operator's own slots lead its indices; anything else and it is left alone
        isop && (s > length(inds) || inds[s] != own[s]) && break
        u = idx.uid
        u in dums || continue
        (isLightcone(idx) || idx isa DummyPatternIndex) && continue
        if !isempty(encl)       #the metric the set builds, as popping it would write
            g = generate_metric(idx => upper, idx => upper)
            (g !== nothing && all(D -> commutes_with_metric(D, g), encl)) || continue
        end
        push!(seen, u)
        push!(sets, idx.index_set)
        push!(symonly, isop)
    end
    #the operand's own occurrences, rather than the copies of them the derivative exposes
    if isop
        inner = (encl..., t)
        for (_, tt) in summands(operand(t)), f in terms(tt)
            _eligible_slots!(seen, sets, symonly, f, dums, inner)
        end
    end
    return seen
end

#the classes of a tensor whose one tableau, a single row or column, holds every slot in order, as a
#SymmetricTensor's does: shared between all such tensors of a rank, and only ever read
const _WHOLE_CLASS_CACHE = Dict{Tuple{Int,Int},Tuple{Vector{Int},Vector{Int}}}()

function _slot_classes(t::TensorSuperType)
    t isa OperatorSuperType && return nothing
    tabs = tableaux(t)
    isempty(tabs) && return nothing           #`nothing` overall when every slot is its own class
    n = length(indices(t))
    if length(tabs) == 1 && n >= 2
        tab1 = first(tabs)
        sl1, rl1 = tab1.index_assignment, tab1.row_lengths
        k1 = length(rl1) == 1 ? 1 : all(==(1), rl1) ? -1 : 0
        if k1 != 0 && length(sl1) == n && all(i -> @inbounds(sl1[i] == i), 1:n)
            return get!(() -> (ones(Int, n), fill(k1, n)), _WHOLE_CLASS_CACHE, (n, k1))
        end
    end
    #plain loops and concrete types: this runs per factor. Built at the first tableau that makes a class
    cls = kind = nothing
    for tab in tabs
        slots, rl = tab.index_assignment, tab.row_lengths
        length(slots) < 2 && continue
        k = length(rl) == 1 ? 1 : all(==(1), rl) ? -1 : 0
        k == 0 && continue
        cls === nothing && (cls = collect(1:n); kind = zeros(Int, n))
        free = true
        for s in slots
            (s <= n && kind[s] == 0) || (free = false; break)
        end
        free || continue
        c = minimum(slots)
        for s in slots
            cls[s] = c
            kind[s] = k
        end
    end
    return cls === nothing ? nothing : (cls, kind)
end
_slot_classes(T::TensorTerm) = Union{Nothing,Tuple{Vector{Int},Vector{Int}}}[_slot_classes(t) for t in terms(T)]

"""The symmetry class of slot `s` of a factor, from its `_slot_classes`."""
_slot_class(c, s) = c === nothing ? s : c[1][s]

"""
Append one factor's contribution to a canonical signature.

Free indices contribute their uid; dummies contribute their set and a sequence number assigned on
first appearance, so two terms differing only by a consistent dummy renaming produce identical keys.
`order` lists the slots in the order to emit them. Operators additionally contribute their operand's
structure, recursively.
"""
function _sig_factor!(key::Vector{UInt64}, t::TensorSuperType,
    dummyno::Dict{UInt64,UInt64}, dums::Set{UInt64}, order=nothing, orient=nothing, flipped=nothing, fsgn=nothing)
    push!(key, _factor_identity(t))
    _sig_slots!(key, indices(t), dummyno, dums, orient, flipped, fsgn, order)
    if t isa OperatorSuperType
        for (c, tt) in summands(operand(t))
            push!(key, hash(c))
            for u in terms(tt)
                _sig_factor!(key, u, dummyno, dums, nothing, orient, flipped, fsgn)
            end
            push!(key, 0xFFFF_FFFF_FFFF_FFFF)               # summand separator
        end
    end
    return key
end

"""
    canonical_signature(T::TensorTerm) -> (sign, key)

A key such that two terms share it exactly when they are equal up to dummy renaming and the
ordering of commuting factors, with `sign * (term with key) == T`.
"""
function canonical_signature(T::TensorTerm)
    ts = terms(T)
    n = length(ts)
    classes = _slot_classes(T)
    dums = dummy_uids(T)
    inds = Vector{IndexPair}[indices(t) for t in ts]      #once, concretely typed, for every pass below
    orient = _orientable_dummies(ts, dums)
    flipped, fsgn = orient === nothing ? (nothing, nothing) : (UInt64[], Ref(1))
    perm = collect(1:n)
    sgn = 1
    flat = nothing
    if n > 1 && _freely_reorderable(ts)
        flat = _flat_slots(inds)
        col = _contraction_colours(T, classes, orient, inds, flat)
        sort!(perm, by=i -> (col[i], i))
        perm == 1:n || (sgn = _graded_perm_sign(ts, perm))
    end
    dummyno = Dict{UInt64,UInt64}()
    key = UInt64[]
    sizehint!(key, n + 2 * sum(t -> length(indices(t)), ts; init=0))
    if all(isnothing, classes)
        @inbounds for p in perm
            _sig_factor!(key, ts[p], dummyno, dums, nothing, orient, flipped, fsgn)
        end
        return (fsgn === nothing ? sgn : sgn * fsgn[], key)
    end
    rank = invperm(perm)
    #every slot's uid and factor, flat, for the partner lookups of the sorts below
    off, fu, ff = flat === nothing ? _flat_slots(inds) : flat
    mem = Int[]                       #a class's slots, their keys, and the order sorting them
    kb = NTuple{5,UInt64}[]
    ob = Int[]
    @inbounds for p in perm
        order = nothing
        if classes[p] !== nothing
            cls, kind = classes[p]
            for c in eachindex(cls)
                (cls[c] == c && kind[c] != 0) || continue    #each class once, from its first slot
                second = findnext(==(c), cls, c + 1)
                if findnext(==(c), cls, second + 1) === nothing
                    #two slots, as in F_{μν}: one comparison, no sort
                    _slot_key(inds, fu, ff, off, classes, rank, dums, p, second, orient, dummyno) <
                        _slot_key(inds, fu, ff, off, classes, rank, dums, p, c, orient, dummyno) || continue
                    order === nothing && (order = collect(eachindex(cls)))
                    order[c], order[second] = second, c
                    kind[c] == -1 && (sgn = -sgn)
                    continue
                end
                empty!(mem); empty!(kb); empty!(ob)
                for s in eachindex(cls)
                    cls[s] == c || continue
                    push!(mem, s)
                    push!(kb, _slot_key(inds, fu, ff, off, classes, rank, dums, p, s, orient, dummyno))
                    push!(ob, length(mem))
                end
                #a stable insertion sort, as sortperm is stable, whose shifts are the inversions it undoes
                inv = 0
                for a in 2:length(ob)
                    x = ob[a]
                    b = a - 1
                    while b >= 1 && kb[ob[b]] > kb[x]
                        ob[b+1] = ob[b]
                        b -= 1
                        inv += 1
                    end
                    ob[b+1] = x
                end
                inv == 0 && continue
                order === nothing && (order = collect(eachindex(cls)))
                kind[c] == -1 && isodd(inv) && (sgn = -sgn)
                for a in eachindex(mem)
                    order[mem[a]] = mem[ob[a]]
                end
            end
        end
        _sig_factor!(key, ts[p], dummyno, dums, order, orient, flipped, fsgn)
    end
    return (fsgn === nothing ? sgn : sgn * fsgn[], key)
end

"""
How `canonical_signature` orders slot `s` of factor `p` within its symmetry class: a free index by its
name, a dummy by where its partner sits, which a renaming keeps. Two dummies running to one class of an
earlier factor are ordered by the codes `dummyno` gave them there, so both ends read them in one order;
their levels cannot, as a pair that may be reoriented is ordered without them.
"""
function _slot_key(inds::Vector{Vector{IndexPair}}, fu::Vector{UInt64}, ff::Vector{Int}, off::Vector{Int},
                   classes, rank::Vector{Int}, dums::Set{UInt64}, p::Int, s::Int, orient=nothing, dummyno=nothing)
    k = off[p] + s
    @inbounds u = fu[k]
    #a pair that may be reoriented is ordered without its levels
    ps = orient !== nothing && u in orient ? UInt64(0) : UInt64(inds[p][s][2].position)
    if u in dums
        code = dummyno === nothing ? UInt64(0) : get(dummyno, u, UInt64(0))   #0 until its pair is met
        @inbounds for j in eachindex(fu)
            (fu[j] == u && j != k) || continue
            q = ff[j]
            return (UInt64(1), UInt64(rank[q]), UInt64(_slot_class(classes[q], j - off[q])), code, ps)
        end
    end
    return (UInt64(0), u, UInt64(0), UInt64(0), ps)
end

"""
Collect summands that are equal up to dummy index naming and the ordering of commuting factors.

Groups by [`canonical_signature`](@ref) and merges on exact key equality, so no pair of terms is
ever compared structurally.
"""
function collect_terms(TEx::TensorExpression)
    s = summands(TEx)
    n = length(s)
    n <= 1 && return TEx

    out = CoeffPair[]
    sizehint!(out, n)
    keys = Vector{Vector{UInt64}}()
    signs = Int[]
    buckets = Dict{UInt64,Vector{Int}}()
    fps = UInt[sum(t -> hash(_factor_identity(t)), terms(tt); init=UInt(0)) for (_, tt) in s]
    seen = Dict{UInt,Int}()
    for f in fps
        seen[f] = get(seen, f, 0) + 1
    end

    for (i, (coef, tt)) in enumerate(s)
        isZero(coef) && continue
        if seen[fps[i]] == 1
            push!(out, CoeffPair(coef, tt))
            push!(keys, UInt64[])
            push!(signs, 1)
            continue
        end
        sgn, key = canonical_signature(tt)
        h = hash(key)
        bk = get!(() -> Int[], buckets, h)
        found = false
        for k in bk
            if keys[k] == key
                # both sit over the same canonical term; re-express this one against that rep
                out[k] = CoeffPair(out[k][1] + coef * sgn * signs[k], out[k][2])
                found = true
                break
            end
        end
        if !found
            push!(out, CoeffPair(coef, tt))
            push!(keys, key)
            push!(signs, sgn)
            push!(bk, length(out))
        end
    end

    filter!(p -> !isZero(p[1]), out)
    isempty(out) && return zero()
    return TensorExpression(_RAW, out)
end

"""Two summands, both coefficient 1: merge if the terms are equal. Cannot produce a zero coefficient."""
@inline function _merge2(v::Vector{CoeffPair})
    @inbounds if v[1][2] == v[2][2]
        return TensorExpression(_RAW, CoeffPair[CoeffPair(v[1][1] + v[2][1], v[1][2])])
    end
    return TensorExpression(_RAW, v)
end

"""Insert one (coef, term) into an already-merged summand vector. Mutates `v`."""
function _insert_merged(v::Vector{CoeffPair}, c, t::TensorTerm)
    @inbounds for i in eachindex(v)
        if v[i][2] == t
            nc = v[i][1] + c
            if isZero(nc)
                deleteat!(v, i)
                isempty(v) && return zero()
            else
                v[i] = CoeffPair(nc, t)
            end
            return TensorExpression(_RAW, v)
        end
    end
    push!(v, CoeffPair(c, t))
    return TensorExpression(_RAW, v)
end

const _EMPTY_IDX = Int[]
"""Above this many summands `merge_sum` hashes terms instead of comparing all pairs."""
const _MERGE_HASH_CUTOFF = 16

function merge_sum(TEx::TensorExpression)
    src = TEx.summands
    n = length(src)
    n == 0 && return zero()
    if n == 1
        @inbounds return isZero(src[1][1]) ? zero() : TEx
    end

    # detection
    firstbad = 0
    buckets = nothing
    if n <= _MERGE_HASH_CUTOFF
        @inbounds for i in 1:n
            isZero(src[i][1]) && (firstbad=i; break)
            ti = src[i][2]
            for j in 1:(i-1)
                src[j][2] == ti && (firstbad=i; break)
            end
            firstbad != 0 && break
        end
    else
        b = Dict{UInt,Vector{Int}}()
        sizehint!(b, n)
        @inbounds for i in 1:n
            if isZero(src[i][1])
                firstbad = i
                break
            end
            ti = src[i][2]
            bk = get!(() -> Int[], b, hash(ti))
            dup = false
            for j in bk
                if src[j][2] == ti
                    dup = true
                    break
                end
            end
            if dup
                firstbad = i
                break
            end
            push!(bk, i)          # only reached for a clean, distinct summand
        end
        buckets = b
    end
    firstbad == 0 && return TEx          # nothing to do
    out = Vector{CoeffPair}(undef, n)
    m = firstbad - 1
    @inbounds for i in 1:m
        out[i] = src[i]                  # clean prefix, already known distinct & nonzero
    end

    if n <= _MERGE_HASH_CUTOFF
        @inbounds for i in firstbad:n
            s = src[i]
            isZero(s[1]) && continue
            t = s[2]
            found = false
            for k in 1:m
                if out[k][2] == t
                    c = out[k][1] + s[1]
                    if isZero(c)
                        for q in k:(m-1)
                            ;
                            out[q] = out[q+1];
                        end
                        m -= 1
                    else
                        out[k] = CoeffPair(c, t)
                    end
                    found = true;
                    break
                end
            end
            found || (m += 1; out[m]=s)
        end
    else
        if buckets === nothing
            buckets = Dict{UInt,Vector{Int}}()
            @inbounds for k in 1:m
                push!(get!(() -> Int[], buckets, hash(out[k][2])), k)
            end
        end
        ncancelled = 0
        @inbounds for i in firstbad:n
            s = src[i]
            isZero(s[1]) && continue
            t = s[2];
            h = hash(t)
            found = false
            for k in get(buckets, h, _EMPTY_IDX)
                if out[k][2] == t
                    c = out[k][1] + s[1]
                    out[k] = CoeffPair(c, t)
                    isZero(c) && (ncancelled += 1)
                    found = true;
                    break
                end
            end
            if !found
                m += 1;
                out[m] = s
                push!(get!(() -> Int[], buckets, h), m)
            end
        end
        if ncancelled > 0                # in-place compaction, no closure, no filter!
            w = 0
            @inbounds for k in 1:m
                isZero(out[k][1]) && continue
                w += 1;
                out[w] = out[k]
            end
            m = w
        end
    end

    m == 0 && return zero()
    m < n && resize!(out, m)
    TEx.summands = out
    return TEx
end

"""
Keying for symmetry equivs
"""
function _sym_bucket_key(T::TensorSuperType)
    T isa OperatorSuperType && return UInt(0)
    h = hash(name(T))
    h = hash(tableaux(T), h)
    inds = indices(T)
    h = hash(length(inds), h)
    x = UInt(0)
    @inbounds for i in eachindex(inds)
        _is_dummy_at(inds, i) && continue
        x ⊻= hash(inds[i][1].uid, hash(inds[i][2].position))
    end
    return hash(x, h)
end

function merge_sum_syms(TEx::TensorExpression)
    out = filter!(x -> !isZero(x[1]), TEx.summands) #removes 0 coefficnets
    if isempty(out)
        return zero()
    end
    TEx.summands = out
    new_merged = CoeffPair[]
    exact = Dict{UInt,Vector{Int}}()
    symb = Dict{UInt,Vector{Int}}()

    for summand in TEx.summands
        t = summand[2]
        found = false

        for i in get(exact, hash(t), _EMPTY_IDX)
            if t == new_merged[i][2]
                new_merged[i] = (new_merged[i][1] + summand[1] => t)
                found = true
                break
            end
        end

        if !found && length(t.terms) == 1
            key = _sym_bucket_key(t.terms[1])
            for i in get(symb, key, _EMPTY_IDX)
                length(new_merged[i][2].terms) == 1 || continue
                rescomp = same_tensor_mod_syms(t.terms[1], new_merged[i][2].terms[1])
                if rescomp[1] == true
                    new_merged[i] = (rescomp[2] * new_merged[i][1] + summand[1] => t)
                    found = true
                    break
                end
            end
        end

        if !found
            push!(new_merged, summand)
            k = length(new_merged)
            push!(get!(() -> Int[], exact, hash(t)), k)
            length(t.terms) == 1 && push!(get!(() -> Int[], symb, _sym_bucket_key(t.terms[1])), k)
        end
    end

    TEx.summands = new_merged
    TEx.summands = filter!(x -> !isZero(x[1]), TEx.summands)
    return isempty(TEx.summands) ? zero() : TEx
end


##############################################################################################################################
##############################################################################################################################

# printing and display

function Base.show(io::IO, T::TensorSuperType)
    print(io, tens2str(T))
end
#for future use with bras and kets to hijack the rendering
function show_factor(io::IO, T::TensorSuperType, pos::Int64, TT::TensorTerm)
    show(io, T)
end

function Base.show(io::IO, T::TensorTerm)
    for tens in 1:length(terms(T))
        show_factor(io, terms(T)[tens], tens, T)
    end
end
_unwrapnum(x) = SymbolicUtils.isconst(x) ? SymbolicUtils.unwrap_const(x) : x

"""
A coefficient split into its numeric factor and its symbolic rest, the rest `nothing` when there is
none: `2q*r` gives `(2, q*r)`, `q` gives `(1, q)`.
"""
function _num_sym(x)
    x isa Number && return (x, nothing)
    SymbolicUtils.ismul(x) || return (1, x)
    args = _unwrapnum.(SymbolicUtils.arguments(x))
    syms = filter(a -> !(a isa Number), args)
    return (prod(filter(a -> a isa Number, args); init=1),
            isempty(syms) ? nothing : length(syms) == 1 ? syms[1] : *(syms...))
end

"""The terms of `x` as a sum: its arguments if it is one, none if it is zero, else `x` alone."""
_addends(x) = SymbolicUtils.isadd(x) ? _unwrapnum.(SymbolicUtils.arguments(x)) : isequal(x, 0) ? () : (x,)

"""The terms of a sum, already rendered, joined with each one's sign as the operator between."""
function _join_terms(parts)
    out = ""
    for (k, s) in enumerate(parts)
        neg = startswith(s, "-")
        out *= k == 1 ? s : (neg ? " - " * s[nextind(s, 1):end] : " + " * s)
    end
    return out
end

#a rational whose denominator is 1 is an integer, and printed as one
_real_string(x::Real) = x isa Rational && isone(denominator(x)) ? string(numerator(x)) : string(x)

"""
`n` times the symbolic `sym` for display, written by juxtaposition, `2q`, `-iq`, with a sum bracketed:
`2(q + r)`, which unbracketed read as `2q + r`.
"""
function _scaled_string(n, sym, T)
    sym === nothing && return coef_string(n, T; drop_unit=false)
    p = isone(n) ? "" : isequal(n, -1) ? "-" : coef_string(n, T; drop_unit=false)
    return p * (SymbolicUtils.isadd(sym) ? coef_string(sym, T; drop_unit=false) : string(sym))
end

function coef_string(c, T::TensorExpression; drop_unit::Bool=true)
    #unwrapped, since a `Num` is a `Number` and would be printed whole, unbracketed, as one
    c = _unwrapnum(Symbolics.unwrap(Symbolics.simplify(Symbolics.expand(c))))
    #a unit coefficient on a tensor is not written, but inside a sum it is a term like any other
    if drop_unit && !is_const(T)
        isequal(c, 1) && return ""
        isequal(c, -1) && return "-"
    end
    if c isa Number
        isreal(c) && return _real_string(real(c))
        re, ip = reim(c)
        mag = isone(abs(ip)) ? "i" : _real_string(abs(ip)) * "i"
        iszero(re) && return (ip < 0 ? "-" : "") * mag
        return "(" * _real_string(re) * (ip < 0 ? " - " : " + ") * mag * ")"
    end
    try
        args = _unwrapnum.(SymbolicUtils.arguments(c))
        op = SymbolicUtils.operation(c)
        op === (+) && return "(" * _join_terms([coef_string(a, T; drop_unit=false) for a in args]) * ")"
        op === (*) && return _scaled_string(_num_sym(c)..., T)
        if op === complex                     #a symbolic real and imaginary part: re + i·im
            re, ip = args
            n, sym = _num_sym(ip)
            parts = String[coef_string(a, T; drop_unit=false) for a in _addends(re)]
            push!(parts, _scaled_string(n * im, sym, T))
            return length(parts) == 1 ? parts[1] : "(" * _join_terms(parts) * ")"
        end
        if op === (/)                         #bracketed whole, or `1 / qA` reads as 1/(qA)
            n, d = args
            ns = coef_string(n, T; drop_unit=false)
            neg = startswith(ns, "-")
            neg && (ns = ns[nextind(ns, 1):end])
            ds = SymbolicUtils.iscall(d) ? "(" * string(d) * ")" : string(d)
            return (neg ? "-" : "") * "(" * ns * " / " * ds * ")"
        end
    catch
    end
    return string(c)
end
function Base.show(io::IO, T::TensorExpression)
    summandss = summands(T)
    n = length(summandss)
    isempty(summandss) && return print(io, "0")
    cs = coef_string(summandss[1][1], T)
    for (i, tens) in enumerate(summandss)
        print(io, cs)
        print(io, tens[2])
        if i < n
            cs = coef_string(summandss[i+1][1], T)
            if !startswith(cs, "-")
                print(io, "+")
            end
        end
    end
end


##############################################################################################################################
##############################################################################################################################

# sorting


Base.sort(T::TensorSuperType) = sort(TensorExpression(T))
function Base.sort(TEx::TensorExpression)::TensorExpression
    isZeroTens(TEx) && return zero()
    out = CoeffPair[]
    sizehint!(out, length(TEx.summands))
    for summand in TEx.summands
        coeff = summand[1]
        (fac, new_terms) = _sort_factor_indices(summand[2])
        fac == 0 && continue
        (fac2, expout) = bubble_sort_with_grading!(new_terms)
        append!(out, summands((coeff * fac) * expout))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end

"""
The term with each factor's indices in canonical order, and the sign that cost
"""
function _sort_factor_indices(tt::TensorTerm)
    ts = terms(tt)
    out = Vector{TensorSuperType}(undef, length(ts))
    fac = 1
    @inbounds for k in eachindex(ts)
        t = ts[k]
        (s, inds) = find_canonical_rep(t)
        s == 0 && return (0, tt)            #the factor vanishes: nothing to build
        if t isa OperatorSuperType || s != 1 || inds != indices(t)
            t = copy(t)                     #then as sort_symmetric!, on the copy
            set_indices!(t, inds)
        end
        fac *= s
        out[k] = t
    end
    return (fac, TensorTerm(out))
end


function Base.isless(a::TensorSuperType, b::TensorSuperType)
    if name(a) < name(b)
        return true
    elseif name(a) == name(b)
        return indices(a) < indices(b)
    end
    return false
end


function Base.sort(term::TensorTerm)::TensorExpression
    (fac, new_terms) = _sort_factor_indices(term)
    fac == 0 && return zero()
    (fac2, expout) = bubble_sort_with_grading!(new_terms)
    return fac * expout
end


"""
Dummy indices of a term as a set of uids, for O(1) membership during sorting.

`bubble_sort_with_grading!` consults this once per comparison, so a `Vector` here costs a linear scan every time
"""
function dummy_uids(T::Tensorial)
    ds = dummies(T)
    s = sizehint!(Set{UInt64}(), length(ds))      #sized once, not rehashed as it fills
    for d in ds
        push!(s, d.uid)
    end
    return s
end

@inline _is_dummy_uid(dummies::Set{UInt64}, idx) = idx.uid in dummies
@inline _is_dummy_uid(dummies, idx) = idx in dummies   # legacy: any container of indices

function is_less_mod_dummies(T1::TensorSuperType, T2::TensorSuperType, dummies)
    if dummies === nothing
        return T1 < T2
    elseif length(dummies) == length(indices(T1)) && length(indices(T1)) == length(indices(T2))
        i1 = indices(T1)
        i2 = indices(T2)

        @inbounds for i = eachindex(i1)
            a, ap = i1[i]
            b, bp = i2[i]
            adef = ap == a.index_set.default_position
            bdef = bp == b.index_set.default_position
            if adef && !bdef
                return true
            elseif !adef && bdef
                return false
            elseif a == b
                continue
            else
                return a < b
            end
        end
        return false
    else
        i1 = indices(T1)
        i2 = indices(T2)
        n = min(length(i1), length(i2))#may not have same lengths!
        @inbounds for i = 1:n
            a, ap = i1[i]
            b, bp = i2[i]
            if _is_dummy_uid(dummies, a)
                if _is_dummy_uid(dummies, b)
                    continue
                else
                    return false

                end
            elseif _is_dummy_uid(dummies, b)
                return true
            else
                adef = ap == a.index_set.default_position
                bdef = bp == b.index_set.default_position
                if adef && !bdef
                    return true
                elseif !adef && bdef
                    return false
                elseif a == b
                    continue
                else
                    return a < b
                end
            end
        end
        # All common entries were equal or dummies
        return length(i1) < length(i2)
    end
end


function custom_order_isless(T1::TensorSuperType, T2::TensorSuperType, dummies=nothing)
    # tensors cloned from a common source share their (immutable) data, so this settles name
    # equality with a pointer comparison instead of hashing and comparing strings
    data(T1) === data(T2) && return is_less_mod_dummies(T1, T2, dummies)
    if !isempty(sort_order) && haskey(sort_order, name(T1))
        if haskey(sort_order, name(T2))
            if sort_order[name(T1)] < sort_order[name(T2)]
                return true
            elseif name(T1) == name(T2)
                return is_less_mod_dummies(T1, T2, dummies) ###################finish fixing
            else
                return false
            end
        end
    end

    if name(T1) == name(T2)
        return is_less_mod_dummies(T1, T2, dummies) ###################finish fixing
    else
        return name(T1) < name(T2)
    end
end

"""
Whether two factors may be exchanged for free, contributing at most a grading sign and never an
extra term.

Anything that is not a pair of matrices reorders freely. Two matrices do so only when they are
known explicitly to commute. See `has_commutation_rel` for the third case, a registered
relation that does not vanish, where the exchange is allowed but costs a remainder.
"""
function can_swap(a::TensorSuperType, b::TensorSuperType)
    (a isa MatrixSuperType && b isa MatrixSuperType) || return true
    c = get(commutators, minmax(_ckey(a), _ckey(b)), nothing)
    c === nothing && return false
    return c.grading == 1 && isZeroTens(c.C)
end

"""
Whether a commutation relation has been registered for `a` and `b`. When `can_swap` is false this
separates a pair that may be exchanged at the cost of a remainder from one that may not be
exchanged at all.
"""
function has_commutation_rel(a::TensorSuperType, b::TensorSuperType)
    (a isa MatrixSuperType && b isa MatrixSuperType) || return false
    return haskey(commutators, minmax(_ckey(a), _ckey(b)))
end

"""
True when no two factors carry a registered commutation relation, so reordering them can only ever
contribute a sign and never an extra term.
"""
function _all_freely_commuting(v)
    mats = Int[]
    @inbounds for i in eachindex(v)
        v[i] isa MatrixSuperType && push!(mats, i)
    end
    length(mats) < 2 && return true
    @inbounds for a in eachindex(mats), b in (a+1):length(mats)
        can_swap(v[mats[a]], v[mats[b]]) || return false
    end
    return true
end

function bubble_sort_with_grading!(T::TensorTerm)::Tuple{Int,TensorExpression}
    dums = dummy_uids(T)      # uid set: the comparator hits this on every comparison
    v = terms(T)
    n = length(v)

    # Adjacent-swap sorting is only forced when a swap can emit a commutator remainder, since an
    # algorithm that moves factors over long distances cannot account for the terms generated at
    # each intermediate transposition. 
    if n > 1 && _all_freely_commuting(v)
        perm = sortperm(1:n; lt=(i, j) -> custom_order_isless(v[i], v[j], dums), alg=MergeSort)
        sgn = _graded_perm_sign(v, perm)
        v .= v[perm]
        update_indices!(T)
        return (sgn, sgn * T)
    end

    swapped = true
    swaps = 0
    out_expression = zero()
    while swapped
        swapped = false
        for i in 1:(n-1)
            # custom_order_isless already falls back to alphabetical when no custom order is registered
            custom_order_isless(v[i+1], v[i], dums) || continue
            if can_swap(v[i], v[i+1])
                v[i], v[i+1] = v[i+1], v[i]
                swaps += grading(v[i]) * grading(v[i+1])
                swapped = true
            elseif !has_commutation_rel(v[i], v[i+1])
                # two matrices with no relation registered do not commute, so this pair is simply
                # left as written. The sort is a partial one: everything else still moves.
                continue
            else #need to apply the commutator and generate a new term
                first_half = v[1:(i-1)] #should be same
                final_half = v[(i+2):end]
                (grade, comm) = get_commutator(v[i], v[i+1])#should get the commutator

                out_expression += (-1)^(swaps)*sort(TensorTerm(first_half) * comm * TensorTerm(final_half))#recursively sort the rest
                v[i], v[i+1] = v[i+1], v[i]
                swaps += grading(v[i]) * grading(v[i+1]) + (grade < 0) #extra grade from operators
                swapped = true
            end
        end
    end
    update_indices!(T)
    # @show (-1)^(swaps)* T + out_expression
    return ((-1)^(swaps), (-1)^(swaps) * T + out_expression)
end

function bubble_sort_index_with_grading!(T::TensorSuperType, mapping::Dict)
    n = length(mapping)
    # println(n)
    # println(mapping)
    swapped = true
    swaps = 1
    while swapped
        swapped = false
        for i in 1:(n-1)
            # println(i)
            if mapping[i+1] < mapping[i]
                #compute grading
                for tab in tableaux(T)
                    #adj_matrix is indexed by slot but reaches only the highest slot the tableau
                    #covers, so a tableau on part of a tensor is smaller than it.
                    i + 1 <= size(tab.adj_matrix, 1) || continue
                    if tab.adj_matrix[i, i+1] != 0
                        mapping[i], mapping[i+1] = mapping[i+1], mapping[i]
                        swaps *= tab.adj_matrix[i, i+1]
                        swapped = true
                        break
                    end
                end
            end
        end
    end
    return swaps, mapping

end


"""
LaTeX for one factor of a product, given where it sits.

The counterpart of `show_factor`: a factor whose markup depends on its neighbours -
a bra followed by a ket, which together are one `\\braket` - needs to see the term
it belongs to. Everything else ignores the position and defers to `to_latex`.
"""
to_latex_factor(T, pos::Int64, TT::TensorTerm) = to_latex(T)

function to_latex(T::TensorTerm)
    out = ""
    fs = terms(T)
    for tens in eachindex(fs)
        out *= to_latex_factor(fs[tens], tens, T)
    end
    return Base.replace(out * " ", "\$" => "")
end

function needs_parentheses(c)
    return (length(c)==0) ? false : c[1]=='('
end

function latex_coeff(c)
    s = latexify(c)
    # @show needs_parentheses(c)
    return needs_parentheses(c) ? "\\left(" * s * "\\right)" : s
end

"""
A coefficient split into its sign and the LaTeX of its magnitude.

The sign has to be pulled out before latexifying for fractions.
"""
function _latex_signed_coeff(cs::AbstractString)
    neg = startswith(cs, "-")
    mag = neg ? String(SubString(cs, 2)) : String(cs)
    return neg, (isempty(mag) ? "" : latex_coeff(mag))
end

"""
LaTeX for a value on its own, with no delimiters. `env=:raw` matters: latexify wraps a number in
`\$…\$` but a Symbolics expression in a whole `equation` environment.
"""
_tex(x) = Base.replace(string(latexify(x; env=:raw)), "\$" => "")

"""
`x` times the LaTeX `extra`, with anything multiplying a fraction written into its numerator as it
would be by hand: `\\frac{2i}{3}` rather than `\\frac{2}{3} i`.
"""
function _latex_times(x::Number, extra::AbstractString)
    if x isa Rational && denominator(x) != 1
        n = numerator(x)
        top = isempty(extra) ? string(n) : (isone(n) ? extra : string(n) * extra)
        return "\\frac{" * top * "}{" * string(denominator(x)) * "}"
    end
    x isa Rational && (x = numerator(x))
    isempty(extra) && return _tex(x)
    return isone(x) ? extra : _tex(x) * extra
end

"""A numeric coefficient times the LaTeX `extra` of its symbolic factors, split from its sign."""
function _latex_number(c::Number, extra::AbstractString="")
    isreal(c) && return (real(c) < 0, _latex_times(abs(real(c)), extra))
    re, ip = real(c), imag(c)
    unit = isempty(extra) ? "i" : "i " * extra
    iszero(re) && return (ip < 0, _latex_times(abs(ip), unit))
    body = "\\left(" * _join_terms(Tuple{Bool,String}[_latex_number(re), (ip < 0, _latex_times(abs(ip), "i"))]) * "\\right)"
    return (false, isempty(extra) ? body : body * " " * extra)
end

#signed LaTeX terms, (negative, magnitude), joined the same way
_join_terms(terms::Vector{Tuple{Bool,String}}) = _join_terms([(neg ? "-" : "") * tex for (neg, tex) in terms])

"""
LaTeX for `n` times the symbolic `sym`, split from its sign. A sum is bracketed unless it lands alone
in a numerator, as in `\\frac{q + r}{2}`: `2q + r` and `i q + r` would misread it.
"""
function _latex_scaled(n, sym, T)
    sym === nothing && return _latex_number(n)
    symtex = _tex(sym)
    if SymbolicUtils.isadd(sym) && !(n isa Rational && denominator(n) != 1)
        symtex = "\\left(" * symtex * "\\right)"
    end
    return _latex_number(n, symtex)
end

"""
A coefficient as LaTeX, split into sign and magnitude, built from its value.

`coef_string` is for display, and attaches factors by juxtaposition - `2//3i`, `1//2m` - which reads
naturally but which Julia's parser binds tighter than `//`. Handing that string to latexify therefore
computed `2//(3i)` and printed `\\frac{2}{3i}`. Working from the value keeps the numeric and symbolic
parts apart until the LaTeX is assembled. Shapes not recognised here fall back to the display string.
"""
function _latex_coef(c, T::TensorExpression; drop_unit::Bool=true)
    c = _unwrapnum(Symbolics.unwrap(Symbolics.simplify(Symbolics.expand(c))))
    #a unit coefficient on a tensor is not written, but inside a sum it is a term like any other
    if drop_unit && !is_const(T)
        isequal(c, 1) && return (false, "")
        isequal(c, -1) && return (true, "")
    end
    c isa Number && return _latex_number(c)
    try
        op = SymbolicUtils.operation(c)
        args = _unwrapnum.(SymbolicUtils.arguments(c))
        if op === (+)
            return (false, "\\left(" * _join_terms(Tuple{Bool,String}[_latex_coef(a, T; drop_unit=false) for a in args]) * "\\right)")
        elseif op === (*)
            return _latex_scaled(_num_sym(c)..., T)
        elseif op === complex                 #a symbolic real and imaginary part: re + i·im
            re, ip = args
            n, sym = _num_sym(ip)
            parts = Tuple{Bool,String}[_latex_coef(a, T; drop_unit=false) for a in _addends(re)]
            push!(parts, _latex_scaled(n * im, sym, T))
            length(parts) == 1 && return parts[1]
            return (false, "\\left(" * _join_terms(parts) * "\\right)")
        elseif op === (/)                     #the sign in front, not inside the numerator
            n, d = args
            neg, ntex = SymbolicUtils.isadd(n) ? (false, _tex(n)) : _latex_coef(n, T; drop_unit=false)
            dneg, dtex = SymbolicUtils.isadd(d) ? (false, _tex(d)) : _latex_coef(d, T; drop_unit=false)
            return (neg ⊻ dneg, "\\frac{" * ntex * "}{" * dtex * "}")
        end
    catch
    end
    return _latex_signed_coeff(coef_string(c, T))
end

function to_latex(T::TensorExpression)
    summandss = summands(T)
    isempty(summandss) && return "0"
    out = ""
    for (i, tens) in enumerate(summandss)
        neg, tex = _latex_coef(tens[1], T)
        if i == 1
            neg && (out *= "-")
        else
            out *= neg ? "-" : "+"
        end
        out *= tex * " "
        out *= to_latex(tens[2])
    end
    return Base.replace(out * " ", "\$" => "")
end

"""
A name as LaTeX.
"""
function uni_to_TeX(charr::String)
    isempty(charr) && return ""
    mark, base = indexmark(charr)
    mark === nothing && return Base.replace(latexify(charr), "#\"L" => "")
    cmd = mark == COMBINING_MACRON ? "overline" : _TEX_MARK_TO_SUFFIX[mark]
    return "\\" * cmd * "{" * uni_to_TeX(base) * "}"
end

"""
An index list as LaTeX, one spacer group per index.

The spacers are what keep the order of raised and lowered slots: `T^{a}{}_{b}` and `T_{b}{}^{a}`
are different objects, and collecting the uppers and lowers into one group each would print them
the same. Tensors and operators both come through here, so neither can lose that order.
"""
function _latex_indices(inds)
    out = ""
    for (ind, pos) in inds
        tex = ind.name == "±" ? "\\pm" : ind.name == "∓" ? "\\mp" : uni_to_TeX(ind.name)
        mark = ind isa DummyPatternIndex ? "?" : ""
        out *= (pos.position == 0 ? "{ }^{" : "{ }_{") * mark * tex * "}"
    end
    return out
end

#allow for outer layer overwrite
function _to_latex(T::TensorSuperType)
    out = ""
    if name(T) != ""
        #γ5 is named "γ5" and also carries its 5 as a label, so only the letter is the name - the
        #same trim `tens2str` makes for display. Rendering the whole name gave γ5^5.
        nm = T isa Gamma5 ? string(first(name(T))) : name(T)
        out *= uni_to_TeX(nm) * " "
    end
    out *= _latex_indices(indices(T))

    lab = labels(T)
    if !isempty(lab)
        for (ind, pos) in lab
            if pos.position == 0
                out *= "{ }^{"
                out *= ind
                out*="}"

            else
                out *= "{ }_{"

                out*=ind
                out*="}"
            end
        end

    end
    if function_of_on
        fo = function_of(T)
        if !isempty(fo)
            names = sort!([c.name for c in fo], lt=_coord_lt)
            out *= "(" * join(uni_to_TeX.(names), ",") * ")"
        end
    end
    return Base.replace(out * " ", "\$" => "")

end

function to_latex(T::TensorSuperType)
    _to_latex(T)
end


function tens2str(T::TensorSuperType)
    out = ""
    if T isa Gamma5
        out *= first(name(T))
    else
        out *= _display_name(name(T))
    end
    tinds = indices(T)
    if !isempty(tinds)
        lowers = ""
        uppers = ""
        lowersempt::Bool = 1
        uppersempt::Bool = 1
        for (ind, pos) in tinds
            if pos.position == 0
                if typeof(ind) == DummyPatternIndex
                    uppers *= "?"
                end
                uppers = uppers * _display_name(ind.name)
                uppersempt = 0
                lowers = lowers * "∘"
            else
                uppers = uppers * "∘"
                lowersempt = 0
                if typeof(ind) == DummyPatternIndex
                    lowers *= "?"
                end
                lowers = lowers * _display_name(ind.name)
            end
        end
        if !uppersempt
            out *= "^{" * uppers * "}"
        end
        if !lowersempt
            out *= "_{" * lowers * "}"
        end
    end
    lab = labels(T)
    if !isempty(lab)
        lowers = ""
        uppers = ""
        lowersempts::Bool = 1
        uppersempts::Bool = 1
        for (ind, pos) in lab
            if pos.position == 0
                uppers = uppers * ind
                uppersempts = 0
                lowers = lowers * "∘"
            else
                uppers = uppers * "∘"
                lowersempts = 0
                lowers = lowers * ind
            end
        end
        if !uppersempts
            out *= "^{" * uppers * "}"
        end
        if !lowersempts
            out *= "_{" * lowers * "}"
        end
    end
    if function_of_on
        fo = function_of(T)
        if !isempty(fo)
            names = sort!([c.name for c in fo], lt=_coord_lt)
            out *= "(" * join(_display_name.(names), ",") * ")"
        end
    end
    return out
end

"""
Order coordinate names: Latin before Greek, lowercase before uppercase,
alphabetical within each group.
"""
function _coord_lt(a::AbstractString, b::AbstractString)
    ka, kb = _coord_key(a), _coord_key(b)
    ka != kb && return ka < kb
    return a < b
end

function _coord_key(s::AbstractString)
    c = first(s)
    greek = c >= '\u0370'          # Greek block starts at U+0370
    upper = isuppercase(c)
    return (greek ? 1 : 0, upper ? 1 : 0)
end


function replace(T::TensorTerm, old_T::TensorSuperType, new_T::TensorSuperType, repeat::Bool=true)
    replace(TensorExpression(T), old_T, new_T, repeat)
end
"""Naive replacements..."""
function replace(TEx::TensorExpression, old_T::TensorSuperType, new_T::TensorSuperType, repeat::Bool=true)
    for k in eachindex(TEx.summands)
        summand = TEx.summands[k]
        matches = []
        for tind in eachindex(summand[2].terms)
            if summand[2].terms[tind] == old_T
                push!(matches, tind)
            end
        end

        if repeat
            for match in matches
                summand[2].terms[match] = new_T
            end
        else
            if !isempty(matches)
                summand[2].terms[1] = new_T
            end
        end
    end
    return update_indices!(TEx)
end

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# YOUNG PROJECTION
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

#complicated logic to convert a youngtab to a projection operator of permuations, and apply to a tensor
#this was quite sensitive to get working for generic tableaux
##############################################################################################

const _PROJ_CAP = typemax(Int) ÷ 4

@inline _sat_mul(a::Int, b::Int) = (b != 0 && a > _PROJ_CAP ÷ b) ? _PROJ_CAP : a * b

function _fact_sat(n::Int)
    n >= 21 && return _PROJ_CAP
    f = 1
    for i in 2:n
        f *= i
    end
    return f
end

"""
    young_projector_size(x) -> Int
"""
function young_projector_size(tab::YoungTableau)
    rows = tab.row_lengths
    isempty(rows) && return 1
    n = 1
    for r in rows
        n = _sat_mul(n, _fact_sat(r))
    end
    for j in 1:maximum(rows)                 # column lengths are the conjugate partition
        n = _sat_mul(n, _fact_sat(count(>=(j), rows)))
    end
    return n
end

function young_projector_size(T::TensorSuperType)
    tabs = tableaux(T)
    isempty(tabs) && return 1
    n = 0
    for tab in tabs
        n = min(_PROJ_CAP, n + young_projector_size(tab))
    end
    return n
end

young_projector_size(T::OperatorSuperType) =
    _sat_mul(invoke(young_projector_size, Tuple{TensorSuperType}, T), young_projector_size(operand(T)))

function young_projector_size(T::TensorTerm)
    n = 1
    for t in terms(T)
        n = _sat_mul(n, young_projector_size(t))
    end
    return n
end

young_projector_size(TEx::TensorExpression) =
    isempty(summands(TEx)) ? 1 : maximum(young_projector_size(t) for (_, t) in summands(TEx))

"""Default ceiling on `young_projector_size` before multi-term reduction is skipped."""
const DEFAULT_MULTITERM_BUDGET = 1024

"""
    multiterm_reduce(TEx; budget=DEFAULT_MULTITERM_BUDGET) -> TensorExpression
"""
function multiterm_reduce(TEx::TensorExpression; budget::Int=DEFAULT_MULTITERM_BUDGET)
    sz = young_projector_size(TEx)
    (sz <= 1 || sz > budget) && return TEx
    proj = hard_merge(sort(collect_terms(young_projection(TEx))))
    isZeroTens(proj) && return zero()
    return length(summands(proj)) < length(summands(TEx)) ? proj : TEx
end

function young_projection(TT::TensorTerm, proj::Union{Vector{String},Nothing}=nothing)
    ts = terms(TT)
    isempty(ts) && return TensorExpression(TT)

    term_acc = nothing
    for tens in ts
        piece = if typeof(tens) <: OperatorSuperType
            newop = _with_operand(tens, young_projection(operand(tens), proj))
            young_projection(newop, proj)
        elseif typeof(tens) <: TensorSuperType
            young_projection(tens, proj)
        else
            continue
        end
        term_acc = term_acc === nothing ? piece : term_acc * piece
    end

    term_acc === nothing && return TensorExpression(TT)
    return term_acc
end

function young_projection(TEx::TensorExpression, proj::Union{Vector{String},Nothing}=nothing)
    s = summands(TEx)
    out = CoeffPair[]
    for (coef, tterm) in s
        piece = coef * young_projection(tterm, proj)
        # isZeroTens(piece) && continue
        append!(out, summands(piece))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end


function factor_weight(factor::String)
    weight = 0
    for subf in split_additive_factors(factor)
        coef, _, _ = parse_factor(subf)
        weight += abs(coef)
    end
    return weight
end

function factor_abs_sum(factor::String)
    sum_coeff = 0
    for subf in split_additive_factors(factor)
        coef, _, _ = parse_factor(subf)
        sum_coeff += abs(coef)
    end
    return sum_coeff
end

function compute_normalization_from_factors(factors::Vector{String})
    total = prod(factor_abs_sum(f) for f in factors)
    return 1 // total
end
function young_projection(T::TensorSuperType, proj::Union{Vector{String},Nothing}=nothing)
    projectors = proj === nothing ? [compute_young_projector(tab) for tab in tableaux(T)] : [proj for tab in tableaux(T)]
    # println(tab for tab in T.tableaux)
    norms = proj === nothing ? [1 // compute_hook_product(tab) for tab in tableaux(T)] : [compute_normalization_from_factors(prr) for prr in projectors]
    tex = zero()
    isempty(projectors) && return TensorExpression(T)
    for (proj, norm) in zip(projectors, norms)
        tex += apply_young_projector(T, proj) * norm
    end
    return tex
end

function apply_young_projector(T::TensorSuperType, projector::Vector{String})
    terms = Vector{Tuple{Rational{Int},Vector{Pair{IndexTypes,IndexPosition}}}}([(1, collect(indices(T)))])
    for proj in reverse(projector)
        terms = apply_string_factor(terms, proj)
    end
    result = Dict{Vector,Number}()
    for (s, t) in terms
        result[t] = get(result, t, 0) + s
    end
    no_dupes = [(s, t) for (t, s) in result if s != 0]
    tex = zero() #new tensor expression to hold result of projection
    for (prefact, inds) in no_dupes
        Tn = copy(T)
        set_indices!(Tn, inds)
        # println("prefact: ", prefact)
        tex += prefact * (Tn)
        # println(tex)
    end
    return tex #normalisation fixed elsewhere above

end

function compute_hook_product(tab::YoungTableau)
    hook_lengths = compute_hook_length(tab)
    return prod(hook_lengths[hook_lengths .> 0]; init=1)
end

function compute_hook_length(tab::YoungTableau)
    rows = tab.row_lengths
    n_rows = length(rows)
    n_cols = maximum(rows)

    hook_lengths = zeros(Int, n_rows, n_cols)

    for r in 1:n_rows
        for c in 1:rows[r]
            right = rows[r] - c
            below = sum(r2 >= c for r2 in rows[(r+1):end]; init=0)
            hook_lengths[r, c] = 1 + right + below
        end
    end

    return hook_lengths
end

# Swap multiple pairs in a tensor
function swap_string_pairs(t::Vector, pairs)
    t = copy(t)
    for (i, j) in pairs
        t[i], t[j] = t[j], t[i]
    end
    return t
end

# Apply a general n-cycle like (2 3 4)
function apply_cycle(t::Vector, cycle)
    t = copy(t)
    temp = t[cycle[end]]
    for i in reverse(2:length(cycle))
        t[cycle[i]] = t[cycle[i-1]]
    end
    t[cycle[1]] = temp
    return t
end


function parse_number(s::AbstractString)
    s = strip(s)
    if occursin('/', s)
        return parse(Rational{Int}, s)
    elseif occursin('.', s)
        return parse(Float64, s)
    else
        return parse(Int, s)
    end
end

function parse_factor(factor::SubString{String})
    factor = strip(factor)

    # Remove outer parentheses
    if startswith(factor, "(") && endswith(factor, ")")
        factor = strip(factor[2:(end-1)])
    end

    coef = one(Rational{Int})

    # Leading sign
    if startswith(factor, "+")
        factor = strip(factor[2:end])
    elseif startswith(factor, "-")
        coef = -one(Rational{Int})
        factor = strip(factor[2:end])
    end

    # Extract coefficient
    parts = split(factor, "(", limit=2)
    coeff = strip(parts[1])
    coeff = strip(replace(coeff, r"[+-]\s*$" => ""))

    if !(coeff == "" || coeff == "1")
        coef *= parse_number(coeff)
    end

    pairs = Tuple{Int,Int}[]
    cycles = Vector{Int}[]

    for m in eachmatch(r"\(([^)]+)\)", factor)
        nums = parse.(Int, split(strip(m.captures[1]), " "))

        if length(nums) == 2
            push!(pairs, (nums[1], nums[2]))
        elseif length(nums) > 2
            push!(cycles, nums)
        end
    end

    return coef, pairs, cycles
end


# Split "(1 - A + B)" → ["1 ", "- A", "+ B"]
function split_additive_factors(factor::String)
    inner = factor[2:(end-1)]
    return split(inner, r"(?=[+-])")
end


function apply_string_factor(terms, factor::String)
    new_terms = []

    for subf in split_additive_factors(factor)
        coef, pairs, cycles = parse_factor(subf)

        for (s, t) in terms
            t1 = swap_string_pairs(t, pairs)
            for cyc in cycles
                t1 = apply_cycle(t1, cyc)
            end
            push!(new_terms, (s * coef, t1))
        end
    end

    return new_terms
end

# permutation to cycle-list string
function perm_str(perm::Vector{<:Integer}, labels::Vector{<:Integer})
    cycles = perm_to_cycles(perm, labels)

    if isempty(cycles)
        return "1"
    end

    parts = ["(" * join(c, " ") * ")" for c in cycles]
    return join(parts, "")
end
# sign of permutation
function perm_sign(p::Vector{Int})
    s = 1
    for i in 1:(length(p)-1)
        for j in (i+1):length(p)
            if p[i] > p[j]
                s *= -1
            end
        end
    end
    return s
end

function perm_to_cycles(perm::Vector{<:Integer}, labels::Vector{<:Integer})
    n = length(perm)
    visited = falses(n)
    cycles = Vector{Vector{Int}}()

    for i in 1:n
        if !visited[i] && perm[i] != i
            cycle = Int[]
            j = i

            while !visited[j]
                visited[j] = true
                push!(cycle, labels[j])   # label at position j
                j = perm[j]               # move via permutation
            end

            push!(cycles, cycle)
        end
    end

    return cycles
end

# full symmetric-group sum (row symmetriser)
function row_sym_sum(labels)
    k = length(labels)
    perms = collect(permutations(1:k))

    terms = String[]
    first = true

    for perm in perms
        term = perm_str(perm, labels)

        if first
            push!(terms, term)
            first = false
        else
            push!(terms, "+ " * term)
        end
    end

    return "(" * join(terms, " ") * ")"
end

# full signed sum (column antisymmetriser)
function col_antisym_sum(labels)
    k = length(labels)
    perms = collect(permutations(1:k))

    terms = String[]
    first = true

    for perm in perms
        s = perm_sign(perm)
        term = perm_str(perm, labels)

        if first
            if s == 1
                push!(terms, term)
            else
                push!(terms, "- " * term)
            end
            first = false
        else
            if s == 1
                push!(terms, "+ " * term)
            else
                push!(terms, "- " * term)
            end
        end
    end

    return "(" * join(terms, " ") * ")"
end

function compute_young_projector(tab::YoungTableau)

    labels = tab.index_assignment
    row_lengths = tab.row_lengths

    # -------------------------
    # build rows
    # -------------------------
    rows = []
    idx = 1
    for rlen in row_lengths
        push!(rows, labels[idx:(idx+rlen-1)])
        idx += rlen
    end

    # -------------------------
    # build columns
    # -------------------------
    n_cols = maximum(row_lengths)

    columns = []
    for c in 1:n_cols
        col = Int[]
        for r in rows
            if c <= length(r)
                push!(col, r[c])
            end
        end
        push!(columns, col)
    end
    factors = String[]

    # columns first (antisym)

    # rows second (sym)
    for row in rows
        if length(row) > 1
            push!(factors, row_sym_sum(row))
        end
    end
    for col in columns
        if length(col) > 1
            push!(factors, col_antisym_sum(col))
        end
    end

    return factors
end


# Function to generate transpositions for a row or column
function transpositions(indices::Vector)
    ts = []
    n = length(indices)
    for i in 1:(n-1)
        for j in (i+1):n
            push!(ts, (indices[i], indices[j]))
        end
    end
    return ts
end


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# COPYING
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
symcopy(x) = error("symcopy not defined on ", typeof(x), ".")
symcopy(x::Pair) = x
symcopy(x::String) = x
symcopy(t::TensorSuperType) = typeof(t)(t.data, t.indices)
symcopy(v::Vector{<:Pair}) = copy(v)                    # indices, labels, pairs are immutable
symcopy(d::Dict) = copy(d)                              # weights
symcopy(s::Vector{Coordinate}) = s
symcopy(d::TensorData) = d           # immutable and shared, nothing to copy
function symcopy(t::TensorExpression)
    s = summands(t)
    out = Vector{CoeffPair}(undef, length(s))
    @inbounds for i in eachindex(s)
        out[i] = s[i][1] => symcopy(s[i][2])
    end
    return TensorExpression(_RAW, out)
end
function symcopy(tt::TensorTerm)
    nt = TensorTerm(TensorSuperType[symcopy(t) for t in terms(tt)])
    return nt
end
#see if below like above?
function copy(T::TensorTerm)
    ts = terms(T)
    isempty(ts) && return TensorTerm(Vector{TensorSuperType}())
    out = Vector{TensorSuperType}(undef, length(ts))
    @inbounds for i in eachindex(ts)
        out[i] = copy(ts[i])
    end
    return TensorTerm(out)
end
copy(T::TensorSuperType) = typeof(T)(T.data, T.indices)
function copy(T::TensorExpression)
    s = summands(T)
    if length(s) == 1 && isZero(s[1][1])
        return TensorExpression()
    end
    out = CoeffPair[]
    sizehint!(out, length(s))
    @inbounds for i in eachindex(s)
        isZero(s[i][1]) && continue
        push!(out, CoeffPair(s[i][1], copy(s[i][2])))
    end
    isempty(out) && return TensorExpression()
    return TensorExpression(_RAW, out)
end
"""
    copy_owned(x)

Like symcopy but the result has its own index vectors.
Safe to mutate in place without affecting the original.
"""
copy_owned(x) = symcopy(x)

copy_owned(t::TensorSuperType) = typeof(t)(t.data, copy(t.indices))
copy_owned(tt::TensorTerm) = TensorTerm(TensorSuperType[copy_owned(t) for t in terms(tt)])
function copy_owned(t::TensorExpression)
    nt = TensorExpression()
    nt.summands = CoeffPair[coef => copy_owned(term) for (coef, term) in summands(t)]
    return nt
end


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# UTILS
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
#more or less generic functions that either are simple or have many uses across the codebase

"""
The zero() tensor expression used globally in Alakazam.

This is a method on `Base.zero` rather than a separate function: a package-owned `zero` would clash
with the one `Base` exports, forcing users to write `Alakazam.zero()`. `Base` has no zero-argument
`zero`, so this signature is unambiguous.
"""
const _zero = TensorExpression()
@inline zero()::TensorExpression = _zero
const _consttensor=Tensor("")
@inline consttens()::Tensor = _consttensor


function first_element(arr, elem1, elem2)
    idx1 = findfirst(x -> x == elem1, arr)
    idx2 = findfirst(x -> x == elem2, arr)
    if idx1 == nothing || (idx2 != nothing && idx2 < idx1)
        return elem2
    else
        return elem1
    end
end

fresh_indices(T::Union{TensorSuperType,TensorTerm}, E::Union{TensorSuperType,TensorTerm}) = fresh_indices(T, TensorExpression(E))
function fresh_indices(T::Union{TensorSuperType,TensorTerm},
    E::TensorExpression)
    T = copy(T)
    # all indices appearing in the expression
    forbidden = OrderedSet{IndexSuperType}()
    collect_indices!(forbidden, E)
    collect_indices!(forbidden, T)
    # available fresh indices
    avail = Dict{IndexSet,OrderedSet{IndexSuperType}}()

    for set in get_index_sets(T)
        avail[set] = OrderedSet(
            filter(x -> !(x in forbidden),
                get_unused_indices(T, set))
        )
    end
    dummie_matchup = Dict{IndexSuperType,IndexSuperType}()
    rep_guide = Dict{Pair{IndexSuperType,IndexPosition},
        Pair{IndexSuperType,IndexPosition}}()

    for ip in indices(T)

        oldind = first(ip)
        oldpos = last(ip)
        isLightcone(oldind) && continue
        if oldind in forbidden

            set = oldind.index_set

            if oldind in dummies(T)
                if haskey(dummie_matchup, oldind)
                    newind = dummie_matchup[oldind]
                else
                    isempty(avail[set]) && begin
                        nid = generate_indices(set)
                        nid === nothing && throw(OutOfIndices(set))
                        push!(avail[set], nid)
                    end

                    newind = popfirst!(avail[set])
                    dummie_matchup[oldind] = newind
                end
            else
                isempty(avail[set]) && begin
                    nid = generate_indices(set)
                    nid === nothing && throw(OutOfIndices(set))
                    push!(avail[set], nid)
                end

                newind = popfirst!(avail[set])
            end

            rep_guide[oldind=>oldpos] = newind => oldpos
        end
    end

    set_indices_deep!(T, rep_guide)

    return T
end


"""
Get all of the used indices, including dummies, in an expression.
"""
function collect_indices!(inds::OrderedSet{IndexSuperType}, T)
    if T isa OperatorSuperType
        # indices carried by the operator itself
        for ind in op_indices(T)
            push!(inds, ind[1])
        end
        # indices hidden inside the operand
        collect_indices!(inds, operand(T))

    elseif T isa TensorExpression

        for summand in summands(T)
            collect_indices!(inds, summand[2])
        end

    elseif T isa TensorTerm
        for term in terms(T)
            collect_indices!(inds, term)
        end
    else
        # ordinary tensor term
        for ind in indices(T)
            push!(inds, ind[1])
        end
    end

    return inds
end
function collect_indices!(inds::Vector{IndexTypes}, T)
    if T isa OperatorSuperType
        for ind in op_indices(T)
            ind[1] in inds || push!(inds, ind[1])
        end
        collect_indices!(inds, operand(T))
    elseif T isa TensorExpression
        for summand in summands(T)
            collect_indices!(inds, summand[2])
        end
    elseif T isa TensorTerm
        for term in terms(T)
            collect_indices!(inds, term)
        end
    else
        for ind in indices(T)
            ind[1] in inds || push!(inds, ind[1])
        end
    end
    return inds
end
