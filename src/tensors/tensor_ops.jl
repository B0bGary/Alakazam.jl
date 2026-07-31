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
function ==(a::TensorSuperType, b::TensorSuperType)
    if (name(a) == name(b) && indices(a) == indices(b) && function_of(a)==function_of(b) && tableaux(a)==tableaux(b))
        # println("apparnetly equl inds")
        return true
    else
        # println("false at tens level")
        return false
    end
end

# Base.hash(d::TensorSuperType, h::UInt) = hash(
#     (
#         name(d), indices(d), function_of(d), tableaux(d)
#     ),
#     h
# )
function Base.hash(d::TensorSuperType, h::UInt)::UInt
    h = hash(data(d).name, h)
    @inbounds for p in indices(d)
        h = hash(p[1].uid, hash(p[2].position, h))
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
==(a::TensorExpression, b::TensorTerm) = singleterm(a) && first_coef(a) == 1 && (first_summand(a) == b)
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
# _hash_tensor(term, h::UInt)::UInt = hash(term, h)#specialises per concrete type, optimises, TTerm hash
# function Base.hash(t::TensorTerm, h::UInt)::UInt
#     for term in terms(t)
#         h = _hash_tensor(term, h)
#     end
#     return h
# end
function ==(a::TensorExpression, b::TensorExpression)
    if length(a.summands) != length(b.summands)
        # @show "wrong leng"
        # @show coeff_equal(first_coef(b),0)

        (length(a.summands)==0 && length(b.summands)==1 && coeff_equal(first_coef(b), 0)) || (length(b.summands)==0 && length(a.summands)==1 && coeff_equal(first_coef(a), 0)) ? (return true) : (return false)
        # return false
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
# #add single term check to make things less expensive
# ==(a::TensorExpression, b::TensorTerm) = a == TensorExpression(b)
# ==(a::TensorTerm, b::TensorExpression) = TensorExpression(a) == b
# ==(a::TensorSuperType, b::TensorExpression) = TensorExpression(a) == b
# ==(a::TensorExpression, b::TensorSuperType) = a == TensorExpression(b)
# ==(a::TensorSuperType, b::TensorTerm) = ==(TensorTerm(a), b)
# ==(a::TensorTerm, b::TensorSuperType) = ==(a, TensorTerm(b))
# ==(a::TensorTerm, b::TensorTerm) = terms(a) == terms(b)



# """Stronger coef checking for symbolics nastiness."""
# function coeff_equal(a::SymorNum, b::SymorNum)
#     isequal(a, b) && return true
#     return Symbolics.simplify(Symbolics.expand(a - b))===0
# end
coeff_equal(a::Number, b::Number) = a == b            # no Symbolics path at all
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

function hard_equals_graded(a::TensorTerm, b::TensorTerm)
    return same_tensorterm_mod_syms_grade(a, b)
end


cong(T1, T2) = T1≅T2

function ≅(a::TensorTerm, b::TensorTerm)
    return hard_equals(a, b)
end

≅(a::TensorSuperType, b::TensorExpression) = (b ≅ a)
≅(a::TensorExpression, b::TensorSuperType) = a ≅ TensorExpression(b)


# function ==(a::Vector{T} where {T<:TensorSuperType}, b::Vector{T} where {T<:TensorSuperType})
#     if length(a) != length(b)
#         return false
#     end
#     for i = eachindex(a)
#         if !(a[i] == b[i])
#             return false
#         end
#     end
#     return true
# end

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
function ≅(a::TensorExpression, b::TensorExpression)
    a = hard_merge(a)
    b = hard_merge(b)
    C=merge_sum_syms(a-b)
    m = hard_merge(symmetrise_partials(C))
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
    new = zero()
    for summand2 in T2.summands
        for summand1 in T1.summands
            new += (summand2[1] * summand1[1]) * (summand2[2] × summand1[2])
        end
    end

    return new
end

function ×(T1::TensorExpression, T2::Union{TensorSuperType,TensorTerm})
    new_terms = TensorExpression()
    T1 = rename_dummies(T1)
    for summandin in eachindex(T1.summands)
        t = (T1.summands[summandin][2]×T2)
        new_terms += (T1.summands[summandin][1] * t)
    end
    return new_terms
end

function ×(T2::Union{TensorSuperType,TensorTerm}, T1::TensorExpression)
    new_terms = TensorExpression()
    T1 = rename_dummies(T1)
    for summandin in eachindex(T1.summands)
        t = (T2 × T1.summands[summandin][2])
        new_terms += (T1.summands[summandin][1] * t)
    end
    return new_terms
end


function ×(T1::Union{TensorTerm,TensorSuperType}, T2::Union{TensorTerm,TensorSuperType})
    T2 = copy(T2)
    t1inds = OrderedSet{IndexSuperType}()
    t2inds = OrderedSet{Pair{Integer,IndexSuperType}}() #all indices, numbered in order
    for inds in indices(T1)
        push!(t1inds, inds[1])
    end
    for inds in eachindex(indices(T2))
        iddd=indices(T2)[inds][1]
        if !isLightcone(iddd) && iddd in t1inds #if theres a repeated indice  
            push!(t2inds, inds => iddd) # which index to replace
        else
            # push!(t2ok,T2.indices[inds][1])# index doesnt need to be replaced
        end
    end
    # we now have all the indices we need to replace and where they are on the tensor
    sets = get_index_sets(T2) #Set{IndexSet}()
    avail = Dict{IndexSet,OrderedSet{IndexSuperType}}() #pairs of index sets, with the unused indices
    for set in sets
        # set.name=="lightcone" && continue
        a = OrderedSet{IndexSuperType}(intersect(get_unused_indices(T1, set), get_unused_indices(T2, set)))
        push!(avail, set => a) #ones we can replace with
    end
    dummie_matchup = Dict{IndexSuperType,IndexSuperType}()
    t2dums = dummies(T2)
    replc=Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
    for ind in t2inds # for each we need to replace (vec pos => index)
        # if isLightcone(ind[2]) 
        #     nind=ind[2]
        # else
        # @show ind "not lc"
        if ind[2] in t2dums

            if ind[2] in keys(dummie_matchup)
                if isLightcone(ind[2])
                    nind=ind[2]
                else
                    nind = pop!(dummie_matchup, ind[2])
                end
            else
                nind=safe_pop!(avail[ind[2].index_set], ind[2].index_set)
                push!(dummie_matchup, ind[2] => nind)
            end
        else
            if isLightcone(ind[2])
                nind=ind[2]
            else
                nind=safe_pop!(avail[ind[2].index_set], ind[2].index_set)
            end
        end
        # end
        replc[ind[2]=>indices(T2)[ind[1]][2]]=nind => indices(T2)[ind[1]][2]
    end
    set_indices_deep!(T2, replc)
    return T1 * T2
end

##############################################################################################################################
##############################################################################################################################
#indices
##############################################################################################################################
#update indices

update_indices!(T::TensorSuperType) = T

"""
Updates external indices of a tensor term, expression, or operator, by looking at the indices of all the terms it contains.
"""
function update_indices!(T::TensorTerm, n::Int=-1)
    n==-1 && (n=length(indices(T)))
    inds = Vector{Pair{IndexTypes,IndexPosition}}(undef, n)
    k=1
    for t in terms(T)
        typeof(t) <: OperatorSuperType && update_indices!(t)
        k=_copy_inds!(inds, k, t)
    end
    # @assert k==n+1
    #new checker 
    @inbounds for i in eachindex(inds)
        idx = inds[i][1]
        isLightcone(idx) && continue
        c = 1
        for j in (i+1):length(inds)
            inds[j][1] == idx && (c += 1)
        end
        c > 2 && (println(T); throw(OverusedIndices(T)))
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
    used = IndexSuperType[]
    collect_indices!(used, T)                  # needs a Vector-accepting collect_indices!
    out = IndexSuperType[]
    for i in IS.indices
        i isa DummyPatternIndex && continue    # inline filter_wildcard_dummies
        i in used && continue                   # inline setdiff
        push!(out, i)
    end
    return out
end
function filter_wildcard_dummies(I::Vector{IndexSuperType})
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

    all_indices = IndexSuperType[]
    for i in indices(T)
        i[1] in all_indices || push!(all_indices, i[1])
    end

    replacements = Pair{IndexSuperType,IndexSuperType}[]
    inset = IndexSuperType[]
    remain = IndexSuperType[]
    ordered = IndexSuperType[]

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
"""Swap the positions of all indices in an expression, i.e. all lower indices are moved to upper, and vice versa."""
function op_pos(tt::Tensorial)
    transform(t -> op_pos(t), tt)
end

op_pos(tt::TensorSuperType) = set_indices_deep!(copy(tt), op_pos(indices(tt)))




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

    all_indices = IndexSuperType[]
    for t in terms(T), i in indices(t)
        i[1] in all_indices || push!(all_indices, i[1])
    end

    inset = IndexSuperType[]
    remain = IndexSuperType[]

    for set in index_sets
        # indices of this set already in use, excluding the dummies being replaced
        empty!(inset)
        for a in all_indices
            a.index_set === set || continue
            a in dumms && continue
            push!(inset, a)
        end

        # candidates: this set's non-wildcard indices not already in use
        empty!(remain)
        for x in set.indices
            x isa DummyPatternIndex && continue
            x in inset && continue
            push!(remain, x)
        end

        replacements = Dict{IndexSuperType,IndexSuperType}()
        for d in dumms
            d.index_set === set || continue
            if isempty(remain)
                nid = generate_indices(set)
                nid === nothing && throw(OutOfIndices(set))
                push!(remain, nid)
            end
            replacements[d] = popfirst!(remain)
        end
        isempty(replacements) || (T = rename_recursive!(T, replacements))
    end
    return T
end
function rename_recursive!(T::TensorExpression, replacements)
    new_exp = zero()
    for term in summands(T)
        new_exp += term[1] * rename_recursive!(term[2], replacements)
    end
    return new_exp
end

function rename_recursive!(T::TensorTerm, replacements)
    for term in terms(T)
        rename_recursive!(term, replacements)
    end
    update_ext_indices!(T)
    return T
end

function rename_recursive!(t::TensorSuperType, replacements)
    if typeof(t) <: OperatorSuperType
        if !isempty(op_indices(t))
            set_op_indices!(t, [Pair((haskey(replacements, x) ? replacements[x] : x), y) for (x, y) in op_indices(t)])
        end
        oprd = operand(t)
        rename_recursive!(oprd, replacements)
        return t
    else
        if !isempty(indices(t))
            set_indices!(t, [Pair((haskey(replacements, x) ? replacements[x] : x), y) for (x, y) in indices(t)])
        end
        return t
    end
end


function rename_dummies(T::TensorExpression)
    new_exp = zero()
    for term in summands(T)
        new_exp += term[1] * rename_dummies(term[2])
    end
    return new_exp
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




"""Look through the symmetries of the tensor, and find the configuration of indices that is closest to the lexo-canonical order, as measured by the rep_distance function.
Returns the sign and the new configuration of indices. If the cost is  zero(), we can return immediately, as we have found the lexo-canonical representation."""
# function find_canonical_rep(T::TensorSuperType)
#     # sort indices into canonical order, serving as a baseline of cost  zero().
#     lexo = indices(T)[sortperm(indices(T))]#this will automatically put in alphabetical
#     best_inds = indices(T)
#     sign = 1
#     best_cost = Inf
#     for tableau in tableaux(T)
#         syms = tableau.symmetries
#         for perm in syms
#             new_inds = apply_permutation_cycle(indices(T), perm)
#             cost = rep_distance(lexo, new_inds)
#             if cost == 0
#                 return (1, new_inds)
#             end
#             if cost < best_cost || (cost == best_cost && new_inds < best_inds)
#                 best_inds = new_inds
#                 best_cost = cost
#             end
#         end
#         asyms = tableau.antisymmetries
#         for perm in asyms
#             new_inds = apply_permutation_cycle(indices(T), perm)
#             cost = rep_distance(lexo, new_inds)
#             if cost == 0
#                 return (-1, new_inds)
#             end
#             if cost < best_cost || (cost == best_cost && new_inds < best_inds)
#                 best_inds = new_inds
#                 best_cost = cost
#                 sign = -1
#             end
#         end

#     end
#     # apply symmetries to the tensor, and measure the cost by the difference from the baseline.
#     return (sign, best_inds)
# end
# function find_canonical_rep(T::TensorSuperType)::Tuple{Int,Vector{IndexPair}}
#     inds=indices(T)
#     lexo = indices(T)[sortperm(indices(T))]
#     if isSymmetric(T)
#         return (1, lexo)
#     end

#     if isAntisymmetric(T)
#         p = sortperm(inds)
#         return (permutation_sign(p), lexo)
#     end
#     start_inds = indices(T)
#     best_inds = start_inds
#     best_cost = rep_distance(lexo, start_inds)
#     best_sign = 1

#     queue = [(start_inds, 1)]
#     head = 1

#     # permutation => signs already seen (+1 and/or -1)
#     seen = Dict{Tuple,Set{Int}}()
#     # seen = Dict{typeof(start_inds),Set{Int}}()
#     # seen[start_inds] = Set([1])
#     seen[Tuple(start_inds)] = Set([1])

#     while head <= length(queue)
#         inds, sign = queue[head]
#         head += 1

#         cost = rep_distance(lexo, inds)

#         if cost == 0
#             return (sign, inds)
#         end

#         if cost < best_cost || (cost == best_cost && inds < best_inds)
#             best_cost = cost
#             best_inds = inds
#             best_sign = sign
#         end

#         for tab in tableaux(T)

#             row_gens = row_generators(tab.tableau)
#             col_gens = column_generators(tab.tableau)

#             # symmetric generators
#             for g in row_gens
#                 cycle = cycles_to_string(dict_to_cycles(g))
#                 new_inds = apply_permutation_cycle(inds, cycle)

#                 signs = get!(seen, Tuple(new_inds), Set{Int}())
#                 if !(sign in signs)
#                     push!(signs, sign)
#                     push!(queue, (new_inds, sign))
#                 end
#             end

#             # antisymmetric generators
#             for g in col_gens
#                 cycle = cycles_to_string(dict_to_cycles(g))
#                 new_inds = apply_permutation_cycle(inds, cycle)

#                 new_sign = -sign
#                 signs = get!(seen, Tuple(new_inds), Set{Int}())
#                 if !(new_sign in signs)
#                     push!(signs, new_sign)
#                     push!(queue, (new_inds, new_sign))
#                 end
#             end
#         end
#     end

#     return (best_sign, best_inds)
# end
@inline _signbit(s::Int) = s > 0 ? 0x01 : 0x02
"""
    find_canonical_rep(T::TensorSuperType) -> (sign, indices)

Find the canonical index ordering for `T` under the symmetries described by its
Young tableaux, together with the sign picked up in reaching it.

Canonical form is the lexicographic ordering. Shortcut for fully
symmetric or fully antisymmetric tensor using perm.
Otherwise the tableaux generate a group of allowed index permutations, and the function searches it
breadth-first for a representative matching the lexicographic target.
"""
function find_canonical_rep(T::TensorSuperType)::Tuple{Int,Vector{IndexPair}}
    inds = indices(T)
    perm = sortperm(inds)
    lexo = inds[perm]

    isSymmetric(T) && return (1, lexo)
    isAntisymmetric(T) && return (permutation_sign(perm), lexo)

    start_inds = inds
    best_inds = start_inds
    best_cost = rep_distance(lexo, start_inds)
    best_sign = 1

    # generators depend only on T's tableaux — build once
    row_cycles = String[]
    col_cycles = String[]
    for tab in tableaux(T)
        for g in row_generators(tab.tableau)
            push!(row_cycles, cycles_to_string(dict_to_cycles(g)))
        end
        for g in column_generators(tab.tableau)
            push!(col_cycles, cycles_to_string(dict_to_cycles(g)))
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
            new_inds = apply_permutation_cycle(curr, cycle)
            b = _signbit(sign)
            m = get(seen, new_inds, 0x00)
            if m & b == 0x00
                seen[new_inds] = m | b
                push!(queue, (new_inds, sign))
            end
        end

        for cycle in col_cycles
            new_inds = apply_permutation_cycle(curr, cycle)
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
# function find_canonical_rep(T::TensorSuperType)
#     lexo = indices(T)[sortperm(indices(T))]
#     inds = indices(T)
#     best_inds = inds
#     sign = 1
#     best_cost = rep_distance(lexo, inds)

#     for tab in tableaux(T)
#         row_inds = sorted_by_blocks(inds, tab.tableau)
#         row_cost = rep_distance(lexo, row_inds)
#         if row_cost == 0
#             return (1, row_inds)
#         elseif row_cost < best_cost || (row_cost == best_cost && row_inds < best_inds)
#             best_inds, best_cost, sign = row_inds, row_cost, 1
#         end

#         col_inds = sorted_by_blocks(inds, get_columns(tab))
#         col_cost = rep_distance(lexo, col_inds)
#         if col_cost == 0
#             return (-1, col_inds)
#         elseif col_cost < best_cost || (col_cost == best_cost && col_inds < best_inds)
#             best_inds, best_cost, sign = col_inds, col_cost, -1
#         end
#     end

#     return (sign, best_inds)
# end
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
    exp=symmetrise_partials(TEx)
    # exp = rename_dummies(exp)
    # exp = contract_metrics(exp, EpsilonTensor)
    # @show exp,"start"
    # @show TEx
    exp = pop_metric(exp)
    # @show exp,"popped"
    exp = sort(exp)
    # @show exp,"sorted"
    exp = eliminate_epsilon(exp)
    # @show exp,"elim"
    # exp = rename_dummies(exp)
    new_exp = simp_ops_n_sym(exp)
    new_exp = contract_metrics(new_exp, Metric)
    new_exp = eliminate_kronecker(new_exp)
    # new_exp = sort(new_exp)
    new_exp = rename_dummies(new_exp)
    new_exp = sort(new_exp)
    return new_exp
end

canonicalise(T::Tensorial) = canonicalise(TensorExpression(T),true)
function canonicalise(TEx::TensorExpression,cp=false)
    cp==false && (TEx=copy(TEx))
    TEx=merge_sum_syms(TEx)
    return sort(TEx)
end
canonicalize(T) = canonicalise(T)
# println(rename_dummies(oo2))
# println(merge_sum(rename_dummies(oo2)))
# outexp = hard_merge(new_exp)
# new_exp=symmetrise_partials(new_exp)
# println("ouput symmetrised partials: ",new_exp)

# println(new_exp)

function simp_ops_n_sym(exp::TensorExpression)
    new_exp = zero()
    for summand in exp.summands
        # println("start ",summand)
        new_term = summand[1]
        for tens in summand[2].terms
            temp = tens
            if typeof(tens) <: OperatorSuperType
                # println("old oprd: ", tens.operand)
                op = simp_ops_n_sym(operand(tens))
                #  println("new oprd: ",op)
                if op == zero() || op === nothing
                    new_term *= zero()
                else
                    set_operand!(tens, op) #operator on nth term
                    temp = tens
                    new_term *= sort_symmetric(temp)
                end
            else
                # println(temp)
                # println(new_term)
                new_term *= sort_symmetric(temp)
            end
        end
        ##println(new_exp)
        ##println(new_term)
        new_exp += new_term
    end
    return new_exp
end



hard_simplify(T::TensorTerm) = hard_simplify(TensorExpression(T))
hard_simplify(T::TensorSuperType) = hard_simplify(TensorExpression(T))

function hard_simplify(TEx::TensorExpression)
    return hard_merge(sort(TEx))
end



####################################################################
#merging terms with same tensor content
# default is merge_sum, naive merging of coeffs
# hard: full aggressive simplification
# syms: use basic syms of tensors to canonicalise
#_merge2,_insert_merge: helpers for Texpr construction to save time
## need to add functionality for symmetric and antisym tensors


# function hard_merge(TEx::TensorExpression)
#     if TEx == zero()
#         return zero()
#     end
#     ##println(" zero() ma:", zero())
#     out = filter!(x -> !isZero(x[1]), TEx.summands) #removes 0 coefficnets
#     if isempty(out)
#         return zero()
#     end
#     TEx.summands = out
#     new_merged = CoeffPair[]
#     for summand in TEx.summands
#         find_weighted_square(construct_tensorterm_graph(summand[2])) && continue #drop if  zero()
#         found = false
#         for i in eachindex(new_merged)
#             isZero(new_merged[i][1]) && continue #skip comparing against zeros
#             (eq, gr) = hard_equals_graded(summand[2], new_merged[i][2])
#             if eq
#                 new_merged[i] = (gr * new_merged[i][1] + summand[1] => summand[2])
#                 found = true
#                 break
#             end
#         end
#         if !found
#             push!(new_merged, summand)
#         end
#     end
#     TEx.summands = new_merged
#     TEx.summands = filter!(x -> !isZero(x[1]), TEx.summands)
#     return isempty(TEx.summands) ? zero() : TEx
# end

function hard_merge(TEx::TensorExpression)
    src = TEx.summands
    n = length(src)
    n == 0 && return zero()

    out = Vector{CoeffPair}(undef, n)
    m = 0

    @inbounds for i in 1:n
        s = src[i]
        isZero(s[1]) && continue
        t = s[2]
        find_weighted_square(construct_tensorterm_graph(t)) && continue   # identically zero

        found = false
        for k in 1:m
            eq, gr = hard_equals_graded(t, out[k][2])
            eq || continue
            c = gr * out[k][1] + s[1]
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

        found || (m += 1; out[m]=s)
    end

    m == 0 && return zero()
    m < n && resize!(out, m)
    TEx.summands = out
    return TEx
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

# function merge_sum(TEx::TensorExpression)
#     out = filter!(x -> !isZero(x[1]), TEx.summands)
#     isempty(out) && return  zero()
#     n = length(out)
#     new_merged = CoeffPair[]
#     if n <= 8 #for small sums, do direclty. for large, use a hashing
#         for summand in out
#             found = false
#             for i in eachindex(new_merged)
#                 if summand[2] == new_merged[i][2]
#                     new_merged[i] = (new_merged[i][1] + summand[1]) => summand[2]
#                     found = true;
#                     break
#                 end
#             end
#             found || push!(new_merged, summand)
#         end
#         TEx.summands = filter!(x -> !isZero(x[1]), new_merged)
#         return isempty(TEx.summands) ?  zero() : TEx
#     else
#         buckets = Dict{UInt,Vector{Int}}()
#         for summand in out
#             h = hash(summand[2])
#             found = false
#             for i in get(buckets, h, Int[])
#                 if summand[2] == new_merged[i][2]
#                     new_merged[i] = (new_merged[i][1] + summand[1]) => summand[2]
#                     found = true;
#                     break
#                 end
#             end
#             found || (push!(new_merged, summand); push!(get!(buckets, h, Int[]), length(new_merged)))
#         end
#     end
#     TEx.summands = filter!(x -> !isZero(x[1]), new_merged)
#     return isempty(TEx.summands) ?  zero() : TEx
# end
# function merge_sum(TEx::TensorExpression)
#     src = TEx.summands
#     n = length(src)
#     n == 0 && return zero()

#     new_merged = CoeffPair[]
#     sizehint!(new_merged, n)

#     if n <= 8
#         ##for small sums, do direclty. for large, use a hashing
#         @inbounds for s in src
#             isZero(s[1]) && continue
#             found = false
#             for i in eachindex(new_merged)
#                 if s[2] == new_merged[i][2]
#                     c = new_merged[i][1] + s[1]
#                     if isZero(c)
#                         deleteat!(new_merged, i)
#                     else
#                         new_merged[i] = CoeffPair(c, s[2])
#                     end
#                     found = true
#                     break
#                 end
#             end
#             found || push!(new_merged, s)
#         end
#     else
#         # large hash, mark cancellations
#         buckets = Dict{UInt,Vector{Int}}()
#         ncancelled = 0
#         @inbounds for s in src
#             isZero(s[1]) && continue
#             h = hash(s[2])
#             found = false
#             for i in get(buckets, h, Int[])
#                 if s[2] == new_merged[i][2]
#                     c = new_merged[i][1] + s[1]
#                     new_merged[i] = CoeffPair(c, s[2])
#                     isZero(c) && (ncancelled += 1)
#                     found = true
#                     break
#                 end
#             end
#             if !found
#                 push!(new_merged, s)
#                 push!(get!(buckets, h, Int[]), length(new_merged))
#             end
#         end
#         ncancelled > 0 && filter!(x -> !isZero(x[1]), new_merged)
#     end

#     isempty(new_merged) && return zero()
#     TEx.summands = new_merged
#     return TEx
# end
const _EMPTY_IDX = Int[]

function merge_sum(TEx::TensorExpression)
    src = TEx.summands
    n = length(src)
    n == 0 && return zero()
    if n == 1
        @inbounds return isZero(src[1][1]) ? zero() : TEx
    end

    # --- detection: first summand that needs anything doing to it ---
    firstbad = 0
    @inbounds for i in 1:n
        isZero(src[i][1]) && (firstbad=i; break)
        ti = src[i][2]
        for j in 1:(i-1)
            src[j][2] == ti && (firstbad=i; break)
        end
        firstbad != 0 && break
    end
    firstbad == 0 && return TEx          # nothing to do — ZERO allocations

    # --- one exact-size buffer, manual length, no push!/deleteat! ---
    out = Vector{CoeffPair}(undef, n)
    m = firstbad - 1
    @inbounds for i in 1:m
        out[i] = src[i]                  # clean prefix, already known distinct & nonzero
    end

    if n <= 16
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
        buckets = Dict{UInt,Vector{Int}}()
        @inbounds for k in 1:m
            push!(get!(() -> Int[], buckets, hash(out[k][2])), k)
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

function merge_sum_syms(TEx::TensorExpression)
    # @show TEx
    out = filter!(x -> !isZero(x[1]), TEx.summands) #removes 0 coefficnets
    if isempty(out)
        return zero()
    end
    TEx.summands = out
    new_merged = CoeffPair[]
    for summand in TEx.summands
        found = false
        for i in eachindex(new_merged)

            if summand[2] == new_merged[i][2]
                # @show summand[2], new_merged[i][2]
                # @show summand[2] , new_merged[i][2] "equal"
                new_merged[i] = (new_merged[i][1] + summand[1] => summand[2])
                found = true
                break
            elseif length(summand[2].terms) == 1 && length(new_merged[i][2].terms) == 1
                rescomp = same_tensor_mod_syms(summand[2].terms[1], new_merged[i][2].terms[1])
                # @show rescomp
                if rescomp[1] == true #i had this block removed before idk why...
                    new_merged[i] = (rescomp[2] * new_merged[i][1] + summand[1] => summand[2])
                    found = true
                end
            elseif length(summand[2].terms) == length(new_merged[i][2].terms)
            end
        end
        if !found
            push!(new_merged, summand)
        end
    end

    TEx.summands = new_merged
    TEx.summands = filter!(x -> !isZero(x[1]), TEx.summands)
    return isempty(TEx.summands) ? zero() : TEx
end

# function merge_sum_syms(TEx::TensorExpression)
#     src = TEx.summands
#     n = length(src)
#     n == 0 && return zero()
#     if n == 1
#         @inbounds return isZero(src[1][1]) ? zero() : TEx
#     end

#     out = Vector{CoeffPair}(undef, n)
#     m = 0

#     @inbounds for i in 1:n
#         s = src[i]
#         isZero(s[1]) && continue
#         t = s[2]
#         found = false

#         # cheap pass: exact term equality
#         for k in 1:m
#             if out[k][2] == t
#                 c = out[k][1] + s[1]
#                 if isZero(c)
#                     for q in k:m-1; out[q] = out[q+1]; end
#                     m -= 1
#                 else
#                     out[k] = CoeffPair(c, t)
#                 end
#                 found = true; break
#             end
#         end

#         # expensive pass: symmetry-equivalence, single-tensor terms only
#         if !found && length(terms(t)) == 1
#             t1 = terms(t)[1]
#             for k in 1:m
#                 ok = out[k][2]
#                 length(terms(ok)) == 1 || continue
#                 same, sgn = same_tensor_mod_syms(terms(ok)[1], t1)
#                 same || continue
#                 c = sgn * out[k][1] + s[1]
#                 if isZero(c)
#                     for q in k:m-1; out[q] = out[q+1]; end
#                     m -= 1
#                 else
#                     out[k] = CoeffPair(c, out[k][2])
#                 end
#                 found = true; break
#             end
#         end

#         found || (m += 1; out[m] = s)
#     end

#     m == 0 && return zero()
#     m < n && resize!(out, m)
#     TEx.summands = out
#     return TEx
# end


##############################################################################################################################
##############################################################################################################################

# printing and display

function Base.show(io::IO, T::TensorSuperType)
    print(io, tens2str(T))
end

function Base.show(io::IO, T::TensorTerm)
    for tens in terms(T)
        print(io, tens)
    end
end
_unwrapnum(x) = SymbolicUtils.isconst(x) ? SymbolicUtils.unwrap_const(x) : x
function coef_string(c, T::TensorExpression)
    c = Symbolics.simplify(Symbolics.expand(c))
    c = _unwrapnum(c)  
    isequal(c, 1) && !isconst(T) && return ""
    isequal(c, -1) && !isconst(T) && return "-"
    if c isa Number

        isreal(c) && return string(real(c))
        r, im_part = real(c), imag(c)
        if iszero(r)
            isequal(im_part, 1) && return "i"
            isequal(im_part, -1) && return "-i"
            s = string(im_part)
            return occursin('+', s) || occursin('-', s) && !startswith(s, "-") ?
                   "(" * s * ")i" :
                   s * "i"
        end
        return "(" * string(c) * ")"
    end
    # symbolic: check for complex numeric factor
    try
        args = _unwrapnum.(SymbolicUtils.arguments(c))     # <- unwrap each argument
        op = SymbolicUtils.operation(c)

        if op === (+)
            return "(" * join(coef_string.(args, Ref(T)), " + ") * ")"
        end

        if op === (*)
            num_args = filter(a -> a isa Number, args)
            sym_args = filter(a -> !(a isa Number), args)

            if !isempty(num_args)
                num = prod(num_args)
                sym = isempty(sym_args) ? 1 :
                      length(sym_args) == 1 ? sym_args[1] :
                      op(sym_args...)

                prefix = coef_string(num, T)
                return prefix == "" ? string(sym) :
                       prefix == "-" ? "-" * string(sym) :
                       prefix * string(sym)
            end
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

# function sort_exp(TEx::TensorExpression, Order::Vector{T} where {T<:TensorSuperType})
#     for summand in TEx.summands



#     end
#     return TEx
# end

# function Base.sort(T::TensorTerm)
#     return sort(TensorExpression(T))
# end

Base.sort(T::TensorSuperType) = sort(TensorExpression(T))
function Base.sort(TEx::TensorExpression)::TensorExpression
    new_exp = zero()
    if isZeroTens(TEx)
        return new_exp
    end
    for summand in TEx.summands
        coeff = summand[1]
        new_terms = copy(summand[2])
        fac = 1
        for term in terms(new_terms)
            fac *= sort_symmetric!(term)[1]#indices to canonical order, edit in place and extract factors.
        end
        (fac2, expout) = bubble_sort_with_grading!(new_terms)
        new_exp += coeff * fac * expout
    end
    return new_exp
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
    new_terms = copy(term)
    fac = 1
    for term in terms(new_terms)
        fac *= sort_symmetric!(term)[1]#indices to canonical order, edit in place and extract factors.
    end
    (fac2, expout) = bubble_sort_with_grading!(new_terms)
    return fac * expout
end




function is_less_mod_dummies(T1::TensorSuperType, T2::TensorSuperType, dummies)
    try
        if dummies === nothing
            return T1 < T2
        elseif length(dummies) == length(indices(T1))
            i1 = indices(T1)
            i2 = indices(T2)

            for i = eachindex(i1)
                if i1[i][1].index_set.default_position == i1[i][2] && i2[i][1].index_set.default_position != i2[i][2]
                    return true
                elseif i1[i][1].index_set.default_position != i1[i][2] && i2[i][1].index_set.default_position == i2[i][2]
                    return false
                elseif i1[i][1] == i2[i][1]
                    continue
                else
                    return i1[i][1] < i2[i][1]
                end
            end
            return false
        else
            # println("in dumms")
            # println(dummies)
            i1 = indices(T1)
            i2 = indices(T2)
            n = min(length(i1), length(i2))#may not have same lengths!
            for i = 1:n
                if i1[i][1] in dummies
                    if i2[i][1] in dummies
                        continue
                    else
                        return false

                    end
                elseif i2[i][1] in dummies
                    return true
                else
                    if i1[i][1].index_set.default_position == i1[i][2] && i2[i][1].index_set.default_position != i2[i][2]
                        return true
                    elseif i1[i][1].index_set.default_position != i1[i][2] && i2[i][1].index_set.default_position == i2[i][2]
                        return false
                    elseif i1[i][1] == i2[i][1]
                        continue
                    else
                        return i1[i][1] < i2[i][1]
                    end
                end
            end
            # All common entries were equal or dummies
            return length(i1) < length(i2)
        end
        return false
    catch e
        @show T1
        @show T2
        @show i1 = indices(T1)
        @show i2 = indices(T2)
        rethrow(e)
    end
end

# function is_less_mod_dummies(T1::TensorSuperType, T2::TensorSuperType, dummies)
#     i1 = indices(T1)
#     i2 = indices(T2)

#     # helper: unified "is default position?"
#     is_default(idx) = idx[2] == idx[1].index_set.default_position

#     for i in eachindex(i1)

#         a = i1[i][1]
#         b = i2[i][1]

#         # unified dummy/default priority
#         pa = (dummies !== nothing && a in dummies) || !is_default(i1[i]) ? 1 : 0
#         pb = (dummies !== nothing && b in dummies) || !is_default(i2[i]) ? 1 : 0

#         # 1. primary rule: default < non-default
#         if pa != pb
#             return pa < pb
#         end

#         # 2. tie-breaker: structural equality
#         if a != b
#             return a < b
#         end
#     end

#     return false
# end

function custom_order_isless(T1::TensorSuperType, T2::TensorSuperType, dummies=nothing)
    if haskey(sort_order, name(T1))
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
Check to see if tensors are matrices, and if they can commute.
"""
function can_swap(a::TensorSuperType, b::TensorSuperType)
    (a isa MatrixSuperType && b isa MatrixSuperType) && return !haskey(commutators, minmax(name(a), name(b)))
    return true
end

function bubble_sort_with_grading!(T::TensorTerm)::Tuple{Int,TensorExpression}
    dums = dummies(T)
    v = terms(T)
    n = length(v)
    swapped = true
    swaps = 0
    out_expression = zero()
    while swapped
        swapped = false
        for i in 1:(n-1)
            if !isempty(sort_order)
                if can_swap(v[i], v[i+1]) && custom_order_isless(v[i+1], v[i], dums) #i want to add an ignore dummy names in sort order
                    v[i], v[i+1] = v[i+1], v[i]
                    swaps += grading(v[i]) * grading(v[i+1])
                    swapped = true
                elseif custom_order_isless(v[i+1], v[i], dums) #need to apply the commutator and generate a new term
                    first_half = v[1:(i-1)] #should be same
                    final_half = v[(i+2):end]
                    (grade, comm) = get_commutator(v[i], v[i+1])#should get the commutator

                    # @show (-1)^(swaps)*TensorTerm(first_half) * comm * TensorTerm(final_half) #need actually to consider the current num swaps
                    out_expression += (-1)^(swaps)*sort(TensorTerm(first_half) * comm * TensorTerm(final_half))#recursively sort the rest
                    # @show out_expression
                    v[i], v[i+1] = v[i+1], v[i]
                    swaps += grading(v[i]) * grading(v[i+1]) + (grade < 0) #extra grade from operators
                    swapped = true
                end
            elseif can_swap(v[i], v[i+1]) && custom_order_isless(v[i+1], v[i], dums)###make always
                v[i], v[i+1] = v[i+1], v[i]
                swaps += grading(v[i]) * grading(v[i+1])
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


function to_latex(T::TensorTerm)
    out = ""
    for tens in terms(T)
        out *= to_latex(tens)
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

function to_latex(T::TensorExpression)
    out = ""

    summandss = summands(T)
    n = length(summandss)
    isempty(summandss) && return "0"
    cs = latex_coeff(coef_string(summandss[1][1], T))
    for (i, tens) in enumerate(summandss)
        out*=cs*" "
        out*=to_latex(tens[2])
        if i < n
            cs = coef_string(summandss[i+1][1], T)
            if !startswith(cs, "-")
                # @show cs
                out *= "+"
            end
            cs=latexify(cs)
        end
    end
    return Base.replace(out * " ", "\$" => "")
end

function uni_to_TeX(charr::String)
    out = ""
    # dot=false
    #   println(charr[1])
    #   println(codepoint(charr[1]))
    if charr[1] == '\u0307'
        out *= "\\dot{"
        out *= Base.replace(latexify(charr[3:end]), "#\"L" => "")
        out *= "}"
        return out
    elseif codepoint(charr[1]) == 772
        out *= "\\overline{"
        out *= Base.replace(latexify(charr[3:end]), "#\"L" => "")
        out *= "}"
        return out
    else
        # println("not dotted")
        return Base.replace(latexify(charr), "#\"L" => "")
    end
end

# function to_latex_basic(T::TensorSuperType)
#     out = ""
#     out *= (T.name) * " "
#     if !isempty(indices(T))
#         lowers = ""
#         uppers = ""
#         lowersempt::Bool = 1
#         uppersempt::Bool = 1
#         for (ind, pos) in indices(T)
#             if ind.name == "±"
#                 ntemp = "\\pm"
#             elseif ind.name == "∓"
#                 ntemp = "\\mp"
#             else
#                 ntemp = (ind.name)
#             end
#             if pos.position == 0
#                 if typeof(ind) == DummyPatternIndex
#                     uppers *= "?"
#                 end

#                 uppers = uppers * (ntemp) * " "
#                 uppersempt = 0

#             else
#                 lowersempt = 0
#                 if typeof(ind) == DummyPatternIndex
#                     lowers *= "?"
#                 end
#                 lowers = lowers * (ntemp) * " "
#             end
#         end
#         if !uppersempt
#             out *= "^{" * uppers * "}"
#         end
#         if !lowersempt
#             out *= "_{" * lowers * "}"
#         end
#     end
#     return Base.replace(out * " ", "\$" => "")
# end

function to_latex(T::TensorSuperType)
    out = ""
    if name(T) != ""
        out *= uni_to_TeX(name(T)) * " "
    end
    if !isempty(indices(T))
        for (ind, pos) in indices(T)
            if ind.name == "±"
                ntemp = "\\pm"
            elseif ind.name == "∓"
                ntemp = "\\mp"
            else
                ntemp = uni_to_TeX(ind.name)
            end
            if pos.position == 0

                out *= "{ }^{"

                if typeof(ind) == DummyPatternIndex
                    out *= "?"
                end

                # uppers = uppers * (ntemp) * " "
                # uppersempt = 0
                out *= ntemp
                out *= "}"
            else
                # lowersempt = 0

                out *= "{ }_{"

                if typeof(ind) == DummyPatternIndex
                    out *= "?"
                end
                # lowers = lowers * (ntemp) * " "

                out *= ntemp
                out *= "}"
            end
        end

    end

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


function tens2str(T::TensorSuperType)
    out = ""
    if typeof(T) == DummyPatternTensor
        # out*="?"
    end

    if T isa Gamma5
        out *= first(name(T))
    else
        out *= name(T)
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
                uppers = uppers * ind.name
                uppersempt = 0
                lowers = lowers * "∘"
            else
                uppers = uppers * "∘"
                lowersempt = 0
                if typeof(ind) == DummyPatternIndex
                    lowers *= "?"
                end
                lowers = lowers * ind.name
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
            out *= "(" * join(names, ",") * ")"
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

function compute_normalization_weighted(factors::Vector{String})
    total_weight = prod(factor_weight(f) for f in factors)
    return 1 // total_weight
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


function apply_factorized_projector(factors::Vector{String}, tensor::String)
    terms = [(one(Rational{Int}), collect(tensor))]

    for f in reverse(factors)
        terms = apply_string_factor(terms, f)
    end

    result = Dict{Vector{Char},Float64}()
    for (s, t) in terms
        result[t] = get(result, t, 0.0) + s
    end

    return [(s, t) for (t, s) in result if abs(s) > 1e-12]
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

    # -------------------------
    # TRUE Young symmetriser:
    # product of:
    #   column antisym sums
    #   row sym sums
    # -------------------------
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
symcopy(t::TensorSuperType) = typeof(t)(symcopy(data(t)))
# the field types that currently fall through to deepcopy:
symcopy(v::Vector{<:Pair}) = copy(v)                    # indices, labels — Pairs are immutable
symcopy(d::Dict) = copy(d)                              # weights
# symcopy(s::OrderedSet{Coordinate}) = copy(s)            # function_of — Coordinates shared via registry
symcopy(s::Set{Coordinate}) = s
function symcopy(d::TensorData)
    return TensorData(
        d.name,
        d.indices,     # only this is mutated ever...
        d.function_of,       # shared
        d.weights,           # shared
        d.tableaux,          # shared
        d.labels,            # shared
    )
end
symcopy_owned(d::TensorData) = begin
    #  NCOPY[]+=1 
    return TensorData(d.name, copy(d.indices), d.function_of, d.weights, d.tableaux, d.labels)
end
#dont need to check over
# symcopy(t::TensorExpression) = TensorExpression(_RAW, CoeffPair[coef => symcopy(term) for (coef, term) in summands(t)])
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
    if isempty(terms(T))
        return TensorTerm(Vector{TensorSuperType}())
    end
    return TensorTerm([copy(t) for t in terms(T)])
end
copy(T::TensorSuperType) = typeof(T)(symcopy(data(T)))
function copy(T::TensorExpression)
    out = TensorExpression()
    if length(terms(T)) == 1 && isZero(first_coef(T))
        return out
    end
    for summand in summands(T)
        out += summand[1] * copy(summand[2])
    end
    return out
end
"""
    copy_owned(x)

Like symcopy but the result has its own index vectors.
Safe to mutate in place without affecting the original.
"""
copy_owned(x) = symcopy(x)

copy_owned(d::TensorData) = begin
    return TensorData(
        d.name,
        copy(d.indices),                # OWNED
        d.function_of,
        d.weights,
        d.tableaux,
        d.labels,
        # d.uid,
    )
end
copy_owned(t::TensorSuperType) = typeof(t)(copy_owned(data(t)))
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

"""The  zero() tensor expression used globally in Alakazam."""
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
function collect_indices!(inds::Vector{IndexSuperType}, T)
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