
"""
    Creates a metric tensor with (-,+) signature. 
    By default, the signature is (1,d-1) where d is the range of the IndexSet of the attached indices, i.e. we use mostly-plus Minkowski convention unless otherwise specified.
    This is used in places where we take a trace of the metric, for example.
"""
mutable struct Metric <: SymMetricSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Metric(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    # signature::Pair{Int,Int} #added metric signature
    function Metric(name::String, args...)
        (t, tinds) = build_tensor(name, args...)
        if length(tinds) != 2
            throw(ArgumentError("Metric must have exactly 2 indices!"))
        end
        if tinds[1][2] != tinds[2][2]
            throw(ArgumentError("Metric indices must have the same level!"))
        end
        t = _add_tableau(t, generate_symmetric_tableau(tinds))
        return new(t, tinds)
    end
    # function Metric(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(), signature=nothing) where {I<:IndexSuperType}
    #     if length(indices) != 2
    #         throw(ArgumentError("Metric must have exactly 2 indices"))
    #     end
    #     new(name, indices, Vector{Coordinate}(), Dict{String,Number}(), generate_symmetric_tableau(indices))
    # end
    # function Metric(d::TensorData, signature::Union{Nothing,Pair{Int,Int}}=nothing)
    #     if signature === nothing
    #         signature = (1 => indices(d)[1][1].index_set.range - 1)
    #     end
    #     return new(d, signature)
    # end
    function Metric(id1::Pair{I1,IndexPosition}, id2::Pair{I2,IndexPosition}) where {I1<:IndexSuperType,I2<:IndexSuperType}
        return _metric_or_error(id1, id2)
    end
    function Metric(id::Vector{Pair{I1,IndexPosition}}) where I1<:IndexSuperType
        return _metric_or_error(id[1], id[2])
    end
end


mutable struct EpsilonTensor <: EpsilonTensorSuperType
    data::TensorData
    indices::Vector{IndexPair}
    EpsilonTensor(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    EpsilonTensor(args...) = EpsilonTensor("ϵ", args...)
    function EpsilonTensor(name::String, args...)
        (t, tinds) = build_tensor(name, args...)
        i = tinds
        #the below enforces only two index epsilon tensor. 
        # if i[1][2] != i[2][2] || i[1][1].index_set != i[2][1].index_set # length(i) != 2 || #removed the length check for now
        #     throw(InvalidEpsilonIndices()) #should be both upper or both lower
        # end
        #any length of indices
        #Levels must be uniform
        if length(unique(idx.index_set for (idx, pos) in i)) != 1
            throw(InvalidEpsilonIndices(i, "all indices must belong to the same index set"))
        elseif length(unique(pos for (idx, pos) in i)) != 1
            throw(InvalidEpsilonIndices(i, "all indices must be in the same position (all raised or all lowered)"))
        elseif length(i) != i[1][1].index_set.range
            throw(InvalidEpsilonIndices(i,
                "an epsilon tensor on index set '$(i[1][1].index_set.name)' requires " *
                "$(i[1][1].index_set.range) indices (the dimension of the set) in order to be a top (volume) form - got $(length(i)) indices."))
        end
        t = _add_tableau(t, generate_antisymmetric_tableau(tinds))
        return new(t, tinds)
    end
end


mutable struct SymplecticForm <: EpsilonTensorSuperType
    data::TensorData
    indices::Vector{IndexPair}
    SymplecticForm(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    SymplecticForm(args...) = SymplecticForm("ω", args...)
    function SymplecticForm(name::String, args...)
        (t, tinds) = build_tensor(name, args...)
        i = tinds
        #ensure all indices are from same index set, all upper, and that even number of indices
        if length(unique(idx.index_set for (idx, pos) in i)) != 1
            throw(InvalidEpsilonIndices(i, "all indices must belong to the same index set."))
        elseif length(unique(pos for (idx, pos) in i)) != 1
            throw(InvalidEpsilonIndices(i, "all indices must be in the same position (all raised or all lowered)."))
        elseif length(i) % 2 != 0
            throw(InvalidEpsilonIndices(i, "a sympletic form requires an even number of indices, got $(length(i))."))
        end
        t = _add_tableau(t, generate_antisymmetric_tableau(tinds))
        return new(t, tinds)
    end
end



mutable struct KroneckerDelta <: KroneckerSuperType
    data::TensorData
    indices::Vector{IndexPair}
    KroneckerDelta(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    KroneckerDelta(args...) = KroneckerDelta("δ", args...)
    function KroneckerDelta(name::String, args...)
        (t, tinds) = build_tensor(name, args...)
        i = tinds
        if i[1][2] == i[2][2] || i[1][1].index_set != i[2][1].index_set #should be one upper and one lower #should be from the same index set
            throw(InvalidKroneckerIndices(i))
        end
        t = _add_tableau(t, generate_symmetric_tableau(tinds))

        return new(t, tinds)
    end
end

const Delta = KroneckerDelta


##################
#methods



##############################################################################################################################
##############################################################################################################################
### metric ops
##############################################################################################################################
###################################################################################################################################


"""Create a metric of the correct type, based on the type of indices. This could probably be moved to the constructor... (maybe - ϵ)"""
function generate_metric(id1::Pair{I1,IndexPosition}, id2::Pair{I2,IndexPosition}) where {I1<:IndexSuperType,I2<:IndexSuperType}
    # @show id1, id2
    #the set's metric_function, so every metric built here agrees with a user's g(x)
    fo = metric_function(id1[1].index_set)
    if id1[1].index_set.metric_type===nothing
        metric=nothing
    elseif id1[1].index_set.metric_type <: SymMetricSuperType
        if !isempty(fo) #varies, so must never be mistaken for a constant η or δ
            metric = id1[1].index_set.metric_type("g", [id1, id2], fo)
        elseif id1[1].index_set.signature[1] == 0 #euclidean signature
            metric = id1[1].index_set.metric_type("δ", [id1, id2])
        else
            metric = id1[1].index_set.metric_type("η", [id1, id2])
        end
    elseif id1[1].index_set.metric_type <: AntisymMetricSuperType ### should be careful about which inds to pop... should perhaps check grading etc!!!
        metric = id1[1].index_set.metric_type([id1, id2], fo)
    else #if hermitian etc, need to be more careful, maybe don't want to pop
        metric = nothing
    end
    return metric
end


"""
Whether indices of `S` can be raised and lowered: exactly when `generate_metric` builds a metric for it.
A set with no metric (SU(N≥3), whose ε is invariant but not a metric) cannot, nor one whose metric is
neither symmetric nor antisymmetric.
"""
_can_raise_lower(S::IndexSet) = S.metric_type !== nothing &&
    (S.metric_type <: SymMetricSuperType || S.metric_type <: AntisymMetricSuperType)

"""
The metric on two indices, for a caller that asked for one by name.
"""
function _metric_or_error(id1, id2)
    m = generate_metric(id1, id2)
    m === nothing || return m
    S = id1[1].index_set
    S.metric_type === nothing && throw(ArgumentError(
        "index set '$(S.name)' has no metric, so its indices cannot be raised or lowered. An upper " *
        "index still contracts with a lower one directly; to raise and lower, declare a metric type " *
        "on the set."))
    throw(ArgumentError("index set '$(S.name)' has a metric of type $(S.metric_type), which cannot " *
        "be built from a pair of indices."))
end

""""Pop" metrics off all tensors, putting indices in default position and returning the popped metric as a prefactor. This is used to put tensors in a canonical form for simplification and comparison.
"""
function pop_metric(TEx::TensorExpression)
    out = CoeffPair[]
    sizehint!(out, length(summands(TEx)))
    for summand in summands(TEx)
        append!(out, summands(TensorExpression(summand[1] * pop_metric(summand[2]))))
    end
    isempty(out) && return zero()
    return update_indices!(TensorExpression(out))
end

""""Pop" metrics off all tensors, putting indices in default position and returning the popped metric as a prefactor. This is used to put tensors in a canonical form for simplification and comparison.
"""
function pop_metric(tt::TensorTerm)
    # Which indices are available depends only on what the term already uses plus what the metrics popped so far have taken.
    claimed = IndexTypes[]
    collect_indices!(claimed, tt)
    new_term = 1
    for term in terms(tt)
        new_term *= _pop_metric(term, claimed)
    end
    return update_indices!(new_term)
end

"""
"Pop" metrics off a tensors, putting indices in default position and returning the popped metric as a prefactor. This is used to put tensors in a canonical form for simplification and comparison.
Currenltly, non-spacetime indices and raised and lowers with eps tensors, which is not necessarily what we want in the long run, but should work for now. 
This is because we have SU(2) flavour in mind. Need to decide the best way of implementing group metrics.
"""
# function pop_metric(T::TensorSuperType, orig_exp::Union{TensorExpression,TensorTerm,Nothing}=nothing, running_term::Union{TensorExpression,TensorTerm,Nothing}=nothing)# makes all indices on the object go to the default position according to the index set
#     if T isa AntisymMetricSuperType || T isa SymMetricSuperType || T isa OperatorSuperType ||T isa KroneckerSuperType
#         return T
#     end
#     T = copy_owned(T)
#     tensor = 1
#     new_index_pattern = indices(T)
#     used_indices = Set{IndexSuperType}()
#     for pos in eachindex(new_index_pattern)
#         index = new_index_pattern[pos]
#         isLightcone(index)&& continue
#         if !in_default_pos(index)
#             ###println("resetting to default positoin")
#             if orig_exp === nothing
#                 spare = get_unused_indices(T, index[1].index_set)
#             else
#                 spare = get_unused_indices(orig_exp, index[1].index_set)
#                 if running_term !== nothing
#                     spare = intersect(spare, get_unused_indices(running_term, index[1].index_set))
#                 end
#             end
#             spare = setdiff(spare, used_indices)
#             # if isempty(spare)
#             #     @show index[1].index_set
#             #     @show objectid(index[1].index_set)
#             #     nid = generate_indices(index[1].index_set)
#             #     @show objectid(nid.index_set)
#             #     @show  index[1].index_set
#             #     @show nid
#             #     println("gend new index")
#             #     @show collect(UInt32, nid.name)
#             #     nid === nothing && throw(OutOfIndices(index[1].index_set))
#             #     push!(spare, nid)
#             # end
#             use_this = safe_pop!(spare,index[1].index_set)
#             dummy = (use_this => (index[1].index_set.default_position == upper ? lower : upper))
#             push!(used_indices, use_this)
#             push!(used_indices, index[1])
#             # @show index,dummy
#             metric = generate_metric(index, dummy)
#             if metric === nothing
#                 continue
#             end
#             # if dummy[1].index_set.metric_type isa SymMetricSuperType
#             #     if dummy[1].index_set.signature[1]==0 #euclidean signature
#             #         metric = Metric("δ", [index, dummy])
#             #     else
#             #         metric = Metric("η", [index, dummy])
#             #     end
#             # elseif dummy[1].index_set.metric_type isa AntisymMetricSuperType ### should be careful about which inds to pop... should perhaps check grading etc!!!
#             #     metric = EpsilonTensor([index, dummy])
#             # else #if hermitian etc, need to be more careful, maybe don't want to pop
#             #     continue
#             # end
#             tensor = tensor * metric
#             new_index_pattern[pos] = dummy[1] => index[1].index_set.default_position
#         end
#     end
#     # set_indices!(T,new_index_pattern) #check if really need this...
#     return tensor * T
# end

function pop_metric(T::TensorSuperType,
    orig_exp::Union{TensorExpression,TensorTerm,Nothing}=nothing,
    running_term::Union{TensorExpression,TensorTerm,Nothing}=nothing)
    claimed = IndexTypes[]
    collect_indices!(claimed, orig_exp === nothing ? T : orig_exp)
    running_term === nothing || collect_indices!(claimed, running_term)
    return _pop_metric(T, claimed)
end

function _fresh_index(iset::IndexSet, claimed::Vector{IndexTypes})
    for i in iset.indices
        i isa DummyPatternIndex && continue
        i in claimed && continue
        return i
    end
    nid = generate_indices(iset)
    nid === nothing && throw(OutOfIndices(iset))
    return nid
end

"""
`target(set)` is the level each index is moved to: a tensor's set's default position unless given. A
derivative asks for its own indices with its own level, opposite that; any other operator is skipped.
"""
function _pop_metric(T::TensorSuperType, claimed::Vector{IndexTypes}, target=nothing)
    # skip types
    (T isa AntisymMetricSuperType || T isa SymMetricSuperType ||
     (T isa OperatorSuperType && target === nothing) || T isa KroneckerSuperType) && return T
    target === nothing && (target = S -> S.default_position)
    # need popping?
    inds = T isa OperatorSuperType ? op_indices(T) : indices(T)
    needs = false
    @inbounds for p in inds
        iset = p[1].index_set          # one load; islightcone and default_position both live here
        iset.islightcone && continue
        if p[2] != target(iset)
            needs = true
            break
        end
    end
    needs || return T                      # no copy needed etc to save
    tensor = 1
    rep = Dict{IndexPair,IndexPair}()
    @inbounds for index in inds
        iset = index[1].index_set
        iset.islightcone && continue
        index[2] == target(iset) && continue
        use_this = _fresh_index(iset, claimed)
        dummy = use_this => op_pos(target(iset))
        push!(claimed, use_this)
        metric = generate_metric(index, dummy)
        metric === nothing && continue
        tensor = tensor * metric
        rep[index] = dummy[1] => target(iset)
    end
    isempty(rep) && return T
    #copy on write: a fresh shell whose indices are replaced, never edited, so whatever shares them
    #is untouched; set_indices_deep! also refreshes an operator's cached full index list
    T = copy(T)
    set_indices_deep!(T, rep)
    return tensor * T
end

"""A metric of the kind popping produces, as opposed to any other factor it leaves behind."""
_is_popped_metric(g) = g isa SymMetricSuperType || g isa AntisymMetricSuperType

"""
Pop the metrics off a derivative and its operand. Each is popped next to its tensor, inside the
derivative, and carried out past it only where `commutes_with_metric` allows.
"""
function _pop_metric(O::DerivativeSuperType, claimed::Vector{IndexTypes})
    s = summands(operand(O))
    outside = TensorSuperType[]
    D = if length(s) != 1
        _pop_metric_per_summand(O, claimed)
    else
        c, tt = s[1]
        inside = TensorSuperType[]
        popped = false
        #a metric written in the operand may come out as one popped there does, unless it is all
        #the operand holds: moving that would leave a derivative of a number
        body = any(t -> !_is_popped_metric(t), terms(tt))
        for f in terms(tt)
            p = _pop_metric(f, claimed)
            if p === f
                if body && _is_popped_metric(f) && commutes_with_metric(O, f)
                    push!(outside, f)
                    popped = true
                else
                    push!(inside, f)
                end
                continue
            end
            popped = true
            #a nested derivative's result carries whatever it let out; those may go further still
            for g in terms(isTensorTerm(p))
                _is_popped_metric(g) && commutes_with_metric(O, g) ? push!(outside, g) : push!(inside, g)
            end
        end
        if popped
            Oc = copy(O)
            set_operand!(Oc, c * TensorTerm(inside))
            Oc
        else
            O
        end
    end
    D isa PartialDerivative || (D = _pop_metric(D, claimed, S -> op_pos(S.default_position)))
    return isempty(outside) ? D : TensorTerm(outside) * D
end

"""
With several summands in the operand no metric is common to all of them, so each summand has its
metrics popped where it stands and none leave the derivative.
"""
function _pop_metric_per_summand(O::DerivativeSuperType, claimed::Vector{IndexTypes})
    n = length(claimed)
    taken = IndexTypes[]
    new_operand = zero()
    for (c, tt) in summands(operand(O))
        #a dummy is private to its summand, so every summand may reuse the same fresh names
        trial = copy(claimed)
        new_term = 1
        for f in terms(tt)
            new_term *= _pop_metric(f, trial)
        end
        new_operand += c * new_term
        for i in @view trial[(n + 1):end]
            i in taken || push!(taken, i)
        end
    end
    isempty(taken) && return O
    append!(claimed, taken)
    Oc = copy(O)
    set_operand!(Oc, new_operand)
    return Oc
end

"""
Returns true when the metric has two indices in the same level.
False for KroneckerDelta.
"""
function is_riemann_metric(metric::MetricSuperType)
    # met_pos1 = indices(metric)[1][2].position
    # met_pos2 = indices(metric)[2][2].position
    sumind = (met_pos(metric, 1) + met_pos(metric, 2)) % 2
    return sumind == 0
end
function is_riemann_metric(metric::KroneckerSuperType)
    return false
end

signature(metric::MetricSuperType) = indices(metric)[1][1].index_set.signature

"""Returns whether index in slot is upper or lower"""
met_pos(metric::MetricSuperType, slot::Int) = indices(metric)[slot][2].position

""""
Trivial contraction. Add nested operator version.
"""
contract_metrics(T::TensorSuperType) = T

""""
Contract all metrics of a particular type in an expression, whereever possible.
"""
# function contract_metrics(T::TensorTerm, typ::Type=Metric)
#     contract_metrics(1 * T, typ)
# end
function contract_metrics(T::TensorTerm, typ::Type=Metric)
    T = _contract_in_operands(T, typ)
    T === nothing && return zero()
    (sign, factor, extra_deltas, newterm) = _contract_metrics_term(T, typ)
    return sign * factor * (extra_deltas * newterm)
end

"""
The term with the metrics inside each operator's operand contracted, recursively. An operator whose
operand has nothing to contract is kept as it is, and so is the term if none has. `nothing` when an
operand contracts to zero, since the term then vanishes.
"""
function _contract_in_operands(tt::TensorTerm, typ::Type)
    ts = terms(tt)
    any(t -> t isa OperatorSuperType, ts) || return tt
    changed = false
    newts = Vector{TensorSuperType}(undef, length(ts))
    for (k, t) in enumerate(ts)
        if t isa OperatorSuperType && !isZeroTens(operand(t))
            op = operand(t)
            newop = contract_metrics(op, typ)
            isZeroTens(newop) && return nothing
            if newop !== op
                #as applying the operator does: set_operand! also refreshes the cached indices
                t = copy(t)
                set_operand!(t, newop)
                changed = true
            end
        end
        newts[k] = t
    end
    return changed ? TensorTerm(newts) : tt
end


""""
Contract all metrics of a particular type in an expression, whereever possible.
"""
# function contract_metrics(T::TensorExpression, typ::Type=Metric)
#     found = false
#     for s in summands(T), t in terms(s[2])
#         if t isa typ
#             found = true;
#             break
#         end
#     end
#     found || return T
#     T = copy_owned(T)
#     resultout = zero()
#     for k in eachindex(summands(T))
#         summand = summands(T)[k]
#         terms = summand[2].terms
#         metrics, non_metrics = let ms = eltype(terms)[], nms = eltype(terms)[]
#             for x in terms
#                 push!(x isa typ ? ms : nms, x)
#             end
#             ms, nms
#         end
#         uncontr_mets = Vector{TensorSuperType}() #anything that isnt contracted will be put here
#         sign = 1 #for antisym metrics (spinor)
#         factor = 1 #(for contracted deltas etc)
#         extra_deltas = 1
#         to_remove = Set{Int}()
#         #what we do here is look at the remaining and see if they can be special contracted
#         i = 1
#         while i <= length(metrics)
#             j = 1
#             while j < i
#                 if (j in to_remove)
#                     j += 1
#                     continue
#                 end
#                 metric = metrics[i]
#                 tens = metrics[j]
#                 #need case by case
#                 if metric isa EpsilonTensorSuperType # see if they have same indices, and then contract
#                     (elim, fct, ed) = reduce_epsilon(metric, tens)
#                     factor *= fct
#                     extra_deltas *= ed
#                     if elim
#                         push!(to_remove, i)
#                         push!(to_remove, j)
#                         j += 1
#                         break
#                     end
#                 elseif metric isa KroneckerSuperType
#                     (elim, fct, ed) = reduce_kronecker(metric, tens)
#                     factor *= fct
#                     append!(metrics, ed)
#                     if elim
#                         push!(to_remove, i)
#                         push!(to_remove, j)
#                         j += 1
#                         break
#                     end
#                 elseif metric isa Metric
#                     (elim, fct, ed) = reduce_metric(metric, tens)
#                     factor *= fct
#                     extra_deltas *= ed
#                     if elim
#                         push!(to_remove, i)
#                         push!(to_remove, j)
#                         j += 1
#                         break
#                     end
#                 end
#                 j += 1
#             end
#             i += 1
#         end
#         # deleteat!(metrics, sort!(collect(to_remove)))
#         if !isempty(to_remove)
#             keep = 0
#             for i in eachindex(metrics)
#                 i in to_remove && continue
#                 keep += 1;
#                 metrics[keep] = metrics[i]
#             end
#             resize!(metrics, keep)
#         end
#         for metric in metrics #for each in the found metrics
#             rie = is_riemann_metric(metric)#indices in same level?
#             flag = false # terminate when metric is contracted
#             index1, index2 = indices(metric)[1][1], indices(metric)[2][1]
#             # println("trying to contract",metric )
#             if metric isa KroneckerSuperType
#                 if isLightcone(index1)
#                     factor *= (index1 == index2) ? 1 : 0
#                     flag = true
#                     break #break because identically  zero()
#                 else
#                     if index1 == index2 # contracts same index kroneckers
#                         factor *= index1.index_set.range
#                         flag = true #something here breaks!!!!! 
#                     end
#                 end
#             elseif metric isa EpsilonTensorSuperType
#                 if length(indices(metric)) > 2 #skip over eps with many inds. this should avoid stuff breaking... 
#                     #actually complicated, since if num>2 no longer nec. bosonic...
#                     push!(uncontr_mets, metric)
#                     continue #skip
#                 end
#                 if index1 == index2 # contracts same index kroneckers
#                     factor *= 0
#                     flag = true
#                     break#break because identically  zero()
#                 end
#             end
#             temptens=nothing
#             for tens in non_metrics #try to contract with each other tensor
#                 temptens=tens
#                 flag && break #stop if already contracted
#                 metric isa EpsilonTensorSuperType && tens isa KroneckerSuperType && continue #dont raise/lower eps or kronecker
#                 ###########################################################################
#                 if tens isa OperatorSuperType
#                     #check outermost first!
#                     opinds=op_indices(tens)
#                     for indpairnum in eachindex(opinds)
#                         flag && break
#                         (flag, sgn) = contraction_replacement!(opinds, indpairnum, metric, rie)
#                         sign *= sgn
#                     end
#                     #need to find the part of the op that has the inds.
#                     op = isTensor(operand(tens))
#                     if !flag && !(op === nothing) #is a single tensor#operand isa TensorSuperType || operand isa TensorTerm || (operand isa TensorExpression && length(summands(operand)) == 1 && length(operand.summands[1][2].terms) == 1)
#                         while isTensor(op) isa OperatorSuperType #gets to the root #progressively checks
#                             oppinds=op_indices(isTensor(op))
#                             flag && break
#                             for indpairnum in eachindex(oppinds)
#                                 if flag #if contracted, escape this operator
#                                     break
#                                 end
#                                 (flag, sgn) = contraction_replacement!(oppinds, indpairnum, metric, rie) #true contract,
#                                 sign *= sgn

#                             end
#                             if flag #if contracted, escape, we're done
#                                 break
#                             end
#                             op = operand(isTensor(op))
#                             # operand = operand(operand.summands[1][2].terms[1])
#                         end
#                         flag && break
#                         #ran out of operators, try to contract with summands if we havent contracted yet
#                         for indpairnum in eachindex(indices(op))
#                             flag && break
#                             (flag, sgn) = contraction_replacement!(indices(op), indpairnum, metric, rie) #true contract,
#                             sign *= sgn
#                         end
#                         #ensure that changes propagate up the chain....
#                         # end

#                     end
#                 else###not an operator
#                     tins = indices(tens)
#                     for indpairnum in eachindex(tins)
#                         flag && break
#                         (flag, sgn) = contraction_replacement!(tins, indpairnum, metric, rie)
#                         sign *= sgn
#                     end
#                 end
#             end
#             flag || push!(uncontr_mets, metric)
#             temptens===nothing || update_indices!(temptens)
#         end
#         trm = extra_deltas * TensorTerm(vcat(uncontr_mets, non_metrics))
#         # update_indices!(summands(T)[k][2])
#         # summands(T)[k] = summands(T)[k][1] * sign * factor => summands(T)[k][2]
#         # @show summands(T)[k][1] , sign, factor, trm
#         resultout += summands(T)[k][1] * sign * factor * trm
#     end
#     return resultout
# end
function contract_metrics(T::TensorExpression, typ::Type=Metric)
    s = summands(T)
    isempty(s) && return T

    # inside operands first, innermost first: a metric left inside a derivative, as popping does
    # where it cannot pass, contracts there, ∇_a(g^{ab}∇_b u) = ∇_a∇^a u
    if any(p -> any(t -> t isa OperatorSuperType, terms(p[2])), s)
        newtts = [_contract_in_operands(tt, typ) for (_, tt) in s]
        #rebuilt only if something changed: the caller tells an untouched operand by identity
        if any(k -> newtts[k] !== s[k][2], eachindex(s))
            kept = CoeffPair[CoeffPair(s[k][1], newtts[k]) for k in eachindex(s) if newtts[k] !== nothing]
            isempty(kept) && return zero()
            T = TensorExpression(kept)
            s = summands(T)
        end
    end

    # cheap scan: is there anything of this type to contract?
    found = false
    @inbounds for k in eachindex(s)
        for t in terms(s[k][2])
            if t isa typ
                found = true;
                break
            end
        end
        found && break
    end
    found || return T

    # resultout = zero()
    # @inbounds for k in eachindex(s)
    #     coef, tterm = s[k]
    #     (sign, factor, extra_deltas, newterm) = _contract_metrics_term(tterm, typ)
    #     resultout =resultout ⨥ coef * sign * factor * (extra_deltas * newterm)
    # end
    out = CoeffPair[]
    sizehint!(out, length(s))
    @inbounds for k in eachindex(s)
        coef, tterm = s[k]
        (sign, factor, extra_deltas, newterm) = _contract_metrics_term(tterm, typ)
        contribution = coef * sign * factor * (extra_deltas * newterm)
        isZeroTens(contribution) && continue
        append!(out, summands(contribution))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
    return resultout
end

"""
Apply `metric`'s contraction to `tens`, descending through whatever operators it wraps.
Returns whether the metric was consumed, and the sign the replacement cost.
"""
function _contract_into_operand!(tens, metric, rie)
    t = isTensor(tens)
    t === nothing && return (false, 1)      #a sum, which this contraction does not descend into
    sign = 1
    if t isa OperatorSuperType
        opinds = op_indices(t)
        for k in eachindex(opinds)
            (flag, sgn) = contraction_replacement!(opinds, k, metric, rie)
            sign *= sgn
            flag && return (true, sign)
        end
        #into the operand only past a derivative the metric passes: g^{ab}∂_c(A_b) is ∂_c(A^a) only if ∂_c g = 0
        t isa DerivativeSuperType && !commutes_with_metric(t, metric) && return (false, sign)
        (flag, sgn) = _contract_into_operand!(operand(t), metric, rie)
        return (flag, sign * sgn)
    end
    tins = indices(t)
    for k in eachindex(tins)
        (flag, sgn) = contraction_replacement!(tins, k, metric, rie)
        sign *= sgn
        flag && return (true, sign)
    end
    return (false, sign)
end

function _contract_metrics_term(T::TensorTerm, typ::Type)
    ts = terms(T)
    any(t -> t isa typ, ts) || return (1, 1, 1, T)

    T = copy_owned(T)
    ts = terms(T)

    nm = 0
    for t in ts
        ;
        t isa typ && (nm += 1);
    end
    metrics = Vector{TensorSuperType}(undef, nm)
    non_metrics = Vector{TensorSuperType}(undef, length(ts) - nm)
    im_ = 1;
    inm = 1
    for t in ts
        if t isa typ
            metrics[im_] = t;
            im_ += 1
        else
            non_metrics[inm] = t;
            inm += 1
        end
    end

    uncontr_mets = TensorSuperType[]
    sign = 1
    factor = 1
    extra_deltas = 1
    to_remove = Set{Int}()

    i = 1
    while i <= length(metrics)
        j = 1
        while j < i
            if j in to_remove
                j += 1;
                continue
            end
            metric = metrics[i]
            tens = metrics[j]
            if metric isa EpsilonTensorSuperType
                (elim, fct, ed) = reduce_epsilon(metric, tens)
                factor *= fct;
                extra_deltas *= ed
                if elim
                    push!(to_remove, i);
                    push!(to_remove, j);
                    j += 1;
                    break
                end
            elseif metric isa KroneckerSuperType
                (elim, fct, ed) = reduce_kronecker(metric, tens)
                factor *= fct
                append!(metrics, ed)
                if elim
                    push!(to_remove, i);
                    push!(to_remove, j);
                    j += 1;
                    break
                end
            elseif metric isa Metric
                (elim, fct, ed) = reduce_metric(metric, tens)
                factor *= fct;
                extra_deltas *= ed
                if elim
                    push!(to_remove, i);
                    push!(to_remove, j);
                    j += 1;
                    break
                end
            end
            j += 1
        end
        i += 1
    end

    if !isempty(to_remove)
        keep = 0
        for i in eachindex(metrics)
            i in to_remove && continue
            keep += 1
            metrics[keep] = metrics[i]
        end
        resize!(metrics, keep)
    end

    for metric in metrics
        rie = is_riemann_metric(metric)
        flag = false
        index1, index2 = indices(metric)[1][1], indices(metric)[2][1]

        if metric isa KroneckerSuperType
            if isLightcone(index1)
                factor *= (index1 == index2) ? 1 : 0
                flag = true
                break
            elseif index1 == index2
                factor *= index1.index_set.range
                flag = true
            end
        elseif metric isa EpsilonTensorSuperType
            if length(indices(metric)) > 2
                push!(uncontr_mets, metric)
                continue
            end
            if index1 == index2
                factor *= 0
                flag = true
                break
            end
        end

        #the tensor the metric went into, whose cached indices then need refreshing: an operator's
        #op_indices are edited in place, but its full index list is a separate cache
        contracted = nothing
        for tens in non_metrics
            flag && break
            metric isa EpsilonTensorSuperType && tens isa KroneckerSuperType && continue
            #Do not let a metric raise or lower an epsilon's indices. Kronecker is exempt
            !(metric isa KroneckerSuperType) && tens isa EpsilonTensorSuperType && continue

            if tens isa OperatorSuperType
                (flag, sgn) = _contract_into_operand!(tens, metric, rie)
                sign *= sgn
            else
                tins = indices(tens)
                for indpairnum in eachindex(tins)
                    flag && break
                    (flag, sgn) = contraction_replacement!(tins, indpairnum, metric, rie)
                    sign *= sgn
                end
            end
            flag && (contracted = tens)
        end
        flag || push!(uncontr_mets, metric)
        contracted === nothing || update_indices!(contracted)
    end

    n = length(uncontr_mets) + length(non_metrics)
    allterms = Vector{TensorSuperType}(undef, n)
    @inbounds for i in eachindex(uncontr_mets)
        ;
        allterms[i] = uncontr_mets[i];
    end
    off = length(uncontr_mets)
    @inbounds for i in eachindex(non_metrics)
        ;
        allterms[off+i] = non_metrics[i];
    end

    return (sign, factor, extra_deltas, TensorTerm(allterms))
end

"""
    antisym_contraction(metric, slot) -> (indpair, sign)

Raising/lowering rule for an antisymmetric form.

    ϵ_{ab} F^{b} =  F_{a}                     (contraction with slot 2, sign 1)
    ϵ_{ab} F^{a} = -ϵ_{ba} F^{a} = -F_{b}     (contraction with slot 1, sign -1)

"""
antisym_contraction(::AntisymMetricSuperType, slot::Integer)= (slot == 2) ? 1 : -1


"""
Repeated logic for replacing tensor inds with ones from a metric. Modifies the tensor in question. Returns flag to remove and the sign accumulated
"""
function contraction_replacement!(tins, indpairnum, metric, isRiemannType::Bool)
    # @show metric indpairnum tins
    indpair = tins[indpairnum]
    sign = 1
    flag = false
    inds = indices(metric)
    m1 = inds[1][1]
    m1pos = inds[1][2]
    m2 = inds[2][1]
    m2pos = inds[2][2]
    issym = metric isa SymMetricSuperType
    isasym = metric isa AntisymMetricSuperType

    if indpair[1] == m2 #tensor has an index matching the 2nd index of the metric
        #some complicated logic to make sure indices go to the right position# this should let you contract Kronckers and metrics in a consistent way
        if isRiemannType # metric, not a kronecker
            if issym
                indpair = m1 => m1pos # raise/lower index #η_ab F^b=F_a
            elseif isasym
                indpair=m1=>m1pos
                sign = antisym_contraction(metric, 2)
            end
        else #TTHIS IS A KroneckerDelta
            indpair = m1 => indpair[2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
        end
        flag = true #this term has been contracted, break
    elseif indpair[1] == m1
        #some complicated logic to make sure indices go to the right position, based on symmetry of metric. if contracting with antisym, pick up a minus
        if isRiemannType
            if issym
                indpair = m2 => m2pos # raise/lower. η_baF^b=η_abF^b=F_a
            elseif isasym
                indpair= m2=>m2pos
                sign = antisym_contraction(metric, 1)
            end
        else #kronecker
            indpair = m2 => indpair[2] # keep index in same pos ( δ_i^j A^i= A^j for example)
        end
        flag = true #this term has been contracted, break
    end
    tins[indpairnum] = indpair
    return flag, sign
end


function reduce_metric(metric1::MetricSuperType, metric2::MetricSuperType)
    factor = 1
    extra_deltas = 1
    eliminated = false

    inds1 = indices(metric1)
    inds2 = indices(metric2)

    # must have opposite positions to contract
    inds1[1][2] == inds2[1][2] && return (eliminated, factor, extra_deltas)

    # find shared index with opposite positions
    upper_metric = inds1[1][2] == upper ? metric1 : metric2
    lower_metric = inds1[1][2] == upper ? metric2 : metric1

    u_inds = indices(upper_metric)
    l_inds = indices(lower_metric)

    shared_u = findfirst(p -> p[1] in free_indices_nopos(lower_metric), u_inds)
    shared_u === nothing && return (eliminated, factor, extra_deltas)

    shared_l = findfirst(p -> p[1] == u_inds[shared_u][1], l_inds)

    # remaining free indices form the new delta
    rem_u = u_inds[shared_u == 1 ? 2 : 1]
    rem_l = l_inds[shared_l == 1 ? 2 : 1]

    extra_deltas *= KroneckerDelta([rem_u[1] => rem_u[2], rem_l[1] => rem_l[2]])
    eliminated = true

    return (eliminated, factor, extra_deltas)
end

##############################################################################################################################
##############################################################################################################################
### KroneckerDelta
##############################################################################################################################
##############################################################################################################################

""""
Eliminate Kronecker deltas by contracting indices or taking traces.
"""
eliminate_kronecker(T::Tensorial) = eliminate_kronecker(TensorExpression(T))

""""
Eliminate Kronecker deltas by contracting indices or taking traces.
"""
eliminate_kronecker(T::TensorExpression) = contract_metrics(T, KroneckerSuperType)

function reduce_kronecker(metric::KroneckerSuperType, tens::KroneckerSuperType)
    factor = 1 #(for contracted deltas etc)
    extra_deltas = Vector{TensorSuperType}()
    eliminated = false
    if free_indices_nopos(metric) == free_indices_nopos(tens)
        factor *= indices(metric)[1][1].index_set.range
        isempty(free_indices_ordered(metric)) && (factor *= indices(metric)[1][1].index_set.range)
        eliminated = true
    elseif !isempty(intersect(free_indices_nopos(metric), free_indices_nopos(tens)))
        #shared index, turn to a kronecker via contraction
        if indices(metric)[1][2].position == indices(tens)[1][2].position  #same positions
            if indices(metric)[1][1] == indices(tens)[1][1] #eg eps^{ij} eps^{ik}
                push!(extra_deltas, KroneckerDelta([indices(metric)[2][1] => indices(metric)[2][2], indices(tens)[2][1] => indices(tens)[1][2]]))
            elseif indices(metric)[2][1] == indices(tens)[2][1] #2nd index same, #eg eps^{ij} eps^{kj}
                push!(extra_deltas, KroneckerDelta([indices(metric)[1][1] => indices(metric)[1][2], indices(tens)[1][1] => indices(tens)[2][2]]))
            elseif indices(metric)[1][1] == indices(tens)[2][1] #alt index same, #eg eps^{ij} eps^{ki}
                push!(extra_deltas, KroneckerDelta([indices(metric)[2][1] => indices(metric)[2][2], indices(tens)[1][1] => indices(tens)[1][2]]))
            elseif indices(metric)[2][1] == indices(tens)[1][1] #alt index same, #eg eps^{ij} eps^{jk}
                push!(extra_deltas, KroneckerDelta([indices(metric)[1][1] => indices(metric)[1][2], indices(tens)[2][1] => indices(tens)[2][2]]))
            end
        elseif indices(metric)[1][2].position != indices(tens)[1][2].position   #opposite positions 
            if indices(metric)[1][1] == indices(tens)[1][1]    #eg eps^{ij} eps_{ik}
                push!(extra_deltas, KroneckerDelta([indices(metric)[2][1] => indices(metric)[2][2], indices(tens)[2][1] => indices(tens)[2][2]]))
            elseif indices(metric)[2][1] == indices(tens)[2][1] #2nd index same, #eg eps^{ij} eps_{kj}
                push!(extra_deltas, KroneckerDelta([indices(metric)[1][1] => indices(metric)[1][2], indices(tens)[1][1] => indices(tens)[1][2]]))
            elseif indices(metric)[1][1] == indices(tens)[2][1] #alt index same, #eg eps^{ij} eps_{ki}
                push!(extra_deltas, KroneckerDelta([indices(metric)[2][1] => indices(metric)[2][2], indices(tens)[1][1] => indices(tens)[2][2]]))
            elseif indices(metric)[2][1] == indices(tens)[1][1] #alt index same, #eg eps^{ij} eps_{jk}
                push!(extra_deltas, KroneckerDelta([indices(metric)[1][1] => indices(metric)[1][2], indices(tens)[2][1] => indices(tens)[2][2]]))
            end
        end
        eliminated = true
    end
    return (eliminated, factor, extra_deltas)
end

##############################################################################################################################
##############################################################################################################################
### Eps Tensor
##############################################################################################################################

""""
Eliminate the Levi-Civita epsilon tensor by contracting indices or taking traces.
"""
eliminate_epsilon(T::TensorTerm) = eliminate_epsilon(TensorExpression(T))

""""
Eliminate the Levi-Civita epsilon tensor by contracting indices or taking traces.
"""
eliminate_epsilon(T::TensorExpression) = eliminate_kronecker(contract_metrics(T, EpsilonTensorSuperType))

"""
Contract a pair of epsilons, expanding only over the slots that are actually free.

    ε_{a₁…a_k c₁…c_{n-k}} ε^{a₁…a_k d₁…d_{n-k}} = sgn(det η) · k! · δ^{[d₁…d_{n-k}]}_{[c₁…c_{n-k}]}

A pair sharing k indices therefore needs a sum over (n-k)! permutations with k! in front,
not over all n!.
"""
function reduce_epsilon(tens::EpsilonTensorSuperType, metric::EpsilonTensorSuperType)
    # @show tens,metric
    inds_t = indices(tens)
    inds_m = indices(metric)
    n = length(inds_m)
    n != length(inds_t) && return (false, 1, 1)
    inds_t[1][1].index_set != inds_m[1][1].index_set && return (false, 1, 1)
   
    if n != inds_m[1][1].index_set.range
        n == 2 || return (false, 1, 1)      # higher-rank non-volume forms: leave alone
        slot_t = findfirst(r -> any(q -> q[1] == r[1], inds_m), inds_t)
        slot_t === nothing && return (false, 1, 1)
        shared = inds_t[slot_t][1]
        slot_m = findfirst(q -> q[1] == shared, inds_m)
        # a contraction needs the shared index raised on one form and lowered, if at the same level there is nothing to cancel
        inds_t[slot_t][2] == inds_m[slot_m][2] && return (false, 1, 1)
        free_t = inds_t[3-slot_t]
        free_m = inds_m[3-slot_m]
        free_t[2] == free_m[2] && return (false, 1, 1)
        sgn = (slot_t == slot_m) ? -1 : 1
        return (true, sgn,
                TensorExpression(KroneckerDelta([free_t[1] => free_t[2],
                                                 free_m[1] => free_m[2]])))
    end

    shared_t = Int[]
    shared_m = Int[]
    for (j, r_t) in enumerate(inds_t)
        i = findfirst(q -> q[1] == r_t[1] && q[2] != r_t[2], inds_m)
        if i !== nothing && !(i in shared_m)
            push!(shared_t, j)
            push!(shared_m, i)
        end
    end
    free_t = Int[j for j in 1:n if !(j in shared_t)]
    free_m = Int[i for i in 1:n if !(i in shared_m)]
    k = length(shared_t)
    #Writing each epsilon with its shared slots last costs the sign of that reordering.
    presign = perm_sign(vcat(free_t, shared_t)) * perm_sign(vcat(free_m, shared_m))

    extra_deltas = zero()
    if k == n
        extra_deltas += 1
    else
        for σ in permutations(1:(n - k))
            sgn = perm_sign(σ)
            deltas = TensorSuperType[]
            for l in 1:(n - k)
                r_t = inds_t[free_t[σ[l]]]
                r_m = inds_m[free_m[l]]
                if r_t[2] == r_m[2]
                    (r_t[1].index_set.metric_type===nothing || r_t[1].index_set.metric_type <: AntisymMetricSuperType) && return (false, 1, 1) #skip if metrics are more epsilons etc
                    # same position: use metric to lower/raise
                    if r_t[2] == upper
                        mt = generate_metric(r_t[1] => upper, r_m[1] => upper)
                        if mt === nothing
                            error("Invalid reduction of epsilons (index set not supported)")
                        end
                        push!(deltas, mt)
                    else
                        mt = generate_metric(r_t[1] => lower, r_m[1] => lower)
                        if mt === nothing
                            error("Invalid reduction of epsilons (index set not supported)")
                        end
                        push!(deltas, mt)
                    end
                else
                    push!(deltas, KroneckerDelta([r_t[1] => r_t[2], r_m[1] => r_m[2]]))
                end
            end
            extra_deltas += sgn * TensorTerm(deltas)
        end
    end
    #det of the metric.
    mt = inds_t[1][1].index_set.metric_type
    dt = det(inds_t[1][1]) * ((mt !== nothing && mt <: AntisymMetricSuperType) ? -1 : 1)
    return (true, dt * presign * factorial(k), extra_deltas)
end


