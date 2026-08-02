
"""
Basic Matrix object. Generically, is non-commuting. Commutation relations must be specified. 
It will (anti)commute with all tensor for which an explicit com-rel is not defined.
    """
mutable struct Matrix <: MatrixSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Matrix(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Matrix(name::String, args...)
        return new(build_tensor(name, args...)...)
    end
end

# mutable struct Ket <: MatrixSuperType
#     data::TensorData
#     function Ket(name::String, args...)
#         return new(build_tensor(name, args...)...)
#     end
#     function Ket(d::TensorData)
#         return new(d)
#     end
# end
# mutable struct Bra <: MatrixSuperType
#     data::TensorData
#     function Bra(name::String, args...)
#         return new(build_tensor(name, args...)...)
#     end
#     function Bra(d::TensorData)
#         return new(d)
#     end
# end


mutable struct Pauli <: PauliSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Pauli(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Pauli(name::String, args...)
        return new(build_tensor(name, args...)...)
    end
    function Pauli(args...)
        return new(build_tensor("σ", args...)...)
    end

end

mutable struct Gamma <: GammaSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Gamma(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Gamma(name::String, args...)
        return new(build_tensor(name, args...)...)
    end
    function Gamma(args...)
        return new(build_tensor("γ", args...)...)
    end

end

mutable struct Gamma5 <: Gamma5SuperType
    data::TensorData
    indices::Vector{IndexPair}
    Gamma5(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Gamma5(name::String, args...)
        return new(build_tensor(name, "5" => upper, args...)...)
    end
    function Gamma5(args...)
        return new(build_tensor("γ5", "5" => upper, args...)...)
    end

end


"""
    @matrix MatrixType

Macro to define a  new MatrixType<:MatrixSuperType, with basic constructors.
For example,
    @matrix GMMatrix
    A = GMMatrix("A", [μ => lower])
"""
macro matrix(name)
    name isa Symbol || error("@matrix expects a type name symbol, got $name")
    esc(quote
        mutable struct $name <: MatrixSuperType
            data::$(TensorData)
            indices::Vector{$(IndexPair)}
            $name(args...) = new($(build_tensor)(args...)...)
            $name(d::$(TensorData), i::Vector{$(IndexPair)}) = new(d, i)
        end
    end)
end

####
#reduce special objects
#####
function contract_sigma(T::Tensorial)
    result = contract_adjacent(T, PauliSuperType)
    return eliminate_kronecker(contract_metrics(eliminate_epsilon(result)))
end

"""
Generic reduce matrix function - can be overloaded to define new identities.
"""
reduce_matrix(m::MatrixSuperType, m2::MatrixSuperType, tterm::TensorTerm) = (false,  zero())

function reduce_matrix(σ1::PauliSuperType, σ2::PauliSuperType, tterm::TensorTerm)
    if length(indices(σ2)) == length(indices(σ1))#same #inds for contraction.
        if length(indices(σ1)) == 1#basic pauli
            refindex = indices(σ1)[1]
            if refindex[2] != indices(σ2)[1][2] #different pos
                return (false,  zero())
            end
            spare = get_unused_indices(tterm, refindex[1].index_set)
            newind = safe_pop!(spare,refindex[1].index_set)
            mt = generate_metric(indices(σ1)[1], indices(σ2)[1])
            if mt === nothing
                error("Reduction of index set not supported.")
            end
            metind=safe_pop!(spare,refindex[1].index_set)
            res = mt + im * EpsilonTensor(vcat(indices(σ1), indices(σ2), newind => refindex[2])) * generate_metric(newind => op_pos(refindex[2]),metind=>op_pos(refindex[2]))*Pauli(metind => (refindex[2]))
            return (true, res)
        elseif length(indices(σ1)) == 1#explicit spinor inds 

        end
    else
        return (false,  zero())
    end
end


"""
Moves Matrices of the given type to be adjacent, then contracts them when possible using simplification rules.
"""
# function contract_adjacent(To::Tensorial, target_type::Type)

#     T = 1 * To
#     T = copy(T)

#     resultout =  zero()
#     any_contraction = false
#     changed = true

#     while changed
#         # T = pop_metric(T)
#         changed = false
#         new_resultout =  zero()

#         for k in eachindex(summands(T))
#             summand = summands(T)[k]
#             all_terms = summand[2].terms
#             contracted = false
#             target_indices = findall(t -> t isa target_type, all_terms)
#             if length(target_indices) >= 2
#                 anchor = target_indices[1]
#                 n = length(target_indices)
#                 dragged_terms = drag_to_adjacent_n(summand[1], all_terms, target_indices)
#                 # @show dragged_terms
#                 for (c, trms) in summands(dragged_terms)
#                     contracted=false
#                     t = terms(trms)
#                     if all(t[anchor+m-1] isa target_type for m in 1:n)
#                         for r in n:-1:2
#                             try
                             
                                
#                                 (elim, result) = reduce_matrix(t[anchor:(anchor+r-1)]..., trms)
#                                 # @show result
#                                 if elim
#                                     new_resultout += c *
#                                                      TensorTerm(t[1:(anchor-1)]) *
#                                                      result *
#                                                      TensorTerm(t[(anchor+r):end])
#                                     contracted = true
#                                     any_contraction = true
#                                     changed = true
#                                     break
#                                 end

#                             catch e
#                                 # @show e
#                                 e isa MethodError || rethrow(e)
#                                 continue
#                             end
#                         end
#                     end

#                     # ONLY fallback if NOTHING in THIS summand contracted
#                     if !contracted
#                         new_resultout += c * trms
#                     end
#                 end

#             else
#                 new_resultout += summand[1] * TensorTerm(all_terms)
#             end
#         end

#         T = new_resultout
#         resultout = new_resultout
#     end
#     if !any_contraction
#         return To
#     end
#     return eliminate_kronecker(
#         contract_metrics(
#             eliminate_epsilon(resultout)
#         )
#     )
# end
"""
Moves matrices of the given type to be adjacent, then contracts them when possible
using simplification rules. When no identity applies to the full set of target
matrices, progressively smaller subsets are dragged adjacent and tried instead.
"""
function contract_adjacent(To::Tensorial, target_type::Type)
    T = copy(1 * To)

    resultout =  zero()
    any_contraction = false
    changed = true

    while changed
        changed = false
        new_resultout =  zero()

        for k in eachindex(summands(T))
            summand = summands(T)[k]
            all_terms = summand[2].terms
            target_indices = findall(t -> t isa target_type, all_terms)

            if length(target_indices) < 2
                new_resultout += summand[1] * TensorTerm(all_terms)
                continue
            end

            summand_contracted = false

            # Largest subsets first; the full set is tried before any fallback.
            for subset in subsets_decreasing(target_indices)
                anchor = subset[1]
                m = length(subset)

                dragged = drag_subset_adjacent(summand[1], all_terms, subset, target_type)

                candidate =  zero()
                subset_hit = false

                for (c, trms) in summands(dragged)
                    t = terms(trms)
                    run_ok = (anchor + m - 1 <= length(t)) &&
                             all(t[anchor + j - 1] isa target_type for j in 1:m)

                    reduced_here = false
                    if run_ok
                        try
                            (elim, result) = reduce_matrix(t[anchor:(anchor + m - 1)]..., trms)
                            if elim
                                candidate += c *
                                             TensorTerm(t[1:(anchor - 1)]) *
                                             result *
                                             TensorTerm(t[(anchor + m):end])
                                reduced_here = true
                                subset_hit = true
                            end
                        catch e
                            e isa MethodError || rethrow(e)
                        end
                    end

                    # commutator side-terms (or a failed reduction on this
                    # particular summand) pass through unchanged; they get
                    # re-scanned on the next `while changed` pass
                    if !reduced_here
                        candidate += c * trms
                    end
                end

                if subset_hit
                    # commit the dragged+reduced form only on success
                    new_resultout += candidate
                    summand_contracted = true
                    any_contraction = true
                    changed = true
                    break
                end
                # on failure: discard `candidate` entirely and try the next
                # subset starting from the ORIGINAL all_terms
            end

            if !summand_contracted
                new_resultout += summand[1] * TensorTerm(all_terms)
            end
        end

        T = new_resultout
        resultout = new_resultout
    end

    any_contraction || return To
    return eliminate_kronecker(
        contract_metrics(
            eliminate_epsilon(resultout)
        )
    )
end

"""
    drag_subset_adjacent(summand_coeff, all_terms, subset, target_type)

Generalization of `drag_to_adjacent_n` to an arbitrary (sorted, ascending) subset
of target positions. Drags each subset element, in order, to sit directly after
the first subset element (the anchor), accumulating commutator terms.

Elements are relocated inside each summand by their occurrence *rank* among
terms of `target_type`. Because the subset is sorted ascending, dragging an
element leftward only shifts the ranks of occurrences to its left, so the ranks
computed on the original `all_terms` remain valid for all later subset elements.
"""
function drag_subset_adjacent(summand_coeff, all_terms, subset::Vector{Int}, target_type::Type)
    anchor = subset[1]
    n = length(subset)

    all_occ = findall(t -> t isa target_type, all_terms)
    ranks = [findfirst(==(i), all_occ) for i in subset]

    term_list = summand_coeff * TensorTerm(all_terms)
    n == 1 && return term_list

    for k in 2:n
        target_pos = anchor + k - 1
        new_term_list =  zero()
        for (c, tterm) in summands(term_list)
            trms = terms(tterm)
            occurrences = findall(t -> t isa target_type, trms)
            # commutator side-terms may have consumed some targets: skip those
            if length(occurrences) < ranks[k]
                new_term_list += c * tterm
                continue
            end
            current_pos = occurrences[ranks[k]]
            if current_pos == target_pos
                new_term_list += c * tterm
                continue
            end
            new_term_list += drag_to_adjacent(c, trms, target_pos - 1, current_pos)
        end
        term_list = new_term_list
    end

    return term_list
end

"""
    subsets_decreasing(v)

All subsets of `v` of size >= 2, largest first. Within a size class, subsets are
ordered by span (fewest interleaved non-subset elements first), so cheaper drags
that generate fewer commutator terms are attempted before long-range ones.
"""
function subsets_decreasing(v::Vector{Int})
    n = length(v)
    out = Vector{Vector{Int}}()
    for r in n:-1:2
        group = Vector{Vector{Int}}()
        _combinations!(group, v, r)
        sort!(group, by = s -> s[end] - s[1])   # prefer tight subsets
        append!(out, group)
    end
    return out
end

function _combinations!(out::Vector{Vector{Int}}, v::Vector{Int}, r::Int)
    n = length(v)
    r > n && return
    idx = collect(1:r)
    while true
        push!(out, v[idx])
        i = r
        while i > 0 && idx[i] == n - r + i
            i -= 1
        end
        i == 0 && return
        idx[i] += 1
        for j in (i + 1):r
            idx[j] = idx[j - 1] + 1
        end
    end
end
#######
"""
In a vector of tensorterms, move j to i+1, keeping track of grading and accumulated tensors from commutation relations
"""
function drag_to_adjacent(summand_coeff, all_terms, i, j)
    # Returns a list of (coeff, terms) pairs representing the result
    # of dragging all_terms[j] left to position i+1

    results = []
    target = all_terms[j]
    current_terms = copy(all_terms)

    # drag target from position j to i+1 one step at a time
    pos = j
    coeff = summand_coeff
    term_list = coeff * TensorTerm(current_terms) #return tensor expression
    if j == i + 1
        return term_list
    end

    while pos > i + 1
        new_term_list =  zero()
        for (c, tterm) in summands(term_list)
            trms = terms(tterm)
            left = trms[pos-1]
            right = trms[pos]  # this is target
            # check if they have non-trivial commutation
            if can_swap(left, right)
                # @show "canswap",left,right
                grade = 0
                extra = nothing
            else
                # @show "need commuator",get_commutator(left, right),left,right
                (grade, extra) = get_commutator(left, right)  # returns (sign*swap, extra commutator terms)
            end
            # simple swap term
            new_terms = copy(trms)
            new_terms[pos-1] = right
            new_terms[pos] = left
            swaps = grading(trms[pos-1]) * grading(trms[pos]) + (grade < 0)
            new_term_list += (c * (-1)^(swaps) * TensorTerm(new_terms))
            # extra commutator terms if any
            if extra !== nothing
                # @show "adding extra",extra
                extra_terms = TensorTerm(trms[1:(pos-2)]) * (extra) * TensorTerm(trms[(pos+1):end])
                new_term_list += (c * extra_terms)
            end
        end
        term_list = new_term_list
        pos -= 1
    end

    return term_list
end

function drag_to_adjacent_n(summand_coeff, all_terms, indices::Vector{Int})
    anchor = indices[1]
    n = length(indices)

    # start with a single summand
    term_list = summand_coeff * TensorTerm(all_terms)

    for k in 2:n
        target_pos = anchor + k - 1
        new_term_list =  zero()
        for (c, tterm) in summands(term_list)
            trms = terms(tterm)
            # find current position of k-th target by searching for target_type
            # it should be the k-th occurrence of target_type in trms
            occurrences = findall(t -> t isa typeof(all_terms[indices[k]]), trms)
            if length(occurrences) < k || occurrences[k] == target_pos
                new_term_list += c * tterm
                continue
            end
            current_pos = occurrences[k]
            new_term_list += drag_to_adjacent(c, trms, target_pos - 1, current_pos)
        end
        term_list = new_term_list
    end
    # @show term_list,"outer"
    return term_list
end


"""
If there is a simplification rule found, manipulate by raising and lowering indicies until the layout matches
"""
function match_and_apply_rule(pattern_terms, actual_terms, rhs_template, full_term)
    # @show actual_terms
    # @show full_term
    if length(free_indices_ordered(TensorTerm(pattern_terms))) != length(free_indices_ordered(TensorTerm(actual_terms)))
        return (false,  zero())
    end
    substitution = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()

    corrected_terms = TensorSuperType[]
    used_dummies = Set{IndexSuperType}()

    for (pat, act) in zip(pattern_terms, actual_terms)
        index_map = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()  # a_idx => dummy at p_pos, for replacing on act
        for (p_idx, a_idx) in zip(indices(pat), indices(act))
            p_name, p_pos = p_idx
            a_name, a_pos = a_idx

            if p_pos != a_pos
                dummies = get_unused_indices(full_term, a_name.index_set)
                filter!(x -> !(x in used_dummies), dummies)
                dummy = safe_pop!(dummies, a_name.index_set)
                push!(used_dummies, dummy)
                substitution[p_idx] = dummy => p_pos
                
                push!(corrected_terms, Metric([a_name => a_pos, dummy => a_pos]))
             
                index_map[a_idx] = (dummy => p_pos)  # act gets dummy at pattern's level
            else
                substitution[p_idx] = a_idx
                index_map[a_idx] = a_idx  # no change
            end
        end
        # replace indices on act with corrected versions
        push!(corrected_terms, set_indices_deep!(act, index_map))
    end
    # @show corrected_terms
    #handle dummies that clash in the pattern
    rhs_dummies = dummies(rhs_template)
    actual_free = Set(i[1] for i in union(indices(full_term), indices(TensorTerm(corrected_terms))))
    clashing = filter(d -> any(x -> x== d, actual_free), rhs_dummies)
    if !isempty(clashing)
        rhs_template = copy(rhs_template)
        for clash_idx in clashing
            fresh = get_unused_indices(full_term, clash_idx.index_set)
            # @show full_term
            # @show indices(full_term)
            # @show fresh
            filter!(x -> !(x in actual_free), fresh)
            # @show fresh
            replacement = safe_pop!(fresh,clash_idx.index_set)
            # @show replacement
            push!(used_dummies, replacement)
            substitution[clash_idx=>upper] = replacement => upper
            substitution[clash_idx=>lower] = replacement => lower
        end
    end
    rhs = set_from_pattern(rhs_template, substitution)
    target_indices = findall(x -> x isa typeof(actual_terms[1]), corrected_terms)
    dragged = drag_to_adjacent_n(1, corrected_terms, target_indices)
    # @show dragged
    # replace adjacent block with rhs, keeping surrounding metrics
    result =  zero()
    n = length(actual_terms)
    for (c, tterm) in summands(dragged)
        t = terms(tterm)
        new_anchor = findfirst(x -> x isa typeof(actual_terms[1]), t)
        # @show new_anchor
        result += c * TensorTerm(t[1:(new_anchor-1)]) * rhs * TensorTerm(t[(new_anchor+n):end])
    end
    return (true, result)
end



"""
Generate a function reduce_matrix that reduces a product of matrices.
"""
macro mrule(expr)
    @assert expr.head == :(=)

    lhs = expr.args[1]
    rhs = expr.args[2]

    function extract_terms(ex)
        if ex isa Expr && ex.head == :call && ex.args[1] == :*
            out = Any[]
            for arg in ex.args[2:end]
                append!(out, extract_terms(arg))
            end
            return out
        else
            return [ex]
        end
    end

    lhs_terms = extract_terms(lhs)
    n = length(lhs_terms)

    arg_names = [gensym() for _ in 1:n]
    type_exprs = [:(typeof($(esc(t)))) for t in lhs_terms]

    args_with_types = [
        :($(arg_names[i])::$(type_exprs[i]))
        for i in 1:n
    ]

    return quote
        function Alakazam.reduce_matrix($(args_with_types...), tterm::Alakazam.TensorTerm)
            pattern_terms = TensorSuperType[$(map(esc, lhs_terms)...)]
            actual_terms = TensorSuperType[$(arg_names...)]
            rhs_template = TensorExpression($(esc(rhs)))
            return match_and_apply_rule(
                pattern_terms,
                actual_terms,
                rhs_template,
                tterm
            )
        end

        # Force this specialization into the pkgimage. The method is reached
        # by runtime dispatch from contract_adjacent, so without an explicit
        # directive its native code is compiled during precompilation but
        # dropped from the image (no static backedge to root it).
        precompile(Alakazam.reduce_matrix,
                   ($(type_exprs...), Alakazam.TensorTerm))
    end
end