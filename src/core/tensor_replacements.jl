####################################################################################################
# replace_pattern: the single entry point.

"""Normalise anything the user may pass as an expression or a replacement."""
_as_expression(v::TensorExpression) = v
_as_expression(v::Tensorial) = TensorExpression(v)
_as_expression(v::String) = TensorExpression(parse_pattern2tensor(v))
# a contraction standing for a scalar invariant is the commonest rule of all, so numbers and
# symbolic scalars are accepted as replacements directly
_as_expression(v::SymorNum) = TensorExpression(v)

"""Normalise anything the user may pass as a pattern into the `TensorTerm` the matcher wants."""
_as_pattern(p::TensorTerm) = p
_as_pattern(p::TensorSuperType) = TensorTerm(p)
_as_pattern(p::String) = parse_pattern2tensor(p)

"""
    replace_pattern(exp, pattern, replace, repeat=true, excluded=Set{TensorSuperType}())

Find `pattern` in `exp` and substitute `replace` for it, repeating until no match remains unless
`repeat` is false.

Each of the three may be given in whatever form is convenient: a bare tensor, a `TensorTerm`, or a
`TensorExpression`, and the pattern and the replacement may also be pattern strings such as
`"?A_{?x}B^{?x}"`. They are normalised here, so every combination is accepted rather than only the
pairings that happen to have a method.

A wildcard is `?` followed by a single character: `?A` as a name matches any tensor, `?x` in a slot
matches any index, and repeated wildcards must bind consistently. In the replacement, `?A` supplies
the matched tensor's name while `??A` supplies the matched object itself.
"""
function replace_pattern(exp::Tensorial, pattern::Union{Tensorial,String},
                         replace::Union{Tensorial,String,SymorNum}, repeat::Bool=true,
                         excluded::Set{TensorSuperType}=Set{TensorSuperType}(),
                         descend::Bool=true)
    pat = _as_pattern(pattern)
    rep = _as_expression(replace)
    new_exp = zero()
    for summand in summands(_as_expression(exp))
        tt = descend ? _substitute_in_operands(term(summand), pat, rep, repeat, excluded) : term(summand)
        new_exp += coef(summand) * substitute_terms(tt, pat, rep, repeat, excluded)
    end
    return new_exp
end

"""
Rewrite inside the operand of every operator factor of a term, leaving the factors themselves to
`substitute_terms`.

Operands are descended into once, on the way in. The `repeat` pass of `substitute_terms` passes
`descend = false`, since by then the only operators left are those the replacement introduced
"""
function _substitute_in_operands(tt::TensorTerm, pattern::TensorTerm, replace::TensorExpression,
                                 repeat::Bool, excluded::Set{T}) where {T<:TensorSuperType}
    fs = terms(tt)
    any(f isa OperatorSuperType for f in fs) || return tt
    return TensorTerm(TensorSuperType[
        f isa OperatorSuperType ?
            _with_operand(f, replace_pattern(operand(f), pattern, replace, repeat, excluded)) : f
        for f in fs])
end



####################################################################################################
#technical backend

## O(n) runtime approx amortised
"""
Constructs a histogram of frequency of tensors in a term, and includes their positions
"""
function term2histogram(term::TensorTerm, excluded::Set{T}=Set{TensorSuperType}()) where T<:TensorSuperType
    pat = Dict{String,Pair{Int,Vector{Int}}}() # tensor, (counts, [positions in expression]) #dont check indices yet
    #some sort of dictionary with repeated elements, constant time access etc after hashing
    for (i, tens) in enumerate(terms(term))
        tens in excluded && continue #skip some
        old_val = get(pat, name(tens), nothing)
        if old_val === nothing
            push!(pat, name(tens) => (1 => [i])) #if not in hist, add it 
        else

            push!(pat, name(tens) => (old_val[1] + 1 => push!(old_val[2], i))) #increment the count

        end
    end
    return invert_histogram(pat)
end


"""
Inverts a histogram of the form [tensorname => (count, [positions in term])] to [count => [tensorname => [positions in term]]]
"""
#change the histogram to make the counts the keys
function invert_histogram(hist::Dict{String,Pair{Int,Vector{Int}}})
    #maybe n^2 time? less i think. n for going through all, Ordered dict is linked list, at worst O(n) insertion
    new_hist = OrderedDict{Int,Vector{Pair{String,Vector{Int}}}}() #counts => [tensorname =>[positions in term]]
    for (key, value) in zip(keys(hist), values(hist))
        old_val = get(new_hist, value[1], nothing)
        if old_val === nothing
            push!(new_hist, value[1] => [key => value[2]]) #if not in hist, add it 
        else

            push!(new_hist, value[1] => push!(old_val, key => value[2])) #increment the count

        end
    end
    return new_hist

end



"""
True when the matched factors cannot be brought together by a graded sign alone: some factor lying
between them has a registered commutation relation, so moving past it emits extra terms rather than
just a sign.
"""
function _needs_drag(allt, poss::Vector{Int})
    length(poss) > 1 || return false
    anchor = poss[1]
    mset = Set(poss)
    for p in poss
        p == anchor && continue
        for x in (anchor+1):(p-1)
            x in mset && continue
            can_swap(allt[p], allt[x]) || return true
        end
    end
    return false
end

"""
Bring the factors at `poss` together at `poss[1]`, honouring the registered commutation relations.

Returns a `TensorExpression`: the summand in which the matched factors are adjacent, plus whatever
commutator remainders the moves generated. Those remainders no longer contain all of the matched
factors, so they are simply left for the caller's next matching pass rather than substituted into.

Factors are located by object identity, because after a drag the positions have shifted and some
summands will have lost the factor entirely to a commutator.
"""
function _gather_matched(allt, poss::Vector{Int})
    targets = [allt[p] for p in poss]
    out = TensorExpression(TensorTerm(copy(allt)))
    for k in 2:length(targets)
        acc = zero()
        for (c, tt) in summands(out)
            trms = terms(tt)
            i = findfirst(t -> t === targets[k-1], trms)
            j = findfirst(t -> t === targets[k], trms)
            if i === nothing || j === nothing || j <= i + 1
                acc += c * tt                      # adjacent already, or eaten by a commutator
            else
                dragged = drag_to_adjacent(c, trms, i, j)
                acc += dragged === nothing ? c * tt : dragged
            end
        end
        out = acc
    end
    return out
end

"""
Find instances of a tensor term inside a larger term, and replace it with another term. 
Uses pattern matching to find the correct replacements, and accounts for grading. 
Repeats until no more matches are found if repeat is true.
Groups tensors into instances of the term before replacing, so that grading is correctly considered.
"""
function substitute_terms(exp::TensorTerm, pattern::TensorTerm, replace::TensorExpression, repeat::Bool=true, excluded::Set{T}=Set{TensorSuperType}()) where T<:TensorSuperType

    # An empty pattern would match vacuously at every position and substitute nothing, so it is
    # rejected rather than silently rewriting the expression.
    isempty(terms(pattern)) && throw(InvalidReplacementPattern(
        "the pattern is empty", "a pattern must name at least one tensor"))

    # operands are rewritten by `_substitute_in_operands` before this is called
    flg=isTensorTerm(exp)
    if (flg === nothing)
        return exp
    else
        exp=flg
    end

    exp_hist = term2histogram(exp, excluded)
    pat_hist = term2histogram(pattern)
    out = patt_match(pat_hist, exp_hist, pattern, exp, Dict{String,String}(), Dict{String,IndexSuperType}(), excluded, repeat)

    if out.matched
        # Non-commuting factors between the matched ones cannot be crossed with a sign: the move emits commutator terms. 
        let allt = terms(exp), poss = sort!([v[2] for v in values(out.all_matches[1].position_matches)])
            if _needs_drag(allt, poss)
                gathered = _gather_matched(allt, poss)
                # A drag blocked by a non-commuting pair hands the term back untouched. The matched factors can then never be brought together, 
                # so the pattern does not apply here and recursing on an unchanged term would not terminate.
                gathered == TensorExpression(exp) && return exp
                return replace_pattern(gathered, pattern, replace, repeat, excluded)
            end
        end

        new_exp = zero()

        # Only the first match is substituted. `all_matches` holds the alternative ways the pattern
        # could bind, which are choices rather than summands: adding them would multiply-count a
        # single occurrence, and for a wildcard name matching several different tensors the
        # alternatives need not even share free indices. Further occurrences are reached by `repeat`.
        let match = first(out.all_matches)
            allt = terms(exp)
            poss = sort!([v[2] for v in values(match.position_matches)])
            matched = Set(poss)
            anchor = poss[1]

            # The replacement is now spliced in where the match began, not pulled to the front of the term. 
            # Only the matched factors move: each one travels left to the anchor, 
            # past exactly the unmatched factors lying between the anchor and itself. 
            #That is what the graded sign counts.
            grade = 0
            for p in poss
                p == anchor && continue
                for x in (anchor+1):(p-1)
                    x in matched && continue
                    grade += grading(allt[p]) * grading(allt[x])
                end
            end
            pre = TensorSuperType[allt[i] for i in 1:(anchor-1)]
            post = TensorSuperType[allt[i] for i in (anchor+1):length(allt) if !(i in matched)]

            # How a wildcard in the replacement resolves:
            #   ?A   takes the matched tensor's NAME, and nothing else. Coordinates, tableaux,
            #        labels and the concrete type all come from the replacement as written, so
            #        `?A` can be used to strip or restate them.
            #   ??A  takes the matched OBJECT: its data and its type, with only the index
            #        placement coming from the replacement. `??A -> ??A` is the identity, and a
            #        matched Fermion stays a Fermion rather than losing its grading.
            # Both refer to the same binding, so `??A` looks up the `?A` recorded by the pattern.
            pat_terms = terms(pattern)
            matched_src = Dict{String,TensorSuperType}()
            for pr in match.position_matches
                matched_src[name(pat_terms[pr[1]])] = allt[pr[2]]
            end

            rep2 = copy(replace)
            resolve_wildcard_names!(rep2, match.tensor_name_matches, matched_src)
            freshen_private_dummies!(rep2, allt, nothing, new_exp)
            for summand in summands(rep2)
                for rep in terms(term(summand))
                    operator_index_subs!(match.index_matches, rep)
                    update_indices!(rep)
                end
            end

            # add replacement tensors to excluded so repeat doesn't re-substitute into them
            for summand in summands(rep2)
                for t in terms(term(summand))
                    push!(excluded, t)
                end
            end


            new_exp += ((-1)^grade) * TensorTerm(pre) * rep2 * TensorTerm(post)
            update_indices!(rep2)
        end

        if repeat
            tex = replace_pattern(new_exp, pattern, replace, repeat, excluded, false)
            typeof(tex) != TensorExpression && (tex = TensorExpression(tex))
            return tex
        end

        return new_exp
    else
        return exp
    end
end

"""
    freshen_private_dummies!(rep2, factors, extra, sofar)

Rename a replacement's own dummy indices away from everything the match is being spliced into.
"""
function freshen_private_dummies!(rep2::TensorExpression, factors, extra, sofar)
    taken = OrderedSet{IndexSuperType}()
    for ctx in (factors, extra)                 # factor lists
        (ctx === nothing || isempty(ctx)) && continue
        collect_indices!(taken, TensorTerm(ctx))
    end
    sofar === nothing || collect_indices!(taken, sofar)   # an expression
    isempty(taken) && return rep2

    blocked = OrderedSet{IndexSuperType}()
    for q in taken
        push!(blocked, q)
    end
    for summand in summands(rep2), (q, _) in indices(term(summand))
        push!(blocked, q)          # nothing the replacement already uses may be a rename target
    end
    renamed = Dict{IndexSuperType,IndexSuperType}()

    for summand in summands(rep2)
        tt = term(summand)
        dums = dummies(tt)
        isempty(dums) && continue

        guide = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
        for (ind, pos) in indices(tt)
            (ind in dums) || continue          # free indices carry the pattern binding: leave them
            (ind in taken) || continue         # no clash with the context, nothing to do
            ind isa DummyPatternIndex && continue
            if !haskey(renamed, ind)
                set = ind.index_set
                pool = OrderedSet{IndexSuperType}(
                    q for q in get_unused_indices(tt, set) if !(q in blocked))
                new = safe_pop!(pool, set)   # generates a further index if the pool runs dry
                renamed[ind] = new
                push!(blocked, new)
                push!(taken, new)
            end
            guide[ind=>pos] = renamed[ind] => pos
        end
        isempty(guide) || set_indices_deep!(tt, guide)
    end
    return rep2
end

"""
Resolve wildcard tensor names in a replacement, descending through operator operands.

`names` maps a pattern wildcard to the concrete name it bound; `srcs` maps it to the matched object
itself. A name left unresolved would print as a literal `?A` in the result, so an unbound wildcard is
reported rather than passed through.
"""
function resolve_wildcard_names!(rep, names::Dict{String,String}, srcs::Dict{String,TensorSuperType})
    if rep isa TensorExpression
        for summand in summands(rep)
            trms = terms(term(summand))
            for i in eachindex(trms)
                trms[i] = resolve_wildcard_names!(trms[i], names, srcs)
            end
        end
        return rep
    elseif rep isa OperatorSuperType
        resolve_wildcard_names!(operand(rep), names, srcs)
        return rep
    end

    pn = name(rep)
    startswith(pn, "?") || return rep
    whole = startswith(pn, "??")                      # ??A adopts the object, ?A only its name
    key = whole ? pn[nextind(pn, 1):end] : pn

    nm = get(names, key, nothing)
    nm === nothing && throw(InvalidReplacementPattern(
        "wildcard tensor '$key' in the replacement has no match in the pattern",
        isempty(names) ? "the pattern names no wildcard tensors" :
                         "available: " * join(sort!(collect(keys(names))), ", ")))

    if whole
        src = get(srcs, key, nothing)
        src !== nothing && applicable(typeof(src), data(src), indices(rep)) &&
            return typeof(src)(data(src), copy(indices(rep)))
    end
    set_name!(rep, nm)
    return rep
end

"""Deep function for replacing dummies in operators"""
function operator_index_subs!(out, rep)
    if rep isa TensorExpression
        for summand in summands(rep)
            for t in terms(term(summand))
                operator_index_subs!(out, t)
            end
        end
    elseif !(rep isa OperatorSuperType)
        ri = copy(indices(rep))
        for index in eachindex(ri)
            if ri[index][1] isa DummyPatternIndex
                if get(out, ri[index][1].name, nothing) === nothing #problem with checking equality since technically different objects #use names instead
                    throw(InvalidReplacementPattern(
                        "pattern index '?$(ri[index][1].name)' has no match in the replacement",
                        "available: " * join(sort!(collect(keys(out))), ", ")))
                else
                    ri[index] = out[name(ri[index][1])] => ri[index][2] # maps index from original expression to replacement
                end #need a deep version etc to handle operators...
            end
        end
        set_indices!(rep, ri)
    else
        #do op inds
        ri = copy(op_indices(rep))
        for index in eachindex(ri)
            if ri[index][1] isa DummyPatternIndex
                if get(out, ri[index][1].name, nothing) === nothing #problem with checking equality since technically different objects #use names instead
                    throw(InvalidReplacementPattern(
                        "pattern index '?$(ri[index][1].name)' has no match in the replacement",
                        "available: " * join(sort!(collect(keys(out))), ", ")))
                else
                    ri[index] = out[name(ri[index][1])] => ri[index][2] # maps index from original expression to replacement
                end #need a deep version etc to handle operators...
            end
        end
        #recursively do the operand.
        set_op_indices!(rep, ri)                # was missing entirely
        operator_index_subs!(out, operand(rep))
    end
end



"""Matches a pattern dummy tensor to an expression tensor"""
function match_dummy_logic(dummy_matching::Dict{String,String}, dummy_name::String, exp_dummy::String)
    dummy = copy(dummy_matching)
    if get(dummy, dummy_name, nothing) === nothing #if dummy isnt already assigned a val in match
        dummy[dummy_name] = exp_dummy
    elseif dummy[dummy_name] != exp_dummy #dummy doesnt match
        return (false, dummy_matching) #dummy tensor name pattern
    end
    return (true, dummy)

end


"""
Assign every pattern slot of one name class to a distinct expression position.

Each slot is tried against every free position rather than scanned forward, because a forward scan
is order-dependent: it cannot match a pattern whose slots occur in the expression in the opposite
order. That case arises as soon as several factors share a name and are told apart by their labels
or index positions (`k_1` and `k_2`, say) which is precisely what labels are for.

Returns `(matched, [pattern position => expression position], bindings)`. The search is bounded by
the size of the name class, which the histogram has already cut down.
"""
function assign_slots(patpos::Vector{Int}, exppos::Vector{Int}, j::Int, used::AbstractVector{Bool},
                      acc::Dict{String,IndexSuperType}, patt::TensorTerm, expr::TensorTerm)
    j > length(patpos) && return (true, Pair{Int,Int}[], acc)
    for t in eachindex(exppos)
        used[t] && continue
        outt = check_index_match(patpos[j], exppos[t], patt, expr, acc)
        outt[1] || continue
        used[t] = true
        done, rest, acc3 = assign_slots(patpos, exppos, j + 1, used, outt[2], patt, expr)
        used[t] = false
        done && return (true, pushfirst!(rest, patpos[j] => exppos[t]), acc3)
    end
    return (false, Pair{Int,Int}[], acc)
end

"""Uses pattern matching to find a pattern in an expression, and returns the match if found. Accounts for dummy tensors and indices.
Uses a recursive dynamic programming approach, peeling off one term at a time from the pattern and expression, and checking for matches."""
function patt_match(pat::OrderedDict{Int,Vector{Pair{String,Vector{Int}}}}, exp::OrderedDict{Int,Vector{Pair{String,Vector{Int}}}},
    tensorpatt::TensorTerm, tensorexp::TensorTerm, dummy_matching::Dict{String,String}=Dict{String,String}(), dummy_index_matching::Dict{String,IndexSuperType}=Dict{String,IndexSuperType}(),
    excluded::Set{T}=Set{TensorSuperType}(), repeat::Bool=false) where T<:TensorSuperType

    matches = NamedTuple[]

    if length(pat) == 1
        (pkeys, pvalues) = first(pat)
        if length(first(pat)[2]) == 1
            if length(exp) == 1
                looper = [first(exp)]
            else
                looper = zip(keys(exp), values(exp))
            end
            for (ekey, evalue) in looper
                if pkeys[1] <= ekey
                    for eval in evalue
                        pdummy = false
                        if pvalues[1][1][1] == '?'
                            pdummy = true
                        end
                        if pdummy || pvalues[1][1] == eval[1]
                            dummy_temp = copy(dummy_matching)
                            dumindtemp = copy(dummy_index_matching)
                            if pdummy
                                out = match_dummy_logic(dummy_temp, pvalues[1][1], eval[1])
                                if out[1] == false
                                    continue
                                else
                                    dummy_temp = out[2]
                                end
                            end
                            # Bindings are threaded through the assignment, so a repeated wildcard
                            # must bind consistently: passing the original `dummy_index_matching`
                            # to every check let one rebind freely, and "q_{?x}q^{?x}" then matched
                            # `q^i q_j`, silently zeroing terms that do not vanish.
                            ok, match, acc = assign_slots(pvalues[1][2], eval[2], 1,
                                                          falses(length(eval[2])),
                                                          dummy_index_matching, tensorpatt, tensorexp)
                            if !ok
                                if !pdummy
                                    return (matched=false, all_matches=NamedTuple[])
                                else
                                    continue
                                end
                            end
                            dumindtemp = acc
                            push!(matches, (position_matches=match, tensor_name_matches=copy(dummy_temp), index_matches=copy(dumindtemp)))
                            if !repeat
                                return (matched=true, all_matches=matches)
                            end
                        end
                    end
                end
                continue
            end
            isempty(matches) && return (matched=false, all_matches=NamedTuple[])
            return (matched=true, all_matches=matches)

        else
            if length(exp) == 1
                looper = [first(exp)]
            else
                looper = zip(keys(exp), values(exp))
            end
            for (ekey, evalue) in looper
                for eval in evalue
                    check = patt_match(OrderedDict(first(pat)[1] => [first(pat)[2][1]]), OrderedDict(ekey => [eval]), tensorpatt, tensorexp, dummy_matching, dummy_index_matching, excluded, repeat)
                    if check.matched
                        for cmatch in check.all_matches
                            exp2 = copy(exp)
                            pat2 = copy(pat)
                            rem = copy(first(pat)[2])
                            popfirst!(rem)
                            push!(pat2, first(pat)[1] => rem)
                            #changed this to non-mutating filter...
                            push!(exp2, ekey => filter(x -> x[1] != eval[1], evalue))
                            if ekey - pkeys[1] >= 1
                                if get(exp2, ekey - pkeys[1], nothing) === nothing
                                    push!(exp2, (ekey - pkeys[1]) => [eval[1] => eval[2][pkeys:end]])
                                else
                                    #used vcat instead of push here...
                                    push!(exp2, (ekey - pkeys[1]) => vcat(exp2[ekey-pkeys[1]], eval[1] => eval[2][pkeys:end]))
                                end
                            end
                            out = patt_match(pat2, exp2, tensorpatt, tensorexp, cmatch.tensor_name_matches, cmatch.index_matches, excluded, repeat)
                            if out.matched
                                for omatch in out.all_matches
                                    push!(matches, (
                                        position_matches=vcat(cmatch.position_matches, omatch.position_matches),
                                        tensor_name_matches=copy(omatch.tensor_name_matches),
                                        index_matches=copy(omatch.index_matches)
                                    ))
                                    if !repeat
                                        return (matched=true, all_matches=matches)
                                    end
                                end
                            end
                        end
                    end
                end
            end
            isempty(matches) && return (matched=false, all_matches=NamedTuple[])
            return (matched=true, all_matches=matches)
        end

    else
        # multiple frequencies in pattern, must all match
        # seed with empty match
        all_current = NamedTuple[(position_matches=Pair{Int,Int}[], tensor_name_matches=copy(dummy_matching), index_matches=copy(dummy_index_matching))]

        for (pkeys, pvalues) in pat
            next = NamedTuple[]
            for current in all_current
                out = patt_match(OrderedDict(pkeys => pvalues), exp, tensorpatt, tensorexp, current.tensor_name_matches, current.index_matches, excluded, repeat)
                if out.matched
                    for omatch in out.all_matches
                        push!(next, (
                            position_matches=vcat(current.position_matches, omatch.position_matches),
                            tensor_name_matches=copy(omatch.tensor_name_matches),
                            index_matches=copy(omatch.index_matches)
                        ))
                        if !repeat
                            return (matched=true, all_matches=next)
                        end
                    end
                end
            end
            isempty(next) && return (matched=false, all_matches=NamedTuple[])
            all_current = next
        end

        return (matched=true, all_matches=all_current)
    end
end

"""Special logic for remapping indices in a replacement when there are operators."""
function construct_index_rep_map(tenpat::OperatorSuperType, tenact::OperatorSuperType)
    map = Dict{String,Index}()
    for pos in eachindex(op_indices(tenpat))
        if op_indices(tenpat)[pos][2] != op_indices(tenact)[pos][2] || op_indices(tenpat)[pos][1].index_set != op_indices(tenact)[pos][1].index_set #if wrong positions or wrong index set
            pi_, ai = op_indices(tenpat)[pos], op_indices(tenact)[pos]
            if pi_[2] != ai[2]
                throw(InvalidReplacementPattern(
                    "index position mismatch at slot $pos: pattern has $(pi_[1]) $(pi_[2] == upper ? "upper" : "lower"), " *
                    "target has $(ai[1]) $(ai[2] == upper ? "upper" : "lower")"))
            elseif pi_[1].index_set != ai[1].index_set
                throw(InvalidReplacementPattern(
                    "index set mismatch at slot $pos: pattern index $(pi_[1]) is from '$(pi_[1].index_set.name)', " *
                    "target index $(ai[1]) is from '$(ai[1].index_set.name)'"))
            end
        end
        map[indices(tenpat)[pos][1].name] = op_indices(tenact)[pos][1]
    end
    return map
end


"""
The labels of a factor, seen through any operators standing over it. An operator does not carry
its operand's labels, so would otherwise match. A chain over something that is not a single tensor
reports nothing, which leaves it a wildcard over labels exactly as before.
"""
_labels_through_operators(t::Tensorial) = labels(t)
function _labels_through_operators(O::OperatorSuperType)
    ls = labels(O)
    isempty(ls) || return ls
    inner = isTensor(operand(O))
    return inner === nothing ? ls : _labels_through_operators(inner)
end

"""
The name of the tensor a factor, seen through any operators over it. 
"""
_operand_name_through_operators(t::Tensorial) = name(t)
function _operand_name_through_operators(O::OperatorSuperType)
    inner = isTensor(operand(O))
    return inner === nothing ? "" : _operand_name_through_operators(inner)
end

"""
Checks if the indices of a pattern tensor and expression tensor match, accounting for dummy indices. Returns true if they match, false otherwise.
    Also returns a mapping of pattern dummy indices to expression indices if they match.
"""
function check_index_match(loc_pat::Int, loc_exp::Int, pat::TensorTerm, exp::TensorTerm, dumind::Dict{String,IndexSuperType})
    dummy = copy(dumind)
    # i think we dont need this if statement...
    ptens = terms(pat)[loc_pat]
    etens = terms(exp)[loc_exp]
    plabs = _labels_through_operators(ptens)
    if !isempty(plabs) && plabs != _labels_through_operators(etens)
        return (false, dumind)
    end
    pname = _operand_name_through_operators(ptens)
    if !isempty(pname) && pname[1] != '?' && pname != _operand_name_through_operators(etens)
        return (false, dumind)
    end

    pii = indices(ptens)
    ei = indices(etens)
    if length(pii) != length(ei)
        return (false, dumind)
    end
    for ind = eachindex(pii)
        if typeof(pii[ind][1]) <: DummyPatternIndex
            if pii[ind][2] == ei[ind][2]#checks if same position
                indexx = pii[ind][1]
                if get(dummy, indexx.name, nothing) === nothing #if dummy isnt already assigned a val in match
                    dummy[indexx.name] = ei[ind][1]
                elseif dummy[indexx.name].name != ei[ind][1].name #dummy doesnt match
                    return (false, dumind) #dummy tensor name pattern
                end
            else
                return (false, dumind) #if different positions, then no match
            end
        else
            if pii[ind] != ei[ind]
                return (false, dumind)
            end
        end

    end
    return (true, dummy)
    # else
    #     return (indices(terms(pat)[loc_pat]) == indices(terms(exp)[loc_exp]), dumind)
    # end
end




"""Given a string to represent a replacement pattern, constructs a dummy tensor term out of it.
The syntax is as follows:
    -Tensors are represented by capital letters, eg A, B, C etc. 
    -Indices are represented by lowercase letters, eg a, b, c etc. 
    -Upper indices are denoted by ^, lower indices by _. 
    -If there are multiple upper or lower indices, they should be enclosed in {}, eg ^{ab}, _{cd} etc. 
    -Dummy tensors can be represented by prefixing the tensor name with the wildcard?, eg ?A, ?B etc. Dummy indices can be represented similarly, eg ?a, ?b etc. 
    -Wildcards ?A, ?B etc will match with any tensor, but the same wildcard will match with the same tensor across the pattern. Similarly for dummy indices.
    -Spaces are ignored."""
function parse_pattern2tensor(pattern::String)
    #can probs unify some of the parsing logic with the other "macro" that constructs tensors from strings
    pattern = strip(pattern)#remove spaces
    it = Iterators.Stateful(pattern)
    #add a way user can do custom patterns, eg \eps for eps tensor, \delta for delta func etc. some sort of vector containg possible presets (or a tree!)
    #else itll always output as a default tensor
    tensor_term = Vector{TensorSuperType}() #array of tensors
    tensor = []
    # need to : make many tensors, put them all in a tensor term
    #Construct the tensor
    for char in it #should only get here at the start of a tensor if valid syntax!!!!!!!!!!
        if string(char) == "?"
            # "?A" is a wildcard name; "??A" additionally adopts the matched object's data
            prefix = "?"
            if peek(it) == '?'
                iterate(it)
                prefix = "??"
            end
            nxt = iterate(it)
            nxt === nothing && throw(InvalidReplacementPattern(
                "pattern ends with '$prefix', which names no tensor",
                "a wildcard name is '?' (or '??') followed by a single character, as in \"?A\""))
            tensor = Tensor(prefix * string(nxt[1]))
        else
            tensor = Tensor(string(char)) #A     name is important #maybe something about coord dependance
        end
        # Parse this tensor's trailing index groups. PEEK before consuming: consuming the
        # terminator swallows the next tensor's name, which is what broke multi-tensor patterns.
        while true
            nxt = peek(it)
            (nxt === nothing || (nxt != '^' && nxt != '_')) && break
            char = iterate(it)[1]
                if string(char) == "^"
                    if string(peek(it)[1]) == "{"
                        iterate(it)
                        ## do indices
                        ind = 1
                        while string(peek(it)[1]) != "}"
                            if length(indices(tensor)) < ind #if not already an index
                                if string(peek(it)[1]) == "?"
                                    iterate(it)
                                    push!(indices(tensor), DummyPatternIndex(string(iterate(it)[1])) => upper)
                                else
                                    push!(indices(tensor), Index(string(iterate(it)[1])) => upper)
                                end

                            elseif string(peek(it)[1]) != "*" #if already, and current value is not a *
                                if string(peek(it)[1]) == "?"
                                    iterate(it)
                                    indices(tensor)[ind] = DummyPatternIndex(string(iterate(it)[1])) => upper
                                else
                                    indices(tensor)[ind] = Index(string(iterate(it)[1])) => upper
                                end
                            else
                                iterate(it)
                            end
                            ind += 1
                        end
                        iterate(it)

                    else
                        push!(indices(tensor), Index(string(iterate(it)[1])) => upper) #single upper index
                    end
                end
                # char = iterate(it)[1]
                if string(char) == "_"
                    if string(peek(it)[1]) == "{"
                        iterate(it)
                        ## do indices
                        ind = 1
                        while string(peek(it)[1]) != "}"
                            if length(indices(tensor)) < ind #if not already an index, stops out ofbounds error. add new index
                                if string(peek(it)[1]) == "?"
                                    iterate(it)
                                    push!(indices(tensor), DummyPatternIndex(string(iterate(it)[1])) => lower)
                                else
                                    push!(indices(tensor), Index(string(iterate(it)[1])) => lower)
                                end

                            elseif string(peek(it)[1]) != "*" #if already, and current value is not a *
                                if string(peek(it)[1]) == "?" #next index is a dummy, won't need exact match in replacement
                                    iterate(it)
                                    indices(tensor)[ind] = DummyPatternIndex(string(iterate(it)[1])) => lower
                                else
                                    indices(tensor)[ind] = Index(string(iterate(it)[1])) => lower
                                end
                            else
                                iterate(it) # removes * 
                            end
                            ind += 1
                        end
                        iterate(it) #remove }
                    else
                        push!(indices(tensor), Index(string(iterate(it)[1])) => lower) #single upper index
                    end
                end
        end
        push!(tensor_term, copy(tensor))#adds a copy
    end

    return TensorTerm(tensor_term) #should construct a term out of our set of patterns

end

"""
    solve_linear(eqn, target) -> (target, value) | (nothing, eqn)

Basic function to olve the linear equation `eqn = 0` for `target`.
`target` may be any single factor, an operator included, so `∂_μ[A]` can be solve_lineard as readily as `A`.

    solve_linear(2*A + B, A)              # (A, -1//2 B)
    solve_linear(∂(A) - B, ∂(A))          # (∂[A], B)
    solve_linear(TensorExpression(B), A)   # (nothing, B)
"""
function solve_linear(eqn::TensorExpression, target::Tensorial)
    coef = 0
    rest = zero()
    for (k, tt) in summands(eqn)
        fs = terms(tt)
        if length(fs) == 1 && fs[1] == target
            coef += k
        elseif any(f -> f == target, fs)
            return (nothing, eqn)        # multiplied by something else, so not linear in `target`
        else
            rest = rest + TensorExpression(CoeffPair[CoeffPair(k, tt)])
        end
    end
    coef == 0 && return (nothing, eqn)   # `target` does not appear at all
    return (target, simplify((-(1 // 1) / coef) * rest))
end
solve_linear(eqn::Tensorial, target::Tensorial) = solve_linear(TensorExpression(eqn), target)
