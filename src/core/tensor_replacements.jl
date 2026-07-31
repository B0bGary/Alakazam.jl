####################################################################################################
#various functions to be called

function replace_pattern(exp::TensorTerm, pattern::String, replace::TensorTerm, repeat::Bool=true)
    return substitute_terms(exp, parse_pattern2tensor(pattern), TensorExpression(replace), repeat)
end
function replace_pattern(exp::TensorTerm, pattern::String, replace::TensorExpression, repeat::Bool=true)
    return substitute_terms(exp, parse_pattern2tensor(pattern), replace, repeat)
end
function replace_pattern(exp::TensorTerm, pattern::String, replace::TensorSuperType, repeat::Bool=true)
    return substitute_terms(exp, parse_pattern2tensor(pattern), TensorExpression(replace), repeat)
end
function replace_pattern(exp::TensorSuperType, pattern::String, replace::TensorSuperType, repeat::Bool=true)
    return substitute_terms(TensorTerm(exp), parse_pattern2tensor(pattern), TensorExpression(replace), repeat)
end
function replace_pattern(exp::TensorExpression, pattern::String, replace::TensorTerm, repeat::Bool=true)
    new_exp = zero()
    for summand in summands(exp)
        new_exp += coef(summand)*substitute_terms(term(summand), parse_pattern2tensor(pattern), TensorExpression(replace), repeat)
    end
    return new_exp
end
function replace_pattern(exp::TensorExpression, pattern::String, replace::TensorSuperType, repeat::Bool=true)
    new_exp = zero()
    for summand in summands(exp)
        new_exp += coef(summand)*substitute_terms(term(summand), parse_pattern2tensor(pattern), TensorExpression(replace), repeat)
    end
    return new_exp
end

function replace_pattern(exp::TensorTerm, pattern::String, replace::String, repeat::Bool=true)
    return substitute_terms(exp, parse_pattern2tensor(pattern), TensorExpression(parse_pattern2tensor(replace)), repeat)
end
function replace_pattern(exp::TensorSuperType, pattern::String, replace::String, repeat::Bool=true)
    return substitute_terms(TensorTerm(exp), parse_pattern2tensor(pattern), TensorExpression(parse_pattern2tensor(replace)), repeat)
end
function replace_pattern(exp::TensorExpression, pattern::String, replace::String, repeat::Bool=true, excluded=nothing)
    new_exp = zero()
    for summand in summands(exp)
        new_exp += coef(summand)*substitute_terms(term(summand), parse_pattern2tensor(pattern), TensorExpression(parse_pattern2tensor(replace)), repeat, excluded)
    end
    return new_exp
end


function replace_pattern(exp::TensorExpression, pattern::String, replace::TensorExpression, repeat::Bool=true)
    new_exp = zero()
    for summand in summands(exp)
        new_exp += coef(summand)*substitute_terms(term(summand), parse_pattern2tensor(pattern), replace, repeat)
    end
    return new_exp
end

function replace_pattern(exp::TensorExpression, pattern::TensorTerm, replace::TensorExpression, repeat::Bool=true, excluded=nothing)
    new_exp = zero()
    for summand in summands(exp)
        new_exp += coef(summand)*substitute_terms(term(summand), pattern, replace, repeat, excluded)
    end
    return new_exp
end

function replace_pattern(exp::TensorSuperType, pattern::String, replace::TensorExpression, repeat::Bool=true)
    return substitute_terms(TensorTerm(exp), parse_pattern2tensor(pattern), replace, repeat)
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
Find instances of a tensor term inside a larger term, and replace it with another term. 
Uses pattern matching to find the correct replacements, and accounts for grading. 
Repeats until no more matches are found if repeat is true.
Groups tensors into instances of the term before replacing, so that grading is correctly considered.
"""
function substitute_terms(exp::TensorTerm, pattern::TensorTerm, replace::TensorExpression, repeat::Bool=true, excluded::Set{T}=Set{TensorSuperType}()) where T<:TensorSuperType

    #call recursviely on the oeprand...
    if isTensor(exp) isa OperatorSuperType
        # println("did this")
        # @show exp
        set_operand!(isTensor(exp), replace_pattern(operand(isTensor(exp)), pattern, replace, repeat, excluded))
        # @show exp
    end
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
        new_exp = zero()

        for match in out.all_matches
            grade = 0
            exp2 = copy(terms(exp))
            vlast = 0

            for v in values(match.position_matches)
                v = v[2]
                if vlast != 0 && vlast < v
                    v -= 1
                end
                if v != 1
                    for i in 1:(v-1)
                        grade += grading(exp2[i]) * grading(exp2[v])
                    end
                end
                splice!(exp2, v)
                vlast = v
            end

            rep2 = copy(replace)
            for summand in summands(rep2)
                for repind in eachindex(terms(term(summand)))
                    rep = terms(term(summand))[repind]
                    rep_dum = typeof(rep) <: DummyPatternTensor
                    if rep_dum
                        set_name!(rep, match.tensor_name_matches[name(rep)])
                    end
                    operator_index_subs!(match.index_matches, rep)
                    update_indices!(rep)
                    if rep isa DummyPatternTensor
                        rep = typeof(terms(exp)[1])(data(rep))
                    end
                end
            end

            # add replacement tensors to excluded so repeat doesn't re-substitute into them
            for summand in summands(rep2)
                for t in terms(term(summand))
                    push!(excluded, t)
                end
            end

            new_exp += ((-1)^grade) * rep2 * TensorTerm(exp2)
            update_indices!(rep2)
        end

        if repeat
            tex = replace_pattern(new_exp, pattern, replace, repeat, excluded)
            typeof(tex) != TensorExpression && (tex = TensorExpression(tex))
            return tex
        end

        return new_exp
    else
        return exp
    end
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
                            match = []
                            j = 1
                            for i in 1:ekey
                                outt = check_index_match(pvalues[1][2][j], eval[2][i], tensorpatt, tensorexp, dummy_index_matching)
                                if outt[1]
                                    push!(match, pvalues[1][2][j] => eval[2][i])
                                    j += 1
                                    if j > pkeys[1]
                                        dumindtemp = outt[2]
                                        break
                                    end
                                end
                            end
                            if j <= pkeys[1]
                                if !pdummy
                                    return (matched=false, all_matches=NamedTuple[])
                                else
                                    continue
                                end
                            end
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
        # multiple frequencies in pattern — must all match
        # seed with empty match
        all_current = NamedTuple[(position_matches=[], tensor_name_matches=copy(dummy_matching), index_matches=copy(dummy_index_matching))]

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
Checks if the indices of a pattern tensor and expression tensor match, accounting for dummy indices. Returns true if they match, false otherwise. 
    Also returns a mapping of pattern dummy indices to expression indices if they match.
"""
function check_index_match(loc_pat::Int, loc_exp::Int, pat::TensorTerm, exp::TensorTerm, dumind::Dict{String,IndexSuperType})
    dummy = copy(dumind)
    # i think we dont need this if statement...
    # if typeof(terms(pat)[loc_pat]) <: DummyPatternTensor
    pii=indices(terms(pat)[loc_pat])
    ei = indices(terms(exp)[loc_exp])
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
            s = string(iterate(it)[1])
            tensor = DummyPatternTensor("?" * s) #?A  ie the name isn't important
        else
            tensor = Tensor(string(char)) #A     name is important #maybe something about coord dependance
        end
        next = peek(it)
        if next !== nothing
            char = iterate(it)[1]
            #now, do indices
            #   while string(char) == "^" ||string(char)=="_"
            while string(char) == "^" || string(char) == "_"
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
                next = peek(it)
                if next !== nothing
                    char = iterate(it)[1]
                else
                    break
                end
            end

        end
        push!(tensor_term, copy(tensor))#adds a copy
    end

    return TensorTerm(tensor_term) #should construct a term out of our set of patterns

end