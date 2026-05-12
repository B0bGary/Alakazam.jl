
###############################################################
###############################################################
###############################################################
###############################################################

#derivation actions and patterns


derivative_actions = Dict{DerivativeSuperType,Dict{TensorSuperType,TensorExpression}}() #where {D<:}
derivative_algebra = Dict{Pair{DerivativeSuperType,DerivativeSuperType},Pair{Integer,TensorExpression}}() #integer is 1 for comm, -1 for anticommut, expression is the result of the commutator/anticommutator (eg for fermionic derivatives)



function hard_equals(a::DerivativeSuperType, b::DerivativeSuperType)#check to see if derivatives commute
    if a.name == b.name && a.operator_indices == b.operator_indices && a.operand == b.operand #basic check
        return true
    end
    opra = isTensor(a.operand)
    oprab = isTensor(b.operand)
    if oprab != false && opra != false && typeof(opra) <: DerivativeSuperType && typeof(oprab) <: DerivativeSuperType #nested ops

        # println("checking nested derivatives")
        # println(opra)        
        # println(oprab)
        #check to see if the operand of one derivative is the other, and if the operators are the same (if they commute then this is true)
        if opra.operator_indices == b.operator_indices && oprab.operator_indices == a.operator_indices && opra.name == b.name && oprab.name == a.name && opra.operand == oprab.operand
            # println("derivatives could commute")
            return derivatives_commute(a::DerivativeSuperType, b::DerivativeSuperType)
        else
            return false
        end
    end
    return false
end


function ==(a::DerivativeSuperType, b::DerivativeSuperType)#check to see if derivatives commute
    if a.name == b.name && a.operator_indices == b.operator_indices && a.operand == b.operand #basic check
        return true
    end
    opra = isTensor(a.operand)
    oprab = isTensor(b.operand)
    if oprab != false && opra != false && typeof(opra) <: DerivativeSuperType && typeof(oprab) <: DerivativeSuperType #nested ops

        # println("checking nested derivatives")
        # println(opra)        
        # println(oprab)
        #check to see if the operand of one derivative is the other, and if the operators are the same (if they commute then this is true)
        if opra.operator_indices == b.operator_indices && oprab.operator_indices == a.operator_indices && opra.name == b.name && oprab.name == a.name && opra.operand == oprab.operand
            # println("derivatives could commute")
            return derivatives_commute(a::DerivativeSuperType, b::DerivativeSuperType)
        else
            return false
        end
    end
    return false
end



function signed_equals(a::DerivativeSuperType, b::DerivativeSuperType)#check to see if derivatives commute
    if a.name == b.name && a.operator_indices == b.operator_indices && a.operand == b.operand #basic check
        return (true, 1)
    end
    opra = isTensor(a.operand)
    oprab = isTensor(b.operand)
    if oprab != false && opra != false && typeof(opra) <: DerivativeSuperType && typeof(oprab) <: DerivativeSuperType #nested ops

        # println("checking nested derivatives")
        # println(opra)        
        # println(oprab)
        #check to see if the operand of one derivative is the other, and if the operators are the same (if they commute then this is true)
        if opra.operator_indices == b.operator_indices && oprab.operator_indices == a.operator_indices && opra.name == b.name && oprab.name == a.name && opra.operand == oprab.operand
            # println("derivatives could commute")
            return derivatives_commute_signed(a::DerivativeSuperType, b::DerivativeSuperType)
        else
            return (false, 0)
        end
    end
    return (false, 0)
end



function register_derivative(name::String, functionof::Union{Coordinate,Vector{Coordinate}}, operator_indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}()) where {I<:IndexSuperType}
    return Derivative(name, zero, functionof, operator_indices)
end


#logic needs to be more complex, since we want not literal matches but symbolic matches
function register_derivative_commutation(D1::DerivativeSuperType, D2::DerivativeSuperType, res::TensorExpression, grading::Integer)
    derivative_algebra[D1=>D2] = grading => res
    derivative_algebra[D2=>D1] = grading => -res
end

function register_derivative_commutation(D1::DerivativeSuperType, D2::DerivativeSuperType, res::TensorExpression)
    derivative_algebra[D1=>D2] = 1 => res
    derivative_algebra[D2=>D1] = 1 => -res
end

function register_derivative_anticommutation(D1::DerivativeSuperType, D2::DerivativeSuperType, res::TensorExpression)
    derivative_algebra[D1=>D2] = -1 => res
    derivative_algebra[D2=>D1] = -1 => res
end

function register_derivative_commutation(D1::DerivativeSuperType, D2::DerivativeSuperType)
    derivative_algebra[D1=>D2] = 1 => zero
    derivative_algebra[D2=>D1] = 1 => zero
end

function register_derivative_anticommutation(D1::DerivativeSuperType, D2::DerivativeSuperType)
    derivative_algebra[D1=>D2] = -1 => zero
    derivative_algebra[D2=>D1] = -1 => zero
end

function derivatives_commute_signed(D1::DerivativeSuperType, D2::DerivativeSuperType)
    if typeof(D1) <: PartialDerivative && typeof(D2) <: PartialDerivative
        # println(D1)
        # println(outer_grading(D1))
        # println(D2)
        # println(outer_grading(D2))
        #use grading of the derivative itself, not the grading of full object since we commute ops not expressions
        if (outer_grading(D1) == 0 && outer_grading(D2) == 0) || (outer_grading(D1) != outer_grading(D2))
            return (true, 1) #assume partial derivatives commute
        else
            return (true, -1) #fermionic derivatives anticommute... add anti behaviour later
        end
    end
    if haskey(derivative_algebra, D1 => D2)
        return derivative_algebra[D1=>D2][1] == (1, zero)
    else
        return (false, 0)
    end
end

function derivatives_commute(D1::DerivativeSuperType, D2::DerivativeSuperType)
    return derivatives_commute_signed(D1, D2) == (true, 1)
end


function operators_commute(D1::DerivativeSuperType, D2::DerivativeSuperType)
    if haskey(derivative_algebra, D1 => D2)
        return derivative_algebra[D1=>D2][1] == 1
    else
        return false
    end
end

function same_operator(D1::T, D2::T2) where {T<:DerivativeSuperType,T2<:DerivativeSuperType}
    ind1 = D1.operator_indices
    ind2 = D2.operator_indices
    doubles = Dict{Index,Index}()

    if D1.name != D2.name
        return (false => doubles)
    end

    #quick check for inequality
    if length(ind1) != length(ind2)
        return (false => doubles)
    end

    # check if indices are the same in each position, modulo names of pairs of dummies
    for ind in eachindex(ind1)
        if ind1[ind][1].index_set == ind2[ind][1].index_set && ind1[ind][2].position == ind2[ind][2].position
            if haskey(doubles, ind1[ind][1])
                val = doubles[ind1[ind][1]]
                if val != ind2[ind][1]
                    #if the patterns don't match ie repeated indices etc
                    return (false => doubles)
                end
            else
                doubles[ind1[ind][1]] = ind2[ind][1]
            end
        else
            #mistmatch
            return (false => doubles)
        end

    end
    return (true => doubles)
end


#define action of a derivative operator on some basis of tensors


function add_derivative_action(D::T, exp::TensorTerm, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}

    if length(exp.terms) != 1
        throw(InvalidDerivativeBasis())
    else
        add_derivative_action(D, exp.terms[1], res)
    end

end

function add_derivative_action(D::T, exp::TensorExpression, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}

    if length(exp.summands) != 1 || length(exp.summands[1][2].terms) != 1
        throw(InvalidDerivativeBasis())
    else
        add_derivative_action(D, exp.summands[1][2].terms[1], res)
    end

end

function add_derivative_action(D::T, exp::TensorSuperType, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}
    key = nothing
    in_dict = false
    pattern_ref = Dict()
    res = TensorExpression(res)
    for ky in keys(derivative_actions)
        out = same_operator(D, ky)
        B = out[1]
        o = out[2]
        if B == true
            key = ky
            in_dict = true
            pattern_ref = o
            break
        end
    end

    #check to make sure action not already defined, else add

    #now we have index pattern if we found a match
    if in_dict
        action_dict = derivative_actions[key] #all registered derivatives to act on
        for k in keys(action_dict)#looks through all tensors
            if same_tensor_mod_pattern(exp, k, pattern_ref)[1]
                return false #already in the dict
            end
        end
        #no matches found.... proceed

        #need to make sure correct indices
        pat_adj_tensor = index_to_pattern(exp, pattern_ref)
        pat_adj_out = index_to_pattern(res, pattern_ref)
        derivative_actions[key][pat_adj_tensor] = pat_adj_out #adds the expression
        return true
    else
        derivative_actions[D] = Dict()
        derivative_actions[D][exp] = res #adds the expression
        return true
    end
end


#given an expression/tensor rename/relabel indices according to a particular pattern

function index_to_pattern(T1::TensorExpression, pattern_ref::Dict{Index,Index})
    new_exp = zero
    for summand in T1.summands
        new_summand = summand[1]
        for term in summand[2].terms
            new_summand *= index_to_pattern(term, pattern_ref)
        end
        new_exp += new_summand
    end
    return new_exp
end


function index_to_pattern(T1::TensorSuperType, pattern_ref::Dict{Index,Index})
    newT = copy(T1)
    ind1 = newT.indices
    for ind in eachindex(ind1)
        if haskey(pattern_ref, ind1[ind][1])
            new_ind = pattern_ref[ind1[ind][1]]
            ind1[ind] = new_ind => ind1[ind][2]
        end
    end
    return newT
end


# checks if tensors are the same in terms of structure of indices, ignoring the specfic name
# eg differentiating between contracted and non contracted

function same_tensor_mod_pattern(T1::T, T2::T_2, index_dict::Dict{Index,Index}) where {T<:TensorSuperType,T_2<:TensorSuperType}

    ind1 = T1.indices
    ind2 = T2.indices
    if length(ind1) != length(ind2)
        return (false, Dict())
    end
    temp_dict = Dict{Index,Index}()#stores non-dummies
    for ind in eachindex(ind1)
        if ind1[ind][2] != ind2[ind][2]
            return (false, Dict()) #wrong positions
        end
        if ind1[ind][1].index_set != ind2[ind][1].index_set
            return (false, Dict()) #wrong index set
        end
        if haskey(index_dict, ind1[ind][1])#check if already in pattern
            temp = index_dict[ind1[ind][1]]
            if temp != ind2[ind][1]
                return (false, Dict())  #after subbing appropritate pattern, the index has a different name
            end
            #not a dummy index
        elseif haskey(temp_dict, ind1[ind][1])
            temp = temp_dict[ind1[ind][1]]
            if temp != ind2[ind][1]
                return (false, Dict()) #after subbing appropritate pattern, the index has a different name
            end
        else
            temp_dict[ind1[ind][1]] = ind2[ind][1] #matches "free" (non derivative) indices up, this ensures same pattern
        end
    end
    return (true, temp_dict)
end



##############################################################################################################################
##############################################################################################################################

#applying derivatives



#leibniz out the operands of all derivatives in an expression

function apply_derivative(TEX::TensorExpression)
    out = zero
    for summand in TEX.summands
        nt = summand[1]
        for term in summand[2].terms
            if typeof(term) <: DerivativeSuperType
                # println("applying derivative to term", term)
                nt *= apply_derivative(term)
            else
                nt *= term
            end

        end
        out = TensorExpression(out)
        # println("nt is ",nt)
        # println(out)
        out += TensorExpression(nt)
    end
    return TensorExpression(out)
end


#apply the derivative and cancel terms that are annihilated

function apply_derivative(D::DerivativeSuperType)
    D = product_rule(D)
    out = zero
    for summandind in eachindex(D.summands) # for each summand (coefficient, tensorterm)
        summand = D.summands[summandind]
        term = summand[1]
        # println("applying derivative to summand ", summand)
        for tensor in summand[2].terms # for each individual tensor in tensorterm
            if typeof(tensor) <: DerivativeSuperType # if the tensor is an operator
                operands = tensor.operand #get the operand (should only be one thanks to prod rule)
                tens = operands.summands[1][2].terms[1] #get the only item in the operand
                funcof = get_function_of(tens)
                if isempty(intersect(funcof, Set(tensor.coord))) #see if tensor is a function of the derivative coords
                    # println("annihilated by derivative since not a function of the coord")
                    # println("tens is ", tens)
                    # println("derivative is ", tensor)
                    # println("function of is ", tens.function_of)
                    # println("coord is ", tensor.coord)
                    term *= zero #annihilate the term   
                else
                    out1 = perform_derivative_action(tensor)
                    term *= out1
                end
            else
                term *= tensor
            end
        end

        if term == 0
            term = zero
        end
        if out == 0
            out = zero
        end
        out += term

    end
    out = TensorExpression(out)
    return merge_sum(out)
end


function get_function_of(T::TensorSuperType)
    if typeof(T) <: DerivativeSuperType
        return get_function_of(T.operand)
    elseif hasproperty(T, :function_of)
        return T.function_of
    else
        return Coordinate[]
    end
end

function get_function_of(T::TensorExpression)
    coords = Coordinate[]
    for (_, term) in T.summands
        for t in term.terms
            append!(coords, get_function_of(t))
        end
    end
    return coords
end

"""
Performs the actual derivative operation on tensors for which it is defined.
If there are rules (as defined by add_derivative_action), then the derivative will be performed, else expressions will be left alone eg ∇B
"""
function perform_derivative_action(D::DerivativeSuperType)
    operand = D.operand.summands[1][2].terms[1] #get the operand
    for ky in keys(derivative_actions)
        A = same_operator(D, ky)
        B = A[1]
        o = A[2]
        # println(B)
        # println(o)
        if B
            action_dict = derivative_actions[ky] # dict of tensors we've defined deriv on
            for k in keys(action_dict)#looks through all tensors
                op = same_tensor_mod_pattern(operand, k, o)
                if op[1]
                    # println("found a match")
                    # println(D)
                    # println(operand)
                    # println(ky)
                    # println(op)
                    # println("same op out")
                    # println(A)
                    #do derivative:
                    replacement = action_dict[k] # need to go through and replace this with correct pattern
                    # println("replacement is")
                    # println(replacement)
                    inverted_dict = Dict(v => k for (k, v) in o)
                    inverted_dict2 = Dict(v => k for (k, v) in op[2])

                    all_invs = merge(inverted_dict, inverted_dict2)
                    # pat_rep = index_to_pattern(replacement, inverted_dict)
                    # println(pat_rep)
                    # pat_rep = index_to_pattern(pat_rep, inverted_dict2)
                    pat_rep = index_to_pattern(replacement, all_invs)

                    # println(pat_rep)
                    return pat_rep
                end
            end

        end

    end
    return D
end


"""
expand_derivative(TEX,D,DExpansion)

Replaces a symbolic derivative with its full form, eg substituting a covariant derivative D with all of its connection terms.
'TEX' is the original expression. 'D' is the derivative to replace, and 'DExpansion' is a tensor expression that will be a sum of tensors/operators.
For each operator in 'DExpansion', the operand will be replaced by the matching operand in 'TEX'.
"""
function expand_derivative(TEX::Union{TensorExpression,<:DerivativeSuperType}, D::DerivativeSuperType, DExpansion::TensorExpression, outertex=nothing)

    #### keeping all the comments here in case something breaks in future

    if (get_free_indices(D) != get_free_indices(DExpansion))
        println(D)
        println(DExpansion)
        println(TEX)
        println("It appears that your expansion is not compatible.")
        throw(FreeIndexMismatch())
    end

    #think about: hellish dummie indice relabeling//////

    new_TEX = zero #construct new expression term by term
    TEX = 1 * TEX #ensures tensor exp
    avail = Dict{IndexSet,OrderedSet{Index}}() #pairs of index sets, with the unused indices
    dummymap = Dict{String,Index}()#dict of what to replace pattern indices with
    for summand in TEX.summands
        new_summand = summand[1]
        for tensor in summand[2].terms #look at each tensor in the tensor term
            # println("checking tensor")
            # println(tensor)
            if typeof(tensor) <: DerivativeSuperType #apply to operands first recursively
                # println("expanding operand")
                # println("expanding operand")
                # println(tensor.operand)
                # println("im gunna recooooooorse")
                # println(tensor)
                # println(typeof(tensor))
                # println(typeof(D))
                if outertex === nothing
                    tensor.operand = expand_derivative(tensor.operand, D, DExpansion, TEX)
                else
                    tensor.operand = expand_derivative(tensor.operand, D, DExpansion, outertex)
                end
                #feed in term... need to make sure we dont reuse outer indices
                #recursively expand the operand if it has the same derivative in it
                update_ext_indices(tensor)
                # println("did a hecking recursion")
                # println(tensor)
                # println(get_inds(tensor))
            end
            if typeof(tensor) == typeof(D) && tensor.name == D.name #operator found
                operand = tensor.operand
                # println(TensorExpression(D))
                # println(TensorExpression(tensor))
                mapping = construct_index_rep_map(D, tensor) #index map so we use correct indices in expansion
                # println("deriv index remap:")
                # println(mapping)
                # println("mapping ", mapping)
                replace = copy(DExpansion) #new copy of expansion to replace in our exp

                #figure out what indices we're allowed to use in the replacement
                sets = get_index_sets(tensor) #Set{IndexSet}()
                if outertex === nothing
                    ref = TEX
                else
                    ref = outertex
                end
                for set in sets

                    a1 = OrderedSet{Index}(get_unused_indices(ref, set))
                    b = OrderedSet{Index}(get_unused_indices(tensor.operand, set)) ##will check for dummies inside...
                    a = intersect(a1, b)
                    # println("unused")
                    # println(a)
                    # println("tex v ref")
                    # println(TEX)
                    # println(ref)
                    push!(avail, set => a) #ones we can replace with
                end
                #construct the replacement
                new_rep_summand = zero
                for repsummand in replace.summands
                    new_rep_term = repsummand[1]


                    ###############fix index logic for multiple tensors in the expressions
                    ## currently, K^{*β}_{α*}∂^{i*}_{*β}[[0*]]   => goes to K^{*β}_{α*}∂^{i*}_{*β}[[F^{ik**}_{**jβ}]], where F^{ik**}_{**jβ}

                    found_op = false
                    for reptensor in repsummand[2].terms ##add somewhere to handle a term with no op (like a YM covar).
                        ##println(reptensor)
                        if typeof(reptensor) <: OperatorSuperType
                            found_op = true
                            # #if we found an operator,
                            # println("operators")
                            # println("Running on operand")
                            # println(operand)    
                            # println("we will try and expand",D)                        
                            # println(expand_derivative(operand,D,DExpansion))
                            # println("^result")
                            reptensor.operand = operand#set to correct operand

                            #dums=get_dummies(reptensor)#gets the placeholder dummy indices from the pattern
                            for pair in eachindex(reptensor.operator_indices) # replace dummies with correct etc
                                #println(reptensor.operator_indices[pair][1])
                                # println("checking")
                                # println(reptensor.operator_indices[pair][1].name)

                                if haskey(mapping, reptensor.operator_indices[pair][1].name) #if the index is fixed already from before
                                    # println("replacing")
                                    # println(reptensor.operator_indices[pair][1].name)

                                    ##println("SYMBOILIC INDEX ALREADY FIXED by pattern:")
                                    ##println(reptensor.operator_indices[pair][1])
                                    ##println(mapping[reptensor.operator_indices[pair][1].name])
                                    ##println(typeof(reptensor.operator_indices[pair]))
                                    # #println(typeof(mapping[reptensor.operator_indices[pair][1].name]))
                                    reptensor.operator_indices[pair] = mapping[reptensor.operator_indices[pair][1].name] => reptensor.operator_indices[pair][2]#should replace with the right thing
                                #println("maps to ",mapping[reptensor.operator_indices[pair][1].name])
                                elseif haskey(dummymap, reptensor.operator_indices[pair][1].name)
                                    ##println("there's already a dummy assigned")
                                    #println(haskey(dummymap,reptensor.operator_indices[pair][1].name))
                                    #println(dummymap[reptensor.operator_indices[pair][1].name])
                                    #old = reptensor.operator_indices[pair][1].name
                                    # #println(, " maps to ",dummymap[reptensor.operator_indices[pair][1].name])
                                    # #println(dummymap)
                                    reptensor.operator_indices[pair] = dummymap[reptensor.operator_indices[pair][1].name] => reptensor.operator_indices[pair][2]
                                    # #println(dummymap)
                                    # #println(dummymap[reptensor.operator_indices[pair][1].name])
                                    #println(old," maps to ",reptensor.operator_indices[pair][1].name)
                                else #if just a placeholder index in the pattern, take from the list
                                    ##println(avail)
                                    if !haskey(avail, reptensor.operator_indices[pair][1].index_set)
                                        #if we're using a dummy index type thats not in the orignal expression
                                        #changed to ref!!
                                        a11 = OrderedSet{Index}(get_unused_indices(ref, reptensor.operator_indices[pair][1].index_set))
                                        b1 = OrderedSet{Index}(get_unused_indices(tensor.operand, reptensor.operator_indices[pair][1].index_set)) ##will check for dummies inside...
                                        aa = intersect(a11, b1)
                                        new_index = pop!(aa)
                                    else
                                        new_index = pop!(avail[reptensor.operator_indices[pair][1].index_set]) # get an unused index
                                    end
                                    dummymap[reptensor.operator_indices[pair][1].name] = new_index #makes sure summation dummies do the right thing
                                    reptensor.operator_indices[pair] = new_index => reptensor.operator_indices[pair][2]
                                    #println("maps to ",new_index)
                                end
                                # for repindex in reptensor.
                            end
                        else
                            ##println("nonoperator")
                            for pair in eachindex(reptensor.indices) # replace dummies with correct etc
                                if haskey(mapping, reptensor.indices[pair][1].name) #if the index is fixed already from before
                                    # #println(reptensor.operator_indices[pair][1])
                                    # #println(mapping[reptensor.operator_indices[pair][1].name])
                                    # #println(typeof(reptensor.operator_indices[pair]))
                                    # #println(typeof(mapping[reptensor.operator_indices[pair][1].name]))
                                    reptensor.indices[pair] = mapping[reptensor.indices[pair][1].name] => reptensor.indices[pair][2]#should replace with the right thing
                                elseif haskey(dummymap, reptensor.indices[pair][1].name)
                                    reptensor.indices[pair] = dummymap[reptensor.indices[pair][1].name] => reptensor.indices[pair][2]

                                else #if just a placeholder index in the pattern, take from the list
                                    if !haskey(avail, reptensor.indices[pair][1].index_set)
                                        #if we're using a dummy index type thats not in the orignal expression
                                        ##println("getta new index for index sets")
                                        #changed to ref!!!
                                        a11 = OrderedSet{Index}(get_unused_indices(ref, reptensor.indices[pair][1].index_set))
                                        b1 = OrderedSet{Index}(get_unused_indices(tensor.operand, reptensor.indices[pair][1].index_set)) ##will check for dummies inside...
                                        aa = intersect(a11, b1)
                                        new_index = pop!(aa)
                                        # new_index = pop!(get_unused_indices(ref, reptensor.indices[pair][1].index_set))
                                        ##println("new index, ", new_index)
                                    else
                                        ##println("none yet taken")
                                        new_index = pop!(avail[reptensor.indices[pair][1].index_set]) # get an unused index
                                        ##println(new_index)
                                    end
                                    #new_index = pop!(avail[reptensor.indices[pair][1].index_set]) # get an unused index
                                    dummymap[reptensor.indices[pair][1].name] = new_index #makes sure summation dummies do the right thing
                                    reptensor.indices[pair] = new_index => reptensor.indices[pair][2]

                                end
                            end
                        end

                        # println("prev new rep: ",new_rep_term) 
                        # println(typeof(new_rep_term))   
                        # println("multiplying by: ",reptensor)
                        #println(typeof(reptensor))
                        new_rep_term *= reptensor
                        # println("new_rep_term is ",new_rep_term)

                    end
                    # println(found_op)
                    # println("new rep summand is ",new_rep_term)
                    #  update_ext_indices(new_rep_term)
                    #println("old summand is ",new_rep_summand)
                    #println("adding ",new_rep_term)
                    # println("new rep term is ",new_rep_term)
                    update_ext_indices(new_rep_term)
                    update_ext_indices(new_rep_summand)
                    if found_op == false
                        new_rep_term *= operand #if there was no operator, we need to add the operand at the end, since the pattern is just for the operator part
                    end
                    # println("updated ",new_rep_term)
                    ##println(new_rep_term)
                    ##println(get_free_indices(new_rep_term))
                    ##println(new_rep_summand)
                    ##println(get_free_indices(new_rep_summand))
                    new_rep_summand += new_rep_term
                end
                # println("timesing here")
                new_summand *= new_rep_summand
                # println(new_summand)
            else#normal tensor
                # println("norsml tenstimesing here22222")
                # println(tensor)
                # println(typeof(tensor))
                # println(D)
                # println(typeof(D))
                new_summand *= tensor
                # println(new_summand)
            end
        end
        #        update_ext_indices(new_summand)
        #println("adding here22222")
        new_TEX += new_summand #*summand[1] * #might be squaring...
        #println(new_TEX)
    end
    # println(D)
    # println("final new tex is")
    # println(new_TEX)
    return new_TEX
end

# graded leibniz

function product_rule(TEX::TensorExpression)
    new_summand = zero
    for summand in TEX.summands
        nterm = summand[1]
        for term in summand[2].terms
            ###println("doing this1")
            if typeof(term) <: LeibnizOperatorSuperType
                ###println("doing this")
                nterm *= product_rule(term)

            else
                nterm *= term
            end
        end
        new_summand += nterm
    end
    return new_summand
end


function product_rule(O::LeibnizOperatorSuperType)
    new_summand = zero
    ##println("ns: ",new_summand)
    ##println(zero)
    ##println("operator is")
    ##println(O)
    for summandind = eachindex(O.operand.summands) #for each sum term in the expression
        summand = O.operand.summands[summandind] #get the summand (tensor term)

        # #this should hopefully handle nested operators :)))
        # if length(summand[2].terms)==1 && typeof(summand[2].terms[1])<:LeibnizOperatorSuperType
        #     summand=product_rule(summand[2].terms[1])
        # end


        # for summ in summand.summands


        # end
        grad = 0 #initialise grading
        #    println("summand is: ",summand)
        ##println("zero is , ", zero)
        for tensorind in eachindex(summand[2].terms)# for each tensor in the summand (tensorterm)
            grd = (-1)^(grad)
            operator_on_tensor = typeof(O)(O.name, summand[2].terms[tensorind], O.coord, O.operator_indices) #operator on nth term
            #     println("this is the op on a tens")
            # println(operator_on_tensor)
            # #println("help")
            new_prod = 1
            #    println(summand[2].terms)
            ##println("zero is now , ", zero)
            for indsss in eachindex(summand[2].terms)
                ###println("looking a thte ", indsss," summand: ", summand[2].terms[indsss])
                if indsss != tensorind
                    ###println(new_prod)
                    ###println(summand[2].terms[indsss])
                    # println(new_prod)
                    # println(summand[2].terms[indsss])
                    new_prod = new_prod * summand[2].terms[indsss]
                    ##println("me when zero is , ", zero)
                else
                    # update_ext_indices(operator_on_tensor)
                    ###println("updated ext inidices to: ",operator_on_tensor.indices)
                    new_prod = new_prod * operator_on_tensor
                    ##println("zero is wwaasfd, ", zero)
                end
            end
            #    println("we am HERE")
            #     # new_prod=[t  where t=summand[2].terms]
            #    println(new_prod)
            # new_prod.terms[tensorind]=operator_on_tensor
            if new_summand == Nothing
                new_summand = zero + grd * summand[1] * new_prod #adding the new term to the leibniz expanded version of the summand
            #    println(new_summand, "lets see if meessing up")
            else
                ##println(new_summand)
                ##println(new_summand.summands)
                ##println(typeof(new_summand))
                ##println("try me")
                ##println(grd)
                ##println( summand[1])
                ##println(new_prod)
                ##println(grd * summand[1] * new_prod)
                ##println(typeof(grd * summand[1] * new_prod))
                new_summand = new_summand + grd * summand[1] * new_prod #adding the new term to the leibniz expanded version of the summand
                #    println(new_summand, "addition problem??")
            end

            grad = (grad + grading(summand[2].terms[tensorind]) * op_grading(O)) % 2
            ##println("where..........")
        end

    end
    return new_summand
end

function symmetrise_partials(TEX::TensorExpression)
    new_summand = zero

    for summand in TEX.summands
        nterm = summand[1]
        for term in summand[2].terms
            if typeof(term) <: PartialSuperType
                symd, remap = symmetrise_partials(term)
                temp = zero
                # temp += (term) / length(symd)
                # println("symd is ", symd)
                for par in symd
                    # println("par is ", par)
                    temp += (nterm * par) / length(symd) #divide by number of terms to symmetrise
                    # println("temp is ", temp)
                end
                # println("nnnnnnterms is ", nterm)
                nterm = temp
                # println("now nterm is ", nterm)
            else
                nterm *= term
            end
        end
        # println("nterms is ", nterm)
        new_summand += nterm
        # println("new summand is ", new_summand)
    end
    # println("symmetrised partials is ", new_summand)
    # println("original was ", TEX)
    # println("averaged is ", (new_summand + TEX) / 2)
    return (new_summand + TEX) / 2
end

function symmetrise_partials(TT::TensorTerm)
    new_term = zero
    for term in TT.terms
        if typeof(term) <: PartialSuperType
            (trm,remap)=symmetrise_partials(term)
            new_term *= trm
        else
            new_term *= term
        end
    end
    # println("symmetrised partials is NT ", new_term)
    # println(remap)
    return new_term
end

function root_oprnd(D::PartialSuperType)
    opr = isTensor(D.operand)
    if typeof(opr) <: PartialSuperType
        return root_oprnd(opr)
    else
        return D.operand
    end
end

function symmetrise_partials(D::PartialSuperType)
    # println("get_op",get_operator_chain(D))
    root = root_oprnd(D)
    # println("root is ", root)
    chn = get_operator_chain(D) #get the chain of ops
    # println("operator chain:")
    # println(chn)
    inds = [op.operator_indices[1][1] for op in chn]
    # println("indices are ", inds)
    swps = graded_swaps(chn)
    # println("graded swaps:")
    # println(swps)
    term = []
    indpermsfor = graded_swaps_opinds(chn)
    indpermsforpass = [(inds[swp[1]], swp[2]) for swp in indpermsfor]
    # println("graded swaps for indices:")
    # println(indpermsforpass)
    for swp in swps
        # println("permuted is ", swp[1], " with sign ", swp[2])
        permuted = swp[1]
        sign = swp[2]
        push!(term, sign * rebuild_nested(permuted, root))
    end
    # println("symmetrised term is ", term)
    return (term,indpermsforpass)
    # println("symmetrising ", D)
    # println("operand is ", D.operand)
    # tyot=isTensor(D.operand)
    # if typeof(tyot) <: PartialSuperType
    #     # println("got a nested one")
    #     new_op = symmetrise_derivs(D.operand)
    #     # println("new op is ", new_op)
    #     return new_op
    # else
    #     return D
    # end
end

function rebuild_nested(permuted, root)
    expr = root

    for op in reverse(permuted)
        expr = typeof(op)(
            op.name,
            expr,
            op.coord,
            op.function_of,
            op.operator_indices
        )
    end

    return expr
end


function graded_sign(perm, ops)
    sign = 1

    for i in 1:length(perm)
        for j in i+1:length(perm)

            # inversion in permutation
            if perm[i] > perm[j]

                gi = grading(ops[perm[i]].operator_indices)
                gj = grading(ops[perm[j]].operator_indices)

                sign *= (-1)^(gi * gj)
            end
        end
    end

    return sign
end

function graded_swap(ops, perm)
    permuted = ops[perm]
    sign = graded_sign(perm, ops)

    return (permuted, sign)
end

function graded_swaps(ops)
    results = []

    n = length(ops)

    for p in permutations(1:n)
        permuted = ops[p]
        sign = graded_sign(p, ops)

        push!(results, (permuted, sign))
    end

    return results
end

function graded_swaps_opinds(ops)
    results = []

    n = length(ops)

    for p in permutations(1:n)
        # println("permutation: ", p)
        permuted = ops[p]
        sign = graded_sign(p, ops)

        push!(results, (p, sign))
    end

    return results
end

function get_operator_chain(T::OperatorSuperType)
    chain = []
    current = T
    while typeof(current) <: OperatorSuperType
        push!(chain, current)
        current = isTensor(current.operand) #get the operand tensor
    end
    return chain
end



"""Creates a new derivative tensor, with multiple indices if the operand is itself a derivative.
Hopefully makes things easier when it comes to simplify algos.
"""
function to_nested(D::PartialSuperType, top=true)
    inds = D.operator_indices
    # println(D.operand)
    # println(isTensor(D.operand))
    # println(typeof(D.operand))
    tyot = isTensor(D.operand)
    if tyot != false && typeof(tyot) <: PartialSuperType
        println("got")
        tmp = to_nested(tyot, false)
        append!(inds, tmp[1])
        if top == false
            return (inds, tmp[2]) #recursive step
        else
            # println(inds)
            return Derivative("D" * string(length(inds)), tmp[2], D.coord, D.function_of, inds) #base case
        end
    elseif top == false
        return (inds, D.operand) ## this is not a deriv
    else
        return D
    end
end


"""
Deep check to see if two operators are the same, by recursively checking if their operands are the same in the case of nested operators.
"""
function nested_opertor_equality(D1::DerivativeSuperType, D2::DerivativeSuperType)
    if D1.name == D2.name && typeof(D1.operand) == typeof(D2.operand) && D1.coord == D2.coord && D1.function_of == D2.function_of && D1.operator_indices == D2.operator_indices
        return true
    else
        return false
    end
end



function update_ext_indices(O::OperatorSuperType) #probably needs to include dummies...
    if typeof(O.operand) <: TensorSuperType || typeof(O.operand) <: TensorTerm
        O.indices = vcat(O.operator_indices, O.operand.indices)
    elseif typeof(O.operand) <: TensorExpression
        O.indices = vcat(O.operator_indices, O.operand.summands[1][2].indices)
    end
end


function outer_grading(T::OperatorSuperType)
    grade = 0
    for (ind, pos) in T.operator_indices
        grade = grade + ind.grading
    end
    return grade % 2
end
function *(TEx::TensorExpression, T::OperatorSuperType)
    # println("doing this sort of multiplcation5")
    return update_indices(TEx * TensorExpression(T))
end

function *(T::OperatorSuperType, TEx::TensorExpression)
    # println("doing this sort of multiplcation6")
    return update_indices(TensorExpression(T) * TEx)
end

function copy(T::DerivativeSuperType) #this needs to actually be a deep copy...
    return typeof(T)(T.name, copy(T.operand), copy(T.coord), copy(T.function_of), copy(T.operator_indices))
end