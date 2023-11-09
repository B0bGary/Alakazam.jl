# include("operator.jl")

##############################################################################################################################
##############################################################################################################################

#derivative rules


derivative_actions = Dict{DerivativeSuperType,Dict{TensorSuperType,TensorExpression}}() #where {D<:}

# function is_pattern_equivalent(D_pat<:DerivativeSuperType,exp_pat::TensorExpression,D_stored<:DerivativeSuperType,exp_stored::TensorExpression)

# end


function same_operator(D1::T, D2::T2) where {T<:DerivativeSuperType,T2<:DerivativeSuperType}
    ind1 = D1.operator_indices
    ind2 = D2.operator_indices
    doubles = Dict{Index,Index}()#####
    if length(ind1)!=length(ind2)
        return (false => doubles)
    end
    
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

function add_derivative_action(D::T, exp::TensorTerm, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}
    
    if length(exp.terms)!=1
        throw(InvalidDerivativeBasis())
    else
        add_derivative_action(D,exp.terms[1],res)
    end

end

function add_derivative_action(D::T, exp::TensorExpression, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}
    
    if length(exp.summands)!=1 || length(exp.summands[1][2].terms)!=1
        throw(InvalidDerivativeBasis())
    else
        add_derivative_action(D,exp.summands[1][2].terms[1],res)
    end

end

function add_derivative_action(D::T, exp::TensorSuperType, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}
    key = nothing
    in_dict = false
    pattern_ref = Dict()
    res = TensorExpression(res)
    for ky in keys(derivative_actions)
        out= same_operator(D, ky)
        B=out[1]
        o=out[2]
        if B == true
            key = ky
            in_dict = true
            pattern_ref = o
            break
        end
    end

    #now we have index pattern if we found a match
    if in_dict
        action_dict = derivative_actions[key] #all registered derivatives to act on
        for k in keys(action_dict)#looks through all tensors
            if same_tensor_mod_pattern(exp, k, pattern_ref)[1]
                return false #already in the dict
            end
        end
        #no matches found....

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
            # println("yeahh hthis might be it")
            # new_ind=ind1[ind][1]
            # println(":O")
            new_ind = pattern_ref[ind1[ind][1]]
            ind1[ind]=new_ind =>ind1[ind][2]
        end
    end
    return newT
end


function same_tensor_mod_pattern(T1::T, T2::T_2, index_dict::Dict{Index,Index}) where {T<:TensorSuperType,T_2<:TensorSuperType}
    # println("this functio nrann")
    ind1 = T1.indices
    ind2 = T2.indices
    if length(ind1)!=length(ind2)
        return (false, Dict())
    end
    temp_dict = Dict{Index,Index}()#stores non-dummies
    for ind in eachindex(ind1)
        # println(ind)
        if ind1[ind][2]!=ind2[ind][2]
            return (false, Dict()) #wrong positions
        end
        if ind1[ind][1].index_set!=ind2[ind][1].index_set
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



# #is the derivative registered by name 
# if haskey(derivative_actions,D.name)
#     deriv_index_options=derivative_actions[D.name]




#     for option in collect(keys(deriv_index_options)) #for all index patterns
#         #check and see if there's a match or nah
#         #have to check: positions and index sets, not the names!
#         for ind in eachindex(option)
#             if option[ind][2].position!=D.operator_indices[ind][2].position || option[ind][1].index_set!=D.operator_indices[ind][1].index_set
#                 in_dict=false
#                 break
#             end
#         end
#         #found a mistatch
#         if !in_dict

#         else
#             #this option is a match
#             key=option

#         end

#         elseif in_dict
#             # if there were no mismatches, then our thing is already there
#             #thus, indices are a key ->>>>>>>>>>>> yeah nah ok yeah i got the key now
#             patterns=deriv_index_options[key]
#             #this should be all the patterns. now just go through and find blah blah pattern
#             #easier said than done θ^i_j is not the same as θ^i_i=θ^j_j but is the same as θ^j_i and θ^k_l
#             #maybe have a think about a better way of doing this

#             return false

#         end
#     end
#     #if we get here, we checked everything and found no matches
#     old_dict=derivative_actions[D.name]
#     old_dict[D.operator_indices]=res

# else
#     derivative_actions[D.name]=...
# end









##############################################################################################################################
##############################################################################################################################

function apply_derivative(TEX::TensorExpression)
    out=zero
    for summand in TEX.summands
        nt=summand[1]
        for term in summand[2].terms
            if typeof(term)<:DerivativeSuperType
              
                nt*=apply_derivative(term)
            else
                nt*=term
            end
            
        end
        out=TensorExpression(out)
        out+=TensorExpression(nt)
    end
        return TensorExpression(out)
end



#apply the derivative and cancel terms that are annihilated

function apply_derivative(D::DerivativeSuperType)
    D = product_rule(D)
    # println("ran product rule")
    out = zero
    # new_summands=D.summands
    for summandind in eachindex(D.summands) # for each summand (coefficient, tensorterm)
        summand = D.summands[summandind]
        term = summand[1]
        #println(term)
        for tensor in summand[2].terms # for each individual tensor in tensorterm
            if typeof(tensor) <: DerivativeSuperType # if the tensor is an operator
                operands = tensor.operand #get the operand (should only be one thanks to prod rule)
                #println(operands, "should be a single term else liebniz wrong")
                tens = operands.summands[1][2].terms[1] #get the only item in the operand
                if isempty(intersect(Set(tens.function_of), Set(tensor.coord))) #see if tensor is a function of the derivative coords
                #    println(tens)
                #    print
                    term*= zero #annihilate the term
                #    println("term:: ",term)
                else
                    #println("old term pre: ", term)
                    # println("performing action")
                    # println(tensor)
                    out1 = perform_derivative_action(tensor)
                    #println("output of perofrming deriv action,", out1)
                    # # do the subs rules
                    # for ky in keys(derivative_actions)
                    #     (B=>o)=same_operator(D,ky) ################ complicated logic....
                    # end
                    #println("old term", term)
                    term *= out1
                    #println("new term is ",term )
                end
            else
                # println("problem zone: term is ", term)
                # println("tensor is, ",tensor)
                term *= tensor
                #println("term is ",term )
            end
        end
        #println(out)
        #println(term)
        #println("adding something")
        if term==0
            term=zero
        end
        if out==0
            out=zero
        end
        out += term
       ##println("scuffma:",out)
    end
    out=TensorExpression(out)
   ##println("out is:", out)

    return merge_sum(out)
end


# performs the actual derivative operation on tensors on which it is defined
# if there are rules, then the derivative will be performed, else expressions will be left as eg ∇B
function perform_derivative_action(D::DerivativeSuperType)
    operand = D.operand.summands[1][2].terms[1] #get the operand
    #println(operand)

    for ky in keys(derivative_actions)
        A = same_operator(D, ky)
        B = A[1]
        o = A[2]
        if B
            action_dict = derivative_actions[ky] # dict of tensors we've defined deriv on
            for k in keys(action_dict)#looks through all tensors
                #println(k)
                #println(o)
                op=same_tensor_mod_pattern(operand, k, o)
                if op[1]
                    #do derivative:
                    replacement = action_dict[k] # need to go through and replace this with correct pattern
                    inverted_dict = Dict(v => k for (k, v) in o)
                    inverted_dict2 = Dict(v => k for (k, v) in op[2])
                    #println("replacement is, ", replacement)
                    pat_rep = index_to_pattern(replacement, inverted_dict)
                    pat_rep=index_to_pattern(pat_rep, inverted_dict2)
                    #println("with correct indices is, ", pat_rep)
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
function expand_derivative(TEX::Union{TensorExpression,<:DerivativeSuperType}, D::DerivativeSuperType, DExpansion::TensorExpression)
    #think about: hellish dummie indice relabeling//////

    new_TEX = zero #construct new expression term by term
    TEX=1*TEX
    avail = Dict{IndexSet,OrderedSet{Index}}() #pairs of index sets, with the unused indices
    dummymap = Dict{String,Index}()#dict of what to replace pattern indices with
    for summand in TEX.summands
        new_summand = summand[1]
        for tensor in summand[2].terms #look at each tensor in the tensor term
            if typeof(tensor) == typeof(D) && tensor.name == D.name #operator found
                operand = tensor.operand
                # println(TensorExpression(D))
                # println(TensorExpression(tensor))
                mapping = construct_index_rep_map(D,tensor) #index map so we use correct indices in expansion
                # println("mapping ", mapping)
                replace = copy(DExpansion) #new copy of expansion to replace in our exp

                #figure out what indices we're allowed to use in the replacement
                sets = get_index_sets(tensor) #Set{IndexSet}()
                
                for set in sets
                    a = OrderedSet{Index}(get_unused_indices(TEX, set))
                    push!(avail, set => a) #ones we can replace with
                end
                #construct the replacement
                new_rep_summand = zero
                for repsummand in replace.summands
                    new_rep_term = repsummand[1]
                    

                    ###############fix index logic for multiple tensors in the expressions
                    ## currently, K^{*β}_{α*}∂^{i*}_{*β}[[0*]]   => goes to K^{*β}_{α*}∂^{i*}_{*β}[[F^{ik**}_{**jβ}]], where F^{ik**}_{**jβ}


                    for reptensor in repsummand[2].terms
                        ##println(reptensor)
                        if typeof(reptensor) <: OperatorSuperType
                            #if we found an operator,
                            #println("operators")
                            reptensor.operand = operand#set to correct operand

                            #dums=get_dummies(reptensor)#gets the placeholder dummy indices from the pattern
                            for pair in eachindex(reptensor.operator_indices) # replace dummies with correct etc
                                #println(reptensor.operator_indices[pair][1])

                                if haskey(mapping, reptensor.operator_indices[pair][1].name) #if the index is fixed already from before
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
                                    if !haskey(avail,reptensor.operator_indices[pair][1].index_set)
                                            #if we're using a dummy index type thats not in the orignal expression
                                            new_index=pop!(get_unused_indices(TEX,reptensor.operator_indices[pair][1].index_set))
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
                                    if !haskey(avail,reptensor.indices[pair][1].index_set)
                                        #if we're using a dummy index type thats not in the orignal expression
                                       ##println("getta new index for index sets")
                                        new_index=pop!(get_unused_indices(TEX,reptensor.indices[pair][1].index_set))
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

                        #println("prev new rep: ",new_rep_term) 
                        #println(typeof(new_rep_term))   
                        #println("multiplying by: ",reptensor)
                        #println(typeof(reptensor))
                        new_rep_term *= reptensor
                        #println("new_rep_term is ",new_rep_term)

                    end
                    #  update_ext_indices(new_rep_term)
                    #println("old summand is ",new_rep_summand)
                    #println("adding ",new_rep_term)
                    update_ext_indices(new_rep_term)
                    update_ext_indices(new_rep_summand)
                   ##println(new_rep_term)
                   ##println(get_free_indices(new_rep_term))
                   ##println(new_rep_summand)
                   ##println(get_free_indices(new_rep_summand))
                    new_rep_summand += new_rep_term
                end
                #println("timesing here")
                new_summand *= new_rep_summand
                #println(new_summand)
            else#normal tensor
                #println("timesing here22222")
                new_summand *= tensor
                #println(new_summand)
            end
        end
        #        update_ext_indices(new_summand)
        #println("adding here22222")
        new_TEX += summand[1] * new_summand
        #println(new_TEX)
    end
   ##println(new_TEX)
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
                new_summand = zero+ grd * summand[1] * new_prod #adding the new term to the leibniz expanded version of the summand
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

            grad = (grad + grading(summand[2].terms[tensorind])*op_grading(O)) % 2
           ##println("where..........")
        end

    end
    return new_summand
end


function update_ext_indices(O::OperatorSuperType)
    if typeof(O.operand)<:TensorSuperType ||typeof(O.operand)<:TensorTerm 
        O.indices=vcat(O.operator_indices, O.operand.indices)
    elseif typeof(O.operand)<:TensorExpression
        O.indices=vcat(O.operator_indices, O.operand.summands[1][2].indices)
    end
end




function *(TEx::TensorExpression, T::OperatorSuperType)
    # println("doing this sort of multiplcation5")
    return update_indices(TEx * TensorExpression(T))
end

function *(T::OperatorSuperType, TEx::TensorExpression)
    # println("doing this sort of multiplcation6")
    return update_indices(TensorExpression(T) * TEx)
end

function copy(T::DerivativeSuperType)
    return typeof(T)(T.name, copy(T.operand), copy(T.coord), copy(T.operator_indices))
end