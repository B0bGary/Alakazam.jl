include("operator.jl")

##############################################################################################################################
##############################################################################################################################

#derivative rules


derivative_actions=Dict{D,Dict{TensorSuperType, TensorExpression}}() where D<:DerivativeSuperType

function is_pattern_equivalent(D_pat<:DerivativeSuperType,exp_pat::TensorExpression,D_stored<:DerivativeSuperType,exp_stored::TensorExpression)

end


function same_operator(D1<:DerivativeSuperType,D2<:DerivativeSuperType)
    ind1=D1.operator_indices
    ind2=D2.operator_indices
    doubles=Dict{String,String}()#####
    for ind in eachindex(ind1)
        if ind1[ind][1].index_set==ind2[ind][1].index_set && ind1[ind][2].position==ind2[ind][2].position
            if haskey(doubles,ind1[ind][1].name)
                val=doubles[ind1[ind][1].name]
                if val!=ind2[ind][1].name
                    #if the patterns don't match ie repeated indices etc
                    return (false =>doubles)
                end
            else
                doubles[ind1[ind][1].name]=ind2[ind][1].name
            end
        else
            #mistmatch
            return (false =>doubles)
        end

    end
    return (true =>doubles)
end


function add_derivative_action(D<:DerivativeSuperType,exp::TensorSuperType,res::TensorExpression)
    key=nothing
    in_dict=false
    pattern_ref={}
    for ky in keys(derivative_actions)
        (B=>o)=same_operator(D,ky)
        if B==true
            key=ky
            in_dict=true
            pattern_reg=o
            break
        end
    end
    
    #now we have index pattern if we found a match
    if in_dict
        #check and see if our pattern matches the expression

    else


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

end










##############################################################################################################################
##############################################################################################################################

#apply the derivative and cancel terms that are annihilated

function apply_derivative(D::DerivativeSuperType)
    D=product_rule(D)
    new_summands=D.summands
    for summandind in eachindex(D.summands) # for each summand (coefficient, tensorterm)
        summand=new_summands[summandind]
        for tensor in summand[2].terms # for each individual tensor in tensorterm
            if typeof(tensor)<:DerivativeSuperType # if the tensor is an operator
                operands=tensor.operand #get the operand (should only be one thanks to prod rule)
                tens=operands.summands[1][2].terms[1] #get the only item in the operand
                if isempty(intersect(Set(tens.function_of),Set(tensor.coord))) #see if tensor is a function of the derivative coords
                    # #println(summand[1])
                    new_summands[summandind]= 0 => summand[2] #annihilate the term
                end
            end
        end
    end
    
        D.summands=new_summands
    return merge_sum(D)
end

# performs the actual derivative operation on tensors on which it is defined
# if there are rules, then the derivative will be performed, else expressions will be left as eg ∇B
function perform_derivative_action(D::DerivativeSuperType)


end


"""
expand_derivative(TEX,D,DExpansion)

Replaces a symbolic derivative with its full form, eg substituting a covariant derivative D with all of its connection terms.
'TEX' is the original expression. 'D' is the derivative to replace, and 'DExpansion' is a tensor expression that will be a sum of tensors/operators.
For each operator in 'DExpansion', the operand will be replaced by the matching operand in 'TEX'.
"""
function expand_derivative(TEX::TensorExpression,D::DerivativeSuperType,DExpansion::TensorExpression)
    #think about: hellish dummie indice relabeling//////

    new_TEX=zero #construct new expression term by term
    for summand in TEX.summands
        new_summand=1
        for tensor in summand[2].terms #look at each tensor in the tensor term
            if typeof(tensor)==typeof(D) && tensor.name==D.name #operator found
                operand=tensor.operand
                mapping=construct_index_rep_map(tensor,D) #index map so we use correct indices in expansion
                replace=copy(DExpansion) #new copy of expansion to replace in our exp

                 #figure out what indices we're allowed to use in the replacement
                 sets=get_index_sets(tensor) #Set{IndexSet}()
                 avail=Dict{IndexSet, OrderedSet{Index}}() #pairs of index sets, with the unused indices
                 for set in sets
                     a=OrderedSet{Index}(get_unused_indices(TEX, set))
                     push!(avail, set=> a) #ones we can replace with
                 end
                #construct the replacement
                 new_rep_summand=zero
                for repsummand in replace.summands
                    new_rep_term=1
                    dummymap=Dict{String,Index}()#dict of what to replace pattern indices with

                    ###############fix index logic for multiple tensors in the expressions
                    ## currently, K^{*β}_{α*}∂^{i*}_{*β}[[0*]]   => goes to K^{*β}_{α*}∂^{i*}_{*β}[[F^{ik**}_{**jβ}]], where F^{ik**}_{**jβ}


                    for reptensor in repsummand[2].terms
                        #println(reptensor)
                        if typeof(reptensor)<:OperatorSuperType
                            #if we found an operator,
                            #println("operators")
                            reptensor.operand=operand#set to correct operand
                        
                            #dums=get_dummies(reptensor)#gets the placeholder dummy indices from the pattern
                            for pair in eachindex(reptensor.operator_indices) # replace dummies with correct etc
                                #println(reptensor.operator_indices[pair][1])

                                if haskey(mapping,reptensor.operator_indices[pair][1].name) #if the index is fixed already from before
                                    # #println(reptensor.operator_indices[pair][1])
                                    # #println(mapping[reptensor.operator_indices[pair][1].name])
                                    # #println(typeof(reptensor.operator_indices[pair]))
                                    # #println(typeof(mapping[reptensor.operator_indices[pair][1].name]))
                                    reptensor.operator_indices[pair]=mapping[reptensor.operator_indices[pair][1].name]=>reptensor.operator_indices[pair][2]#should replace with the right thing
                                    #println("maps to ",mapping[reptensor.operator_indices[pair][1].name])
                                elseif haskey(dummymap,reptensor.operator_indices[pair][1].name)
                                    #println(haskey(dummymap,reptensor.operator_indices[pair][1].name))
                                    #println(dummymap[reptensor.operator_indices[pair][1].name])
                                    old=reptensor.operator_indices[pair][1].name
                                    # #println(, " maps to ",dummymap[reptensor.operator_indices[pair][1].name])
                                    # #println(dummymap)
                                    reptensor.operator_indices[pair]=dummymap[reptensor.operator_indices[pair][1].name]=>reptensor.operator_indices[pair][2]
                                    # #println(dummymap)
                                    # #println(dummymap[reptensor.operator_indices[pair][1].name])
                                    #println(old," maps to ",reptensor.operator_indices[pair][1].name)
                                else #if just a placeholder index in the pattern, take from the list
                                    new_index=pop!(avail[reptensor.operator_indices[pair][1].index_set]) # get an unused index
                                    dummymap[reptensor.operator_indices[pair][1].name]=new_index #makes sure summation dummies do the right thing
                                    reptensor.operator_indices[pair]=new_index=>reptensor.operator_indices[pair][2]
                                    #println("maps to ",new_index)
                                end
                                # for repindex in reptensor.
                            end
                        else
                            #println("nonoperator")
                            for pair in eachindex(reptensor.indices) # replace dummies with correct etc
                                if haskey(mapping,reptensor.indices[pair][1].name) #if the index is fixed already from before
                                    # #println(reptensor.operator_indices[pair][1])
                                    # #println(mapping[reptensor.operator_indices[pair][1].name])
                                    # #println(typeof(reptensor.operator_indices[pair]))
                                    # #println(typeof(mapping[reptensor.operator_indices[pair][1].name]))
                                    reptensor.indices[pair]=mapping[reptensor.indices[pair][1].name]=>reptensor.indices[pair][2]#should replace with the right thing
                                elseif haskey(dummymap,reptensor.indices[pair][1].name)
                                    reptensor.indices[pair]=dummymap[reptensor.indices[pair][1].name]=>reptensor.indices[pair][2]

                                else #if just a placeholder index in the pattern, take from the list
                                    new_index=pop!(avail[reptensor.indices[pair][1].index_set]) # get an unused index
                                    dummymap[reptensor.indices[pair][1].name]=new_index #makes sure summation dummies do the right thing
                                    reptensor.indices[pair]=new_index=>reptensor.indices[pair][2]
                                    
                                end
                            end
                        end
                            
                    #println("prev new rep: ",new_rep_term) 
                    #println(typeof(new_rep_term))   
                    #println("multiplying by: ",reptensor)
                    #println(typeof(reptensor))
                    new_rep_term*=reptensor
                    #println("new_rep_term is ",new_rep_term)
                       
                    end
                  #  update_ext_indices(new_rep_term)
                    #println("old summand is ",new_rep_summand)
                    #println("adding ",new_rep_term)
                    update_ext_indices(new_rep_term)
                    update_ext_indices(new_rep_summand)
                    new_rep_summand+=new_rep_term
                end
                #println("timesing here")
                new_summand*=new_rep_summand
                #println(new_summand)
            else#normal tensor
                #println("timesing here22222")
                new_summand*=tensor
                #println(new_summand)
            end
        end
#        update_ext_indices(new_summand)
        #println("adding here22222")
        new_TEX+=summand[1]*new_summand
        #println(new_TEX)
    end
    return new_TEX
end

# graded leibniz

function product_rule(TEX::TensorExpression)
    new_summand=zero
    for summand in TEX.summands
        nterm=1
        for term in summand[2].terms
            # println("doing this1")
            if typeof(term)<:LeibnizOperatorSuperType
                # println("doing this")
                nterm*=product_rule(term)
            else
                nterm*=term
            end
        end
        new_summand+=summand[1]*nterm
    end
    return new_summand
end


function product_rule(O::LeibnizOperatorSuperType)
    new_summand=zero
    # #println("operator is")
    # #println(O)
    for summandind= eachindex(O.operand.summands) #for each sum term in the expression
        summand=O.operand.summands[summandind] #get the summand (tensor term)
        grad=0 #initialise grading
        for tensorind in eachindex(summand[2].terms)# for each tensor in the summand (tensorterm)
            grd=(-1)^(grad*grading(O))
            operator_on_tensor=typeof(O)(O.name,summand[2].terms[tensorind],O.coord, O.operator_indices) #operator on nth term
            # #println("this is the op on a tens")
            # #println(operator_on_tensor)
            # #println("help")
            new_prod=1
            for indsss in eachindex(summand[2].terms)
                if indsss!= tensorind
               new_prod=new_prod*summand[2].terms[indsss]
                else
                    new_prod=new_prod*operator_on_tensor
                end
            end
            # new_prod=[t  where t=summand[2].terms]
            # #println(new_prod)
            # new_prod.terms[tensorind]=operator_on_tensor
            if new_summand==Nothing
                new_summand=grd*summand[1]*new_prod #adding the new term to the leibniz expanded version of the summand
            else
                new_summand=new_summand+grd*summand[1]*new_prod #adding the new term to the leibniz expanded version of the summand
            end
            
            grad=(grad+grading(summand[2].terms[tensorind]))%2
        end
        
    end
    return new_summand
end


function *(TEx::TensorExpression, T::OperatorSuperType)
    # println("doing this sort of multiplcation5")
    return update_indices(TEx*TensorExpression(T))
end

function *( T::OperatorSuperType,TEx::TensorExpression)
    # println("doing this sort of multiplcation6")
    return update_indices(TensorExpression(T)*TEx)
end

function copy(T::DerivativeSuperType)
    return typeof(T)(T.name,copy(T.operand),copy(T.coord),copy(T.operator_indices))
end