include("operator.jl")

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
                    # println(summand[1])
                    new_summands[summandind]= 0 => summand[2] #annihilate the term
                end
            end
        end
    end
    
        D.summands=new_summands
    return merge_sum(D)
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
                        println(reptensor)
                        if typeof(reptensor)<:OperatorSuperType
                            #if we found an operator,
                            println("operators")
                            reptensor.operand=operand#set to correct operand
                        
                            #dums=get_dummies(reptensor)#gets the placeholder dummy indices from the pattern
                            for pair in eachindex(reptensor.operator_indices) # replace dummies with correct etc
                                if mapping[reptensor.operator_indices[pair][1].name]!==nothing #if the index is fixed already from before
                                    # println(reptensor.operator_indices[pair][1])
                                    # println(mapping[reptensor.operator_indices[pair][1].name])
                                    # println(typeof(reptensor.operator_indices[pair]))
                                    # println(typeof(mapping[reptensor.operator_indices[pair][1].name]))
                                    reptensor.operator_indices[pair]=mapping[reptensor.operator_indices[pair][1].name]=>reptensor.operator_indices[pair][2]#should replace with the right thing
                                elseif dummymap[reptensor.operator_indices[pair][1].name]!==nothing
                                    reptensor.operator_indices[pair]=dummymap[reptensor.operator_indices[pair][1].name]=>reptensor.operator_indices[pair][2]

                                else #if just a placeholder index in the pattern, take from the list
                                    new_index=pop!(avail[reptensor.operator_indices[pair][1].index_set]) # get an unused index
                                    dummymap[reptensor.operator_indices[pair][1].name]=new_index #makes sure summation dummies do the right thing
                                    reptensor.operator_indices[pair]=new_index=>reptensor.operator_indices[pair][2]
                                    
                                end
                                # for repindex in reptensor.
                            end
                        else
                            println("nonoperator")
                            for pair in eachindex(reptensor.indices) # replace dummies with correct etc
                                if mapping[reptensor.indices[pair][1].name]!==nothing #if the index is fixed already from before
                                    # println(reptensor.operator_indices[pair][1])
                                    # println(mapping[reptensor.operator_indices[pair][1].name])
                                    # println(typeof(reptensor.operator_indices[pair]))
                                    # println(typeof(mapping[reptensor.operator_indices[pair][1].name]))
                                    reptensor.indices[pair]=mapping[reptensor.indices[pair][1].name]=>reptensor.indices[pair][2]#should replace with the right thing
                                elseif dummymap[reptensor.indices[pair][1].name]!==nothing
                                    reptensor.indices[pair]=dummymap[reptensor.indices[pair][1].name]=>reptensor.indices[pair][2]

                                else #if just a placeholder index in the pattern, take from the list
                                    new_index=pop!(avail[reptensor.indices[pair][1].index_set]) # get an unused index
                                    dummymap[reptensor.indices[pair][1].name]=new_index #makes sure summation dummies do the right thing
                                    reptensor.indices[pair]=new_index=>reptensor.indices[pair][2]
                                    
                                end
                            end
                        end
                            
                    println("prev new rep: ",new_rep_term) 
                    println(typeof(new_rep_term))   
                    println("multiplying by: ",reptensor)
                    println(typeof(reptensor))
                    new_rep_term*=reptensor
                    println("new_rep_term is ",new_rep_term)
                       
                    end
                  #  update_ext_indices(new_rep_term)
                    println("old summand is ",new_rep_summand)
                    println("adding ",new_rep_term)
                    new_rep_summand+=new_rep_term
                end
                println("timesing here")
                new_summand*=new_rep_summand
                println(new_summand)
            else#normal tensor
                println("timesing here22222")
                new_summand*=tensor
                println(new_summand)
            end
        end
#        update_ext_indices(new_summand)
        println("adding here22222")
        new_TEX+=summand[1]*new_summand
        println(new_TEX)
    end
    return new_TEX
end

# graded leibniz

function product_rule(O::DerivativeSuperType)
    new_summand=Nothing
    # println("operator is")
    # println(O)
    for summandind= eachindex(O.operand.summands) #for each sum term in the expression
        summand=O.operand.summands[summandind] #get the summand (tensor term)
        grad=0 #initialise grading
        for tensorind in eachindex(summand[2].terms)# for each tensor in the summand (tensorterm)
            grd=(-1)^(grad*grading(O))
            operator_on_tensor=typeof(O)(O.name,summand[2].terms[tensorind],O.coord, O.operator_indices) #operator on nth term
            # println("this is the op on a tens")
            # println(operator_on_tensor)
            # println("help")
            new_prod=1
            for indsss in eachindex(summand[2].terms)
                if indsss!= tensorind
               new_prod=new_prod*summand[2].terms[indsss]
                else
                    new_prod=new_prod*operator_on_tensor
                end
            end
            # new_prod=[t  where t=summand[2].terms]
            # println(new_prod)
            # new_prod.terms[tensorind]=operator_on_tensor
            if new_summand==Nothing
                new_summand=grd*summand[1]*new_prod #adding the new term to the leibniz expanded version of the summand
            else
                new_summand=new_summand+grd*summand[1]*new_prod #adding the new term to the leibniz expanded version of the summand
            end
            
            grad=(grad+grading(summand[2].terms[tensorind]))%2
        end
        
    end
    new_summand
end
