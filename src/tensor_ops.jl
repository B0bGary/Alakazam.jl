using Latexify
set_default(convert_unicode=true)

##############################################################################################################################
##############################################################################################################################
### metric ops
##############################################################################################################################
###################################################################################################################################

""""Pop" metrics off all tensors, putting indices in default position and returning the popped metric as a prefactor. This is used to put tensors in a canonical form for simplification and comparison.
"""
function pop_metric(TEx::TensorExpression)
    new_exp = zero

    for summand in TEx.summands
        # println("ORIG EEPXRESSION is: ",summand)

        q = summand[1] * pop_metric(summand[2])
        # println(new_exp)
        # println(q)
        new_exp += q
        # println("OUTP?UT TERM term is: ",new_exp)
    end
    return new_exp
end

""""Pop" metrics off all tensors, putting indices in default position and returning the popped metric as a prefactor. This is used to put tensors in a canonical form for simplification and comparison.
"""
function pop_metric(terms::TensorTerm)
    new_term = 1
    # println("popping",terms)
    for term in terms.terms
        ###println("popping metric off term ",term)
        if new_term != 1
            # println("doing this")
            # println(term)
            # # println(term.indices)
            # q = pop_metric(term, terms, new_term)
            # println("for the love of god")
            # println(new_term)
            # println(new_term.indices)
            # println(q)
            # if typeof(q) == TensorExpression
            #     println(q.summands[end][2].indices)
            # else
            #     println(q.indices)
            # end

            new_term *= pop_metric(term, terms, new_term)
        else
            # println("doing this one")
            new_term *= pop_metric(term, terms)
        end
        ###println(new_term)
    end
    return new_term
end

""""Pop" metrics off a tensors, putting indices in default position and returning the popped metric as a prefactor. This is used to put tensors in a canonical form for simplification and comparison.
"""
function pop_metric(T::TensorSuperType, orig_exp::Union{TensorExpression,TensorTerm,Nothing}=nothing, running_term::Union{TensorExpression,TensorTerm,Nothing}=nothing)# makes all indices on the object go to the default position according to the index set
    # println("doing pop of", T)
    if typeof(T) <: AntisymMetricSuperType || typeof(T) <: SymMetricSuperType || typeof(T) <: OperatorSuperType || typeof(T) <: KroneckerSuperType
        return T
    end
    T = copy(T)
    tensor = 1
    new_index_pattern = T.indices
    used_indices = Set()
    # println(T)
    for pos in eachindex(T.indices)
        index = T.indices[pos]
        if !in_default_pos(index)
            ###println("resetting to default positoin")
            if orig_exp === nothing
                spare = get_unused_indices(T, index[1].index_set)
            else
                spare = get_unused_indices(orig_exp, index[1].index_set)
                if running_term !== nothing
                    spare = intersect(spare, get_unused_indices(running_term, index[1].index_set))
                end
                spare = setdiff(spare, used_indices)

                ###println("using spares from original exp")
                ###println(spare)
            end

            ###println(spare)
            if isempty(spare)
                println(index[1].index_set)
                throw(OutOfIndices())
            end
            use_this = pop!(spare)
            dummy = use_this => index[1].index_set.default_position == upper ? lower : upper
            push!(used_indices, use_this)
            if dummy[1].index_set.name == "spacetime"
                # println("is space index")
                metric = Metric("η", [index, dummy])
            else
                metric = EpsilonTensor([index, dummy])
            end

            tensor = tensor * metric
            ###println(new_index_pattern[pos])
            ###println(dummy => index[1].index_set.default_position)
            new_index_pattern[pos] = dummy[1] => index[1].index_set.default_position
        end
    end
    T.indices = new_index_pattern
    # println("popping metric")
    # println(tensor)
    # println("from tensor")
    # println(T)
    return tensor * T

end

""""
Contract all metrics of a particular type in an expression, whereever possible.
"""
function contract_metrics(T::TensorTerm, typ::Type=Metric)
    contract_metrics(1 * T, typ)
end

""""
Contract all metrics of a particular type in an expression, whereever possible.
"""
function contract_metrics(T::TensorExpression, typ::Type=Metric)
    ###println("running function")
    T = copy(T)
    for k in eachindex(T.summands)
        ###println(k)
        ###println(T.summands)
        summand = T.summands[k]
        # println(summand)
        metrics = filter(y -> typeof(y) <: typ, summand[2].terms)
        # if typ<:EpsilonTensorSuperType
        #     println(metrics)
        # end
        non_metrics = filter(y -> !(typeof(y) <: typ), summand[2].terms)
        uncontr_mets = Vector{TensorSuperType}()
        ###println(metrics) 
        sign = 1 #for antisym metrics (spinor)
        factor = 1 #(for contracted deltas etc)
        extra_deltas = Vector{TensorSuperType}()
        for metric in metrics #for each in the found metrics
            flag = false # terminate when metric is contracted
            index1 = metric.indices[1][1]
            index2 = metric.indices[2][1]
            # println(metric)
            if typeof(metric) <: KroneckerSuperType
                if index1.index_set.name == "lightcone"
                    if index1 != index2
                        factor *= 0
                    end
                    flag = true
                    break
                else
                    if index1 == index2 # contracts same index kroneckers
                        factor *= index1.index_set.range
                        flag = true #something here breaks!!!!! 
                        #       break #idk why but having the break screws things up, make it drop terms in front of the delta?? weird behaviour
                    end
                end
            end
            if typeof(metric) <: EpsilonTensorSuperType
                if index1 == index2 # contracts same index kroneckers
                    factor *= 0
                    flag = true
                    break
                end
            end

            for tens in non_metrics #try to contract with each other tensor
                ###println(tens)
                if typeof(metric) <: EpsilonTensorSuperType && typeof(tens) <: KroneckerSuperType
                    # println("cant contract kronecker")
                    continue
                end
                # println(tens,"not a kronecker")
                if flag # terminat if contracted
                    break
                end

                ###########################################################################
                if typeof(tens) <: OperatorSuperType
                    # if typeof(metric) <: KroneckerSuperType
                    #     println("kron metric,,", tens, metric)
                    # else
                    #     println("non-kron metric,,", tens, metric)
                    # end
                    operand = tens.operand
                    # println(operand)
                    # println(typeof(operand))
                    if typeof(operand) <: TensorSuperType || typeof(operand) <: TensorTerm || (typeof(operand) <: TensorExpression && length(operand.summands) == 1 && length(operand.summands[1][2].terms) == 1)
                        if typeof(operand) <: TensorExpression || typeof(operand) <: OperatorSuperType
                            op_chain = Any[tens]
                            # note: operand is always a tensor exp
                            # println(typeof(operand))
                            # while typeof(operand) <:OperatorSuperType 
                            while length(operand.summands) == 1 && length(operand.summands[1][2].terms) == 1 && typeof(operand.summands[1][2].terms[1]) <: OperatorSuperType
                                #if operand has one summand, and if the object is an operator
                                push!(op_chain, operand.summands[1][2].terms[1])
                                # println("pushed:", operand.summands[1][2].terms[1])
                                operand = operand.summands[1][2].terms[1].operand
                                # println("next operand:", operand)
                                # if typeof(operand)<:TensorExpression
                                #     if typeof(operand.summands[1][2])
                                # end 
                            end
                            # println(typeof(operand))
                            # println(operand)
                            push!(op_chain, operand.summands[1][2].terms[1])
                            # println(op_chain)
                            #operand is expression
                            # println("got hersssse")
                            for indpairnum in eachindex(operand.summands[1][2].terms[1].indices)
                                if flag
                                    break
                                end
                                # println("got hersssssse2222")
                                if operand.summands[1][2].terms[1].indices[indpairnum][1] == index2 #tensor has an index matching the 2nd index of the metric
                                    #some complicated logic to make sure indices go to the right position

                                    # this should let you contract Kronckers and metrics in a consistent way
                                    met_pos1 = metric.indices[1][2].position # identifies when the metric is metric or kronecker etc
                                    met_pos2 = metric.indices[2][2].position
                                    ###println(met_pos1)
                                    ###println(met_pos2)
                                    sumind = (met_pos1 + met_pos2) % 2 # if = 0, both upper/lower, 1 if 1 upper 1 lower
                                    if sumind == 0

                                        if !(typeof(metric) <: SymMetricSuperType)

                                            if (operand.summands[1][2].terms[1].indices[indpairnum][2].position + met_pos1) % 2 == 0 # 2nd metric and tensor pos same, raise and lower with sign change
                                                sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use 
                                                ###println("flipping")
                                                ###println(operand.summands[1][2].terms[1].indices[indpairnum])
                                                if operand.summands[1][2].terms[1].indices[indpairnum][2] == upper
                                                    ###println("sending upper to lower")

                                                    operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[1][1] => lower
                                                    ###println(operand.summands[1][2].terms[1].indices)
                                                else
                                                    ###println("sending lower to upper")
                                                    operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[1][1] => upper
                                                    ###println(operand.summands[1][2].terms[1].indices)
                                                end
                                            else
                                                operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                            end
                                        else

                                            operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                        end
                                    else
                                        ###println("doing this one")
                                        ###println(operand.summands[1][2].terms[1].indices[indpairnum][2])
                                        if (operand.summands[1][2].terms[1].indices[indpairnum][2].position + met_pos2) % 2 == 0 # if tensor index position is same as 2nd metric pos 
                                            operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[1][1] => operand.summands[1][2].terms[1].indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)
                                        else
                                            operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[1][1] => operand.summands[1][2].terms[1].indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                        end
                                    end
                                    #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                                    flag = true #this term has been contracted, break
                                # println("contracted successfully on 2nd index of metric")
                                elseif operand.summands[1][2].terms[1].indices[indpairnum][1] == index1
                                    #some complicated logic to make sure indices go to the right position
                                    #if contracting with antisym, pick up a minus
                                    ###println("symmetric?")
                                    ###println(typeof(metric)<:SymMetricSuperType)
                                    ###println(sign)
                                    met_pos1 = metric.indices[1][2].position
                                    met_pos2 = metric.indices[2][2].position
                                    ###println(met_pos1)
                                    ###println(met_pos2)
                                    sumind = (met_pos1 + met_pos2) % 2
                                    if sumind == 0
                                        ###println(" yo")
                                        if !(typeof(metric) <: SymMetricSuperType)
                                            if (operand.summands[1][2].terms[1].indices[indpairnum][2].position + met_pos1) % 2 == 0 # 1st metric and tensor pos same, raise and lower without sign
                                                if operand.summands[1][2].terms[1].indices[indpairnum][2] == upper
                                                    operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[2][1] => lower
                                                else
                                                    operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[2][1] => upper
                                                end
                                            else
                                                sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use sign
                                                operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower index
                                            end
                                        else

                                            operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower
                                        end
                                    else # shoulder never get here for anti sym
                                        ###println("doing this one")
                                        ###println(operand.summands[1][2].terms[1].indices[indpairnum][2])
                                        if (operand.summands[1][2].terms[1].indices[indpairnum][2].position + met_pos1) % 2 == 0
                                            ###println(" here")
                                            #for anti sym, we need to keep sign but change pos 
                                            #sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                            # if !(typeof(metric)<:SymMetricSuperType)
                                            #    ##println("we antisym")
                                            #     if operand.summands[1][2].terms[1].indices[indpairnum][2].position ==upper
                                            #         operand.summands[1][2].terms[1].indices[indpairnum]=metric.indices[2][1] => lower
                                            #     else
                                            #         operand.summands[1][2].terms[1].indices[indpairnum]=metric.indices[2][1] => upper
                                            #     end
                                            # else 
                                            #       println(operand.summands[1][2].terms[1].indices[indpairnum])
                                            operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[2][1] => operand.summands[1][2].terms[1].indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)
                                        #      println(operand.summands[1][2].terms[1].indices[indpairnum])

                                        else #opposite pos, grading etc
                                            # sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                            ###println("got here")
                                            # if !(typeof(metric)<:SymMetricSuperType) #for anti symm, change positions
                                            #    ##println("we antisym too")
                                            #     if operand.summands[1][2].terms[1].indices[indpairnum][2].position ==upper
                                            #         operand.summands[1][2].terms[1].indices[indpairnum]=metric.indices[2][1] => lower
                                            #     else
                                            #         operand.summands[1][2].terms[1].indices[indpairnum]=metric.indices[2][1] => upper
                                            #     end
                                            # else
                                            #        println(operand.summands[1][2].terms[1].indices[indpairnum])
                                            operand.summands[1][2].terms[1].indices[indpairnum] = metric.indices[2][1] => operand.summands[1][2].terms[1].indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                            # end
                                            # println(operand.summands[1][2].terms[1].indices[indpairnum])
                                            # println(operand.summands[1][2].terms[1].indices)
                                            # println(operand)
                                        end
                                    end
                                    #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                                    flag = true #this term has been contracted, break
                                    # println("contracted successfully on 1st index of metric")
                                end
                                # println("new operand is",operand.summands[1])
                            end
                            for indpairnum in eachindex(tens.operator_indices)

                                if flag
                                    break
                                end
                                if tens.operator_indices[indpairnum][1] == index2 #tensor has an index matching the 2nd index of the metric
                                    #some complicated logic to make sure indices go to the right position

                                    # this should let you contract Kronckers and metrics in a consistent way
                                    met_pos1 = metric.indices[1][2].position # identifies when the metric is metric or kronecker etc
                                    met_pos2 = metric.indices[2][2].position
                                    ###println(met_pos1)
                                    ###println(met_pos2)
                                    sumind = (met_pos1 + met_pos2) % 2 # if = 0, both upper/lower, 1 if 1 upper 1 lower
                                    if sumind == 0

                                        if !(typeof(metric) <: SymMetricSuperType)

                                            if (tens.operator_indices[indpairnum][2].position + met_pos1) % 2 == 0 # 2nd metric and tensor pos same, raise and lower with sign change
                                                sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use 
                                                ###println("flipping")
                                                ###println(tens.operator_indices[indpairnum])
                                                if tens.operator_indices[indpairnum][2] == upper
                                                    ###println("sending upper to lower")

                                                    tens.operator_indices[indpairnum] = metric.indices[1][1] => lower
                                                    ###println(tens.operator_indices)
                                                else
                                                    ###println("sending lower to upper")
                                                    tens.operator_indices[indpairnum] = metric.indices[1][1] => upper
                                                    ###println(tens.operator_indices)
                                                end
                                            else
                                                tens.operator_indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                            end
                                        else

                                            tens.operator_indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                        end
                                    else
                                        ###println("doing this one")
                                        ###println(tens.operator_indices[indpairnum][2])
                                        if (tens.operator_indices[indpairnum][2].position + met_pos2) % 2 == 0 # if tensor index position is same as 2nd metric pos 
                                            tens.operator_indices[indpairnum] = metric.indices[1][1] => tens.operator_indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)
                                        else
                                            tens.operator_indices[indpairnum] = metric.indices[1][1] => tens.operator_indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                        end
                                    end
                                    #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                                    flag = true #this term has been contracted, break
                                elseif tens.operator_indices[indpairnum][1] == index1
                                    #some complicated logic to make sure indices go to the right position
                                    #if contracting with antisym, pick up a minus
                                    ###println("symmetric?")
                                    ###println(typeof(metric)<:SymMetricSuperType)
                                    ###println(sign)
                                    met_pos1 = metric.indices[1][2].position
                                    met_pos2 = metric.indices[2][2].position
                                    ###println(met_pos1)
                                    ###println(met_pos2)
                                    sumind = (met_pos1 + met_pos2) % 2
                                    if sumind == 0
                                        ###println(" yo")
                                        if !(typeof(metric) <: SymMetricSuperType)
                                            if (tens.operator_indices[indpairnum][2].position + met_pos1) % 2 == 0 # 1st metric and tensor pos same, raise and lower without sign
                                                if tens.operator_indices[indpairnum][2] == upper
                                                    tens.operator_indices[indpairnum] = metric.indices[2][1] => lower
                                                else
                                                    tens.operator_indices[indpairnum] = metric.indices[2][1] => upper
                                                end
                                            else
                                                sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use sign
                                                tens.operator_indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower index
                                            end
                                        else

                                            tens.operator_indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower
                                        end
                                    else # shoulder never get here for anti sym
                                        ###println("doing this one")
                                        ###println(tens.operator_indices[indpairnum][2])
                                        if (tens.operator_indices[indpairnum][2].position + met_pos1) % 2 == 0
                                            ###println(" here")
                                            #for anti sym, we need to keep sign but change pos 
                                            #sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                            # if !(typeof(metric)<:SymMetricSuperType)
                                            #    ##println("we antisym")
                                            #     if tens.operator_indices[indpairnum][2].position ==upper
                                            #         tens.operator_indices[indpairnum]=metric.indices[2][1] => lower
                                            #     else
                                            #         tens.operator_indices[indpairnum]=metric.indices[2][1] => upper
                                            #     end
                                            # else 
                                            tens.operator_indices[indpairnum] = metric.indices[2][1] => tens.operator_indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)


                                        else #opposite pos, grading etc
                                            # sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                            ###println("got here")
                                            # if !(typeof(metric)<:SymMetricSuperType) #for anti symm, change positions
                                            #    ##println("we antisym too")
                                            #     if tens.operator_indices[indpairnum][2].position ==upper
                                            #         tens.operator_indices[indpairnum]=metric.indices[2][1] => lower
                                            #     else
                                            #         tens.operator_indices[indpairnum]=metric.indices[2][1] => upper
                                            #     end
                                            # else
                                            tens.operator_indices[indpairnum] = metric.indices[2][1] => tens.operator_indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                            # end
                                        end
                                    end
                                    #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                                    flag = true #this term has been contracted, break
                                end
                            end

                            #operand.summands[1][2].terms[1].operand

                            # println(op_chain)
                            for oind in reverse(eachindex(op_chain))
                                # println("doing this")
                                if oind > 1
                                    # println( op_chain[oind-1],", ",op_chain[oind])
                                    op_chain[oind-1].operand_indices = op_chain[oind].indices
                                    op_chain[oind-1].indices = vcat(op_chain[oind-1].operator_indices, op_chain[oind-1].operand_indices)
                                    # println(op_chain[oind-1])
                                    # println(op_chain[oind-1].indices)
                                end
                            end
                            # tens.operand_indices = operand.summands[1][2].terms[1].indices
                            # tens.indices = vcat(tens.operator_indices, tens.operand_indices)
                            # println("WEIRDDDDDD")
                            # println( tens.operand_indices)
                            # println( operand.summands[1][2].terms[1].indices)
                            update_indices(operand)
                            # println( operand.summands[1][2].terms[1].indices)
                            # println(tens)
                            #####################################################
                            ###################################################3#
                            ########################################################
                        else
                            ########################################
                            ########################################
                            #operand is always a tensor expression!
                            #not an exp
                            # println("111111got here")
                            for indpairnum in eachindex(tens.operand.indices)
                                if flag
                                    break
                                end
                                # println("111111got here2222")
                                if tens.operand.indices[indpairnum][1] == index2 #tensor has an index matching the 2nd index of the metric
                                    #some complicated logic to make sure indices go to the right position

                                    # this should let you contract Kronckers and metrics in a consistent way
                                    met_pos1 = metric.indices[1][2].position # identifies when the metric is metric or kronecker etc
                                    met_pos2 = metric.indices[2][2].position
                                    ###println(met_pos1)
                                    ###println(met_pos2)
                                    sumind = (met_pos1 + met_pos2) % 2 # if = 0, both upper/lower, 1 if 1 upper 1 lower
                                    if sumind == 0

                                        if !(typeof(metric) <: SymMetricSuperType)

                                            if (tens.operand.indices[indpairnum][2].position + met_pos1) % 2 == 0 # 2nd metric and tensor pos same, raise and lower with sign change
                                                sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use 
                                                ###println("flipping")
                                                ###println(tens.operand.indices[indpairnum])
                                                if tens.operand.indices[indpairnum][2] == upper
                                                    ###println("sending upper to lower")

                                                    tens.operand.indices[indpairnum] = metric.indices[1][1] => lower
                                                    ###println(tens.operand.indices)
                                                else
                                                    ###println("sending lower to upper")
                                                    tens.operand.indices[indpairnum] = metric.indices[1][1] => upper
                                                    ###println(tens.operand.indices)
                                                end
                                            else
                                                tens.operand.indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                            end
                                        else

                                            tens.operand.indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                        end
                                    else
                                        ###println("doing this one")
                                        ###println(tens.operand.indices[indpairnum][2])
                                        if (tens.operand.indices[indpairnum][2].position + met_pos2) % 2 == 0 # if tensor index position is same as 2nd metric pos 
                                            tens.operand.indices[indpairnum] = metric.indices[1][1] => tens.operand.indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)
                                        else
                                            tens.operand.indices[indpairnum] = metric.indices[1][1] => tens.operand.indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                        end
                                    end
                                    #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                                    flag = true #this term has been contracted, break
                                elseif tens.operand.indices[indpairnum][1] == index1
                                    #some complicated logic to make sure indices go to the right position
                                    #if contracting with antisym, pick up a minus
                                    ###println("symmetric?")
                                    ###println(typeof(metric)<:SymMetricSuperType)
                                    ###println(sign)
                                    met_pos1 = metric.indices[1][2].position
                                    met_pos2 = metric.indices[2][2].position
                                    ###println(met_pos1)
                                    ###println(met_pos2)
                                    sumind = (met_pos1 + met_pos2) % 2
                                    if sumind == 0
                                        ###println(" yo")
                                        if !(typeof(metric) <: SymMetricSuperType)
                                            if (tens.operand.indices[indpairnum][2].position + met_pos1) % 2 == 0 # 1st metric and tensor pos same, raise and lower without sign
                                                if tens.operand.indices[indpairnum][2] == upper
                                                    tens.operand.indices[indpairnum] = metric.indices[2][1] => lower
                                                else
                                                    tens.operand.indices[indpairnum] = metric.indices[2][1] => upper
                                                end
                                            else
                                                sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use sign
                                                tens.operand.indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower index
                                            end
                                        else

                                            tens.operand.indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower
                                        end
                                    else # shoulder never get here for anti sym
                                        ###println("doing this one")
                                        ###println(tens.operand.indices[indpairnum][2])
                                        if (tens.operand.indices[indpairnum][2].position + met_pos1) % 2 == 0
                                            ###println(" here")
                                            #for anti sym, we need to keep sign but change pos 
                                            #sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                            # if !(typeof(metric)<:SymMetricSuperType)
                                            #    ##println("we antisym")
                                            #     if tens.operand.indices[indpairnum][2].position ==upper
                                            #         tens.operand.indices[indpairnum]=metric.indices[2][1] => lower
                                            #     else
                                            #         tens.operand.indices[indpairnum]=metric.indices[2][1] => upper
                                            #     end
                                            # else 
                                            tens.operand.indices[indpairnum] = metric.indices[2][1] => tens.operand.indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)


                                        else #opposite pos, grading etc
                                            # sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                            ###println("got here")
                                            # if !(typeof(metric)<:SymMetricSuperType) #for anti symm, change positions
                                            #    ##println("we antisym too")
                                            #     if tens.operand.indices[indpairnum][2].position ==upper
                                            #         tens.operand.indices[indpairnum]=metric.indices[2][1] => lower
                                            #     else
                                            #         tens.operand.indices[indpairnum]=metric.indices[2][1] => upper
                                            #     end
                                            # else
                                            tens.operand.indices[indpairnum] = metric.indices[2][1] => tens.operand.indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                            # end
                                        end
                                    end
                                    #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                                    flag = true #this term has been contracted, break
                                end
                            end

                            for indpairnum in eachindex(tens.operator_indices)
                                if flag
                                    break
                                end
                                if tens.operator_indices[indpairnum][1] == index2 #tensor has an index matching the 2nd index of the metric
                                    #some complicated logic to make sure indices go to the right position

                                    # this should let you contract Kronckers and metrics in a consistent way
                                    met_pos1 = metric.indices[1][2].position # identifies when the metric is metric or kronecker etc
                                    met_pos2 = metric.indices[2][2].position
                                    ###println(met_pos1)
                                    ###println(met_pos2)
                                    sumind = (met_pos1 + met_pos2) % 2 # if = 0, both upper/lower, 1 if 1 upper 1 lower
                                    if sumind == 0

                                        if !(typeof(metric) <: SymMetricSuperType)

                                            if (tens.operator_indices[indpairnum][2].position + met_pos1) % 2 == 0 # 2nd metric and tensor pos same, raise and lower with sign change
                                                sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use 
                                                ###println("flipping")
                                                ###println(tens.operator_indices[indpairnum])
                                                if tens.operator_indices[indpairnum][2] == upper
                                                    ###println("sending upper to lower")

                                                    tens.operator_indices[indpairnum] = metric.indices[1][1] => lower
                                                    ###println(tens.operator_indices)
                                                else
                                                    ###println("sending lower to upper")
                                                    tens.operator_indices[indpairnum] = metric.indices[1][1] => upper
                                                    ###println(tens.operator_indices)
                                                end
                                            else
                                                tens.operator_indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                            end
                                        else

                                            tens.operator_indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                        end
                                    else
                                        ###println("doing this one")
                                        ###println(tens.operator_indices[indpairnum][2])
                                        if (tens.operator_indices[indpairnum][2].position + met_pos2) % 2 == 0 # if tensor index position is same as 2nd metric pos 
                                            tens.operator_indices[indpairnum] = metric.indices[1][1] => tens.operator_indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)
                                        else
                                            tens.operator_indices[indpairnum] = metric.indices[1][1] => tens.operator_indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                        end
                                    end
                                    #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                                    flag = true #this term has been contracted, break
                                elseif tens.operator_indices[indpairnum][1] == index1
                                    #some complicated logic to make sure indices go to the right position
                                    #if contracting with antisym, pick up a minus
                                    ###println("symmetric?")
                                    ###println(typeof(metric)<:SymMetricSuperType)
                                    ###println(sign)
                                    met_pos1 = metric.indices[1][2].position
                                    met_pos2 = metric.indices[2][2].position
                                    ###println(met_pos1)
                                    ###println(met_pos2)
                                    sumind = (met_pos1 + met_pos2) % 2
                                    if sumind == 0
                                        ###println(" yo")
                                        if !(typeof(metric) <: SymMetricSuperType)
                                            if (tens.operator_indices[indpairnum][2].position + met_pos1) % 2 == 0 # 1st metric and tensor pos same, raise and lower without sign
                                                if tens.operator_indices[indpairnum][2] == upper
                                                    tens.operator_indices[indpairnum] = metric.indices[2][1] => lower
                                                else
                                                    tens.operator_indices[indpairnum] = metric.indices[2][1] => upper
                                                end
                                            else
                                                sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use sign
                                                tens.operator_indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower index
                                            end
                                        else

                                            tens.operator_indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower
                                        end
                                    else # shoulder never get here for anti sym
                                        ###println("doing this one")
                                        ###println(tens.operator_indices[indpairnum][2])
                                        if (tens.operator_indices[indpairnum][2].position + met_pos1) % 2 == 0
                                            ###println(" here")
                                            #for anti sym, we need to keep sign but change pos 
                                            #sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                            # if !(typeof(metric)<:SymMetricSuperType)
                                            #    ##println("we antisym")
                                            #     if tens.operator_indices[indpairnum][2].position ==upper
                                            #         tens.operator_indices[indpairnum]=metric.indices[2][1] => lower
                                            #     else
                                            #         tens.operator_indices[indpairnum]=metric.indices[2][1] => upper
                                            #     end
                                            # else 
                                            tens.operator_indices[indpairnum] = metric.indices[2][1] => tens.operator_indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)


                                        else #opposite pos, grading etc
                                            # sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                            ###println("got here")
                                            # if !(typeof(metric)<:SymMetricSuperType) #for anti symm, change positions
                                            #    ##println("we antisym too")
                                            #     if tens.operator_indices[indpairnum][2].position ==upper
                                            #         tens.operator_indices[indpairnum]=metric.indices[2][1] => lower
                                            #     else
                                            #         tens.operator_indices[indpairnum]=metric.indices[2][1] => upper
                                            #     end
                                            # else
                                            tens.operator_indices[indpairnum] = metric.indices[2][1] => tens.operator_indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                            # end
                                        end
                                    end
                                    #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                                    flag = true #this term has been contracted, break
                                end
                            end

                            tens.operand_indices = tens.operand.indices
                            tens.indices = vcat(tens.operator_indices, tens.operand_indices)
                        end
                    end
                else
                    ###not an operator
                    #########################################################################################
                    for indpairnum in eachindex(tens.indices)
                        if flag
                            break
                        end
                        if tens.indices[indpairnum][1] == index2 #tensor has an index matching the 2nd index of the metric
                            #some complicated logic to make sure indices go to the right position

                            # this should let you contract Kronckers and metrics in a consistent way
                            met_pos1 = metric.indices[1][2].position # identifies when the metric is metric or kronecker etc
                            met_pos2 = metric.indices[2][2].position
                            ###println(met_pos1)
                            ###println(met_pos2)
                            sumind = (met_pos1 + met_pos2) % 2 # if = 0, both upper/lower, 1 if 1 upper 1 lower
                            if sumind == 0 # metric, not a kronecker

                                if !(typeof(metric) <: SymMetricSuperType) #eps tensor

                                    if (tens.indices[indpairnum][2].position + met_pos1) % 2 == 0 # 2nd metric and tensor pos same, raise and lower with sign change
                                        sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use 
                                        ###println("flipping")
                                        ###println(tens.indices[indpairnum])
                                        if tens.indices[indpairnum][2] == upper
                                            ###println("sending upper to lower")

                                            tens.indices[indpairnum] = metric.indices[1][1] => lower
                                            ###println(tens.indices)
                                        else
                                            ###println("sending lower to upper")
                                            tens.indices[indpairnum] = metric.indices[1][1] => upper
                                            ###println(tens.indices)
                                        end
                                    else
                                        tens.indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                    end
                                else

                                    tens.indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                end
                            else #TTHIS IS A KroneckerDelta
                                ###println("doing this one")
                                # println("is a krnoecker, ",metric)
                                # println("contracting, ",tens)
                                ###println(tens.indices[indpairnum][2])

                                if (tens.indices[indpairnum][2].position + met_pos2) % 2 == 0 # if tensor index position is same as 2nd metric pos 
                                    tens.indices[indpairnum] = metric.indices[1][1] => tens.indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)
                                else
                                    tens.indices[indpairnum] = metric.indices[1][1] => tens.indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                end
                                # println("yields, ",tens)
                            end
                            #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                            flag = true #this term has been contracted, break
                        elseif tens.indices[indpairnum][1] == index1
                            #some complicated logic to make sure indices go to the right position
                            #if contracting with antisym, pick up a minus
                            ###println("symmetric?")
                            ###println(typeof(metric)<:SymMetricSuperType)
                            ###println(sign)
                            met_pos1 = metric.indices[1][2].position
                            met_pos2 = metric.indices[2][2].position
                            ###println(met_pos1)
                            ###println(met_pos2)
                            sumind = (met_pos1 + met_pos2) % 2
                            if sumind == 0
                                ###println(" yo")
                                if !(typeof(metric) <: SymMetricSuperType)
                                    if (tens.indices[indpairnum][2].position + met_pos1) % 2 == 0 # 1st metric and tensor pos same, raise and lower without sign
                                        if tens.indices[indpairnum][2] == upper
                                            tens.indices[indpairnum] = metric.indices[2][1] => lower
                                        else
                                            tens.indices[indpairnum] = metric.indices[2][1] => upper
                                        end
                                    else
                                        sign = sign * ((-1)^(!(typeof(metric) <: SymMetricSuperType))) #use sign
                                        tens.indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower index
                                    end
                                else

                                    tens.indices[indpairnum] = metric.indices[2][1] => metric.indices[2][2] # raise/lower
                                end
                            else # shoulder never get here for anti sym
                                ###println("doing this one")
                                ###println(tens.indices[indpairnum][2])
                                # println("1111is a krnoecker, ",metric)
                                # println("111contracting, ",tens)
                                if (tens.indices[indpairnum][2].position + met_pos1) % 2 == 0
                                    ###println(" here")
                                    #for anti sym, we need to keep sign but change pos 
                                    #sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                    # if !(typeof(metric)<:SymMetricSuperType)
                                    #    ##println("we antisym")
                                    #     if tens.indices[indpairnum][2].position ==upper
                                    #         tens.indices[indpairnum]=metric.indices[2][1] => lower
                                    #     else
                                    #         tens.indices[indpairnum]=metric.indices[2][1] => upper
                                    #     end
                                    # else 
                                    tens.indices[indpairnum] = metric.indices[2][1] => tens.indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)


                                else #opposite pos, grading etc
                                    # sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                    ###println("got here")
                                    # if !(typeof(metric)<:SymMetricSuperType) #for anti symm, change positions
                                    #    ##println("we antisym too")
                                    #     if tens.indices[indpairnum][2].position ==upper
                                    #         tens.indices[indpairnum]=metric.indices[2][1] => lower
                                    #     else
                                    #         tens.indices[indpairnum]=metric.indices[2][1] => upper
                                    #     end
                                    # else
                                    tens.indices[indpairnum] = metric.indices[2][1] => tens.indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                    # end
                                end
                                # println("yields, ",tens)
                            end
                            #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                            flag = true #this term has been contracted, break
                        end
                    end
                end

            end
            if !flag
                # println("didn't contract metric: ", metric)
                push!(uncontr_mets, metric)
            end
        end
        to_remove = []
        for i in eachindex(uncontr_mets)
            if i in to_remove
                continue
            end
            for j in eachindex(uncontr_mets)
                if j in to_remove
                    continue
                end
                if i < j
                    metric = uncontr_mets[i]
                    tens = uncontr_mets[j]
                    # println(metric)
                    # println(tens)
                    if typeof(metric) <: EpsilonTensorSuperType && typeof(tens) <: EpsilonTensorSuperType
                        # println(metric, " and ",tens)
                        # see if they have same indices, and then contract
                        # println(get_free_indices_nopos(metric))
                        # println(get_free_indices_nopos(tens))
                        if get_free_indices_nopos(metric) == get_free_indices_nopos(tens)
                            # println("same indes")
                            #same indices, so do a contraction
                            if metric.indices[1][2] == tens.indices[1][2] #same positions
                                if metric.indices[1][1] == tens.indices[1][1]
                                    factor *= metric.indices[1][1].index_set.range
                                else
                                    factor *= -metric.indices[1][1].index_set.range
                                end
                            else
                                if metric.indices[1][1] == tens.indices[1][1]
                                    factor *= -metric.indices[1][1].index_set.range
                                else
                                    factor *= metric.indices[1][1].index_set.range
                                end
                            end
                            # println(factor)
                            push!(to_remove, i)
                            push!(to_remove, j)
                            break
                        elseif !isempty(intersect(get_free_indices_nopos(metric), get_free_indices_nopos(tens)))
                            #shared index, turn to a kronecker

                            # #check if both upper/lower
                            # println(metric.indices[1][2].position)
                            # println( tens.indices[1][2].position  )
                            if metric.indices[1][2].position == tens.indices[1][2].position  #same positions
                                #    println("same pos")

                                if metric.indices[1][1] == tens.indices[1][1] #eg eps^{ij} eps^{ik}

                                    if metric.indices[1][2].position == upper.position
                                        # println("c")
                                        push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[2][1] => lower]))
                                    else
                                        # println("d")
                                        push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[2][1] => upper]))
                                    end

                                    # factor *= metric.indices[1][1].index_set.range
                                elseif metric.indices[2][1] == tens.indices[2][1] #2nd index same, #eg eps^{ij} eps^{kj}
                                    if metric.indices[1][2].position == upper.position
                                        push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[1][1] => lower]))
                                        # println("e")
                                    else
                                        push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[1][1] => upper]))
                                        # println("f")
                                    end
                                elseif metric.indices[1][1] == tens.indices[2][1] #alt index same, #eg eps^{ij} eps^{ki}
                                    if metric.indices[1][2].position == upper.position
                                        push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[1][1] => lower]))
                                        # println("g")
                                    else
                                        push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[1][1] => upper]))
                                        # println("h")
                                    end

                                    factor *= (-1)
                                    # println(factor)

                                elseif metric.indices[2][1] == tens.indices[1][1] #alt index same, #eg eps^{ij} eps^{jk}
                                    if metric.indices[1][2].position == upper.position
                                        push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[2][1] => lower]))
                                        # println("i")
                                    else
                                        push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[2][1] => upper]))
                                        # println("j")
                                    end
                                    factor *= (-1)
                                    # println(factor)
                                    #  factor *= -metric.indices[1][1].index_set.range
                                end
                            elseif metric.indices[1][2].position != tens.indices[1][2].position   #opposite positions 
                                # println("opp pos")

                                if metric.indices[1][1] == tens.indices[1][1]    #eg eps^{ij} eps_{ik}
                                    push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[2][1] => tens.indices[2][2]]))
                                    factor *= (-1)
                                    # println(factor)
                                    # println("k")
                                elseif metric.indices[2][1] == tens.indices[2][1] #2nd index same, #eg eps^{ij} eps_{kj}
                                    push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[1][1] => tens.indices[1][2]]))
                                    factor *= (-1)
                                    # println("l")
                                elseif metric.indices[1][1] == tens.indices[2][1] #alt index same, #eg eps^{ij} eps_{ki}

                                    push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[1][1] => tens.indices[1][2]]))
                                    # println("m")

                                elseif metric.indices[2][1] == tens.indices[1][1] #alt index same, #eg eps^{ij} eps_{jk}
                                    push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[2][1] => tens.indices[2][2]]))
                                    # println("n")
                                end
                            end
                            # println(factor)
                            push!(to_remove, i)
                            push!(to_remove, j)
                            # println(metric," and ",tens," were contracted to make a delta")
                            break

                        end



                    elseif typeof(metric) <: KroneckerSuperType && typeof(tens) <: KroneckerSuperType
                        # println(metric, " and ",tens)
                        # see if they have same indices, and then contract
                        # println(get_free_indices_nopos(metric))
                        # println(get_free_indices_nopos(tens))
                        if get_free_indices_nopos(metric) == get_free_indices_nopos(tens)
                            factor *= metric.indices[1][1].index_set.range
                            push!(to_remove, i)
                            push!(to_remove, j)
                            break
                        elseif !isempty(intersect(get_free_indices_nopos(metric), get_free_indices_nopos(tens)))
                            #shared index, turn to a kronecker via contraction

                            if metric.indices[1][2].position == tens.indices[1][2].position  #same positions
                                #    println("same pos")

                                if metric.indices[1][1] == tens.indices[1][1] #eg eps^{ij} eps^{ik}
                                    push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[2][1] => tens.indices[1][2]]))

                                    # factor *= metric.indices[1][1].index_set.range
                                elseif metric.indices[2][1] == tens.indices[2][1] #2nd index same, #eg eps^{ij} eps^{kj}

                                    push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[1][1] => tens.indices[2][2]]))
                                    # println("e")

                                elseif metric.indices[1][1] == tens.indices[2][1] #alt index same, #eg eps^{ij} eps^{ki}

                                    push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[1][1] => tens.indices[1][2]]))

                                elseif metric.indices[2][1] == tens.indices[1][1] #alt index same, #eg eps^{ij} eps^{jk}

                                    push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[2][1] => tens.indices[2][2]]))

                                    # println(factor)
                                    #  factor *= -metric.indices[1][1].index_set.range
                                end
                            elseif metric.indices[1][2].position != tens.indices[1][2].position   #opposite positions 
                                # println("opp pos")

                                if metric.indices[1][1] == tens.indices[1][1]    #eg eps^{ij} eps_{ik}
                                    push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[2][1] => tens.indices[2][2]]))

                                    # println(factor)
                                    # println("k")
                                elseif metric.indices[2][1] == tens.indices[2][1] #2nd index same, #eg eps^{ij} eps_{kj}
                                    push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[1][1] => tens.indices[1][2]]))

                                elseif metric.indices[1][1] == tens.indices[2][1] #alt index same, #eg eps^{ij} eps_{ki}

                                    push!(extra_deltas, KroneckerDelta([metric.indices[2][1] => metric.indices[2][2], tens.indices[1][1] => tens.indices[2][2]]))
                                    # println("m")

                                elseif metric.indices[2][1] == tens.indices[1][1] #alt index same, #eg eps^{ij} eps_{jk}
                                    push!(extra_deltas, KroneckerDelta([metric.indices[1][1] => metric.indices[1][2], tens.indices[2][1] => tens.indices[1][2]]))
                                    # println("n")
                                end
                            end
                            # println(factor)
                            push!(to_remove, i)
                            push!(to_remove, j)
                            # println(metric," and ",tens," were contracted to make a delta")
                            break

                        end


                    end

                end
            end
        end
        sorted_idx = sort(to_remove, rev=true)
        for indexx in sorted_idx
            deleteat!(uncontr_mets, indexx)
        end

        uncontr_mets = vcat(extra_deltas, uncontr_mets)
        T.summands[k][2].terms = vcat(uncontr_mets, non_metrics)
        ###println(T.summands[k][1])
        ###println(T.summands[k][1]*sign)
        T.summands[k] = T.summands[k][1] * sign * factor => T.summands[k][2]
        ###println( T.summands[k])
    end
    ###println(T)
    return merge_sum(T)

end



##############################################################################################################################
##############################################################################################################################
### KroneckerDelta
##############################################################################################################################
##############################################################################################################################

""""
Eliminate Kronecker deltas by contracting indices or taking traces.
"""
function eliminate_kronecker(T::TensorTerm)
    eliminate_kronecker(TensorExpression(T))
end

""""
Eliminate Kronecker deltas by contracting indices or taking traces.
"""
function eliminate_kronecker(T::TensorExpression)
    contract_metrics(T, KroneckerSuperType)

end




##############################################################################################################################
##############################################################################################################################

### Eps Tensor


##############################################################################################################################

""""
Eliminate the Levi-Civita epsilon tensor by contracting indices or taking traces.
"""
function eliminate_epsilon(T::TensorTerm)
    eliminate_epsilon(TensorExpression(T))
end

""""
Eliminate the Levi-Civita epsilon tensor by contracting indices or taking traces.
"""
function eliminate_epsilon(T::TensorExpression)
    eliminate_kronecker(contract_metrics(T, EpsilonTensorSuperType))

end


##############################################################################################################################
##############################################################################################################################
### Mathematical ops, arithmetic 
##############################################################################################################################
##############################################################################################################################


########################
#addition
##

function +(A::SymorNum, B::TensorExpression)
    const_ten = Tensor("")
    if get_free_indices(const_ten) == get_free_indices(B)
        return A * const_ten + B
    elseif isempty(Set([x for x in get_free_indices(B) if x[1].index_set.name != "lightcone"]))
        return A * const_ten + B
    else
        throw(FreeIndexMismatch())
    end

end

+(A::TensorExpression, B::SymorNum) = +(B, A)
+(A::TensorTerm, B::SymorNum) = +(B, TensorExpression(A))
+(A::TensorSuperType, B::SymorNum) = +(B, TensorExpression(A))


function +(A::TensorSuperType, B::TensorSuperType)
    if get_free_indices(A) == get_free_indices(B)
        return TensorExpression([A, B])
    elseif isempty(Set([x for x in setdiff(get_free_indices(A), get_free_indices(B)) if x[1].index_set.name != "lightcone"]))
        return TensorExpression([A, B])
    else
        throw(FreeIndexMismatch())
    end

end

function +(A::TensorSuperType, B::TensorExpression)
    if isZeroTens(B) || get_free_indices(A) == get_free_indices(B)
        return TensorExpression(TensorExpression(A), B)
    elseif isempty(Set([x for x in setdiff(get_free_indices(A), get_free_indices(B)) if x[1].index_set.name != "lightcone"]))
        return TensorExpression(TensorExpression(A), B)
    else
        throw(FreeIndexMismatch())
    end
end
+(A::TensorExpression, B::TensorSuperType) = +(B, A)

# function +(A::TensorExpression, B::TensorSuperType)
#     if isZeroTens(A) || get_free_indices_nopos(A) == get_free_indices_nopos(B)
#         return TensorExpression(A, TensorExpression(B))
#     else
#         throw(FreeIndexMismatch())
#     end
# end

function +(A::TensorExpression, B::TensorExpression)
    if isZeroTens(A) || isZeroTens(B) || get_free_indices(A) == get_free_indices(B)
        return TensorExpression(A, B)

    elseif isempty(Set([x for x in setdiff(get_free_indices(A), get_free_indices(B)) if x[1].index_set.name != "lightcone"]))
        return TensorExpression(A, B)
    else
        println(get_free_indices(A), ":", get_free_indices(B))
        println(A)
        println(B)
        println(setdiff(get_free_indices(A), get_free_indices(B)))
        for summand in B.summands
            for tens in summand[2].terms
                println(tens)
                println(get_free_indices(tens))
            end
        end

        throw(FreeIndexMismatch())
    end
end


function +(A::TensorTerm, B::TensorExpression)

    if isZeroTens(B) || get_free_indices(A) == get_free_indices(B) || isempty(Set([x for x in setdiff(get_free_indices(A), get_free_indices(B)) if x[1].index_set.name != "lightcone"]))
        return TensorExpression(TensorExpression(A), B)
    else

        throw(FreeIndexMismatch())
    end
end

function +(B::TensorExpression, A::TensorTerm)

    if isZeroTens(B) || get_free_indices(A) == get_free_indices(B) || isempty(Set([x for x in setdiff(get_free_indices(A), get_free_indices(B)) if x[1].index_set.name != "lightcone"]))
        return TensorExpression(B, TensorExpression(A))
    else

        throw(FreeIndexMismatch())
    end
end

function +(A::TensorTerm, B::TensorTerm)

    if get_free_indices(A) == get_free_indices(B) || isempty(Set([x for x in setdiff(get_free_indices(A), get_free_indices(B)) if x[1].index_set.name != "lightcone"]))
        return TensorExpression(A) + TensorExpression(B)
    else
        throw(FreeIndexMismatch())
    end
end

function +(B::TensorSuperType, A::TensorTerm)

    if get_free_indices(A) == get_free_indices(B) || isempty(Set([x for x in setdiff(get_free_indices(A), get_free_indices(B)) if x[1].index_set.name != "lightcone"]))
        return TensorExpression(B) + TensorExpression(A)
    else
        ###println(get_free_indices(A))
        ###println(get_free_indices(B))
        throw(FreeIndexMismatch())
    end
end

function +(A::TensorTerm, B::TensorSuperType)

    if get_free_indices(A) == get_free_indices(B) || isempty(Set([x for x in setdiff(get_free_indices(A), get_free_indices(B)) if x[1].index_set.name != "lightcone"]))
        return TensorExpression(A) + TensorExpression(B)
    else
        ###println(get_free_indices(A))
        ###println(get_free_indices(B))
        throw(FreeIndexMismatch())
    end
end


########################
#multiplication
##

#helper funcs to ensure indices are valid
function valid_index_placement(T1, T2)
    if check_common_dummies(T1, T2)
        throw(OverusedIndices())
    elseif !check_common_are_lightcone(T1, T2)
        throw(SameLevelDummieError())
    else
        return true
    end
end
function check_common_are_lightcone(T1, T2)
    common = intersect(T1.indices, T2.indices)
    return isempty(common) || all(ind -> ind[1].index_set.name == "lightcone", common)
end
function check_common_are_lightcone(T1::TensorExpression, T2)
    common = intersect(get_inds(T1), get_inds(T2))
    return isempty(common) || all(ind -> ind[1].index_set.name == "lightcone", common)
end
check_common_are_lightcone(T1, T2::TensorExpression) = check_common_are_lightcone(T2, T1)
check_common_are_lightcone(T1::TensorExpression, T2::TensorExpression) = check_common_are_lightcone(T1, T2)


function get_inds(T)
    if typeof(T) <: TensorExpression
        return T.summands[1][2].indices
    else
        return T.indices
    end
end


"""Check if free indices of one tensor are dummies of the other. Return tue if so."""
function check_common_dummies(T1, T2)
    d1 = get_dummies(T1)
    d2 = get_dummies(T2)
    f1 = get_free_indices_nopos(T1)
    f2 = get_free_indices_nopos(T2)
    dp1 = get_dummies_positional(T1)
    dp2 = get_dummies_positional(T2)
    return !(isempty(intersect(d1, f2)) && isempty(intersect(d2, f1)) && isempty(intersect(dp1, dp2)))
    # return isempty(common) || all(ind -> ind[1].index_set.name == "lightcone", common)
end


#basic mult
*(N::SymorNum, B::TensorSuperType) = TensorExpression(N, copy(B))
*(N::SymorNum, B::TensorTerm) = TensorExpression(N => B)
*(B::TensorSuperType, N::SymorNum) = TensorExpression(N, copy(B)) #i dont remember why i put copy...
*(B::TensorTerm, N::SymorNum) = TensorExpression(N => B)

function *(B::TensorExpression, N::SymorNum)
    if B == zero
        return zero
    end
    B.summands = [convert(SymorNum, N) * summand[1] => summand[2] for summand in B.summands]
    return update_indices(B)
end
*(N::SymorNum, B::TensorExpression) = *(B, convert(SymorNum, N)) #fix for grading....?

function *(T::TensorSuperType, TEx::TensorExpression)
    ###println("doing this sort of multiplcation1")
    new_exp = zero
    for sum in TEx.summands
        new_exp += sum[1] * (T * sum[2])#should stop overflow?
    end
    return update_indices(new_exp)
end
function *(TEx::TensorExpression, T::TensorSuperType)
    if valid_index_placement(TEx, T)
        ###println("doing this sort of multiplcation1")
        new_exp = zero
        for sum in TEx.summands
            new_exp += sum[1] * (sum[2] * T)#should stop overflow?
        end
        return update_indices(new_exp)
    end
end

function *(T::TensorTerm, TEx::TensorExpression)
    ###println("doing this sort of multiplcation2")
    return update_indices(TensorExpression([sum[1] => T * sum[2] for sum in TEx.summands]))
end
# *(TEx::TensorExpression, T::TensorTerm) = *(T, TEx)

function *(TEx::TensorExpression, T::TensorTerm)
    ###println("doing this sort of multiplcation3")
    return update_indices(TensorExpression([sum[1] => sum[2] * T for sum in TEx.summands]))
end

function *(T1::TensorSuperType, T2::TensorSuperType)
    if valid_index_placement(T1, T2)
        return update_indices(TensorTerm([copy(T1), copy(T2)]))
    end
end


function *(TEx::TensorTerm, T::TensorSuperType)
    if valid_index_placement(TEx, T)
        return update_indices(TensorTerm(vcat(TEx.terms, copy(T))))
    end
end
function *(T::TensorSuperType, TEx::TensorTerm)
    if valid_index_placement(T, TEx)
        return update_indices(TensorTerm(vcat(copy(T), TEx.terms)))
    end
end
# *(T::TensorSuperType, TEx::TensorTerm) = *(TEx, T)

function *(TEx1::TensorExpression, TEx2::TensorExpression)
    new_tex = zero
    if TEx1 == zero || TEx2 == zero
        return zero
    end
    for T in TEx1.summands
        for T2 in TEx2.summands
            new_tex += TensorExpression(T[1] * T2[1] => T[2] * T2[2])
        end
    end
    return update_indices(new_tex)
end
# *(TEx::TensorExpression, T::TensorSuperType)=*(TEx,TensorExpression(T))
# *(T::TensorSuperType,TEx::TensorExpression, )=*(TEx,TensorExpression(T))

# function *(TEx::TensorExpression, T::TensorSuperType)
#     ###println("doing this sort of multiplcation4")
#     if TEx == zero
#         return zero
#     end
#     for summandind in eachindex(TEx.summands)
#         ###println(TEx.summands[summandind])
#         # println("thuis")
#         # println(T)
#         TEx.summands[summandind] = TEx.summands[summandind][1] => TEx.summands[summandind][2] * T
#         ###println(TEx.summands[summandind])
#         ###println("here")
#     end
#     #ex=(sum(sum[1] => sum[2] * T for sum in TEx.summands))
#     ###println(TEx)
#     return update_indices(TEx)
# end

# *(TEx::TensorExpression,T::TensorSuperType)=*(T,TEx)

# function *(T::TensorSuperType, TEx::TensorTerm)
#     return update_indices(TensorTerm(vcat(copy(T), TEx.terms)))
# end
# *(TEx::TensorTerm, T::TensorSuperType) = *(T, TEx)

function *(T::TensorTerm, TEx::TensorTerm)
    if valid_index_placement(T, TEx)
        return update_indices(TensorTerm(vcat(T.terms, TEx.terms)))
    end
end

########################
#subtraction (defined via mult and add)
##

-(A::TensorSuperType, B::TensorSuperType) = TensorExpression(-1 * B + A)
-(A::TensorSuperType, B::TensorExpression) = TensorExpression(A + (-1 * B))
-(A::TensorExpression, B::TensorSuperType) = TensorExpression(A + (-1 * B))
-(A::TensorExpression, B::TensorExpression) = TensorExpression(A + (-1 * B))
-(A::TensorTerm, B::TensorExpression) = 1 * A - B
-(A::TensorTerm, B::TensorTerm) = 1 * A - B
-(A::TensorExpression, B::TensorTerm) = A - (1 * B)
-(A::TensorTerm, B::TensorSuperType) = (1 * A) - (1 * B)
-(A::TensorSuperType, B::TensorTerm) = (1 * A) - (1 * B)

#negation
-(A::TensorSuperType) = -1 * A
-(A::TensorTerm) = -1 * A
-(A::TensorExpression) = -1 * A



########################
#division (by numerical factors)
##

/(B::TensorSuperType, N::SymorNum) = *(1 / N, B)
/(B::TensorExpression, N::SymorNum) = *(1 / N, B)
/(B::TensorTerm, N::SymorNum) = *(1 / N, B)

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
    for term in T.terms
        if typeof(term) <: OperatorSuperType
            weight += operator_weight(term, w)
        elseif haskey(term.weights, w)
            weight += term.weights[w]
        end
    end
    return weight
end

"""
Check if a TensorExpression is a single TensorTerm.
"""
function singleterm(TEx::TensorExpression)
    if length(TEx.summands) == 1 && TEx.summands[1][1] != 0
        return true
    end
    return false
end
"""
Check if a TensorTerm is a single Tensor.
"""
function singleterm(TT::TensorTerm)
    if length(TT.terms) == 1
        return true
    end
    return false
end
"""
Get all summands in a TensorExpression of a particular weight. For example, this may be terms that carry a particular polynomial weighting.
"""
function get_weighted_terms(T::TensorExpression, w::String, value::Number)
    new_exp = zero
    for summand in T.summands
        if term_weight(summand[2], w) == value
            new_exp += summand[1] * summand[2]
        end
    end
    return new_exp
end

"""
Get all summands in a TensorExpression of below particular weight. For example, this may get leading order terms in a Taylor expansion.
Returns only weights<value. Assumes weights are non-negative.
"""
function drop_weight_above(T::TensorExpression, w::String, value::Number)
    new_exp = zero
    i = 0
    while i < value
        new_exp += get_weighted_terms(T, w, i)
        i += 1
    end
    return new_exp
end

"""
Get all summands in a TensorExpression of above particular weight. For example, this may get terms with dominant powers.
Returns only weights>value. Assumes weights are non-negative.
"""
function drop_weight_below(T::TensorExpression, w::String, value::Number)
    new_exp = zero
    i = 0
    while i <= value
        new_exp += get_weighted_terms(T, w, i)
        i += 1
    end
    return T - new_exp
end

# function +(N::SymorNum, T::TensorExpression)



#     if get_free_indices(A) == get_free_indices(B)
#         return TensorExpression(TensorExpression(A), B)
#     else
#         throw(FreeIndexMismatch())
#     end
# end

##############################################################################################################################


# function -(A::TensorSuperType, B::TensorExpression)
#     if get_free_indices(A)==get_free_indices(B)
#         return TensorExpression(TensorExpression(A),B)
#     else
#         throw(FreeIndexMismatch())
#     end
# end

# function -(A::TensorExpression, B::TensorSuperType)
#     if get_free_indices(A)==get_free_indices(B)
#         return TensorExpression(A,TensorExpression(B))
#     else
#         throw(FreeIndexMismatch())
#     end
# end

# function -(A::TensorExpression, B::TensorExpression)

#     if get_free_indices(A)==get_free_indices(B)
#         return TensorExpression(A,B)
#     else
#         throw(FreeIndexMismatch())
#     end
# end


##############################################################################################################################
##############################################################################################################################
### Advanced Mathematical Operations
##############################################################################################################################
##############################################################################################################################

# juxtaposition/appension/smart times, basically appends a tensor to a term so that all of its indices are free. 
# useful if you want to mutliply something without doing contractions, or want T^3 without manually having to create objects with new indices

# function ×(T1::TensorTerm, T2::TensorSuperType)


# end
function ×(T2::SymorNum, N::Union{TensorExpression,TensorTerm,TensorSuperType})
    N × T2
end

function ×(T2::Union{TensorExpression,TensorTerm,TensorSuperType}, N::SymorNum)
    # T2*N 
    throw(OperatorPrecedenceError())
end


function ×(T2::TensorExpression, T1::TensorExpression)
    new = zero
    used_inds = union(get_dummies(T2), get_free_indices_nopos(T1))
    #need to make sure all different terms have the same free inds
    # unused_flav = Vector(setdiff(flavour_indices.indices, used_inds))
    # unused_spin = Vector(setdiff(spinor_indices.indices, used_inds))
    # # unused_conj = Vector(setdiff(spinor_conj_indices.indices, used_inds))
    ## something like this for all index sets in the the expression

    for summand2 in T2.summands
        for summand1 in T1.summands
            new += (summand2[1] * summand1[1]) * (summand2[2] × summand1[2])
        end
    end

    return new
end


function ×(T2::TensorExpression, T1::Union{TensorSuperType,TensorTerm})

    (-1)^((grading(T1) * grading(T2)) % 2) * (T1 × T2)
end

function ×(T2::Union{TensorSuperType,TensorTerm}, T1::TensorExpression)
    new_terms = TensorExpression()
    T1 = rename_dummies(T1)
    #   ##println(T1)
    #   ##println(T2)
    #   ##println("^in func")
    for summandin in eachindex(T1.summands)
        t = (T2 × T1.summands[summandin][2])
        new_terms += (T1.summands[summandin][1] * t)
        ###println(t)
        ###println(T1.summands[summandin])
        # T1.summands[summandin]=T1.summands[summandin][1] => t
    end

    return new_terms

end

# function ×(T1::TensorTerm, T2::TensorTerm)

#     for term in T2.terms
#         T1 = T1 × term
#     end
#     return T1
# end



function ×(T2::TensorSuperType, T1::TensorTerm)
    return (-1)^((grading(T1) * grading(T2)) % 2) * TensorExpression(T1 × T2)

end

function ×(T1::Union{TensorTerm,TensorSuperType}, T2::T) where {T<:TensorSuperType}
    ###println("crash func")
    ###println(T1)
    ###println(typeof(T1)) #some different behaviour between tensor and tensor term. tensor works..
    ###println(T2)
    T2 = copy(T2)
    t1inds = OrderedSet{IndexSuperType}()
    t2inds = OrderedSet{Pair{Integer,IndexSuperType}}() #all indices, numbered in order
    # t2ok=Set{Index}()
    for inds in T1.indices
        push!(t1inds, inds[1])
    end
    for inds in eachindex(T2.indices)
        if T2.indices[inds][1] in t1inds #if theres a repeated indice
            push!(t2inds, inds => T2.indices[inds][1]) # which index to replace
        else
            # push!(t2ok,T2.indices[inds][1])# index doesnt need to be replaced
        end
    end
    # we now have all the indices we need to replace and where they are on the tensor
    sets = get_index_sets(T2) #Set{IndexSet}()
    avail = Dict{IndexSet,OrderedSet{IndexSuperType}}() #pairs of index sets, with the unused indices
    for set in sets
        a = OrderedSet{IndexSuperType}(intersect(get_unused_indices(T1, set), get_unused_indices(T2, set)))
        push!(avail, set => a) #ones we can replace with
    end
    ###println(avail)
    ###println(T2.indices)
    dummie_matchup = Dict{IndexSuperType,IndexSuperType}()

    t2dums = get_dummies(T2)

    for ind in t2inds # for each we need to replace (vec pos => index)
        ###println(avail[ind[2].index_set])
        ###println(ind)

        if ind[2] in t2dums

            if ind[2] in keys(dummie_matchup)
                nind = pop!(dummie_matchup, ind[2])
                ###println(nind)
                T2.indices[ind[1]] = nind => T2.indices[ind[1]][2]
            else
                if isempty(avail[ind[2].index_set])
                    throw(OutOfIndices)
                end
                nind = popfirst!(avail[ind[2].index_set])
                push!(dummie_matchup, ind[2] => nind)
                T2.indices[ind[1]] = nind => T2.indices[ind[1]][2]
            end
        else
            if isempty(avail[ind[2].index_set])
                throw(OutOfIndices)
            end
            nind = popfirst!(avail[ind[2].index_set])
            ###println(nind)
            T2.indices[ind[1]] = nind => T2.indices[ind[1]][2] #should replace
        end

    end
    # println("before crashhhhh") #same dummies
    # println(T1)
    # println(T2)
    # println(T2.indices)
    return T1 * T2
end


function ×(T1::Union{TensorTerm,TensorSuperType}, T2::TensorTerm)
    ###println("crash func")
    ###println(T1)
    ###println(typeof(T1)) #some different behaviour between tensor and tensor term. tensor works..
    ###println(T2)
    T2 = copy(T2)
    t1inds = OrderedSet{IndexSuperType}()
    t2inds = OrderedSet{Pair{Integer,IndexSuperType}}() #all indices, numbered in order
    # t2ok=Set{Index}()
    for inds in T1.indices
        push!(t1inds, inds[1])
    end

    sets = get_index_sets(T2) #Set{IndexSet}()
    avail = Dict{IndexSet,OrderedSet{IndexSuperType}}() #pairs of index sets, with the unused indices
    for set in sets
        a = OrderedSet{IndexSuperType}(intersect(get_unused_indices(T1, set), get_unused_indices(T2, set)))
        push!(avail, set => a) #ones we can replace with
    end
    ###println(avail)
    ###println(T2.indices)
    dummie_matchup = Dict{IndexSuperType,IndexSuperType}()

    t2dums = get_dummies(T2)

    for tens in T2.terms
        for inds in eachindex(tens.indices)
            if tens.indices[inds][1] in t1inds #if theres a repeated indice
                # push!(t2inds, tens =>(inds => tens.indices[inds][1])) # which index to replace
                ind = inds => tens.indices[inds][1]
                #     else
                #         # push!(t2ok,T2.indices[inds][1])# index doesnt need to be replaced
                #     end
                # end 
                # for ind in t2inds # for each we need to replace (vec pos => index)
                #     ###println(avail[ind[2].index_set])
                #     ###println(ind)

                if ind[2] in t2dums

                    if ind[2] in keys(dummie_matchup)
                        nind = pop!(dummie_matchup, ind[2])
                        ###println(nind)
                        tens.indices[ind[1]] = nind => tens.indices[ind[1]][2]
                    else
                        if isempty(avail[ind[2].index_set])
                            throw(OutOfIndices)
                        end
                        nind = popfirst!(avail[ind[2].index_set])
                        push!(dummie_matchup, ind[2] => nind)
                        tens.indices[ind[1]] = nind => tens.indices[ind[1]][2]
                    end
                else
                    if isempty(avail[ind[2].index_set])
                        throw(OutOfIndices)
                    end
                    nind = popfirst!(avail[ind[2].index_set])
                    ###println(nind)
                    tens.indices[ind[1]] = nind => tens.indices[ind[1]][2] #should replace
                end

            end
        end

    end
    # we now have all the indices we need to replace and where they are on the tensor
    update_indices(T2)
    # println("before") #same dummies
    # println(T1)
    # println(T2)
    # println(T2.indices)
    return T1 * T2
end

##############################################################################################################################
##############################################################################################################################
#equality checking




##############################################################################################################################

##############################################################################################################################

#equality checking

function ==(a::TensorTerm, b::TensorTerm)
    return a.terms == b.terms
end

function hard_equals(a::TensorTerm, b::TensorTerm)
    return same_tensorterm_mod_syms(a, b)
end

function ≅(a::TensorTerm, b::TensorTerm)
    return hard_equals(a, b)
end


function same_summands(a, b)
    eq(x, y) = is_equal(x.first, y.first) && (singleterm(x.second) && singleterm(y.second) ? same_tensor_mod_syms(x.second.terms[1], y.second.terms[1])[1] :  same_tensorterm_mod_syms(x.second, y.second))
    length(a) == length(b) || return false

    used = falses(length(b))

    for x in a
        found = false
        for j in eachindex(b)
            # println(x.second)
            # println(b[j].second)
            # println(same_tensorterm_mod_syms(x.second, b[j].second))
            if !used[j] && eq(x, b[j])
                used[j] = true
                found = true
                break
            end
        end
        # println(x)
        # println(found)
        found || return false
    end

    return true
end


function ==(a::Vector{T} where {T<:TensorSuperType}, b::Vector{T} where {T<:TensorSuperType})
    if length(a) != length(b)
        return false
    end
    for i = eachindex(a)
        if a[i] != b[i]
            return false
        end
    end
    return true
end

function ≅(a::TensorSuperType, b::TensorSuperType)#change to use our tabelau code!
    (same, grweight) = same_tensor_mod_syms(a, b)
    return (same && grweight == 1)
    # if (a.name == b.name && a.indices == b.indices )
    #     return true
    # else
    #     return false
    # end
end


function ==(a::TensorSuperType, b::TensorSuperType)#change to use our tabelau code!
    if (a.name == b.name && a.indices == b.indices)
        return true
    else
        return false
    end
end

# function Base.isequal(a::Pair{<:SymorNum,TensorTerm},b::Pair{SymorNum,TensorTerm})
#     println("aprichecl")
#     return a[1]==b[1] && a[2]==b[2]
# end

# function Base.hash(a::Pair{<:SymorNum,TensorTerm},h::UInt)
#     return Base.hash(a.first, Base.hash(a.second, h))
# end

function ==(a::TensorExpression, b::TensorExpression)
    if length(a.summands) != length(b.summands)
        # println("stop here")
        # println(a)
        # println(b)
        return false
    end
    # println("not stop here")
    # println(a)
    # println(b)
    aa = a.summands
    bb = b.summands
    if same_summands(aa, bb)
        return true
    end
    # println("summands")
    # println(aa)
    # println(bb)
    # # println(typeof(first(aa)))
    # # println(typeof(first(bb)))
    # println("symdiff")
    # for summanda in aa
    #     println(summanda)
    #     println(Base.hash(summanda[2]))
    # end
    # for summandb in bb
    #     println(summandb)
    #     println(Base.hash(summandb[2]))
    # end
    # # println(bb)
    # # println(first(aa))
    # # println(bb[3])
    # # println(Base.hash(aa[1]))
    # # println(Base.hash(bb[3]))
    # println(symdiff(aa,bb))
    # if isempty(symdiff(a.summands,b.summands))
    #     return true
    # end
    # for inde in eachindex(aa)
    #     if xor(typeof(aa[inde][1]) <: Number, typeof(bb[inde][1]) <: Number) || aa[inde][1] != bb[inde][1] || aa[inde][2] != bb[inde][2]
    #         return false
    #     end
    # end

    return false

end


function is_equal(a::SymorNum, b::SymorNum)
    !(xor(typeof(a) <: Number, typeof(b) <: Number) || a != b)
end



##############################################################################################################################
##############################################################################################################################

#indices

##############################################################################################################################


#update indices

function update_indices(T::TensorTerm)
    inds = []
    for t in T.terms
        for i in t.indices
            push!(inds, i)
        end
    end
    # println(inds)
    T.indices = inds
    ordered_counts = OrderedDict((pair[1], 0) for pair in T.indices)
    ###println(ordered_counts)
    for pair in T.indices
        ordered_counts[pair[1]] += 1
    end
    ###println(ordered_counts)
    duplicates = filter(kv -> kv[2] > 2, ordered_counts)
    ###println(duplicates)
    if !isempty(duplicates)
        for d in duplicates
            if d[1].index_set.name != "lightcone"
                println(duplicates)
                throw(OverusedIndices())
            end
        end
    end
    return T
end

function update_indices(T::TensorExpression)
    for summand in T.summands
        update_indices(summand[2])
    end
    return T
end



#get unused indices

function get_unused_indices(T::TensorSuperType, IS::IndexSet)
    ###println(Set(i[1] for i in T.indices))
    ###println(IS.indices)
    return OrderedSet(ordered_setdiff(IS.indices, OrderedSet(i[1] for i in T.indices)))
end

function get_unused_indices(T::TensorTerm, IS::IndexSet)
    ###println("indices of tensor term are: ",T.indices)
    ###println(Set(i[1] for i in T.indices))
    ###println(IS.indices)
    return OrderedSet(ordered_setdiff(IS.indices, OrderedSet(i[1] for i in T.indices)))
end


function get_unused_indices(T::TensorExpression, IS::IndexSet)
    ###println(Set(i[1] for i in T.indices))
    ###println(IS.indices)
    ###println("doing this for some reason")
    temp = IS.indices
    for summand in T.summands
        temp = OrderedSet(ordered_setdiff(temp, OrderedSet(i[1] for i in summand[2].indices)))
    end

    return temp
end


function get_index_index(T::TensorSuperType)
    return OrderedSet([T.indices[i][1] => i for i in 1:size(T.indices)[1]])
end

function get_index_index_rev(T::TensorSuperType)
    return Dict(i => T.indices[i][1] for i in 1:size(T.indices)[1])
end

function index_index_to_dict(s::OrderedSet)
    return Dict(p[2] => p[1] for p in s)
end

#finding free indices

function get_free_indices(T::TensorSuperType)
    return OrderedSet(ordered_setdiff(T.indices, get_dummies_positional(T)))
end

function get_free_indices(T::TensorTerm)
    ###println("doing this one")
    ###println(T.indices)
    d = get_dummies_positional(T)
    # ###println(d)
    # ###println("free")
    ###println(Set(intersect(T.indices,d)))
    return OrderedSet(ordered_setdiff(T.indices, d))
end

function get_free_indices(T::TensorExpression)
    if T == zero
        return Set()
    end
    return OrderedSet(ordered_setdiff(T.summands[1][2].indices, get_dummies_positional(T.summands[1][2])))#all tensors should have matching free indices in a sum already
end

function get_free_indices_nopos(T::Union{TensorExpression,TensorSuperType,TensorTerm})
    inds = get_free_indices(T)
    ind = OrderedSet()
    for d in inds
        push!(ind, d[1])
    end
    return ind
end

##############################################################################################################################

#finding dummie indices


function get_dummies_with_pos_index(T::TensorSuperType)

    return [i => pair[1] for (i, pair) in enumerate(T.indices) if pair[1] in get_dummies(T)]

end

function get_dummies_with_pos_index_and_pos_ul(T::TensorSuperType)

    return [i => pair for (i, pair) in enumerate(T.indices) if pair[1] in get_dummies(T)]

end

function get_dummies_with_pos_index_and_pos_ul(T::TensorTerm)

    return [i => pair for (i, pair) in enumerate(T.indices) if pair[1] in get_dummies(T)]

end


function get_dummies_positional(T::TensorSuperType)
    dum::OrderedSet{Pair{IndexSuperType,IndexPosition}} = OrderedSet()
    for d in get_dummies(T)
        push!(dum, d => upper)
        push!(dum, d => lower)
    end
    return dum
end

function get_dummies_positional(T::Union{TensorTerm,TensorExpression})
    dum::OrderedSet{Pair{IndexSuperType,IndexPosition}} = OrderedSet()
    for d in get_dummies(T)
        push!(dum, d => upper)
        push!(dum, d => lower)
    end
    return dum
end


function get_index_sets(T::Union{TensorSuperType,TensorTerm})
    temp = [ind[1] for ind in T.indices]
    out = OrderedSet{IndexSet}()
    for t in temp
        push!(out, t.index_set)
    end
    return out
end


function get_dummies(T::TensorExpression)
    out = OrderedSet()
    for summand in T.summands
        out = union(out, get_dummies(summand[2]))
    end
    return out
end


#returns a set of dummie indices that are used in summation
function get_dummies(T::Union{TensorSuperType,TensorTerm})

    temp = [ind[1] for ind in T.indices]
    ###println(temp)
    ordered_counts = OrderedDict((pair, 0) for pair in temp)
    for pair in temp
        ordered_counts[pair] += 1
    end
    duplicates = filter(kv -> kv[2] == 2, ordered_counts)
    poscheck = Dict()
    for ind in T.indices
        if haskey(duplicates, ind[1])
            if ind[1].index_set.name != "lightcone"
                if haskey(poscheck, ind[1])
                    if poscheck[ind[1]] == ind[2]
                        # println(T)
                        # println(ind)
                        throw(SameLevelDummieError())
                    end
                else
                    poscheck[ind[1]] = ind[2]
                end
            else
                pop!(duplicates, ind[1])
            end
        end
    end
    return collect(keys(duplicates))
end

# function get_dummies_indpos(T::Union{TensorSuperType,TensorTerm})

#     temp = [(T.indices[ind][1],ind) for ind in 1:length(T.indices)]
#     ###println(temp)
#     ordered_counts = OrderedDict((pair[1], (0, pair[2])) for pair in temp)
#     for pair in temp
#         ordered_counts[pair[1]] = (ordered_counts[pair[1]][1]+1,ordered_counts[pair[1]][2])
#     end
#     println(ordered_counts)
#     duplicates = filter(kv -> kv[2][1] == 2, ordered_counts)
#     poscheck = Dict()
#     for ind in T.indices
#         if haskey(duplicates, ind[1])
#             if ind[1].index_set.name != "lightcone"
#                 if haskey(poscheck, ind[1])
#                     if poscheck[ind[1]] == ind[2]
#                         println(T)
#                         println(ind)
#                         throw(SameLevelDummieError())
#                     end
#                 else
#                     poscheck[ind[1]] = ind[2]
#                 end
#             else
#                 pop!(duplicates,ind[1])
#             end
#         end
#     end
#     return collect(keys(duplicates))
# end

# function get_dummies(T::TensorTerm)
#     temp = [ind[1] for ind in T.indices]
#     ###println(temp)
#     ordered_counts = OrderedDict((pair, 0) for pair in temp)
#     for pair in temp
#         ordered_counts[pair] += 1
#     end
#     duplicates = filter(kv -> kv[2] == 2, ordered_counts)
#     # for ind in T.indices
#     #     if cind 
#     return collect(keys(duplicates))
#     # temp = [ind[1] for ind in T.indices]
#     # ###println(temp)
#     # ###println(size(temp))
#     # # if size(temp)[1] == 8
#     # #    ##println(temp[2])
#     # #    ##println(temp[8])
#     # ###println(temp[2]==temp[8])
#     # # end
#     # # name_pairs = OrderedDict((pair.name, pair) for pair in temp)
#     # ordered_counts = OrderedDict((pair.name, 0) for pair in temp)
#     # ###println("name pairs")
#     # ###println(name_pairs)
#     # for pair in temp
#     #     ordered_counts[pair.name] += 1
#     # end
#     # ###println(ordered_counts)
#     # duplicates = filter(kv -> kv[2] == 2, ordered_counts)
#     # # # ##println(duplicates)
#     # # # for ind in T.indices
#     # # #     if cind 
#     # # ###println("has dummies")
#     # # names=collect(keys(duplicates))
#     # ###println(names) # these are the dummies
#     # # inds=[]
#     # # for name in names # had to do some complicated stuff since sometimes dummies arent recognised, i think due to the way structs are compared?
#     # #     push!(inds, name_pairs[name])
#     # # end
#     # # # filter(kv -> names[kv[1]] != Nothing, name_pairs) #?????????? want to get all indices from the name
#     # ###println(collect(values(inds))) #dummie indices
#     # return collect(keys(duplicates))
# end





function rename_dummies(T::TensorSuperType)
    T = copy(T)
    dummies = get_dummies(T) #returns all the dummies attached to the tensor
    # println(dummies)
    all_indices::OrderedSet{IndexSuperType} = OrderedSet([i[1] for i in T.indices]) #set of all the indices attached to the tensor
    index_sets::OrderedSet{IndexSet} = OrderedSet([d.index_set for d in dummies]) #set of all the index_sets associated with the dummies
    for set in index_sets #for each type of dummie index

        replacements::Vector{Pair} = [] #dummie replacements to make
        inset = OrderedSet(filter(a -> a.index_set == set, all_indices)) # all indices on our tensor that are of the specified index type

        remain = OrderedSet(filter(x -> !(x in inset), OrderedSet(set.indices)))
        # remain=setdiff(OrderedSet(set.indices),inset) # possible indices to replace our dummies that aren't already used on the tensor
        dumset = filter(d -> d.index_set == set, dummies)  #dummie indices to replace of the appropriate type
        for d in dumset #attempt to replace each dummie index
            for r in remain #check candidates
                if first_element(set.indices, d, r) != d # if the dummie doesn't come first lexographically (or in your specified order)
                    push!(replacements, d => r)
                    delete!(remain, r)
                    break
                end
            end
        end
        ###println(replacements)
        for p in replacements

            ###println("here at least")
            if typeof(T) <: OperatorSuperType
                ###println("doing op replacement")
                # println("oper is: ",T)
                T.operator_indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in T.operator_indices]
                T.operand_indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in T.operand_indices]
                # println("new oper is: ", T)
                # println("with indeix", T.indices)
                T.indices = vcat(T.operator_indices, T.operand_indices)
                update_dummies(T.operand, p)
            else
                T.indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in T.indices]
            end
        end
        #     for (x, y) in T.indices
        #         (x,y)=(x == p[1] ? p[2] : x), y for p in replacements 
        #    # T.indices = [Pair((x == p[1] ? p[2] : x), y) for p in replacements for (x, y) in T.indices ]
    end
    return T
end



##############################################################################################################################





function update_ext_indices(T::TensorExpression)
    new_exp = zero
    for tens in T.summands
        if tens[1] == 0
            continue
        end
        update_ext_indices(tens[2])
        new_exp += tens[1] * tens[2]
    end
    T = new_exp
end


function update_ext_indices(T::TensorTerm)
    inds = []
    for tens in T.terms
        if typeof(tens) <: OperatorSuperType
            for summand in tens.operand.summands
                for term in summand[2].terms
                    if typeof(term) <: OperatorSuperType
                        update_ext_indices(term)
                    end
                end
                tens.operand_indices = tens.operand.summands[1][2].indices
                tens.indices = vcat(tens.operand_indices, tens.operator_indices)

            end
        end
        for i in tens.indices
            push!(inds, i)

        end
    end
    T.indices = inds

end




#rename dummies function

function rename_dummies(T::TensorTerm)
    T = copy(T)
    # println(T)
    dummies = get_dummies(T) #returns all the dummies attached to the tensor
    # println(dummies)
    # println("dummies are,",dummies)
    ###println("dummies")
    ###println(dummies)
    all_indices::OrderedSet{IndexSuperType} = OrderedSet([i[1] for t in T.terms for i in t.indices]) #set of all the indices attached to the tensor
    index_sets::OrderedSet{IndexSet} = OrderedSet([d.index_set for d in dummies]) #set of all the index_sets associated with the dummies
    for set in index_sets #for each type of dummie index

        replacements::Vector{Pair} = [] #dummie replacements to make
        inset = OrderedSet(filter(a -> a.index_set == set, all_indices)) # all indices on our tensor that are of the specified index type
        remain = OrderedSet(filter(x -> !(x in inset), OrderedSet(set.indices)))
        #remain=setdiff(OrderedSet(set.indices),inset) # possible indices to replace our dummies that aren't already used on the tensor
        dumset = filter(d -> d.index_set == set, dummies)  #dummie indices to replace of the appropriate type
        for d in dumset #attempt to replace each dummie index
            for r in remain #check candidates
                if first_element(set.indices, d, r) != d # if the dummie doesn't come first lexographically (or in your specified order)
                    push!(replacements, d => r)
                    delete!(remain, r)
                    break
                end
            end
        end
        ###println(replacements)
        for p in replacements
            for t in T.terms
                ###println("got here")
                if !isempty(t.indices)
                    ###println("here")

                    ###println(t)
                    if typeof(t) <: OperatorSuperType
                        ###println("doing op replacement")
                        # println("oper is: ",T)
                        t.operator_indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in t.operator_indices]
                        t.operand_indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in t.operand_indices]
                        # println("new oper is: ",T)

                        t.indices = vcat(t.operator_indices, t.operand_indices)

                        update_dummies(t.operand, p)
                        # println("new oper is: ", t)
                        # println("with indeix", t.indices)

                    else
                        t.indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in t.indices]
                    end
                end
            end
        end
        #     for (x, y) in T.indices
        #         (x,y)=(x == p[1] ? p[2] : x), y for p in replacements 
        #    # T.indices = [Pair((x == p[1] ? p[2] : x), y) for p in replacements for (x, y) in T.indices ]
    end
    update_indices(T)
    # println(T)
    return T
end

function rename_dummies(T::TensorExpression)
    # T=copy(T)
    new_exp = zero
    for term in T.summands
        new_exp += term[1] * rename_dummies(term[2])
    end
    # println("updated version...",new_exp)
    # update_indices(T)
    return new_exp
end


function update_dummies(T::TensorSuperType, reps::Pair{I,I}) where {I<:IndexSuperType}
    ###println(T)
    ###println(reps)
    if typeof(T) <: OperatorSuperType
        update_dummies(T.operand, reps)
        T.operator_indices = [Pair((x == reps[1] ? reps[2] : x), y) for (x, y) in T.operator_indices]
        T.operand_indices = [Pair((x == reps[1] ? reps[2] : x), y) for (x, y) in T.operand_indices]
    end

    if !isempty(T.indices)
        T.indices = [Pair((x == reps[1] ? reps[2] : x), y) for (x, y) in T.indices]
    end
end

function update_dummies(T::TensorTerm, reps::Pair{I,I}) where {I<:IndexSuperType}
    for t in T.terms
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
        for j in i+1:n
            if config[i] > config[j]
                inv += 1
            end
        end
    end

    return inv
end

"""Look through the symmetries of the tensor, and find the configuration of indices that is closest to the lexo-canonical order, as measured by the rep_distance function.
Returns the sign and the new configuration of indices. If the cost is zero, we can return immediately, as we have found the lexo-canonical representation."""
function find_canonical_rep(T::TensorSuperType)
    # sort indices into canonical order, serving as a baseline of cost zero.
    lexo = T.indices[sortperm(T.indices)]#this will automatically put in alphabetical
    best_inds = T.indices
    sign = 1
    best_cost = Inf
    for tableau in T.tableaux
        syms = tableau.symmetries
        for perm in syms
            new_inds = apply_permutation_cycle(T.indices, perm)
            cost = rep_distance(lexo, new_inds)
            if cost == 0
                return (1,new_inds)
            end
            if cost < best_cost || (cost == best_cost && new_inds < best_inds)
                best_inds = new_inds
                best_cost = cost
            end
        end
        asyms = tableau.antisymmetries
        for perm in asyms
            new_inds = apply_permutation_cycle(T.indices, perm)
            cost = rep_distance(lexo, new_inds)
            if cost == 0
                return (-1,new_inds)
            end
            if cost < best_cost || (cost == best_cost && new_inds < best_inds)
                best_inds = new_inds
                best_cost = cost
                sign = -1
            end
        end

    end
    # apply symmetries to the tensor, and measure the cost by the difference from the baseline.
    return (sign, best_inds)
end

"""Given a permutation cycle like "(1 5 2)", applies the corresponding permutation to the vector."""
function apply_permutation_cycle(v::Vector, cycle::String)
    result = copy(v)
    if cycle=="()"
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
    T_new=copy(T)
    T_new.indices = new_inds
    return sign * T_new
end
    # if typeof(T) <: TotallySymmetricSuperType
    #     if T.indices[1][1] == T.indices[2][1] && T.indices[1][2] == op_pos(T.indices[2][2])
    #         # println("deleting?????????")
    #         return zero
    #     end
    # indsorted = sortperm(T.indices)
    #     newinds = T.indices[indsorted]
    #     T.indices = newinds
    #     return T
    # elseif typeof(T) <: TotallyAntisymmetricSuperType
    #     indsorted = sortperm(T.indices)
    #     ###println(T.indices)
    #     ###println(indsorted)
    #     ###println(swaps)
    #     newinds = T.indices[indsorted]
    #     ###println(newinds)
    #     T.indices = newinds
    #     invcount = is_even_perm(indsorted)
    #     return (-1)^(!invcount) * T
    # end
#     return T
# end

function simplify(TEx::TensorTerm)
    simplify(TensorExpression(TEx))
end


function simplify(TEx::TensorExpression)

    exp = contract_metrics(TEx, EpsilonTensor)

    # for summand in exp.summands
    #     println(summand[2])
    #     for t in summand[2].terms
    #         println(t)
    #         println(t.indices)
    #         if typeof(t) <: OperatorSuperType

    #             println(t.operator_indices)
    #             println(t.operand_indices)
    #         end
    #     end
    # end
    # println("eiminating delta")
    exp = pop_metric(exp)

    # for summand in exp.summands
    #     println(summand[2])
    #     for t in summand[2].terms
    #         println(t)
    #         println(t.indices)
    #         if typeof(t) <: OperatorSuperType

    #             println(t.operator_indices)
    #             println(t.operand_indices)
    #         end
    #     end
    # end
    sort(exp)

    exp = eliminate_epsilon(exp)

    exp = rename_dummies(exp)

    # println("renamed:",exp)
    new_exp = simp_ops_n_sym(exp)

    new_exp = eliminate_kronecker(new_exp)
    # rename_dummies()
    # ##println("renamed: ",new_exp)
    # eliminate_epsilon(new_exp)
    new_exp = rename_dummies(new_exp)
    contract_metrics(new_exp, Metric)
    # hard_simplify(new_exp)
    return new_exp
end

function simp_ops_n_sym(exp::TensorExpression)
    new_exp = zero
    for summand in exp.summands
        # println("start ",summand)
        new_term = summand[1]
        for tens in summand[2].terms
            temp = tens
            if typeof(tens) <: OperatorSuperType
                # println("old oprd: ", tens.operand)
                op = simp_ops_n_sym(tens.operand)
                #  println("new oprd: ",op)
                if op == zero || op === nothing
                    new_term *= zero
                else
                    tens.operand = op #operator on nth term
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



function hard_simplify(TEx::TensorExpression)
    hard_merge(TEx)
    return TEx
end




#merging terms with same tensor content

## need to add functionality for symmetric and antisym tensors


function hard_merge(TEx::TensorExpression)
    if TEx == zero
        return zero
    end
    ##println("zero ma:",zero)
    out = filter!(x -> x[1] != 0, TEx.summands) #removes 0 coefficnets
    if isempty(out)
        ##println("killing")
        ##println("zero??:",zero)
        return zero
    end
    ##println("not killing")
    ###println(out)
    TEx.summands = out
    new_merged = []
    for summand in TEx.summands
        ###println("doing")
        found = false

        for i in eachindex(new_merged)
            if hard_equals(summand[2], new_merged[i][2])
                # println("got here")
                # println(summand[1])
                # println(new_merged[i][1])
                # println(summand[2])
                # println(new_merged[i][2])
                new_merged[i] = (new_merged[i][1] + summand[1] => summand[2])
                # ##println("sum is")
                # ##println(new_merged[i])
                # ##println(new_merged) #not editting in place apparaently....
                found = true
                break
                # elseif length(summand[2].terms)==1 && length(new_merged[i][2].terms)==1

                #     rescomp=same_tensor_mod_syms(summand[2].terms[1],new_merged[i][2].terms[1])
                #     if rescomp[1]==true
                #         # println(rescomp)
                #         new_merged[i]=(rescomp[2]*new_merged[i][1] +summand[1] =>summand[2])
                #         found=true
                #     end
            end

            # (res,val)=combine_same_tensors(summand[2],new_merged[i][2])
            # if res
            #     new_merged[i] = (new_merged[i][1] + val*summand[1] => new_merged[i][2])
            #     found = true
            #     break
            # end
        end
        # if found
        #    ##println("breaking....")
        #     break
        # end
        if !found
            push!(new_merged, summand)
        end
    end

    TEx.summands = new_merged
    TEx.summands = filter!(x -> x[1] != 0, TEx.summands)
    return isempty(TEx.summands) ? zero : TEx
end


function merge_sum(TEx::TensorExpression)
    if TEx == zero
        return zero
    end
    ##println("zero ma:",zero)
    out = filter!(x -> !isZero(x[1]), TEx.summands) #removes 0 coefficnets
    if isempty(out)
        ##println("killing")
        ##println("zero??:",zero)
        return zero
    end
    ##println("not killing")
    ###println(out)
    TEx.summands = out
    new_merged = []
    for summand in TEx.summands
        ###println("doing")
        found = false

        for i in eachindex(new_merged)
            if summand[2] == new_merged[i][2]
                # println("got here")
                # println(summand[1])
                # println(new_merged[i][1])
                # println(summand[2])
                # println(new_merged[i][2])
                new_merged[i] = (new_merged[i][1] + summand[1] => summand[2])
                # ##println("sum is")
                # ##println(new_merged[i])
                # ##println(new_merged) #not editting in place apparaently....
                found = true
                break
            elseif length(summand[2].terms) == 1 && length(new_merged[i][2].terms) == 1

                rescomp = same_tensor_mod_syms(summand[2].terms[1], new_merged[i][2].terms[1])
                if rescomp[1] == true
                    # println(rescomp)
                    new_merged[i] = (rescomp[2] * new_merged[i][1] + summand[1] => summand[2])
                    found = true
                end
            end

            # (res,val)=combine_same_tensors(summand[2],new_merged[i][2])
            # if res
            #     new_merged[i] = (new_merged[i][1] + val*summand[1] => new_merged[i][2])
            #     found = true
            #     break
            # end
        end
        # if found
        #    ##println("breaking....")
        #     break
        # end
        if !found
            push!(new_merged, summand)
        end
    end

    TEx.summands = new_merged
    TEx.summands = filter!(x -> !isZero(x[1]), TEx.summands)
    return isempty(TEx.summands) ? zero : TEx
end



##############################################################################################################################
##############################################################################################################################

#gradings


function grading(T::TensorSuperType)
    grade = 0
    for (ind, pos) in T.indices
        grade = grade + ind.grading
    end
    return grade % 2
end

function grading(T::TensorTerm)
    grade = 0
    for tens in T.terms
        grade = grade + grading(tens)
    end
    return grade % 2
end

function grading(T::TensorExpression)
    grade = 0
    for tens in T.summands[1][2].terms # should be same grading for everything in sum since they have matching free indices
        grade = grade + grading(tens)
    end
    return grade % 2
end


##############################################################################################################################
##############################################################################################################################

# printing and display

function Base.show(io::IO, T::TensorSuperType)
    print(io, tens2str(T))
end

function Base.show(io::IO, T::TensorTerm)
    for tens in T.terms
        print(io, tens)
    end
end

function Base.show(io::IO, T::TensorExpression)
    it = Iterators.Stateful(T.summands)
    for tens in it
        # @show (tens[1] == 1)
        print(io, !(typeof(tens[1]) <: Number) ? tens[1] : tens[1] == 1 ? "" : tens[1] == -1 ? "-" : !isreal(tens[1]) ? "(" * (!isZero(real(tens[1])) ? string(tens[1]) : !isZero(imag(tens[1])) ? string(imag(tens[1])) * "Im" : "") * ")" : string(real(tens[1])))
        print(io, tens[2])
        next = peek(it)
        if next !== nothing
            # println(next[1])
            if !isreal(next[1]) || real(next[1]) >= 0
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


function Base.sort(TEx::TensorExpression)
    new_exp = zero
    for summand in TEx.summands
        new_term = summand[1]
        # println("orig coef",summand[1])
        # println("orgin term",summand[2])
        new_terms = copy(summand[2])
        fac = bubble_sort_with_grading!(new_terms.terms)
        new_term *= (fac * new_terms)
        # println(fac)
        # ###println(summand[2])
        # println(new_term)
        new_exp += new_term
    end
    return new_exp
end


function Base.isless(a::TensorSuperType, b::TensorSuperType)
    if a.name < b.name
        return true
    elseif a.name == b.name
        return a.indices < b.indices
    end
    return false
end


function Base.sort(term::TensorTerm)
    new_terms = copy(term)
    fac = bubble_sort_with_grading!(new_terms.terms)
    return fac * new_terms
end


function custom_order_isless(T1::TensorSuperType, T2::TensorSuperType)
    if haskey(sort_order, T1.name)
        if haskey(sort_order, T2.name)
            if sort_order[T1.name] < sort_order[T2.name]
                return true
            else
                return false
            end
        else
            # println(T1.name," is less than ",T2.name)
            return true
        end
    elseif haskey(sort_order, T2.name)
        return false
    else
        return T1.name < T2.name
    end
end

function bubble_sort_with_grading!(v::Vector{T}) where {T<:TensorSuperType}
    n = length(v)
    swapped = true
    swaps = 0
    while swapped
        swapped = false
        for i in 1:n-1
            if !isempty(sort_order)
                if custom_order_isless(v[i+1], v[i])
                    v[i], v[i+1] = v[i+1], v[i]
                    swaps += grading(v[i]) * grading(v[i+1])
                    swapped = true
                end
            elseif isless(v[i+1], v[i])
                v[i], v[i+1] = v[i+1], v[i]
                swaps += grading(v[i]) * grading(v[i+1])
                swapped = true
            end
        end
    end
    return (-1)^(swaps)
end

function bubble_sort_index_with_grading!(T::TensorSuperType, mapping::Dict)

    n = length(mapping)
    # println(n)
    # println(mapping)
    swapped = true
    swaps = 1
    while swapped
        swapped = false
        for i in 1:n-1
            # println(i)
            if mapping[i+1] < mapping[i]

                #compute grading
                for tab in T.tableaux
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


##############################################################################################################################
##############################################################################################################################

# replacements and substitutions


function substitute(TEx::TensorExpression, pattern::String, repeat::Bool=true)
    pattern = strip(pattern) #remove spaces

    #first step: parse terms to a string, that points back to position in term
    #convert pattern to regex

    for summand in TEx.summands
        string = term2strpos(summand[2])

    end



end

function pattern2regex(pattern::String)
    # example "?D^{?i*}_{*?\alpha}[[\theta^{?j*}_{*?\beta}]]: \delta^{?j*}_{*?i} \delta^{?\alpha *}_{*?\beta}"
    reg = ""

    return reg
end

function term2strpos(Term::TensorTerm)
    out = Vector{Pair{String,Int64}}()
    for tind in eachindex(Term.terms)
        push!(out, tens2str(Term.terms[tind]) => tind)
    end
    return out
end


function to_latex(T::TensorTerm)
    out = ""
    for tens in T.terms
        out *= to_latex(tens)
    end
    return Base.replace(out * " ", "\$" => "")
end

function to_latex(T::TensorExpression)
    out = ""
    it = Iterators.Stateful(T.summands)
    for tens in it
        #  println(denominator(rationalize(float(imag(tens[1])))))
        out *= string(tens[1] == 1 ? " " : tens[1] == -1 ? "-" : !isreal(tens[1]) ? "(" * (real(tens[1]) != 0 ? string(tens[1]) : denominator(rationalize(float(imag(tens[1])))) == 1 ? string(float(imag(tens[1]))) * "i" : imag(tens[1]) < 0 ? "-\\frac{" * (numerator(rationalize(float(-1 * imag(tens[1])))) == 1 ? "" : string(numerator(rationalize(float(-1 * imag(tens[1])))))) * "i}{" * string(denominator(rationalize(float(imag(tens[1]))))) * "}" : "\\frac{" * (numerator(rationalize(float(imag(tens[1])))) == 1 ? "" : string(numerator(rationalize(float(imag(tens[1])))))) * "i}{" * string(denominator(rationalize(float(imag(tens[1]))))) * "}") * ")" : denominator(rationalize(float(Real(tens[1])))) == 1 ? string(float(real(tens[1]))) : real(tens[1]) < 0 ? "-\\frac{" * string(numerator(rationalize(float(-1 * Real(tens[1]))))) * "}{" * string(denominator(rationalize(float(Real(tens[1]))))) * "}" : "\\frac{" * string(numerator(rationalize(float(Real(tens[1]))))) * "}{" * string(denominator(rationalize(float(Real(tens[1]))))) * "}")
        # println(string(tens[1] == 1 ? " " : tens[1] == -1 ? " " : "*"))

        out *= string(tens[1] == 1 ? " " : tens[1] == -1 ? " " : "")
        out *= to_latex(tens[2])
        next = peek(it)
        if next !== nothing
            # println(next[1])
            if !isreal(next[1]) || real(next[1]) >= 0
                out *= "+"
            end
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

function to_latex_basic(T::TensorSuperType)
    out = ""
    out *= (T.name) * " "
    if !isempty(T.indices)
        lowers = ""
        uppers = ""
        lowersempt::Bool = 1
        uppersempt::Bool = 1
        for (ind, pos) in T.indices
            if ind.name == "±"
                ntemp = "\\pm"
            elseif ind.name == "∓"
                ntemp = "\\mp"
            else
                ntemp = (ind.name)
            end
            if pos.position == 0
                if typeof(ind) == DummyPatternIndex
                    uppers *= "?"
                end

                uppers = uppers * (ntemp) * " "
                uppersempt = 0

            else
                lowersempt = 0
                if typeof(ind) == DummyPatternIndex
                    lowers *= "?"
                end
                lowers = lowers * (ntemp) * " "
            end
        end
        if !uppersempt
            out *= "^{" * uppers * "}"
        end
        if !lowersempt
            out *= "_{" * lowers * "}"
        end
    end
    return Base.replace(out * " ", "\$" => "")
end

function to_latex(T::TensorSuperType)
    out = ""
    if T.name != ""
        out *= uni_to_TeX(T.name) * " "
    end
    if !isempty(T.indices)
        lowers = ""
        uppers = ""
        lowersempt::Bool = 1
        uppersempt::Bool = 1
        for (ind, pos) in T.indices
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
        # if !uppersempt
        #     out *= "^{" * uppers * "}"
        # end
        # if !lowersempt
        #     out *= "_{" * lowers * "}"
        # end
    end
    return Base.replace(out * " ", "\$" => "")
end


function tens2str(T::TensorSuperType)
    out = ""
    if typeof(T) == DummyPatternTensor
        # out*="?"
    end

    out *= T.name
    if !isempty(T.indices)
        lowers = ""
        uppers = ""
        lowersempt::Bool = 1
        uppersempt::Bool = 1
        for (ind, pos) in T.indices
            if pos.position == 0
                if typeof(ind) == DummyPatternIndex
                    uppers *= "?"
                end
                uppers = uppers * ind.name
                uppersempt = 0
                lowers = lowers * "*"
            else
                uppers = uppers * "*"
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
    return out

end




function replace(T::TensorTerm, old_T::TensorSuperType, new_T::TensorSuperType, repeat::Bool=true)
    replace(TensorExpression(T), old_T, new_T, repeat)
end

function replace(TEx::TensorExpression, old_T::TensorSuperType, new_T::TensorSuperType, repeat::Bool=true)
    for k in eachindex(TEx.summands)
        ###println(k)
        ###println(T.summands)
        summand = TEx.summands[k]
        matches = []
        for tind in eachindex(summand[2].terms)
            if summand[2].terms[tind] == old_T
                push!(matches, tind)
            end
        end
        # matches = filter(y -> y==old_T, summand[2].terms)
        # non_metrics = filter(y -> !(typeof(y) <: typ), summand[2].terms)
        # uncontr_mets = Vector{TensorSuperType}()
        # ###println(metrics) 
        # sign = 1 #for antisym metrics (spinor)
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
    return update_indices(TEx)
end
# """
#     young_projector_block(shape::Vector{Int}, labels::Vector{Int})

# Compute the true Young projector for an arbitrary Young tableau.

# - `shape`: row lengths of the tableau, e.g., [2,2]
# - `labels`: label vector, e.g., [1,2,3,4]

# Returns: Vector of `(sign, permuted_labels)` representing the projected tensor.
# """
# function young_projector_block(shape::Vector{Int}, labels::Vector{Int})
#     n = sum(shape)
#     @assert length(labels) == n "Labels length must match total boxes"

#     # --- 1. Build row and column positions ---
#     rows = []
#     idx = 1
#     for rlen in shape
#         push!(rows, collect(idx:idx+rlen-1))
#         idx += rlen
#     end

#     ncols = maximum(shape)
#     cols = []
#     for c = 1:ncols
#         col = Int[]
#         for r = 1:length(shape)
#             if shape[r] >= c
#                 push!(col, sum(shape[1:r-1])+c)
#             end
#         end
#         push!(cols, col)
#     end

#     # --- 2. Generate all permutations of a vector with sign ---
#     function perms_signed(v::Vector{Int})
#         results = []
#         function recurse!(a::Vector{Int}, i::Int, s::Int)
#             if i == length(a)
#                 push!(results, (s, copy(a)))
#             else
#                 for j = i:length(a)
#                     a[i], a[j] = a[j], a[i]
#                     new_s = s * (-1)^((j-i)%2)
#                     recurse!(a, i+1, new_s)
#                     a[i], a[j] = a[j], a[i]
#                 end
#             end
#         end
#         recurse!(copy(v), 1, 1)
#         return results
#     end

#     # --- 3. Column antisymmetrization (block-based) ---
#     col_terms = [(1, copy(labels))]
#     for col in cols
#         new_terms = []
#         perms = perms_signed(col)
#         for (sign0, lab0) in col_terms
#             for (sign1, perm) in perms
#                 lab_copy = copy(lab0)
#                 for (i, pos) in enumerate(col)
#                     lab_copy[pos] = lab0[perm[i]]
#                 end
#                 push!(new_terms, (sign0*sign1, lab_copy))
#             end
#         end
#         col_terms = new_terms
#     end

#     # --- 4. Row symmetrization (block-based) ---
#     row_terms = [(s, copy(lab)) for (s, lab) in col_terms]
#     for row in rows
#         new_terms = []
#         # generate all row permutations (sign = +1)
#         function row_perms(v)
#             if length(v) == 1
#                 return [v]
#             else
#                 out = []
#                 for i=1:length(v)
#                     rest = vcat(v[1:i-1], v[i+1:end])
#                     for p in row_perms(rest)
#                         push!(out, [v[i]; p])
#                     end
#                 end
#                 return out
#             end
#         end
#         perms = row_perms(row)
#         for (s0, lab0) in row_terms
#             for perm in perms
#                 lab_copy = copy(lab0)
#                 for (i, pos) in enumerate(row)
#                     lab_copy[pos] = lab0[perm[i]]
#                 end
#                 push!(new_terms, (s0, lab_copy))
#             end
#         end
#         row_terms = new_terms
#     end

#     # --- 5. Combine duplicates ---
#     term_dict = Dict{Vector{Int}, Int}()
#     for (s, lab) in row_terms
#         term_dict[lab] = get(term_dict, lab, 0) + s
#     end

#     # Only nonzero terms
#     return [(s, lab) for (lab,s) in term_dict if s != 0]
# end





##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# YOUNG PROJECTION
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

function young_projection(TEx::TensorExpression, proj::Union{Vector{String},Nothing}=nothing)
    tex = zero
    for summand in TEx.summands
        for tens in summand[2].terms
            if typeof(tens) <: TensorSuperType
                tex += summand[1] * young_projection(tens, proj)
            end
        end
    end
    return tex
end

# function count_terms_in_factor(factor::String)
#     # Use your existing split_additive_factors
#     return length(split_additive_factors(factor))
# end

# function compute_normalization_from_counts(factors::Vector{String})
#     # Multiply number of terms in each factor
#     total_terms = prod(count_terms_in_factor(f) for f in factors)
#     return 1 / total_terms
# end

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
    return 1 / total_weight
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
    projectors = proj === nothing ? [compute_young_projector(tab) for tab in T.tableaux] : [proj for tab in T.tableaux]
    # println(tab for tab in T.tableaux)
    norms = proj === nothing ? [1 // compute_hook_product(tab) for tab in T.tableaux] : [compute_normalization_from_factors(prr) for prr in projectors]
    tex = zero
    for (proj, norm) in zip(projectors, norms)
        tex += apply_young_projector(T, proj) * norm
    end
    return tex
end

function apply_young_projector(T::TensorExpression, projector::Vector{String})
    tex = zero
    for summand in T.summands
        for tens in summand[2].terms
            if typeof(tens) <: TensorSuperType
                tex += summand[1] * apply_single_projector!(tens, projector)
            end
        end
    end
    return tex
end

function apply_young_projector(T::TensorTerm, projector::Vector{String})
    tex = zero
    for tens in T.terms
        if typeof(tens) <: TensorSuperType
            tex += apply_young_projector!(tens, projector)
        end
    end
    return tex
end

function apply_young_projector(T::TensorSuperType, projector::Vector{String})
    terms = Vector{Tuple{Number,Vector{Pair{Alakazam.IndexSuperType,Alakazam.IndexPosition}}}}([(1, collect(T.indices))])
    for proj in reverse(projector)
        terms = apply_string_factor(terms, proj)
    end
    result = Dict{Vector,Number}()
    for (s, t) in terms
        result[t] = get(result, t, 0) + s
    end
    no_dupes = [(s, t) for (t, s) in result if s != 0]
    tex = zero #new tensor expression to hold result of projection
    for (prefact, inds) in no_dupes
        Tn = copy(T)
        Tn.indices = inds
        # println("prefact: ", prefact)
        tex += prefact * sort_symmetric(Tn)
        # println(tex)
    end
    return tex #normalisation fixed elsewhere above

end

function compute_hook_product(tab::YoungTableau)
    hook_lengths = compute_hook_length(tab)
    return prod(hook_lengths[hook_lengths.>0]; init=1)
end

function compute_hook_length(tab::YoungTableau)
    rows = tab.row_lengths
    n_rows = length(rows)
    n_cols = maximum(rows)

    hook_lengths = zeros(Int, n_rows, n_cols)

    for r in 1:n_rows
        for c in 1:rows[r]
            right = rows[r] - c
            below = sum(r2 >= c for r2 in rows[r+1:end]; init=0)
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
        factor = strip(factor[2:end-1])
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
    inner = factor[2:end-1]
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

# # Apply a single factor (could be product of swaps)
# function apply_string_factor(terms::Vector{Tuple{Int64, Vector{Pair{Alakazam.IndexSuperType, Alakazam.IndexPosition}}}}, factor::String)
#     # Determine sign
#     plus = occursin("+", factor)

#     # Find all (i j) pairs
#     pairs = [(parse(Int, m.captures[1]), parse(Int, m.captures[2])) 
#              for m in eachmatch(r"\((\d+) (\d+)\)", factor)]

#     new_terms = []
#     for (s, t) in terms
#         # Identity term
#         push!(new_terms, (s, t))
#         # Swapped term
#         push!(new_terms, (s * (plus ? 1 : -1), swap_string_pairs(t, pairs)))
#     end
#     return new_terms
# end


# # Apply a single factor "(1 ± (i j))" to current list of (sign, tensor)
# function apply_string_factor(terms::Vector{Tuple{Int64, Vector{Pair{Alakazam.IndexSuperType, Alakazam.IndexPosition}}}}, factor::String)
#     # Parse factor
#     plus = occursin("+", factor)
#     m = match(r"\((\d+) (\d+)\)", factor)
#     i = parse(Int, m.captures[1])
#     j = parse(Int, m.captures[2])

#     new_terms =Vector{Tuple{Int64, Vector{Pair{Alakazam.IndexSuperType, Alakazam.IndexPosition}}}}([])
#     for (s, t) in terms
#         # Identity term
#         push!(new_terms, (s, t))
#         # Swapped term
#         push!(new_terms, (s * (plus ? 1 : -1), swap_inds(t, i, j)))
#     end
#     # println(new_terms)
#     return new_terms
# end

# # permutation to cycle-list string in your style
# function perm_str(p::Vector{Int})
#     cycles = perm_to_cycles(p)

#     if isempty(cycles)
#         return "1"
#     end

#     parts = ["(" * join(c, " ") * ")" for c in cycles]
#     return join(parts, "")
# end

# # sign of permutation
# function perm_sign(p::Vector{Int})
#     n = length(p)
#     seen = falses(n)
#     s = 1

#     for i in 1:n
#         if !seen[i]
#             j = i
#             cycle_len = 0
#             while !seen[j]
#                 seen[j] = true
#                 j = p[j]
#                 cycle_len += 1
#             end
#             s *= (-1)^(cycle_len - 1)
#         end
#     end

#     return s
# end

# function perm_to_cycles(p::Vector{Int})
#     n = length(p)
#     visited = falses(n)
#     cycles = Vector{Vector{Int}}()

#     for i in 1:n
#         if !visited[i] && p[i] != i
#             cycle = Int[]
#             j = i
#             while !visited[j]
#                 push!(cycle, j)
#                 visited[j] = true
#                 j = p[j]
#             end
#             push!(cycles, cycle)
#         end
#     end

#     return cycles
# end

# # full symmetric-group sum (row symmetriser)
# function row_sym_sum(labels)
#     perms = collect(permutations(labels))
#     terms = [perm_str(p) for p in perms]
#     return "(" * join(terms, " + ") * ")"
# end

# # full signed sum (column antisymmetriser)
# function col_antisym_sum(labels)
#     perms = collect(permutations(labels))

#     terms = String[]
#     for p in perms
#         if perm_sign(p) == 1
#             push!(terms, "+ " * perm_str(p))
#         else
#             push!(terms, "- " * perm_str(p))
#         end
#     end

#     return "(" * join(terms, " ") * ")"
# end

# function compute_young_projector(tab::YoungTableau)

#     labels = tab.index_assignment
#     row_lengths = tab.row_lengths

#     # -------------------------
#     # build rows
#     # -------------------------
#     rows = []
#     idx = 1
#     for rlen in row_lengths
#         push!(rows, labels[idx:idx+rlen-1])
#         idx += rlen
#     end

#     # -------------------------
#     # build columns
#     # -------------------------
#     n_cols = maximum(row_lengths)

#     columns = []
#     for c in 1:n_cols
#         col = Int[]
#         for r in rows
#             if c <= length(r)
#                 push!(col, r[c])
#             end
#         end
#         push!(columns, col)
#     end

#     # -------------------------
#     # TRUE Young symmetriser:
#     # product of:
#     #   column antisym sums
#     #   row sym sums
#     # -------------------------
#     factors = String[]

#     # columns first (antisym)
#     for col in columns
#         if length(col) > 1
#             push!(factors, col_antisym_sum(col))
#         end
#     end

#     # rows second (sym)
#     for row in rows
#         if length(row) > 1
#             push!(factors, row_sym_sum(row))
#         end
#     end

#     return factors
# end


# permutation to cycle-list string in your style
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
    for i in 1:length(p)-1
        for j in i+1:length(p)
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
    cycles = Vector{Vector{Integer}}()

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
        push!(rows, labels[idx:idx+rlen-1])
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


# function compute_young_projector(tab::YoungTableau)
#     labels=tab.index_assignment
#     row_lengths=tab.row_lengths
#     # Split labels into rows
#     rows = []
#     idx = 1
#     for rlen in row_lengths
#         push!(rows, labels[idx:idx+rlen-1])
#         idx += rlen
#     end
#     n_rows = length(rows)
#     n_cols = maximum(row_lengths)

#     # Split labels into columns
#     columns = []
#     for c in 1:n_cols
#         col = Int[]
#         for r in 1:n_rows
#             if c <= length(rows[r])
#                 push!(col, rows[r][c])
#             end
#         end
#         push!(columns, col)
#     end

#     # Column antisymmetrizer factors: (1 - (i j))
#     col_factors = Vector{String}()
#     for col in columns
#         for (i,j) in transpositions(col)
#             push!(col_factors, "(1 - ($i $j))")
#         end
#     end

#     # Row symmetrizer factors: (1 + (i j))
#     row_factors = Vector{String}()
#     for row in rows
#         for (i,j) in transpositions(row)
#             push!(row_factors, "(1 + ($i $j))")
#         end
#     end

#     # Combine column antisymmetrizers first, then row symmetrizers
#     projector = vcat(row_factors,col_factors) #i just swapped THEM!
#     return projector
# end


# Function to generate transpositions for a row or column
function transpositions(indices::Vector)
    ts = []
    n = length(indices)
    for i in 1:n-1
        for j in i+1:n
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


function copy(T::TensorTerm)
    return typeof(T)([copy(t) for t in T.terms])
end

function copy(T::TensorSuperType)
    return typeof(T)(T.name, copy(T.indices), T.function_of, T.weights, T.tableaux)
end

function copy(T::TensorExpression)
    out = TensorExpression()
    if T == zero
        return out
    end

    for summand in T.summands
        ###println(summand)
        ###println(summand[1]*copy(summand[2]))
        out += summand[1] * copy(summand[2])
    end
    return out
end

# name::String
#     operand::TensorExpression
#     coord::Vector{T} where {T<:Coordinate}
#     operator_indices::Vector{Pair{Index, IndexPosition}}
#     operand_indices::Vector{Pair{Index, IndexPosition}}
#     indices::Vector{Pair{Index, IndexPosition}}
#     function_of::Vector{T} where {T<:Coordinate}


function copy(V::Vector{Pair{I,IndexPosition}}) where {I<:IndexSuperType}
    return Vector{Pair{I,IndexPosition}}([pair[1] => pair[2] for pair in V])
end
function copy(T::SymMetricSuperType)
    return typeof(T)(T.name, copy(T.indices))
end
function copy(T::KroneckerSuperType)
    return typeof(T)(copy(T.indices))
end
function copy(T::EpsilonTensorSuperType)
    return typeof(T)(copy(T.indices))
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

function hasIndices(T::TensorSuperType)
    return length(T.indices) != 0
end

function hasIndices(T::TensorTerm)
    return length(T.indices) != 0
end


const zero = TensorExpression()


function isZeroTens(T::TensorExpression)
    if length(T.summands) == 1 && isZero(T.summands[1][1])
        return true
    else
        return false
    end
end


function isZero(x)
    return ((typeof(x) <: Number) && x == 0)
end


# function isZeroTens(T::TEx)
#     if length(T.summands)==1 && T.summands[1][1]==0
#         return true
#     else
#         return false
#     end
# end

# Swap indices at positions i,j in a vector
function swap_inds(s::Vector, i::Int, j::Int)
    t = copy(s)
    t[i], t[j] = t[j], t[i]
    return t
end


function ordered_setdiff(set1::OrderedSet{T1}, set2::OrderedSet{T2}) where {T1,T2}
    result = OrderedSet{Union{T1,T2}}()
    for item in set1
        if !(item in set2)
            push!(result, item)
        end
    end
    return result
end

# Method for vectors
function ordered_setdiff(arr1::Vector{T1}, arr2::Vector{T2}) where {T1,T2}
    result = OrderedSet{Union{T1,T2}}()
    for item in arr1
        if !(item in arr2)
            push!(result, item)
        end
    end
    return result
end

# Method to handle ordered set and vector combination
function ordered_setdiff(set1::OrderedSet{T1}, arr2::Vector{T2}) where {T1,T2}
    result = OrderedSet{Union{T1,T2}}()
    for item in set1
        if !(item in arr2)
            push!(result, item)
        end
    end
    return result
end

# Method to handle vector and ordered set combination
function ordered_setdiff(arr1::Vector{T1}, set2::OrderedSet{T2}) where {T1,T2}
    result = OrderedSet{Union{T1,T2}}()
    for item in arr1
        if !(item in set2)
            push!(result, item)
        end
    end
    return result
end

function first_element(arr, elem1, elem2)
    idx1 = findfirst(x -> x == elem1, arr)
    idx2 = findfirst(x -> x == elem2, arr)
    if idx1 == nothing || (idx2 != nothing && idx2 < idx1)
        return elem2
    else
        return elem1
    end
end

function is_even_perm(p::AbstractVector{T}) where {T<:Real}
    n = length(p)
    if n <= 1
        return true
    end
    # Compute the number of inversions in the left and right halves of the input vector
    mid = div(n, 2)
    left = p[1:mid]
    right = p[mid+1:end]
    inv_left = is_even_perm(left)
    inv_right = is_even_perm(right)
    # Merge the sorted left and right halves while counting inversions
    merged = Vector{T}(undef, n)
    i, j, k = 1, 1, 0
    inv_merge = 0
    while i <= length(left) && j <= length(right)
        if left[i] <= right[j]
            merged[k+1] = left[i]
            i += 1
        else
            merged[k+1] = right[j]
            j += 1
            inv_merge += length(left) - i + 1
        end
        k += 1
    end
    while i <= length(left)
        merged[k+1] = left[i]
        i += 1
        k += 1
    end
    while j <= length(right)
        merged[k+1] = right[j]
        j += 1
        k += 1
    end
    # Return whether the number of inversions is even or odd
    return (inv_left + inv_right + inv_merge) % 2 == 0
end



##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# LETTERING
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
##############################################################################################################################
# unused atm, but potentially wanted for autogenerating more variable names when the user hasnt defined enough 
# (idealy, new indices would be generated from the same character set, eg latin, greek, hebrew etc)

# letters = OrderedSet{Char}()

# # Latin uppercase letters
# for i in 65:90
#    ##println(push!(letters, Char(i)))
# end

# # Latin lowercase letters
# for i in 97:122
#     push!(letters, Char(i))
# end

# # Latin dotted uppercase letters
# for i in 0x1E00:0x1E5F
#     push!(letters, Char(i))
# end

# # Latin dotted lowercase letters
# for i in 0x1E60:0x1E9B
#     push!(letters, Char(i))
# end

# # Latin barred uppercase letters
# for i in 0x1E80:0x1EF9
#     push!(letters, Char(i))
# end

# # Latin barred lowercase letters
# for i in 0x1EFB:0x1F00
#     push!(letters, Char(i))
# end

# # Greek uppercase letters
# for i in 0x391:0x3A9
#     push!(letters, Char(i))
# end

# # Greek lowercase letters
# for i in 0x3B1:0x3C9
#     push!(letters, Char(i))
# end

# # Greek dotted uppercase letters
# for i in 0x1F00:0x1F08
#     push!(letters, Char(i))
# end

# # Greek dotted lowercase letters
# for i in 0x1F10:0x1F16
#     push!(letters, Char(i))
# end

# # Greek barred uppercase letters
# for i in 0x1F08:0x1F0F
#     push!(letters, Char(i))
# end

# # Greek barred lowercase letters
# for i in 0x1F18:0x1F1E
#     push!(letters, Char(i))
# end

# # # Hebrew letters
# # for i in 0x5D0:0x5EA
# #     push!(letters, Char(i))
# # end

# # # Hebrew dotted letters
# # for i in 0xFB1D:0xFB36
# #     push!(letters, Char(i))
# # end

# # # Hebrew barred letters
# # for i in 0xFB4D:0xFB4F
# #     push!(letters, Char(i))
# # end

# # # Hebrew final letters
# # for i in 0x5DA:0x5EA
# #     push!(letters, Char(i))
# # end

# # # Hebrew final dotted letters
# # for i in 0xFB3E:0xFB3F
# #     push!(letters, Char(i))
# # end

# # # Hebrew final barred letters
# # push!(letters, Char(0xFB40))
# # push!(letters, Char(0xFB41))

# # Display the set of letters
###println(contains(letters,̅θ)
###println(letters)