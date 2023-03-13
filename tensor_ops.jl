include("tensor_term.jl")
include("kronecker.jl")
include("epsilon_tensor.jl")
include("metric.jl")



##############################################################################################################################
##############################################################################################################################

### metric ops

function pop_metric(T::TensorSuperType) # makes all indices on the object go to the default position according to the index set
    tensor = 1
    new_index_pattern = T.indices

    for pos in eachindex(T.indices)
        index = T.indices[pos]
        if !in_default_pos(index)
            spare = get_unused_indices(T, index[1].index_set)
            # println(spare)
            dummy = pop!(spare) => index[1].index_set.default_position == upper ? lower : upper
            if dummy[1].grading==0
                metric = Metric("η", [index, dummy])
            else
                metric = EpsilonTensor([index, dummy])
            end
            
            tensor = tensor * metric
            # println(new_index_pattern[pos])
            # println(dummy => index[1].index_set.default_position)
            new_index_pattern[pos] = dummy[1] => index[1].index_set.default_position
        end
    end
    T.indices = new_index_pattern
    # println(tensor)
    return tensor * T

end


function contract_metrics(T::TensorExpression, typ::Type=Metric)
    # println("running function")

    for k in eachindex(T.summands)
        # println(k)
        # println(T.summands)
        summand = T.summands[k]
        metrics = filter(y -> typeof(y) <: typ, summand[2].terms)
        non_metrics = filter(y -> !(typeof(y) <: typ), summand[2].terms)
        uncontr_mets = Vector{TensorSuperType}()
        # println(metrics) 
        sign = 1 #for antisym metrics (spinor)

        for metric in metrics #for each in the found metrics
            flag = false # terminate when metric is contracted
            index1 = metric.indices[1][1]
            index2 = metric.indices[2][1]
            for tens in non_metrics #try to contract with each other tensor
                # println(tens)
                if typeof(metric)<:EpsilonTensorSuperType && typeof(tens)<:MetricSuperType
                    # println("cant contract kronecker")
                    continue
                end
                # println("not a kronecker")
                if flag # terminat if contracted
                    break
                end
                for indpairnum in eachindex(tens.indices)
                    if flag
                        break
                    end
                    if tens.indices[indpairnum][1] == index2 #tensor has an index matching the 2nd index of the metric
                        #some complicated logic to make sure indices go to the right position

                        # this should let you contract Kronckers and metrics in a consistent way
                        met_pos1 = metric.indices[1][2].position # identifies when the metric is metric or kronecker etc
                        met_pos2 = metric.indices[2][2].position
                        # println(met_pos1)
                        # println(met_pos2)
                        sumind = (met_pos1 + met_pos2) % 2 # if = 0, both upper/lower, 1 if 1 upper 1 lower
                        if sumind == 0

                            if !(typeof(metric) <: SymMetricSuperType)

                                if (tens.indices[indpairnum][2].position + met_pos1) % 2 == 0 # 2nd metric and tensor pos same, raise and lower with sign change
                                    sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) )) #use 
                                    # println("flipping")
                                    # println(tens.indices[indpairnum])
                                    if tens.indices[indpairnum][2] == upper
                                        # println("sending upper to lower")
                                        tens.indices[indpairnum] = metric.indices[1][1] => lower
                                        # println(tens.indices)
                                    else
                                        # println("sending lower to upper")
                                        tens.indices[indpairnum] = metric.indices[1][1] => upper
                                        # println(tens.indices)
                                    end
                                else
                                tens.indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                                end
                            else

                                tens.indices[indpairnum] = metric.indices[1][1] => metric.indices[1][2] # raise/lower index
                            end
                        else
                            # println("doing this one")
                            # println(tens.indices[indpairnum][2])
                            if (tens.indices[indpairnum][2].position + met_pos2) % 2 == 0 # if tensor index position is same as 2nd metric pos 
                                tens.indices[indpairnum] = metric.indices[1][1] => tens.indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)
                            else
                                tens.indices[indpairnum] = metric.indices[1][1] => tens.indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                            end
                        end
                        #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                        flag = true #this term has been contracted, break
                    elseif tens.indices[indpairnum][1] == index1
                        #some complicated logic to make sure indices go to the right position
                        #if contracting with antisym, pick up a minus
                        # println("symmetric?")
                        # println(typeof(metric)<:SymMetricSuperType)
                        # println(sign)
                        met_pos1 = metric.indices[1][2].position
                        met_pos2 = metric.indices[2][2].position
                        # println(met_pos1)
                        # println(met_pos2)
                        sumind = (met_pos1 + met_pos2) % 2
                        if sumind == 0
                            # println(" yo")
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
                            # println("doing this one")
                            # println(tens.indices[indpairnum][2])
                            if (tens.indices[indpairnum][2].position + met_pos1) % 2 == 0
                                # println(" here")
                                #for anti sym, we need to keep sign but change pos 
                                #sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                # if !(typeof(metric)<:SymMetricSuperType)
                                #     println("we antisym")
                                #     if tens.indices[indpairnum][2].position ==upper
                                #         tens.indices[indpairnum]=metric.indices[2][1] => lower
                                #     else
                                #         tens.indices[indpairnum]=metric.indices[2][1] => upper
                                #     end
                                # else 
                                tens.indices[indpairnum] = metric.indices[2][1] => tens.indices[indpairnum][2] # keep index in same pos ( δ_j^i A^i= A^j for example)


                            else #opposite pos, grading etc
                                # sign=sign*((-1)^(!(typeof(metric)<:SymMetricSuperType) ))
                                # println("got here")
                                # if !(typeof(metric)<:SymMetricSuperType) #for anti symm, change positions
                                #     println("we antisym too")
                                #     if tens.indices[indpairnum][2].position ==upper
                                #         tens.indices[indpairnum]=metric.indices[2][1] => lower
                                #     else
                                #         tens.indices[indpairnum]=metric.indices[2][1] => upper
                                #     end
                                # else
                                tens.indices[indpairnum] = metric.indices[2][1] => tens.indices[indpairnum][2] # if opposite positions, keep tensor pos ( δ^j_i A^i= A^j for example)
                                # end
                            end
                        end
                        #.indices[indpairnum]=metric.indices[1][1] => metric.indices[1][2]
                        flag = true #this term has been contracted, break
                    end
                end

            end
            if !flag
                push!(uncontr_mets, metric)
            end

        end
        T.summands[k][2].terms = vcat(uncontr_mets, non_metrics)
        # println(T.summands[k][1])
        # println(T.summands[k][1]*sign)
        T.summands[k] = T.summands[k][1] * sign => T.summands[k][2]
    end
    return merge_sum(T)

end

# function contract_metrics(T::Tensor) term
#     contract

# end


##############################################################################################################################
##############################################################################################################################

### KroneckerDelta


##############################################################################################################################

#elim kronecker
function eliminate_kronecker(T::TensorTerm)
    eliminate_kronecker(TensorExpression(T))
end


function eliminate_kronecker(T::TensorExpression)
    contract_metrics(T, KroneckerSuperType)

end




##############################################################################################################################
##############################################################################################################################

### Eps Tensor


##############################################################################################################################

#elim eps
function eliminate_epsilon(T::TensorTerm)
    eliminate_epsilon(TensorExpression(T))
end


function eliminate_epsilon(T::TensorExpression)
    contact_eps(contract_metrics(T, EpsilonTensorSuperType))

end

function contact_eps(T::TensorExpression)
    return T
end



##############################################################################################################################
##############################################################################################################################

### Mathematical ops


##############################################################################################################################

#addition

function +(A::TensorSuperType, B::TensorSuperType)
    if get_free_indices(A) == get_free_indices(B)
        return TensorExpression([A, B])
    else
        throw(FreeIndexMismatch())
    end

end

function +(A::TensorSuperType, B::TensorExpression)
    if  isZeroTens(B) || get_free_indices(A) == get_free_indices(B)
        return TensorExpression(TensorExpression(A), B)
    else
        throw(FreeIndexMismatch())
    end
end

function +(A::TensorExpression, B::TensorSuperType)
    if  isZeroTens(A) || get_free_indices(A) == get_free_indices(B)
        return TensorExpression(A, TensorExpression(B))
    else
        throw(FreeIndexMismatch())
    end
end

function +(A::TensorExpression, B::TensorExpression)
    # println(A)
    # println(B)
    # println(isZeroTens(A))
    # println(isZeroTens(B))
    if isZeroTens(A) || isZeroTens(B) || get_free_indices(A) == get_free_indices(B)
        return TensorExpression(A, B)
    else
        throw(FreeIndexMismatch())
    end
end


function +(A::TensorTerm, B::TensorExpression)

    if  isZeroTens(B) || get_free_indices(A) == get_free_indices(B)
        return TensorExpression(TensorExpression(A), B)
    else
        throw(FreeIndexMismatch())
    end
end

+(A::TensorExpression, B::TensorTerm)= +(B,A)

function +(A::TensorTerm, B::TensorTerm)

    if get_free_indices(A) == get_free_indices(B)
        return TensorExpression(A)+TensorExpression(B)
    else
        throw(FreeIndexMismatch())
    end
end

+(A::TensorSuperType, B::TensorTerm)= +(B,A)

function +(A::TensorTerm, B::TensorSuperType)

    if get_free_indices(A) == get_free_indices(B)
        return TensorExpression(A)+TensorExpression(B)
    else
        # println(get_free_indices(A))
        # println(get_free_indices(B))
        throw(FreeIndexMismatch())
    end
end


# function +(N::SymorNum, T::TensorExpression)



#     if get_free_indices(A) == get_free_indices(B)
#         return TensorExpression(TensorExpression(A), B)
#     else
#         throw(FreeIndexMismatch())
#     end
# end

##############################################################################################################################

#subtraction

-(A::TensorSuperType, B::TensorSuperType) = -1 * B + A
-(A::TensorSuperType, B::TensorExpression) = A + (-1 * B)
-(A::TensorExpression, B::TensorSuperType) = A + (-1 * B)
-(A::TensorExpression, B::TensorExpression) = A + (-1 * B)

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

#multiplication

*(N::SymorNum, B::TensorSuperType) = TensorExpression(N, copy(B))
*(N::SymorNum, B::TensorTerm) = TensorExpression(N => B)
*(B::TensorSuperType, N::SymorNum) = TensorExpression(N, copy(B))
*(B::TensorTerm, N::SymorNum) = TensorExpression(N => B)

function *(B::TensorExpression, N::SymorNum)
    return update_indices(TensorExpression([convert(SymorNum, N) * summand[1] => summand[2] for summand in B.summands]))
end

*(N::SymorNum, B::TensorExpression) = *(B, convert(SymorNum, N))


function *(T::TensorSuperType, TEx::TensorExpression)
    return update_indices(TensorExpression([sum[1] => T * sum[2] for sum in TEx.summands]))
end


function *(T::TensorTerm, TEx::TensorExpression)
    return update_indices(TensorExpression([sum[1] => T * sum[2] for sum in TEx.summands]))
end

function *(TEx::TensorExpression,T::TensorTerm)
    return update_indices(TensorExpression([sum[1] =>  sum[2] *T for sum in TEx.summands]))
end


function *(TEx::TensorExpression, T::TensorSuperType)

    for summandind in eachindex(TEx.summands)
        # println(TEx.summands[summandind])
        TEx.summands[summandind]=TEx.summands[summandind][1] => TEx.summands[summandind][2]*T
        # println(TEx.summands[summandind])
    end
    #ex=(sum(sum[1] => sum[2] * T for sum in TEx.summands))
    # println(TEx)
    return update_indices(TEx)
end


# *(TEx::TensorExpression,T::TensorSuperType)=*(T,TEx)

function *(T::TensorSuperType, TEx::TensorTerm)
    return update_indices(TensorTerm(vcat(copy(T), TEx.terms)))
end

function *(T::TensorTerm, TEx::TensorTerm)
    return update_indices(TensorTerm(vcat(T.terms, TEx.terms)))
end

function *(T1::TensorSuperType, T2::TensorSuperType)
    return update_indices(TensorTerm([copy(T1), copy(T2)]))
end

function *(TEx::TensorTerm, T::TensorSuperType)
    return update_indices(TensorTerm(vcat(TEx.terms, copy(T))))
end

##############################################################################################################################

#smart times, basically appends a tensor so that all of its indices are free. useful if you want to mutliply something without doing contractions

# function ×(T1::TensorTerm, T2::TensorSuperType)


# end
function ×( T2::SymorNum, N::Union{TensorExpression,TensorTerm,TensorSuperType} ) 
    N × T2
end

function ×( T2::Union{TensorExpression,TensorTerm,TensorSuperType}, N::SymorNum) 
    # T2*N 
    throw(OperatorPrecedenceError())
end


function ×( T2::TensorExpression, T1::TensorExpression) 

    vec=[]
    for summand2 in T2.summands
        for summand1 in T1.summands
            push!(vec, summand2[1]*summand1[1] => summand2[2] × summand1[2] )
        end
    end
    T2.summands=vec
    return T2
end


function ×( T2::TensorExpression, T1::Union{TensorSuperType,TensorTerm}) 

    (-1)^((grading(T1)*grading(T2))%2)*(T1 × T2)
end

function ×( T2::Union{TensorSuperType,TensorTerm}, T1::TensorExpression) 
   new_terms=TensorExpression()
   T1=rename_dummies(T1)
#    println(T1)
#    println(T2)
#    println("^in func")
    for summandin in eachindex(T1.summands)
        t=(T2 × T1.summands[summandin][2])
        new_terms+= (T1.summands[summandin][1]* t)
        # println(t)
        # println("poo poo")
        # println(T1.summands[summandin])
        # T1.summands[summandin]=T1.summands[summandin][1] => t
    end

    return new_terms

end

function ×(T1::TensorTerm, T2::TensorTerm) 
    
    for term in T2.terms
        T1=T1 × term
    end
    return T1
end



function ×( T2::TensorSuperType, T1::TensorTerm) 
    return (-1)^((grading(T1)*grading(T2))%2)*TensorExpression(T1 × T2)

end

function ×(T1::Union{TensorTerm, TensorSuperType}, T2::TensorSuperType) 
    # println("crash func")
    # println(T1)
    # println(typeof(T1)) #some different behaviour between tensor and tensor term. tensor works..
    # println(T2)
    T2=copy(T2)
    t1inds=OrderedSet{Index}()
    t2inds=OrderedSet{Pair{Integer,Index}}() #all indices, numbered in order
    # t2ok=Set{Index}()
    for inds in T1.indices
        push!(t1inds, inds[1])
    end
    for inds in eachindex(T2.indices)
        if  T2.indices[inds][1] in t1inds #if theres a repeated indice
            push!(t2inds, inds => T2.indices[inds][1] ) # which index to replace
        else
            # push!(t2ok,T2.indices[inds][1])# index doesnt need to be replaced
        end
    end
    # we now have all the indices we need to replace and where they are on the tensor
    sets=get_index_sets(T2) #Set{IndexSet}()
    avail=Dict{IndexSet, OrderedSet{Index}}() #pairs of index sets, with the unused indices
    for set in sets
        a=OrderedSet{Index}(intersect(get_unused_indices(T1, set), get_unused_indices(T2, set)))
        push!(avail, set=> a) #ones we can replace with
    end
    # println(avail)
    # println(T2.indices)
    dummie_matchup=Dict{Index,Index}()

    t2dums=get_dummies(T2)
    
    for ind in t2inds # for each we need to replace (vec pos => index)
        # println(avail[ind[2].index_set])
        # println(ind)

        if ind[2] in t2dums
           
            if ind[2] in keys(dummie_matchup)
                nind=pop!(dummie_matchup,ind[2])
            # println(nind)
            T2.indices[ind[1]]= nind=> T2.indices[ind[1]][2]
            else
                if isempty(avail[ind[2].index_set])
                    throw(OutOfIndices) 
                 end
                nind=popfirst!(avail[ind[2].index_set])
                push!(dummie_matchup, ind[2] => nind)
                T2.indices[ind[1]]= nind=> T2.indices[ind[1]][2]
            end
        else
            if isempty(avail[ind[2].index_set])
               throw(OutOfIndices) 
            end
            nind=popfirst!(avail[ind[2].index_set])
            # println(nind)
            T2.indices[ind[1]]= nind=> T2.indices[ind[1]][2] #should replace
        end
        
    end
    # println("before crash") #same dummies
    # println(T1)
    # println(T2)
    return T1*T2
end


##############################################################################################################################

#division

/(B::TensorSuperType, N::SymorNum) = *(1 / N, B)
/(B::TensorExpression, N::SymorNum) = *(1 / N, B)

##############################################################################################################################

#equality checking

function ==(a::TensorTerm, b::TensorTerm)
    return a.terms == b.terms
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

function ==(a::TensorSuperType, b::TensorSuperType)
    if a.name == b.name && a.indices == b.indices
        return true
    else
        return false
    end
end



##############################################################################################################################
##############################################################################################################################

#indices

##############################################################################################################################


#update indices

function update_indices(T::TensorTerm)
    inds=[]
    for t in T.terms
        for i in t.indices
            push!(inds, i)
        end
    end
    # println(inds)
    T.indices=inds
    ordered_counts = OrderedDict((pair[1], 0) for pair in T.indices)
    # println(ordered_counts)
    for pair in T.indices
        ordered_counts[pair[1]] += 1
    end
    # println(ordered_counts)
    duplicates = filter(kv -> kv[2] > 2, ordered_counts)
    # println(duplicates)
    if !isempty(duplicates)
        throw(OverusedIndices())
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
    # println(Set(i[1] for i in T.indices))
    # println(IS.indices)
    return setdiff(IS.indices, Set(i[1] for i in T.indices))
end

function get_unused_indices(T::TensorTerm, IS::IndexSet)
    # println(Set(i[1] for i in T.indices))
    # println(IS.indices)
    return setdiff(IS.indices, Set(i[1] for i in T.indices))
end


function get_unused_indices(T::TensorExpression, IS::IndexSet)
    # println(Set(i[1] for i in T.indices))
    # println(IS.indices)
    return setdiff(IS.indices, Set(i[1] for i in T.terms[1].indices))
end

#finding free indices

function get_free_indices(T::TensorSuperType)
    return Set(setdiff(T.indices, get_dummies_positional(T)))
end

function get_free_indices(T::TensorTerm)
    # println("doing this one")
    # println(T.indices)
    d=get_dummies_positional(T)
    # # println(d)
    # # println("free")
    # println(Set(intersect(T.indices,d)))
    return Set(setdiff(T.indices,d))
end

function get_free_indices(T::TensorExpression)
    return Set(setdiff(T.summands[1][2].indices, get_dummies_positional(T.summands[1][2])))#all tensors should have matching free indices in a sum already
end

##############################################################################################################################

#finding dummie indices

function get_dummies_positional(T::TensorSuperType)
    dum::OrderedSet{Pair{Index,IndexPosition}} = Set()
    for d in get_dummies(T)
        push!(dum, d => upper)
        push!(dum, d => lower)
    end
    return dum
end

function get_dummies_positional(T::TensorTerm)
    dum::OrderedSet{Pair{Index,IndexPosition}} = Set()
    for d in get_dummies(T)
        push!(dum, d => upper)
        push!(dum, d => lower)
    end
    return dum
end


function get_index_sets(T::Union{TensorSuperType,TensorTerm})
    temp = [ind[1] for ind in T.indices]
    out=OrderedSet{IndexSet}()
    for t in temp
        push!(out,t.index_set)
    end
    return out
end


#returns a set of dummie indices that are used in summation
function get_dummies(T::Union{TensorSuperType,TensorTerm})

    temp = [ind[1] for ind in T.indices]
    # println(temp)
    ordered_counts = OrderedDict((pair, 0) for pair in temp)
    for pair in temp
        ordered_counts[pair] += 1
    end
    duplicates = filter(kv -> kv[2] == 2, ordered_counts)
    # for ind in T.indices
    #     if cind 
    return collect(keys(duplicates))
end

# function get_dummies(T::TensorTerm)
#     temp = [ind[1] for ind in T.indices]
#     # println(temp)
#     ordered_counts = OrderedDict((pair, 0) for pair in temp)
#     for pair in temp
#         ordered_counts[pair] += 1
#     end
#     duplicates = filter(kv -> kv[2] == 2, ordered_counts)
#     # for ind in T.indices
#     #     if cind 
#     return collect(keys(duplicates))
#     # temp = [ind[1] for ind in T.indices]
#     # # println(temp)
#     # # println(size(temp))
#     # # if size(temp)[1] == 8
#     # #     println(temp[2])
#     # #     println(temp[8])
#     # # println(temp[2]==temp[8])
#     # # end
#     # # name_pairs = OrderedDict((pair.name, pair) for pair in temp)
#     # ordered_counts = OrderedDict((pair.name, 0) for pair in temp)
#     # # println("name pairs")
#     # # println(name_pairs)
#     # for pair in temp
#     #     ordered_counts[pair.name] += 1
#     # end
#     # # println(ordered_counts)
#     # duplicates = filter(kv -> kv[2] == 2, ordered_counts)
#     # # #  println(duplicates)
#     # # # for ind in T.indices
#     # # #     if cind 
#     # # # println("has dummies")
#     # # names=collect(keys(duplicates))
#     # # println(names) # these are the dummies
#     # # inds=[]
#     # # for name in names # had to do some complicated stuff since sometimes dummies arent recognised, i think due to the way structs are compared?
#     # #     push!(inds, name_pairs[name])
#     # # end
#     # # # filter(kv -> names[kv[1]] != Nothing, name_pairs) #?????????? want to get all indices from the name
#     # # println(collect(values(inds))) #dummie indices
#     # return collect(keys(duplicates))
# end





function rename_dummies(T::TensorSuperType)
    dummies = get_dummies(T) #returns all the dummies attached to the tensor
    all_indices::OrderedSet{Index} = OrderedSet([i[1] for i in T.indices]) #set of all the indices attached to the tensor
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
        # println(replacements)
        for p in replacements
            T.indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in T.indices]
        end
        #     for (x, y) in T.indices
        #         (x,y)=(x == p[1] ? p[2] : x), y for p in replacements 
        #    # T.indices = [Pair((x == p[1] ? p[2] : x), y) for p in replacements for (x, y) in T.indices ]
    end
end



##############################################################################################################################


#rename dummies function

function rename_dummies(T::TensorTerm)
    dummies = get_dummies(T) #returns all the dummies attached to the tensor
    # println("dummies")
    # println(dummies)
    all_indices::OrderedSet{Index} = OrderedSet([i[1] for t in T.terms for i in t.indices]) #set of all the indices attached to the tensor
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
        # println(replacements)
        for p in replacements
            for t in T.terms
                t.indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in t.indices]
            end
        end
        #     for (x, y) in T.indices
        #         (x,y)=(x == p[1] ? p[2] : x), y for p in replacements 
        #    # T.indices = [Pair((x == p[1] ? p[2] : x), y) for p in replacements for (x, y) in T.indices ]
    end
    update_indices(T)
    return T
end

function rename_dummies(T::TensorExpression)
    for term in T.summands
        term= term[1] => rename_dummies(term[2])
    end
    update_indices(T)
    return T
end



##############################################################################################################################
##############################################################################################################################

#merging terms with same tensor content

function merge_sum(TEx::TensorExpression)
    out = filter!(x -> x[1] != 0, TEx.summands)
    TEx.summands = out
    new_merged = []
    for summand in TEx.summands
        # println("doing")
        found = false

        for i in eachindex(new_merged)
            if summand[2] == new_merged[i][2]
                # println("got here")
                # println(summand[1])
                # println(new_merged[i][1])
                # println(summand[2])
                new_merged[i] = (new_merged[i][1] + summand[1] => summand[2])
                #  println("sum is")
                #  println(new_merged[i])
                #  println(new_merged) #not editting in place apparaently....
                found = true
                break
            end
        end
        if found
            break
        end
        push!(new_merged, summand)
    end

    TEx.summands = new_merged
    TEx.summands = filter!(x -> x[1] != 0, TEx.summands)
    return isempty(TEx.summands) ? 0 : TEx
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
    print(io,tens2string(T))
end

function Base.show(io::IO, T::TensorTerm)
    for tens in T.terms
        print(io, tens)
    end
end

function Base.show(io::IO, T::TensorExpression)
    it = Iterators.Stateful(T.summands)
    for tens in it
        print(io, tens[1] == 1 ? "" : tens[1] == -1 ? "-" : tens[1])
        print(io, tens[1] == 1 ? "" : tens[1] == -1 ? "" : "*")
        print(io, tens[2])
        next = peek(it)
        if next !== nothing
            if next[1] >= 0
                print(io, "+")
            end
        end
    end
end



##############################################################################################################################
##############################################################################################################################

# sorting

function sort_exp(TEx::TensorExpression, Order::Vector{T} where {T<:TensorSuperType})
    for summand in TEx.summands

        

    end
    return TEx
end



##############################################################################################################################
##############################################################################################################################

# replacements and substitutions


function substitute(TEx::TensorExpression, pattern::String, repeat::Bool=true)

    #first step: parse terms to a string, that points back to position in term



end


function term2stringpos(Term::TensorTerm)
    out=Vector{Pair{String, Integer}}
    for tind in eachindex(Term.terms)
        push!(out, tens2string(Term.terms[tind]) => tind)
    end
   return out
end


function tens2string(T::TensorSuperType)
    out=T.name
    if !isempty(T.indices)
        lowers = ""
        uppers = ""
        lowersempt::Bool = 1
        uppersempt::Bool = 1
        for (ind, pos) in T.indices
            if pos.position == 0
                uppers = uppers * ind.name
                uppersempt = 0
                lowers = lowers * "*"
            else
                uppers = uppers * "*"
                lowersempt = 0
                lowers = lowers * ind.name
            end
        end
        if !uppersempt
            out*="^{" * uppers * "}"
        end
        if !lowersempt
            out*="_{" * lowers * "}"
        end
    end


end




function replace(T::TensorTerm, old_T::TensorSuperType, new_T::TensorSuperType, repeat::Bool=true)
    replace(TensorExpression(T), old_T, new_T, repeat)
end

function replace(TEx::TensorExpression, old_T::TensorSuperType, new_T::TensorSuperType, repeat::Bool=true)
    for k in eachindex(TEx.summands)
        # println(k)
        # println(T.summands)
        summand = TEx.summands[k]
        matches = []
        for tind in eachindex(summand[2].terms)
            if summand[2].terms[tind]== old_T
            push!(matches, tind)
            end
        end
        # matches = filter(y -> y==old_T, summand[2].terms)
        # non_metrics = filter(y -> !(typeof(y) <: typ), summand[2].terms)
        # uncontr_mets = Vector{TensorSuperType}()
        # # println(metrics) 
        # sign = 1 #for antisym metrics (spinor)
        if repeat
            for match in matches
                summand[2].terms[match]=new_T
            end
        else
            if !isempty(matches)
                summand[2].terms[1]=new_T
            end
        end 
    end
    return update_indices(TEx)
end






##############################################################################################################################
##############################################################################################################################

# copying


function copy(T::TensorTerm)
    return typeof(T)([copy(t) for t in T.terms])
end

function copy(T::TensorSuperType)
    return typeof(T)(T.name,copy(T.indices),T.function_of)
end
function copy(V::Vector{Pair{Index,IndexPosition}})
    return [pair[1] => pair[2] for pair in V]
end
function copy(T::MetricSuperType)
    return typeof(T)(T.name,copy(T.indices))
end
function copy(T::KroneckerSuperType)
    return typeof(T)(copy(T.indices))
end
function copy(T::EpsilonTensorSuperType)
    return typeof(T)(copy(T.indices))
end

##############################################################################################################################
##############################################################################################################################

# utils


function isZeroTens(T::TensorExpression)
    if length(T.summands)==1 && T.summands[1][1]==0
        return true
    else
        return false
    end
end

# function isZeroTens(T::TEx)
#     if length(T.summands)==1 && T.summands[1][1]==0
#         return true
#     else
#         return false
#     end
# end



function first_element(arr, elem1, elem2)
    idx1 = findfirst(x -> x == elem1, arr)
    idx2 = findfirst(x -> x == elem2, arr)
    if idx1 == nothing || (idx2 != nothing && idx2 < idx1)
        return elem2
    else
        return elem1
    end
end
