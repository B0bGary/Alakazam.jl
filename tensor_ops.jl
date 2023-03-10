include("tensor_term.jl")


##############################################################################################################################
##############################################################################################################################

### Mathematical ops


##############################################################################################################################

#addition

function +(A::Tensor, B::Tensor)
    if get_free_indices(A)==get_free_indices(B)
        return TensorExpression([A,B])
    else
        throw(FreeIndexMismatch())
    end
  
end

function+(A::Tensor, B::TensorExpression)
    if get_free_indices(A)==get_free_indices(B)
        return TensorExpression(TensorExpression(A),B)
    else
        throw(FreeIndexMismatch())
    end
end

function +(A::TensorExpression, B::Tensor)
    if get_free_indices(A)==get_free_indices(B)
        return TensorExpression(A,TensorExpression(B))
    else
        throw(FreeIndexMismatch())
    end
end

function +(A::TensorExpression, B::TensorExpression)
    
    if get_free_indices(A)==get_free_indices(B)
        return TensorExpression(A,B)
    else
        throw(FreeIndexMismatch())
    end
end

##############################################################################################################################

#multiplication

*(N::SymorNum, B::Tensor)= TensorExpression([N => B])
*(B::Tensor, N::SymorNum)= TensorExpression([N => B])



##############################################################################################################################

#equality checking

function ==(a::TensorTerm, b::TensorTerm)
    return a.terms==b.terms
end

function ==(a::Vector{Tensor}, b::Vector{Tensor})
    if length(a)!=length(b)
        return false
    end
    for i=eachindex(a)
        if a[i]!=b[i]
            return false
        end
    end
    return true
end

function ==(a::Tensor, b::Tensor)
    if a.name==b.name && a.indices==b.indices
        return true
    else 
        return false
    end
end



##############################################################################################################################
##############################################################################################################################

#indices

##############################################################################################################################


#finding free indices

function get_free_indices(T::Tensor)
    return setdiff(T.indices,get_dummies_positional(T))
end

function get_free_indices(T::TensorTerm)
    return setdiff(T.indices,get_dummies_positional(T))
end

function get_free_indices(T::TensorExpression)
    return setdiff(T.summands[1].indices,get_dummies_positional(T.summands[1]))#all tensors should have matching free indices in a sum already
end

##############################################################################################################################

#finding dummie indices

function get_dummies_positional(T::Tensor)
    dum::Set{Pair{Index, IndexPosition}}=Set()
    for d in get_dummies(T)
        push!(dum, d=>upper)
        push!(dum, d=>lower)
    end
    return dum
end

function get_dummies_positional(T::TensorTerm)
    dum::Set{Pair{Index, IndexPosition}}=Set()
    for d in get_dummies(T)
        push!(dum, d=>upper)
        push!(dum, d=>lower)
    end
    return dum
end



 #returns a set of dummie indices that are used in summation
 function get_dummies(T::Tensor)

    temp=[ind[1] for ind in T.indices]
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

function get_dummies(T::TensorTerm)

    temp=[ind[1] for ind in T.indices]
   # println(temp)
   ordered_counts = OrderedDict((pair, 0) for pair in temp)
   for pair in temp
       ordered_counts[pair] += 1
   end
    duplicates = filter(kv -> kv[2] == 2, ordered_counts)
    println(duplicates)
    # for ind in T.indices
    #     if cind 
    return collect(keys(duplicates))
end





function rename_dummies(T::Tensor)
    dummies=get_dummies(T) #returns all the dummies attached to the tensor
    all_indices::OrderedSet{Index}=OrderedSet([i[1] for i in T.indices]) #set of all the indices attached to the tensor
    index_sets::OrderedSet{IndexSet}=OrderedSet([d.index_set for d in dummies]) #set of all the index_sets associated with the dummies
    for set in index_sets #for each type of dummie index

        replacements::Vector{Pair}=[] #dummie replacements to make
        inset=OrderedSet(filter(a -> a.index_set==set ,all_indices)) # all indices on our tensor that are of the specified index type
        
        remain = OrderedSet(filter(x -> !(x in inset), OrderedSet(set.indices)))
       # remain=setdiff(OrderedSet(set.indices),inset) # possible indices to replace our dummies that aren't already used on the tensor
        dumset=filter(d -> d.index_set==set ,dummies)  #dummie indices to replace of the appropriate type
        for d in dumset #attempt to replace each dummie index
            for r in remain #check candidates
                if first_element(set.indices, d, r)!=d # if the dummie doesn't come first lexographically (or in your specified order)
                    push!(replacements, d => r)
                    delete!(remain, r)
                    break
                end
            end
        end
        println(replacements)
        for p in replacements
            T.indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in T.indices ]
        end
    #     for (x, y) in T.indices
    #         (x,y)=(x == p[1] ? p[2] : x), y for p in replacements 
    #    # T.indices = [Pair((x == p[1] ? p[2] : x), y) for p in replacements for (x, y) in T.indices ]
    end
end



##############################################################################################################################


#rename dummies function

function rename_dummies(T::TensorTerm)
    dummies=get_dummies(T) #returns all the dummies attached to the tensor
    # println("dummies")
    # println(dummies)
    all_indices::OrderedSet{Index}=OrderedSet([i[1] for t in T.terms for i in t.indices ]) #set of all the indices attached to the tensor
    index_sets::OrderedSet{IndexSet}=OrderedSet([d.index_set for d in dummies]) #set of all the index_sets associated with the dummies
    for set in index_sets #for each type of dummie index

        replacements::Vector{Pair}=[] #dummie replacements to make
        inset=OrderedSet(filter(a -> a.index_set==set ,all_indices)) # all indices on our tensor that are of the specified index type
        remain = OrderedSet(filter(x -> !(x in inset), OrderedSet(set.indices)))
        #remain=setdiff(OrderedSet(set.indices),inset) # possible indices to replace our dummies that aren't already used on the tensor
        dumset=filter(d -> d.index_set==set ,dummies)  #dummie indices to replace of the appropriate type
        for d in dumset #attempt to replace each dummie index
            for r in remain #check candidates
                if first_element(set.indices, d, r)!=d # if the dummie doesn't come first lexographically (or in your specified order)
                    push!(replacements, d => r)
                    delete!(remain, r)
                    break
                end
            end
        end
        println(replacements)
        for p in replacements
            for t in T.terms
            t.indices = [Pair((x == p[1] ? p[2] : x), y) for (x, y) in t.indices ]
            end
        end
    #     for (x, y) in T.indices
    #         (x,y)=(x == p[1] ? p[2] : x), y for p in replacements 
    #    # T.indices = [Pair((x == p[1] ? p[2] : x), y) for p in replacements for (x, y) in T.indices ]
    end
end



##############################################################################################################################
##############################################################################################################################

#merging terms with same tensor content

function merge_sum(TEx::TensorExpression)
    out = filter!(x -> x[1] != 0, TEx.summands)
    TEx.summands=out
    new_merged=[]
    for summand in TEx.summands
        # println("doing")
        found=false
        
        for i in eachindex(new_merged)
            if summand[2]==new_merged[i][2]
                # println("got here")
                # println(summand[1])
                # println(new_merged[i][1])
                # println(summand[2])
                new_merged[i]=(new_merged[i][1]+summand[1] => summand[2])
                # println("sum is")
                # println(new_merged[i])
                # println(new_merged) #not editting in place apparaently....
                found=true
                break
            end
        end
        if found
            break
        end
        push!(new_merged,summand)    
    end

    TEx.summands=new_merged
    #print(TEx.summands)
    return TEx
end


##############################################################################################################################
##############################################################################################################################

#gradings


function grading(T::Tensor)
    grade=0
    for  (ind, pos) in T.indices
        grade=grade+ind.grading
    end
    return grade%2
 end

 function grading(T::TensorTerm)
    grade=0
    for  tens in T.terms
        grade=grade+grading(tens)
    end
    return grade%2
end

function grading(T::TensorExpression)
    grade=0
    for  tens in T.summands[1].terms # should be same grading for everything in sum since they have matching free indices
        grade=grade+grading(tens)
    end
    return grade%2
end


##############################################################################################################################
##############################################################################################################################

# printing and display

function Base.show(io::IO, T::Tensor)
    print(io, T.name)
    if !isempty(T.indices)
        lowers=""
        uppers=""
        lowersempt::Bool=1
        uppersempt::Bool=1
        for (ind, pos) in T.indices
            if pos.position==0
                uppers=uppers*ind.name
                uppersempt=0
                lowers=lowers * "*"
            else
                uppers=uppers * "*"
                lowersempt=0
                lowers=lowers * ind.name 
            end
        end
        if !uppersempt
            print(io, "^{"*uppers*"}")
        end
        if !lowersempt
            print(io, "_{"*lowers*"}")
        end
    end
end

 function Base.show(io::IO, T::TensorTerm)
    for tens in T.terms
        print(io, tens)
    end
 end

 function Base.show(io::IO, T::TensorExpression)
    it= Iterators.Stateful(T.summands)
    for tens in it
        print(io, tens[1]==1 ? "" : tens[1]==-1 ? "-" : tens[1])
        print(io, tens[1]==1 ? "" : tens[1]==-1 ? "" : "*")
        print(io,tens[2])
        if peek(it)!=nothing
            print(io, "+")
        end
    end
 end



##############################################################################################################################
##############################################################################################################################

# utils

 
function first_element(arr, elem1, elem2)
    idx1 = findfirst(x -> x == elem1, arr)
    idx2 = findfirst(x -> x == elem2, arr)
    if idx1 == nothing || (idx2 != nothing && idx2 < idx1)
        return elem2
    else
        return elem1
    end
end
