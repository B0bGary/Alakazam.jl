include("indices.jl")




mutable struct Tensor
    name::String
    indices::Vector{Pair{Index, IndexPosition}}

    function Tensor(name::String, indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}())
        new(name, indices)
    end
 end


 
mutable struct TensorTerm
    terms::Vector{Tensor} 
    indices::Vector{Pair{Index, IndexPosition}}
    function TensorTerm(terms::Vector{Tensor})
        inds=[]
        for t in terms
            for i in t.indices
                push!(inds, i)

            end
        end
        new(terms,inds)
    end
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

function first_element(arr, elem1, elem2)
    idx1 = findfirst(x -> x == elem1, arr)
    idx2 = findfirst(x -> x == elem2, arr)
    if idx1 == nothing || (idx2 != nothing && idx2 < idx1)
        return elem2
    else
        return elem1
    end
end




# function Base.show(io::IO, P::Vector{Pair{Index, IndexPosition}})
# print(io, p for p in P)
# end



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










 #t.upper_indices for t in terms

#  function update_indices(T::TensorTerm)
#     for t in T.terms

#  end


#  *(A::Tensor, B::Tensor)= TensorTerm([A,B])

#  function remove_element(vector::Vector{T}, element::T) where T
#     index = findfirst(vector .== element)
#     if index === nothing
#         return false  # element not found
#     else
#         splice!(vector, index)
#         return true  # element removed successfully
#     end
# end

#  function getFreeIndices(Ts::TensorTerm)
#     freeup::OrderedSet{Index}=OrderedSet([])
#     freelow::OrderedSet{Index}=OrderedSet([])
#     for term in Ts.terms
#         for ind in term.upper_indices
#         if !(ind in freeup)
#             push!(freeup, ind)
#         else
#             delete!(freeup, ind)
#         end


#         end

#         for ind in term.lower_indices
#             if !(ind in freeup) && !(ind in freelow)
#                 push!(freelow, ind)
#             elseif (ind in freeup)
#                 delete!(freeup, ind)
#             else 
#                 delete!(freelow, ind)
#             end
#             end
#     end
#     return [freeup, freelow]
#  end

#  function getDummieIndices(Ts::TensorTerm)
#     dumup::OrderedSet{Index}=OrderedSet([])
#     dumlow::OrderedSet{Index}=OrderedSet([])
#     tempup::OrderedSet{Index}=OrderedSet([])
#     templow::OrderedSet{Index}=OrderedSet([])
    
#     for term in Ts.terms
#         for ind in term.upper_indices
#         if !(tempup in freeup)
#             push!(freeup, ind)
#         else
#             delete!(freeup, ind)
#         end


#         end

#         for ind in term.lower_indices
#             if !(ind in freeup) && !(ind in freelow)
#                 push!(freelow, ind)
#             elseif (ind in freeup)
#                 delete!(freeup, ind)
#             else 
#                 delete!(freelow, ind)
#             end
#             end
#     end
#     return [freeup, freelow]
#  end


 

# #  function Base.show(io::IO, T::TensorTerm)
# #     for t in T
# #         print(io, "MyType(")
# #     end
    
# #     show(io, x.x)
# #     print(io, ", ")
# #     show(io, x.y)
# #     print(io, ")")
# # end

