include("indices.jl")
include("custom_exceptions.jl")
include("tensor_ops.jl")


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
    function TensorTerm(term::Tensor)
        new([term], term.indices)
    end
 end

 
 mutable struct TensorExpression
    summands::Vector{Pair{SymorNum, TensorTerm}}
    function TensorExpression(summand::TensorTerm)
        merge_sum(new(Vector((1,summand))))
    end
    function TensorExpression(summands::Vector{TensorTerm})
        merge_sum(new(Vector((1,summand) for summand in summands)))
    end
    function TensorExpression(summands::Vector{Tensor})
        merge_sum(new([1=> TensorTerm(summand) for summand in summands]))
    end
    function TensorExpression(exp::TensorExpression, T::Tensor)
        merge_sum(new([exp.summands; [1=>TensorExpression(T)]]))
    end
    function TensorExpression(exp::TensorExpression, T::TensorExpression)
        merge_sum(new([exp.summands; T.summands]))
    end
    function TensorExpression(num::SymorNum, T::Tensor)
        merge_sum(new([num=>TensorTerm(T)]))
    end
    function TensorExpression(vec::Vector{Pair{SymorNum,Tensor}})
        merge_sum(new([v[1]=>TensorTerm(v[2]) for v in vec]))
    end
    function TensorExpression(T::Tensor)
        merge_sum(new([1=> TensorTerm(T)]))
    end
    function TensorExpression(T::TensorTerm)
        merge_sum(new([1 => T]))
    end
end









# function Base.show(io::IO, P::Vector{Pair{Index, IndexPosition}})
# print(io, p for p in P)
# end



 




# function get_dummies_positional(T::TensorExpression)
#     dum::Set{Pair{Index, IndexPosition}}=Set()
#     for d in get_dummies(T.summands[1]) 
#         push!(dum, d=>upper)
#         push!(dum, d=>lower)
#     end
#     return dum
# end









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

