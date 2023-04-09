include("indices.jl")
include("custom_exceptions.jl")

SymorNum = Union{SymbolicUtils.Symbolic,Number}
abstract type TensorSuperType end
abstract type TotallyAntisymmetricSuperType <: TensorSuperType end
abstract type TotallySymmetricSuperType <: TensorSuperType end

mutable struct DummyPatternTensor<:TensorSuperType
    name::String
    indices::Vector{Pair{I, IndexPosition}}  where {I<:IndexSuperType}
    function_of::Vector{Coordinate}

    function DummyPatternTensor(name::String="?", indices::Vector{Pair{I, IndexPosition}} where {I<:IndexSuperType} =Vector{Pair{IndexSuperType, IndexPosition}}()  ,
        function_of::Vector{Coordinate}=Vector{Coordinate}() )  #lets you have scalars with no indices
        new(name, indices,function_of)
    end

    function DummyPatternTensor(name::String="?", indices::Vector{Pair{I, IndexPosition}} where {I<:IndexSuperType} =Vector{Pair{IndexSuperType, IndexPosition}}() )  #lets you have scalars with no indices
        new(name, indices,[])
    end
    function DummyPatternTensor(name::String="?", indices::Vector{Pair{I, IndexPosition}} where {I<:IndexSuperType} =Vector{Pair{IndexSuperType, IndexPosition}}(),function_of::Coordinate=Coordinate("",Vector{IndexSet}()) )  #lets you have scalars with no indices
        new(name, indices,[function_of])
    end

 end



mutable struct Tensor<:TensorSuperType
    name::String
    indices::Vector{Pair{I, IndexPosition}} where {I<:IndexSuperType}
    function_of::Vector{Coordinate}

    function Tensor(name::String, indices::Vector{Pair{I, IndexPosition}} =Vector{Pair{IndexSuperType, IndexPosition}}()  ,
        function_of::Vector{Coordinate}=Vector{Coordinate}())  where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices,function_of)
    end
    

    function Tensor(name::String, indices::Vector{Pair{I, IndexPosition}} =Vector{Pair{IndexSuperType, IndexPosition}}()  ,
        function_of::Coordinate=nothing) where {I<:IndexSuperType} #lets you have scalars with no indices
        new(name, indices,[function_of])
    end

    function Tensor(name::String,
        function_of::Coordinate=nothing) #lets you have scalars with no indices
        new(name, Vector{Pair{IndexSuperType, IndexPosition}}(),[function_of])
    end

    function Tensor(name::String, indices::Vector{Pair{I, IndexPosition}}  =Vector{Pair{IndexSuperType, IndexPosition}}()  ) where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices,[])
    end
    function Tensor(T::DummyPatternTensor) #
        new(T.name,T.indices,T.function_of)
    end
 end




# struct SymmetricTensor <:TotallySymmetricSuperType
# end

mutable struct SymmetricTensor<:TotallySymmetricSuperType
    name::String
    indices::Vector{Pair{I, IndexPosition}} where {I<:IndexSuperType}
    function_of::Vector{Coordinate}

    function SymmetricTensor(name::String, indices::Vector{Pair{I, IndexPosition}} =Vector{Pair{IndexSuperType, IndexPosition}}()  ,
        function_of::Vector{Coordinate}=Vector{Coordinate}())  where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices,function_of)
    end
    

    function SymmetricTensor(name::String, indices::Vector{Pair{I, IndexPosition}} =Vector{Pair{IndexSuperType, IndexPosition}}()  ,
        function_of::Coordinate=nothing) where {I<:IndexSuperType} #lets you have scalars with no indices
        new(name, indices,[function_of])
    end

    function SymmetricTensor(name::String,
        function_of::Coordinate=nothing) #lets you have scalars with no indices
        new(name, Vector{Pair{IndexSuperType, IndexPosition}}(),[function_of])
    end

    function SymmetricTensor(name::String, indices::Vector{Pair{I, IndexPosition}}  =Vector{Pair{IndexSuperType, IndexPosition}}()  ) where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices,[])
    end
    function SymmetricTensor(T::DummyPatternTensor) #
        new(T.name,T.indices,T.function_of)
    end
 end


# struct AntisymmetricTensor <:TotallyAntisymmetricSuperType
# end

mutable struct AntisymmetricTensor<:TotallyAntisymmetricSuperType
    name::String
    indices::Vector{Pair{I, IndexPosition}} where {I<:IndexSuperType}
    function_of::Vector{Coordinate}

    function AntisymmetricTensor(name::String, indices::Vector{Pair{I, IndexPosition}} =Vector{Pair{IndexSuperType, IndexPosition}}()  ,
        function_of::Vector{Coordinate}=Vector{Coordinate}())  where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices,function_of)
    end
    

    function AntisymmetricTensor(name::String, indices::Vector{Pair{I, IndexPosition}} =Vector{Pair{IndexSuperType, IndexPosition}}()  ,
        function_of::Coordinate=nothing) where {I<:IndexSuperType} #lets you have scalars with no indices
        new(name, indices,[function_of])
    end

    function AntisymmetricTensor(name::String,
        function_of::Coordinate=nothing) #lets you have scalars with no indices
        new(name, Vector{Pair{IndexSuperType, IndexPosition}}(),[function_of])
    end

    function AntisymmetricTensor(name::String, indices::Vector{Pair{I, IndexPosition}}  =Vector{Pair{IndexSuperType, IndexPosition}}()  ) where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices,[])
    end
    function AntisymmetricTensor(T::DummyPatternTensor) #
        new(T.name,T.indices,T.function_of)
    end
 end
    # function Tensor(name::String, indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}(),
    #     function_of::Coordinate=nothing) #lets you have scalars with no indices
    #     new(name, indices,[function_of])
    # end

    # function Tensor(name::String, indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) #lets you have scalars with no indices
    #     new(name, indices,[])
    # end


 

mutable struct TensorTerm
    terms::Vector{T} where {T<:TensorSuperType}
    indices::Vector{Pair{IndexSuperType, IndexPosition}}
    function TensorTerm(terms::Vector{T} where {T<:TensorSuperType})
        inds=[]
        for t in terms
            for i in t.indices
                push!(inds, i)

            end
        end
        new(terms,inds)
    end
    function TensorTerm(term::TensorSuperType)
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
    function TensorExpression(summands::Vector{T} where {T<:TensorSuperType})
        merge_sum(new([1=> TensorTerm(summand) for summand in summands]))
    end
    function TensorExpression(exp::TensorExpression, T::TensorSuperType)
        merge_sum(new([exp.summands; [1=>TensorExpression(T)]]))
    end
    function TensorExpression(exp::TensorExpression, T::TensorExpression)
        merge_sum(new([exp.summands; T.summands]))
    end
    function TensorExpression(num::SymorNum, T::TensorSuperType)
        merge_sum(new([num=>TensorTerm(T)]))
    end
    function TensorExpression(vec::Vector{<:Pair{<:SymorNum,<:TensorSuperType}}) 
        merge_sum(new([v[1]=>TensorTerm(v[2]) for v in vec]))
    end
    function TensorExpression(T::TensorSuperType)
        merge_sum(new([1=> TensorTerm(T)]))
    end
    function TensorExpression(T::TensorTerm)
        merge_sum(new([1 => T]))
    end
    function TensorExpression(p::Pair{T ,TensorTerm}) where {T<:SymorNum}
        merge_sum(new([p]))
    end
    function TensorExpression(v::Vector{Pair{T ,TensorTerm}} where {T<:SymorNum})
        merge_sum(new(v))
    end
    function TensorExpression(p::Pair{N ,T}) where {N<:SymorNum,T<:TensorSuperType} 
        merge_sum(new([p[1]=>TensorTerm(p[2])]))
    end
    function TensorExpression(tex::TensorExpression)
        tex
    end
    function TensorExpression()
        new([0=> TensorTerm(Tensor(""))])
    end
    function TensorExpression(n::SymorNum)
        new([n=> TensorTerm(Tensor(""))])
    end
end

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


#  *(A::Tensor, B::TensorSuperType)= TensorTerm([A,B])

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

