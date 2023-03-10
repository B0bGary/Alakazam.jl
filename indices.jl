

mutable struct IndexSet
   name::String
   indices::Vector
   function IndexSet(name::String, indices::Vector=[])
      return new(name, indices)
   end
end

mutable struct Index
   name::String
   index_set::IndexSet
   grading::Int

 function Index(name::String,index_set::IndexSet, grading::Int=0)
   temp=new(name,index_set,grading)
   push!(index_set.indices,temp)
   return temp
 end

 function Index(name::String, grading::Int=0)
   temp=new(name,IndexSet("",[]),grading)
 end
end


struct Coordinate
   name::String
   index_set::IndexSet
end

struct IndexPosition
   position::Bool
   #0 for upper index, 1 for lower index
end

upper=IndexPosition(0)
lower=IndexPosition(1)


function Base.show(io::IO, I::Index)
       print(io, I.name)
end

function Base.show(io::IO, I::IndexPosition)
   print(io, I.position==0 ? "upper" : "lower")
end


function Base.show(io::IO, I::IndexSet)
   print(io,"[")
   it= Iterators.Stateful(I.indices)
   for i in it
      print(io, i.name)
   if peek(it)!=nothing
      print(io, ", ")
   end
   end
   print(io,"]")
end



# function rename_dummies(T::Tensor)









# struct SpinorIndex<:Index
#    indexset::IndexSet
#  end


#  struct IndexSet
#     indices::AbstractVector(Index)
#  end

#  grading(i::SpinorIndex) =1
#  grading(i::FlavourIndex) =0
#  grading(i::SpaceIndex) =0

#  struct FlavourIndex<:Index
    
#  end

#  struct SpaceIndex<:Index
   
#  end

