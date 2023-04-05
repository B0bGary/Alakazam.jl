
abstract type IndexSuperType end


struct IndexPosition
   position::Bool
   #0 for upper index, 1 for lower index
end

mutable struct IndexSet
   name::String
   range::Int
   indices::Vector
   default_position::IndexPosition
   function IndexSet(name::String, range::Int=1, indices::Vector=[], default_position::IndexPosition=upper)
      return new(name, range, indices, default_position)
   end
   function IndexSet(name::String, indices::Vector=[], default_position::IndexPosition=upper)
      return new(name, 0, indices, default_position)
   end
end

mutable struct Index<:IndexSuperType
   name::String
   index_set::IndexSet
   grading::Int

 function Index(name::String,index_set::IndexSet, grading::Int=0)
   temp=new(name,index_set,grading)
   push!(index_set.indices,temp)
   return temp
 end

 function Index(name::String, grading::Int=0)
   new(name,IndexSet(""),grading)
 end
end

mutable struct DummyPatternIndex<:IndexSuperType
   name::String
   index_set::IndexSet
   grading::Int

 function DummyPatternIndex(name::String,index_set::IndexSet, grading::Int=0)
   temp=new(name,index_set,grading)
   push!(index_set.indices,temp)
   return temp
 end

 function DummyPatternIndex(name::String, grading::Int=0)
   new(name,IndexSet("",[]),grading)
 end
end

mutable struct Coordinate
   name::String
   index_set::Vector{IndexSet}
   function Coordinate(name::String, index_set2::IndexSet)
      new(name,[index_set2])
    end
    function Coordinate(name::String, index_set::Vector{IndexSet})
      new(name,index_set)
    end
end








function in_default_pos(i::Pair{Index,IndexPosition})
   return i[2]==i[1].index_set.default_position

end





upper=IndexPosition(0)
lower=IndexPosition(1)


function ==(a::Index, b::Index)
   return a.name==b.name && a.grading==b.grading # index set ....  would cause issues since index set is a vecot. would have to use name but good for now
end

function ==(a::DummyPatternIndex, b::DummyPatternIndex)
   return a.name==b.name && a.grading==b.grading # index set ....  would cause issues since index set is a vecot. would have to use name but good for now
end



function Base.show(io::IO, C::Coordinate)
   print(io, C.name)
   for ind in C.index_set
      print(io,"[")
   print(io,ind.indices[1].name)
   print(io,"]")
   end

end


function Base.show(io::IO, C::Vector{Coordinate})
   print(io,"{")
   it= Iterators.Stateful(C)
   for i in it
      print(io, i)
   if peek(it)!=nothing
      print(io, ", ")
   end
   end
   print(io,"}")
end


function Base.show(io::IO, I::IndexSuperType)
      if typeof(I)==DummyPatternIndex
         print(io, "?")
      end
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

function Base.copy(index::Index)
      return Index(index.name,index.index_set,index.grading)
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

