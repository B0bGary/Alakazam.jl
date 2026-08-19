mutable struct LieGenerator <: GeneratorDerivSuperType
    data::OperatorData
    coord::Vector{Coordinate}   # with respect to; shared between copies, never mutated in place

    LieGenerator(name::String, operand::Tensorial, coord::Coordinate, args...) =
        new(OperatorData(name, operand, args...), Coordinate[coord])

    LieGenerator(name::String, operand::Tensorial, coord::Vector{Coordinate}, args...) =
        new(OperatorData(name, operand, args...), _uniq(coord))

    LieGenerator(name::String, coord::Coordinate, args...) =
        new(OperatorData(name, zero(), args...), Coordinate[coord])

    LieGenerator(name::String, coord::Vector{Coordinate}, args...) =
        new(OperatorData(name, zero(), args...), _uniq(coord))

    LieGenerator(d::OperatorData, coord::Vector{Coordinate}) = new(d, coord)
end



function Base.show(io::IO, T::GeneratorDerivSuperType)
      show_op_base(io,T)
    if operand(T)!=zero()
    print(io,"∘")
    print(io, operand(T))
    end
end



"""
Generates the algebra so(d,2)
"""
function generate_conformal_algebra(d::Int)
   if "spacetime" in keys(SET_REGISTRY)
      iSet=SET_REGISTRY["spacetime"]
      x=COORD_REGISTRY["x"]
   else
      iSet=IndexSet("spacetime", d, lower)
      μ = Index("μ", iSet)#generate a single index so that indexgen has a charset to draw from
      x=Coordinate("x",iSet) 
   end
   
   P=LieGenerator("P",x,[i for i in safe_pop!(iSet,1)])
   M=LieGenerator("M",x,[i for i in safe_pop!(iSet,2)],generate_antisymmetric_tableau([1,2]))
   D=LieGenerator("D",x)
   K=LieGenerator("K",x,[i for i in safe_pop!(iSet,1)])

   I=safe_pop!(iSet,4)

   rhslorentz=im*(Metric(I[1],I[3])*LieGenerator("M",x,[I[2],I[4]])-Metric(I[1],I[4])*LieGenerator("M",x,[I[2],I[3]])-Metric(I[2],I[3])*LieGenerator("M",x,[I[1],I[4]])+Metric(I[2],I[4])*LieGenerator("M",x,[I[1],I[3]]))
   @commutator [M,LieGenerator("M",x,[i for i in I[3:end]])]=rhslorentz
   @commutator [M,LieGenerator("P",x,[I[3]])]=im*(Metric(I[1],I[3])*LieGenerator("P",x,[I[2]])-Metric(I[2],I[3])*LieGenerator("P",x,[I[1]]))
   @commutator [M,LieGenerator("K",x,[I[3]])]=im*(Metric(I[1],I[3])*LieGenerator("K",x,[I[2]])-Metric(I[2],I[3])*LieGenerator("K",x,[I[1]]))

   @commutator [D,P]=im*P
   @commutator [D,K]=-im*K
   @commutator [D,M]=0

   @commutator [P,LieGenerator("P",x,[I[2]])]=0
   @commutator [K,LieGenerator("K",x,[I[2]])]=0
   @commutator [K,LieGenerator("P",x,[I[2]])]= 2*im*(Metric(I[1],I[2])*D-M)

   return [P,M,D,K]
end
