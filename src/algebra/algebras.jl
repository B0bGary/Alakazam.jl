

struct RealNumberSet <: NumberSet end
struct ImNumberSet <: NumberSet end
struct RationalNumberSet <: NumberSet end
struct ComplexNumberSet <: NumberSet end
struct QuaternionNumberSet <: NumberSet end
struct OctonionNumberSet <: NumberSet end

#could make one for each type and not store name...
Base.show(io::IO, ::NumberSet)=print(io, "𝔽")
Base.show(io::IO, ::RationalNumberSet)=print(io, "ℚ")
Base.show(io::IO, ::RealNumberSet)=print(io, "ℝ")
Base.show(io::IO, ::ImNumberSet)=print(io, "iℝ")
Base.show(io::IO, ::ComplexNumberSet)=print(io, "ℂ")
Base.show(io::IO, ::QuaternionNumberSet)=print(io, "ℍ")
Base.show(io::IO, ::OctonionNumberSet)=print(io, "𝕆")

const ℝ=RealNumberSet
const iℝ=ImNumberSet
const ℚ=RationalNumberSet
const ℂ=ComplexNumberSet
const ℍ=QuaternionNumberSet
const 𝕆=OctonionNumberSet


#Lie algebras. I think each should store (somehow) the appropriate type of invariant tensor. Upon construction, register the algebra.

# mutable struct SOAlgebra<:LieAlgebraSuperType
#    d::Integer
#    field::NumberSet
# end

struct SOAlgebra{F<:NumberSet}<:LieAlgebraSuperType
    N::Int
    function SOAlgebra{F}(d::T) where {T<:Integer,F<:NumberSet}
        d > 0 || throw(ArgumentError("dimension must be positive"))
        new{F}(d)
    end
    SOAlgebra(d::Integer) = SOAlgebra{ℝ}(d)
end

struct SLAlgebra{F<:NumberSet}<:LieAlgebraSuperType
    N::Int
    function SLAlgebra{F}(d::T) where {T<:Integer,F<:NumberSet}
        d > 0 || throw(ArgumentError("dimension must be positive"))
        new{F}(d)
    end
    SLAlgebra(d::Integer) = SLAlgebra{ℂ}(d)
end

struct SUAlgebra{F<:NumberSet}<:LieAlgebraSuperType
    N::Int
    function SUAlgebra{F}(d::T) where {T<:Integer,F<:NumberSet}
        d > 0 || throw(ArgumentError("dimension must be positive"))
        new{F}(d)
    end
    SUAlgebra(d::Integer) = SUAlgebra{ℝ}(d)
end

struct SpAlgebra{F<:NumberSet}<:LieAlgebraSuperType
    N::Int
    function SpAlgebra{F}(d::T) where {T<:Integer,F<:NumberSet}
        d > 0 || throw(ArgumentError("dimension must be positive"))
        new{F}(d)
    end
    SpAlgebra(d::Integer) = SpAlgebra{ℝ}(d)
end


"""Returns the metric-like invariant tensor for a given algebra."""
metric(algebra::SUAlgebra)=algebra.N==2 ? EpsilonTensor : nothing
metric(algebra::SLAlgebra)=algebra.N==2 ? EpsilonTensor : nothing
metric(::SOAlgebra)=Metric
metric(::SpAlgebra)=SymplecticForm
epsilon(::L) where L<:LieAlgebraSuperType = EpsilonTensor
epsilon(::SpAlgebra)=SymplecticForm

metric(algebra::Type{SUAlgebra})=algebra.N==2 ? EpsilonTensor : nothing
metric(algebra::Type{SLAlgebra})=algebra.N==2 ? EpsilonTensor : nothing
metric(::Type{SOAlgebra})=Metric
metric(::Type{SpAlgebra})=SymplecticForm
epsilon(::Type{L}) where L<:LieAlgebraSuperType = EpsilonTensor
epsilon(::Type{SpAlgebra})=SymplecticForm

dimension(algebra::L)where L<:LieAlgebraSuperType=algebra.N
rank(algebra::L)where L<:LieAlgebraSuperType=algebra.N-1
rank(algebra::SOAlgebra)=floor(algebra.N/2)
rank(algebra::SpAlgebra)=algebra.N/2


# so(d::Integer, field::NumberSet) = SOAlgebra(d, field)
# su(d::Integer, field::NumberSet) = SUAlgebra(d, field)

