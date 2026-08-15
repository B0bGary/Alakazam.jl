

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


#############################################################################
# Spin algebras, and the division algebra behind them
#############################################################################

"""
    division_algebra(algebra) -> Type{<:NumberSet} | Nothing

The endomorphism algebra of the defining representation, `ℝ`, `ℂ` or `ℍ`.

Frobenius-Schur uses this: a representation is real, complex or quaternionic according to whether its endomorphism
 algebra is ℝ, ℂ or ℍ. For `SL(2,𝔽)` the parameter is the division algebra, so the answer is read straight off the type.
For the others it is a property of the family.
"""
division_algebra(::SLAlgebra{F}) where {F<:NumberSet} = F
division_algebra(::SLAlgebra{OctonionNumberSet}) = ℝ
division_algebra(a::SUAlgebra) = a.N == 2 ? ℍ : ℂ
division_algebra(::SOAlgebra) = ℝ
division_algebra(::SpAlgebra{F}) where {F<:NumberSet} = F <: QuaternionNumberSet ? ℍ : ℝ
division_algebra(::LieAlgebraSuperType) = nothing

"""
    Spin(p, q) -> LieAlgebraSuperType

The spin algebra of a space with `p` timelike and `q` spacelike directions, in the low dimensions
where it coincides with a classical algebra.
"""
struct SpinAlgebra{P,Q} <: LieAlgebraSuperType
    N::Int                                   # index range: components of the minimal spinor
    function SpinAlgebra{P,Q}(n::Integer) where {P,Q}
        (P isa Int && Q isa Int) || throw(ArgumentError("signature must be integers"))
        (P >= 0 && Q >= 0) || throw(ArgumentError("signature entries must be non-negative"))
        new{P,Q}(Int(n))
    end
    SpinAlgebra{P,Q}() where {P,Q} = SpinAlgebra{P,Q}(_minimal_spinor_dim(P + Q))
end

"""
Components of the minimal spinor in `d` dimensions, counted over ℂ: the Dirac spinor has
`2^(d÷2)`, and in even `d` the Weyl projection halves it. A Majorana condition does not change the
count, only relates components to their conjugates.
"""
_minimal_spinor_dim(d::Integer) = 2^(div(d, 2) - (iseven(d) ? 1 : 0))

spacetime_dimension(::SpinAlgebra{P,Q}) where {P,Q} = P + Q
signature(::SpinAlgebra{P,Q}) where {P,Q} = (P, Q)
"""
Invariant form, but only where naming one is safe.
"""
metric(a::SpinAlgebra) = a.N == 2 ? EpsilonTensor : nothing
Base.show(io::IO, ::SpinAlgebra{P,Q}) where {P,Q} = print(io, "spin(", P, ",", Q, ")")

"""
Reality of the minimal chiral spinor.

`Spin(p,q)` sits inside the even subalgebra `Cl⁰(p,q) ≅ Cl(p, q-1)`, whose type is fixed by
`(p - q + 1) mod 8`.

Checked against the standard Lorentzian table for `d = 2..11`: ℝ at 3, 9, 11 for Majorana and ℍ at
5, 7 for symplectic Majorana in odd dimensions, and in even dimensions the type predicts
Majorana-Weyl exactly.

What it does not report is whether a Majorana condition exists, which is a question about the full
`Cl(p,q)` and is genuinely convention-dependent. Working in a Majorana basis means declaring a 
self-conjugate index set directly.
"""
function division_algebra(::SpinAlgebra{P,Q}) where {P,Q}
    m = mod(P - Q + 1, 8)
    m in (0, 1, 2) ? ℝ :
    m in (3, 7)    ? ℂ :
                     ℍ            # 4, 5, 6
end

"""
    Spin(p, q) -> SpinAlgebra{p,q}
    Spin(d)    -> SpinAlgebra{1,d-1}

The spin algebra of a space with `p` timelike and `q` spacelike directions, for any signature.

In low dimensions this coincides with a classical algebra, and `classical(a)` returns that
realisation where one exists:

    Spin(2,1) ≅ SL(2,ℝ)    d = 3     Majorana
    Spin(3,1) ≅ SL(2,ℂ)    d = 4     Weyl
    Spin(5,1) ≅ SL(2,ℍ)    d = 6     symplectic Majorana-Weyl
    Spin(9,1) ≅ SL(2,𝕆)    d = 10    Majorana-Weyl
"""
Spin(p::Integer, q::Integer) = SpinAlgebra{Int(p),Int(q)}()
Spin(d::Integer) = Spin(1, Int(d) - 1)

"""
    classical(a::SpinAlgebra) -> LieAlgebraSuperType | Nothing

The classical algebra `a` coincides with, where an exceptional isomorphism exists, and `nothing`
otherwise. Both orderings of the signature are accepted, `Spin(p,q) ≅ Spin(q,p)`.
"""
function classical(::SpinAlgebra{P,Q}) where {P,Q}
    s = minmax(P, Q)
    s == (1, 2) ? SLAlgebra{ℝ}(2) :
    s == (1, 3) ? SLAlgebra{ℂ}(2) :
    s == (1, 5) ? SLAlgebra{ℍ}(2) :
    s == (1, 9) ? SLAlgebra{𝕆}(2) :
    s == (0, 3) ? SUAlgebra{ℂ}(2) :          # Spin(3) ≅ SU(2)
    s == (0, 5) ? SpAlgebra{ℍ}(2) :          # Spin(5) ≅ Sp(2)
    s == (0, 6) ? SUAlgebra{ℂ}(4) :          # Spin(6) ≅ SU(4)
    nothing
end
classical(::LieAlgebraSuperType) = nothing

"""Spacetime dimension a spin algebra belongs to, inverting [`Spin`](@ref) for the `SL(2,𝔽)` family:
`d = 2 + dim_ℝ(𝔽)`."""
spacetime_dimension(a::SLAlgebra) = a.N == 2 ? 2 + _real_dim(division_algebra(a)) : nothing
spacetime_dimension(::LieAlgebraSuperType) = nothing

_real_dim(::Type{RealNumberSet}) = 1
_real_dim(::Type{ComplexNumberSet}) = 2
_real_dim(::Type{QuaternionNumberSet}) = 4
_real_dim(::Type{OctonionNumberSet}) = 8
_real_dim(::Nothing) = nothing


# so(d::Integer, field::NumberSet) = SOAlgebra(d, field)
# su(d::Integer, field::NumberSet) = SUAlgebra(d, field)

