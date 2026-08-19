###############################################################################################
# Complex conjugation and Hermitian conjugation.
###############################################################################################

#################################################
# which number fields are fixed by conjugation
#################################################

"""
    isReal(T) / isImaginary(T) / isComplex(T)
ℚ counts as real: it is fixed by conjugation, which is the only property these predicates are
for. ℍ and 𝕆 answer `false` to all three, behaviour not yet implemented
"""
isReal(T::TensorSuperType) = _isReal(numberfield(T))
isImaginary(T::TensorSuperType) = _isImaginary(numberfield(T))
isComplex(T::TensorSuperType) = _isComplex(numberfield(T))

# dispatch on the singleton instances, which is what TensorData stores; the bare aliases (ℝ, ℂ...)
# are types, so accept those too for callers that pass them straight through
_isReal(::NumberSet) = false
_isReal(::RealNumberSet) = true
_isReal(::RationalNumberSet) = true

_isImaginary(::NumberSet) = false
_isImaginary(::ImNumberSet) = true

_isComplex(::NumberSet) = false
_isComplex(::ComplexNumberSet) = true

_isReal(F::Type{<:NumberSet}) = _isReal(F())
_isImaginary(F::Type{<:NumberSet}) = _isImaginary(F())
_isComplex(F::Type{<:NumberSet}) = _isComplex(F())

"""True when conjugation is the identity on this leaf, so it need not descend further."""
conj_is_identity(T::TensorSuperType) = _isReal(numberfield(T))

#################################################
# barred names
#################################################

"""
    isbarred(name) -> Bool
    togglebar(name) -> String

Conjugation maps a barred name to its unbarred partner and back, using the same combining macron
as the barred index alphabets (`barred` in indices.jl).
"""
isbarred(nm::AbstractString) = startswith(nm, COMBINING_MACRON)
togglebar(nm::AbstractString) =
    isbarred(nm) ? String(nm[nextind(nm, 1):end]) : string(COMBINING_MACRON) * nm

"""
    isdotted(name) -> Bool
    toggledot(name) -> String

The index-name counterpart of [`togglebar`](@ref), using the combining dot of the dotted index
alphabets (`dotted` in indices.jl), so that conjugating `psi_{alpha}` gives `psibar_{alphadot}`.

Names carry the dot; index sets carry the representation. Both move together under conjugation,
and this handles the name half.
"""
isdotted(nm::AbstractString) = startswith(nm, COMBINING_DOT)
toggledot(nm::AbstractString) =
    isdotted(nm) ? String(nm[nextind(nm, 1):end]) : string(COMBINING_DOT) * nm

"""
The decorations an index name may carry to mark it as belonging to a conjugate representation.
The dot is conventional for spinors, but nothing forces it. A bar, hat or tilde are all valid
and the alphabets in indices.jl provide all four.
"""
const INDEX_MARKS = (COMBINING_DOT, COMBINING_MACRON, COMBINING_HAT, COMBINING_TILDE)

"""
    indexmark(name) -> (mark, base)

Split a leading decoration off an index name: `indexmark("alphadot") == (COMBINING_DOT, "alpha")`.
Returns `(nothing, name)` when the name carries none of [`INDEX_MARKS`](@ref).
"""
function indexmark(nm::AbstractString)
    isempty(nm) && return (nothing, String(nm))
    c = first(nm)
    c in INDEX_MARKS && return (c, String(nm[nextind(nm, 1):end]))
    return (nothing, String(nm))
end

"""
    setmark(S, fallback) -> Char | Nothing

Infer which decoration an index set uses by looking at the indices already declared in it, rather
than assuming the dot. Returns `fallback` when the set is empty (nothing to learn from) or when its
indices disagree, so a set the user marked with a bar produces barred partners.
"""
function setmark(S, fallback)
    isempty(S.indices) && return fallback
    m = indexmark(first(S.indices).name)[1]
    for i in S.indices
        indexmark(i.name)[1] === m || return fallback
    end
    return m
end

#################################################
# the Dagger operator
#################################################

mutable struct Dagger <: DaggerSuperType
    data::OperatorData
    Dagger(operand::Tensorial) = new(OperatorData("✝", operand))
    Dagger() = new(OperatorData("✝", zero()))
end

"""
    T^✝

Postfix conjugation. `✝` is U+271D LATIN CROSS, standing in for `†` because **`†` (U+2020) is not a
legal Julia identifier**

Modifier letters (`ᴴ`,`ᵀ`,`ᵈ`) parse standalone but not here: Julia treats them as identifier *continuation* characters, so
in `T^ᴴ` the mark glues onto the `^` and yields an incomplete operator token. Superscripts (`⁺`)
fail the same way. Only a standalone symbol character works in the exponent position, and `✝` is
the one that looks like a dagger. This has the unfortunate issue that ✝ doesn't have a "\\..." command to enter...

This is why we also use the ` ' ` symbol to denote the adjoint below.

Applies the conjugation rather than parking an unevaluated operator on the expression.
"""
const ✝ = Dagger()

Base.:^(a::Tensorial, ::DaggerSuperType) = _dagger(a)

"""
    T'

Julia's postfix adjoint, same meaning as `T^✝`.
"""
Base.adjoint(a::Tensorial) = _dagger(a)

#################################################
# applying it
#################################################

"""
    apply_dagger(d) -> Tensorial

Carry out the conjugation `d` represents, rather than leaving it as an unevaluated operator.

Kept separate from the `Dagger` object itself so that an expression can hold a symbolic `†` when
the transformation of its operand is not known, and be resolved later once it is.
"""
apply_dagger(dagger::DaggerSuperType) = _dagger(operand(dagger))

"""
    _dagger(x) -> Tensorial

The recursion. Sums distribute (with conjugated coefficients), products reverse, leaves consult
their number field, and a nested dagger collapses by involution.
"""
_dagger(e::TensorExpression) = apply_antimorphism(_dagger, e, conj)

_dagger(tt::TensorTerm) = apply_antimorphism(_dagger, tt)

_dagger(d::DaggerSuperType) = apply_involution(_dagger, d)

function _dagger(T::TensorSuperType)
    # The field and the indices are independent and compose. Field decides the name and any
    # overall factor, each slot decides its own sign.
    f = numberfield(T)
    # A representation that is not isomorphic to its conjugate admits no reality condition, so a
    # tensor carrying such an index cannot be real or imaginary.
    if (_isReal(f) || _isImaginary(f))
        for (idx, _) in indices(T)
            conj_class(idx.index_set) === :complex && error(
                "tensor \"$(name(T))\" is declared over $f but carries index " *
                "\"$(name(idx))\" of set \"$(idx.index_set.name)\", whose conjugate is a " *
                "different representation; such a tensor is complex, so declare it over ℂ")
        end
    end
    _isReal(f)      && return _conj_leaf(T, false)
    _isImaginary(f) && return -1 * _conj_leaf(T, false)
    # only a complex field has a distinct conjugate to name, so only it toggles the bar
    _isComplex(f)   && return _conj_leaf(T, true)
    throw(ArgumentError(
        "conjugation is not defined for a tensor over $(numberfield(T)); " *
        "declare the tensor over ℝ, iℝ or ℂ, or add a rule for this field"))
end

#################################################
# special rules for the invariant forms: epsilon, Kronecker, metric
#################################################

"""
Conjugation of the invariant forms, dispatched on their supertypes rather than sent through the
generic leaf rule.

These are not barred/unbarred pairs, so toggling a name would be meaningless. The rule depends on the indexset conj indices, metric type

    (eps_{αβ})*  = eps_{α̇β̇}      (eps^{αβ})*  = eps^{α̇β̇}      spinor: NO sign
    (eps^{ij})*  = -eps_{ij}      (eps_{ij})*  = -eps^{ij}      SU(2):  sign, AND positions flip
    (δ_α^β)*     = δ_α̇^β̇          (δ_i^j)*     = δ_j^i          Kronecker

The spinor/SU(2) split is exactly `conj_class`, so the epsilon rule is the generic slot rule plus
the observation that an epsilon's own indices are always at their reference position: a quaternionic
epsilon therefore picks up one sign per index, i.e. `(-1)^2 = +1` from the pair... except that the
two indices flip together, which is the single lowering of the pair, giving the stated `-1`. It is
written explicitly for that reason.
"""
function _dagger(T::EpsilonTensorSuperType)
    inds = indices(T)
    isempty(inds) && return T
    S = inds[1][1].index_set
    cls = conj_class(S)
    cls === :fixed && return T
    s = conj_convention(S.metric_type, S)[2] 
    nT = copy(T)
    if cls === :quaternionic
        # SU(2): both indices flip together and the pair costs one sign
        set_indices!(nT, IndexPair[_conj_index(i) => _flip(p) for (i, p) in inds])
    else
        # spinor (complex): partner set, positions kept
        set_indices!(nT, IndexPair[_conj_index(i) => p for (i, p) in inds])
    end
    return s == 1 ? nT : s * nT
end

"""
Bras and kets exchange under the dagger: `⟨ψ|' = |ψ⟩` and `|ψ⟩' = ⟨ψ|`.
A bra is the adjoint of its ket, so the dagger has to return the other type.
The name is carried across unchanged.
"""
function _dagger(T::Union{BraSuperType,KetSuperType})
    # indices are conjugated exactly as for any other tensor (for now), then the type flips
    sgn = 1
    inds = IndexPair[]
    for (idx, pos) in indices(T)
        (ni, np, s) = _conj_slot(idx, pos, false)
        sgn *= s
        push!(inds, ni => np)
    end
    flipped = T isa BraSuperType ? Ket(data(T), inds) : Bra(data(T), inds)
    return sgn == 1 ? flipped : sgn * flipped
end

"""
Kronecker delta. `(δ_α^β)* = δ_α̇^β̇` keeps its positions; `(δ_i^j)* = δ_j^i` is the same object with
its slots read the other way round, which flipping both positions reproduces. No sign either way.
"""
function _dagger(T::KroneckerSuperType)
    inds = indices(T)
    isempty(inds) && return T
    S = inds[1][1].index_set
    cls = conj_class(S)
    cls === :fixed && return T
    nT = copy(T)
    flip = cls === :quaternionic
    set_indices!(nT, IndexPair[_conj_index(i) => (flip ? _flip(p) : p) for (i, p) in inds])
    return nT
end

"""A symmetric metric is real: `(η_{ab})* = η_{ab}`, indices and positions untouched."""
function _dagger(T::SymMetricSuperType)
    inds = indices(T)
    isempty(inds) && return T
    conj_class(inds[1][1].index_set) === :fixed && return T
    nT = copy(T)
    set_indices!(nT, IndexPair[_conj_index(i) => p for (i, p) in inds])
    return nT
end

#################################################
# how conjugation acts on an index
#################################################

"""
    conj_class(S::IndexSet) -> :fixed | :complex | :quaternionic

Frobenius-Schur classification of an index set, read off fields that already exist rather than
declared separately.

    conjugate === nothing         :fixed         no partner; the index is untouched
    conjugate === itself, sym     :fixed         SO(N) vector, spacetime, also untouched
    conjugate is a DIFFERENT set  :complex       spinor <-> dotted spinor
    conjugate === itself, antisym :quaternionic  SU(2) doublet, gets a minus

"""
function conj_class(S)
    c = S.conjugate
    c === nothing && return :fixed
    (c === S || c.name == S.name) || return :complex
    # Frobenius-Schur classifies by the endomorphism algebra, which `division_algebra` reads off
    # the set's algebra, or reports as nothing when it declares none.
    d = division_algebra(S)
    if d !== nothing
        d === ℍ && return :quaternionic
        d === ℝ && return :fixed
        d === ℂ && error("index set $(S.name) carries an algebra whose defining representation "*
                         "is complex, so its conjugate is a different representation, but the "*
                         "set names itself as its own conjugate; declare the partner set and "*
                         "point the two at each other. If a real form is meant, such as a "*
                         "Majorana spinor in a dimension whose chiral spinor is complex, leave "*
                         "the algebra off: a self-conjugate set is real by default.")
        # O has no endomorphism algebra among the three: fall through
    end
    # No algebra: a set that is its own conjugate is taken to be real.
    return :fixed
end

_flip(p::IndexPosition) = IndexPosition(!p.position)

"""Map one index onto its conjugate set, leaving the name alone. `nothing` or self-conjugate index
sets return the index unchanged, since there is no partner set to move to."""
function _conj_index(i)
    S = i.index_set
    c = S.conjugate
    (c === nothing || c === S || c.name == S.name) && return i
    return Index(i.name, c)
end


"""
    conj_convention(metric_type, S::IndexSet) -> (reference, sign)

The conjugation convention for indices of `S`. Dispatch on metric type.
Can override complex conjugation behaviour.

    reference   the level at which conjugation costs nothing
    sign        what a level away from the reference costs

Both the fields and the invariant form read this one pair.

The default takes the reference from the set's `default_position` and the sign from the class, which
reproduces Alakazam's own convention: `(eps_{αβ})* = eps_{αdot βdot}` with no sign, so covariant
spinors conjugate unsigned. Gates' "Superspace" instead sets `Cbar_{αβ} = -C_{αdot βdot}` from which 
`(psi_a)^dag = -psibar_{adot}` follows. This is encoded by eg:

    Alakazam.conj_convention(::Type{MyC}, S::IndexSet) = (upper, -1)
"""
conj_convention(::Type{<:MetricSuperType}, S::IndexSet) =(S.default_position, conj_class(S) === :quaternionic ? -1 : 1)

conj_convention(::Nothing, S::IndexSet) = (S.default_position,1)

"""
The convention for `S`, with the reference flipped for a barred field.
Derived from conj_convention.
θ carries its SU(2) index down where θ̄ carries it up, and the two share one index set.
This determines which picks up a sign.
"""
function _conj_rule(S, barred::Bool)
    (ref, sgn) = conj_convention(S.metric_type, S)
    return (barred ? _flip(ref) : ref, sgn)
end

"""
    _conj_partner(idx) -> Index

The index in the conjugate set standing for `idx`, matched on the undecorated name so that
whichever mark the user chose is respected: `alpha` finds `alphadot`, `alphabar` or a plain `alpha`
in the partner set, whichever is actually declared there.

Only when the partner set holds no match is a name built, and the mark for it is inferred from that
set's existing indices rather than assumed to be a dot. The dot is the last-resort default, used
for a fermionic set that is still empty and so offers nothing to infer from.
"""
function _conj_partner(idx)
    S = idx.index_set
    T = S.conjugate
    srcmark, base = indexmark(idx.name)
    for j in T.indices
        indexmark(j.name)[2] == base && return j
    end
    # nothing declared: undecorated if we came from a marked name, else this set's own convention
    fallback = (srcmark !== nothing || S.grading != fermionic.grading) ? nothing : COMBINING_DOT
    m = setmark(T, fallback)
    return Index(m === nothing ? base : string(m) * base, T)
end

"""
Conjugate one index slot: the new index, its new position, and the sign it contributes.

Position is preserved for a complex set (the partner set carries the variance) and flipped for a
quaternionic one (conjugation exchanges the representation with its dual on the same indices).
"""
function _conj_slot(idx, pos::IndexPosition, barred::Bool)
    S = idx.index_set
    cls = conj_class(S)
    # The default is 1, so a spacetime index is unaffected, while a Majorana set declaring (upper, -1) makes psi_a imaginary.
    # The barred flip belongs to the quaternionic case alone. It exists because θ and θ̄ carry their
    # SU(2) index on opposite levels, which is the same fact as the position flipping. Gor a complex
    # index the position is preserved and ψ and ψ̄ sit at the same level, so flipping the reference
    # would make the two applications asymmetric and conj∘conj would return -1 rather than the
    # identity.
    (ref, s) = _conj_rule(S, barred && cls === :quaternionic)
    sgn = pos == ref ? 1 : s
    cls === :fixed   && return (idx, pos, sgn)
    cls === :complex && return (_conj_partner(idx), pos, sgn)
    return (idx, _flip(pos), sgn)
end

"""
Conjugate a single complex tensor: toggle the bar on its name, send every index through
[`_conj_slot`](@ref), and multiply by the accumulated sign.

Reproduces the 4D 𝒩=2 reality conditions exactly:

    (θ^α_i)*     = θ̄^{α̇i}        SU(2) at reference        -> +
    (θ^{αi})*    = -θ̄^{α̇}_i      SU(2) displaced           -> -
    (θ̄_α̇^i)*     = θ_{αi}         barred, at ITS reference  -> +
    (θ̄_{α̇i})*    = -θ^i_α        barred, displaced         -> -
"""
function _conj_leaf(T::TensorSuperType, bar::Bool)
    barred = isbarred(name(T))
    sgn = 1
    moved = false
    inds = IndexPair[]
    for (idx, pos) in indices(T)
        (ni, np, s) = _conj_slot(idx, pos, barred)
        sgn *= s
        (ni === idx && np == pos) || (moved = true)
        push!(inds, ni => np)
    end
    # nothing to rebuild: no slot moved, no sign, and the name stays put
    (!bar && !moved && sgn == 1) && return T
    nT = copy(T)
    bar && set_name!(nT, togglebar(name(T)))
    (moved && !isempty(inds)) && set_indices!(nT, inds)
    return sgn == 1 ? nT : -1 * nT
end

"""
    apply_antimorphism(op, tt::TensorTerm)               -> Tensorial
    apply_antimorphism(op, e::TensorExpression, scal)    -> Tensorial

`O(AB) = O(B)O(A)`: reverse the factors and apply `op` to each.

The leaf action is a parameter rather than a fixed call to `_dagger`, because reversal is the whole
content of `AntimorphismSuperType.

`scal` is what the antimorphism does to a scalar coefficient.
dagger conjugates (`O(zA) = z* O(A)`), a plain transpose does not. It defaults to `conj` because
that is the *-algebra convention, and `conj` is the identity on reals anyway.

No Koszul signs arise from the reversal itself: the reversed product is built in its new order and
`*` supplies whatever grading signs that ordering carries.
"""
function apply_antimorphism(op, tt::TensorTerm)
    ts = terms(tt)
    length(ts) == 1 && return op(ts[1])
    out = op(ts[end])
    for k in (length(ts) - 1):-1:1
        out = out * op(ts[k])
    end
    return out
end

function apply_antimorphism(op, e::TensorExpression, scal = conj)
    acc = zero()
    for (c, t) in summands(e)
        acc = acc + scal(c) * op(t)
    end
    return acc
end

"""
    apply_involution(op, inv::InvolutionSuperType) -> Tensorial

`O(O(A)) = A`. An operator whose operand is another of the same type collapses to that operand;
otherwise the inner expression is transformed by `op` normally.

Like [`apply_antimorphism`](@ref) the leaf action is passed in, so this serves any
`InvolutionSuperType` rather than dagger alone. The cancellation already keys on `typeof`, so an
involution never collapses against a different operator that happens to sit under it.
"""
function apply_involution(op, inv::InvolutionSuperType)
    inner = operand(inv)
    typeof(inner)==typeof(inv) && return operand(inner) #dagger of a dagger cancels, make sure theyre the same type
    # an expression that is a single unevaluated dagger collapses the same way
    if inner isa TensorExpression && singleterm(inner) && first_coef(inner) == 1
        t = first_summand(inner)
        if t isa TensorTerm && length(terms(t)) == 1 && typeof(terms(t)[1])==typeof(inv)
            return operand(terms(t)[1])
        end
    end
    return op(inner)
end

function Base.show(io::IO, T::DaggerSuperType)
    print(io, "(")#maybe or i keep the [[ brackets
    print(io, operand(T))
    print(io, ")^")
    print(io, name(T))
end
