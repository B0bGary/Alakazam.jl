"""
(Super) Trace operator
"""
mutable struct Trace <: TraceSuperType
    data::OperatorData
    Trace(operand::Tensorial, args...) =
        new(OperatorData("Tr", TensorExpression(operand), args...))
    Trace(args...) =
        new(OperatorData("Tr", zero(), args...))
    Trace(d::OperatorData) = new(d)
end

#"Tr" is a word, and set in math italic it reads as the product T·r
_latex_op_name(::TraceSuperType) = "\\mathrm{Tr}"

"""Eliminate traces by expanding applying matrix reduction rules."""
eliminate_trace(T::Tensorial)=reduce_trace(factor_trace(T))



"""Position of the `n`th factor of `fs` that is (or is not) a matrix."""
function _findnth(fs, want_matrix::Bool, n::Int)
    seen = 0
    @inbounds for k in eachindex(fs)
        (fs[k] isa MatrixSuperType) == want_matrix || continue
        (seen += 1) == n && return k
    end
    return 0
end

"""
Factors out scalars from a trace, leaving only <:MatrixSuperType inside. Returns a TensorExpression.
Uses linearity, `str(aM1 + bM2) = a*str(M1) + b*str(M2)`, respects commutators and grading. 
A summand with no matrices at all leaves an explicit `IdMatrix` inside.
"""
function factor_trace(Tr::TraceSuperType)
    src = summands(operand(Tr))
    out = Vector{CoeffPair}(undef, length(src))
    m = 0

    @inbounds for (c, tt) in src
        fs = terms(tt)
        nmat = 0
        for f in fs
            f isa MatrixSuperType && (nmat += 1)
        end
        m += 1

        if nmat == length(fs)                        # nothing to factor out
            out[m] = CoeffPair(c, TensorTerm(TensorSuperType[_with_operand(Tr, TensorExpression(tt))]))
            continue
        end

        nrest = length(fs) - nmat
        mats = Vector{TensorSuperType}(undef, nmat)
        rest = Vector{TensorSuperType}(undef, nrest + 1)   # +1 slot for the trace itself
        i = j = 0
        for k in eachindex(fs)
            if fs[k] isa MatrixSuperType
                mats[i+=1] = fs[k]
            else
                rest[j+=1] = fs[k]
            end
        end

        perm = Vector{Int}(undef, length(fs))        # the order the factors move into: rest, then mats
        for k in eachindex(fs)
            perm[k] = k <= nrest ? _findnth(fs, false, k) : _findnth(fs, true, k - nrest)
        end
        sgn = _graded_perm_sign(fs, perm)

        inner = nmat == 0 ? TensorTerm(TensorSuperType[IdMatrix()]) : TensorTerm(mats)
        rest[end] = _with_operand(Tr, TensorExpression(inner))
        out[m] = CoeffPair(sgn * c, TensorTerm(rest))
    end

    m == 0 && return zero()
    m < length(out) && resize!(out, m)
    return TensorExpression(out)
end

# traces sitting inside a larger expression, or inside another operator's operand
function factor_trace(exp::TensorExpression)
    src = summands(exp)
    isempty(src) && return zero()
    acc = CoeffPair[]
    for (c, tt) in src
        for (c2, tt2) in summands(factor_trace(tt))
            push!(acc, CoeffPair(c * c2, tt2))
        end
    end
    return isempty(acc) ? zero() : TensorExpression(acc)
end
function factor_trace(tt::TensorTerm)
    out = 1
    for t in terms(tt)
        out *= factor_trace(t)
    end
    return out isa TensorExpression ? out : TensorExpression(out)
end
factor_trace(O::OperatorSuperType) = _with_operand(O, factor_trace(operand(O)))
factor_trace(t::Tensorial) = t


"""
The cyclic rotations of a term, each paired with the graded sign it costs.
Obeys the supertrace: `str(XY) = (-1)^{|X||Y|} str(YX)`.
"""
function cyclic_rotations(tt::TensorTerm)
    fs = terms(tt)
    n = length(fs)
    out = Tuple{Int,TensorTerm}[]
    n == 0 && return push!(out, (1, tt))     # str(1): the only rotation is the identity
    sgn = 1
    for r in 0:(n-1)
        push!(out, (sgn, TensorTerm(TensorSuperType[fs[mod1(k + r, n)] for k in 1:n])))
        if grading(fs[mod1(r + 1, n)]) == 1           # the factor about to move to the back
            crossed = 0
            for k in 2:n
                crossed += grading(fs[mod1(k + r, n)])
            end
            isodd(crossed) && (sgn = -sgn)
        end
    end
    return out
end

"""
    signed_cyclic_equals(a, b) -> (Bool, Int)

Whether two cyclic operators agree up to a rotation of their operands, and the graded sign relating
them.
"""
function signed_cyclic_equals(a::CyclicSuperType, b::CyclicSuperType)
    (typeof(a) === typeof(b) && name(a) == name(b) && op_indices(a) == op_indices(b)) ||
        return (false, 0)
    A, B = operand(a), operand(b)
    ta, tb = isTensorTerm(A), isTensorTerm(B)
    (ta === nothing || tb === nothing) && return (A == B ? (true, 1) : (false, 0))
    length(terms(ta)) == length(terms(tb)) || return (false, 0)
    for (s, rot) in cyclic_rotations(ta)      # rotation 0 is the identity, so equality is covered
        rot == tb && return (true, s)
    end
    return (false, 0)
end

# str(XY) and str(YX) are the same object only when the rotation is sign-free; when it is not they
# differ by a sign and so are not equal.
==(a::CyclicSuperType, b::CyclicSuperType) = signed_cyclic_equals(a, b) === (true, 1)


##############################################################################################
# Trace relations
#  counterpart of the `reduce_matrix` rules for a product of adjacent matrices
##############################################################################################
"""
Generic trace relation, overloaded by `@mrule` when its left-hand side is a `Trace`.
"""
reduce_trace(args...) = (false, zero())

"""
    _reduce_trace_once(T::TraceSuperType) -> (Bool, TensorExpression)

Try the registered trace relations against the operand of `T`, and against each of its cyclic
rotations.

Identities are dropped before matching, when other matrices are there
"""
function _reduce_trace_once(T::TraceSuperType)
    tt = isTensorTerm(operand(T))
    tt === nothing && return (false, zero())
    fs = terms(tt)
    (isempty(fs) || !all(f isa MatrixSuperType for f in fs)) && return (false, zero())

    if any(f isa IdMatrix for f in fs)
        st = nothing
        keep = TensorSuperType[]
        for f in fs
            f isa IdMatrix ? (f.index_set === nothing || (st = f.index_set)) : push!(keep, f)
        end
        # nothing but identities: the trace is the size of the space, if any of them knew it
        isempty(keep) && return st === nothing ? (false, zero()) :
                         (true, TensorExpression(st.range))
        tt = TensorTerm(keep)
    end

    for (sgn, rot) in cyclic_rotations(tt)
        rfs = terms(rot)
        (ok, res) = reduce_trace(rfs..., rot)
        ok && return (true, sgn * res)
    end
    return (false, zero())
end

"""
    reduce_trace(x)

Apply the registered trace relations to every trace in `x`, returning a `TensorExpression`.

Rules are registered with `@mrule` using a `Trace` on the left, for instance
@mrule Trace(Gamma(μ => lower)*Gamma(ν => lower)) = 4*Metric([μ => lower, ν => lower])
These are matched modulo cyclic perms, with approriate grading signs.
"""
function reduce_trace(T::TraceSuperType)
    (ok, res) = _reduce_trace_once(T)
    return ok ? res : TensorExpression(T)
end

function reduce_trace(exp::TensorExpression)
    src = summands(exp)
    isempty(src) && return zero()
    acc = CoeffPair[]
    for (c, tt) in src
        for (c2, tt2) in summands(reduce_trace(tt))
            push!(acc, CoeffPair(c * c2, tt2))
        end
    end
    return isempty(acc) ? zero() : TensorExpression(acc)
end
function reduce_trace(tt::TensorTerm)
    out = 1
    for t in terms(tt)
        out *= reduce_trace(t)
    end
    return out isa TensorExpression ? out : TensorExpression(out)
end
reduce_trace(O::OperatorSuperType) = _with_operand(O, reduce_trace(operand(O)))
reduce_trace(t::Tensorial) = TensorExpression(t)
signed_equals(a::CyclicSuperType, b::CyclicSuperType) = signed_cyclic_equals(a, b)