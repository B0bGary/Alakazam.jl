##############################################################################################
# Integration by parts
#
# Reduce expressions by moving partials, "integrating by parts"
# and discarding "total derivative" boundary terms
#
# Quadratic is easy. Cubics and higher need more sophisticated machinery.
##############################################################################################

"""
    strip_derivatives(t) -> (base, ops)
Similar to root_operand, but returns the operator tree and fails
if not deriv of a single tensor. (Use product rule first)

"""
function strip_derivatives(t::TensorSuperType)
    ops = PartialSuperType[]
    cur = t
    while cur isa PartialSuperType
        nxt = isTensor(operand(cur))
        nxt === nothing && break
        push!(ops, cur)
        cur = nxt
    end
    return (cur, ops)
end
strip_derivatives(t::Tensorial) = (t, PartialSuperType[])

"""Decides which of two base fields collects the derivatives; the larger key wins."""
_ibp_key(t::TensorSuperType) = string(t)

"""
    _shift_all(bases, ops, onto) -> (sign, TensorTerm)

Drive every derivative onto factor `onto`, returning the sign it costs. One move is
`(∂X)Y ≡ -(-1)^{|∂||X|} X(∂Y)`: a minus from integrating by parts, and a graded sign from the derivative
"""
function _shift_all(bases, ops, onto::Int)
    from = 3 - onto
    moving = ops[from]
    n = length(moving)
    # One move is `(∂X)Y ≡ -(-1)^{|∂||X|} X(∂Y)`: the minus from integrating by parts, and a graded sign when the derivative is odd, opposite of product rule basically
    sgn = (-1)^n
    if from == 1
        # moving rightwards, the derivative passes what is left of factor 1 once it and the ones
        # already moved have been peeled off
        g = grading(bases[1]) + sum(op_grading(o) for o in moving; init=0)
        for op in moving
            g -= op_grading(op)
            isodd(op_grading(op) * g) && (sgn = -sgn)
        end
    else
        # moving leftwards, it passes the whole of factor 1
        g = grading(bases[1]) + sum(op_grading(o) for o in ops[1]; init=0)
        for op in moving
            isodd(op_grading(op) * g) && (sgn = -sgn)
            g += op_grading(op)
        end
    end
    f = Vector{TensorSuperType}(undef, 2)
    f[from] = bases[from]
    f[onto] = rebuild_nested(vcat(reverse(moving), ops[onto]), bases[onto])
    return (sgn, TensorTerm(f))
end

"""
The normal form of a single term, as the coefficient pairs it contributes. Empty when the term
is itself a total derivative.
"""
function _ibp_term(c, tt::TensorTerm, onto::Symbol)
    fs = terms(tt)
    n = length(fs)
    n == 0 && return CoeffPair[CoeffPair(c, tt)]
    # a lone ∂[anything] is a total derivative outright
    n == 1 && return fs[1] isa PartialSuperType ? CoeffPair[] : CoeffPair[CoeffPair(c, tt)]
    b1, o1 = strip_derivatives(fs[1])
    b2, o2 = strip_derivatives(fs[2])
    bases, ops = (b1, b2), (o1, o2)
    k1, k2 = _ibp_key(b1), _ibp_key(b2)
    target = onto === :first ? 1 : onto === :last ? 2 : (k1 < k2 ? 2 : 1)
    if k1 == k2
        # Both factors are the same field. The two placements are equal modulo total derivatives, so if they are also negatives of each other as expressions then X ≡ -X
        # and the term drops out. This is what separates φ∂φ = ½∂(φ²)="0", from ψ∂ψ with ψ odd, which does not.
        (sa, ta) = _shift_all(bases, ops, 1)
        (sb, tb) = _shift_all(bases, ops, 2)
        X = simplify(TensorExpression(CoeffPair[CoeffPair(sa, ta)]))
        Y = simplify(TensorExpression(CoeffPair[CoeffPair(sb, tb)]))
        isZeroTens(simplify(X + Y)) && return CoeffPair[]
    end

    (s, t) = _shift_all(bases, ops, target)
    return CoeffPair[CoeffPair(c * s, t)]
end

"""
    drop_total_derivatives(exp; onto = :key) -> TensorExpression

Reduce `exp` modulo total derivatives. Every term quadratic in the
fields is driven to a normal form carrying all of its derivatives on one factor, after which a
total derivative cancels against itself and disappears. The result is zero exactly when `exp`
was a total derivative.

Terms with more than two factors fall through to `drop_total_derivatives_general`, which builds
the total-derivative subspace explicitly instead of relying on a normal form.
"""
function drop_total_derivatives(exp::TensorExpression; onto::Symbol=:key)
    onto in (:key, :first, :last) || throw(ArgumentError("onto must be :key, :first or :last"))
    e = simplify(derive(exp))
    acc = CoeffPair[]
    higher = CoeffPair[]
    for (c, tt) in summands(e)
        if length(terms(tt)) > 2
            push!(higher, CoeffPair(c, tt))
        else
            append!(acc, _ibp_term(c, tt, onto))
        end
    end
    out = isempty(acc) ? zero() : simplify(TensorExpression(acc))
    isempty(higher) && return out
    return simplify(out + drop_total_derivatives_general(TensorExpression(higher)))
end
drop_total_derivatives(t::Tensorial; kw...) = drop_total_derivatives(TensorExpression(t); kw...)

"""
    reduce_bp_with(exp, rules; passes = 12) -> TensorExpression

Reduce `exp` modulo total derivatives while applying replacement pattern `rules`, iterating until neither step
changes anything. Intended use is for rules relating a derivative of one field to another (eg conservation laws).
"""
function reduce_bp_with(exp::TensorExpression, rules; passes::Int=12)
    _sub(e) = foldl((a, r) -> replace_pattern(a, TensorTerm(r[1]), r[2]), rules;
                    init=TensorExpression(e))
    cur = simplify(_sub(exp))
    for _ in 1:passes
        acc = zero()
        for (c, tt) in summands(cur)
            t = TensorExpression(CoeffPair[CoeffPair(c, tt)])
            best, bestfired, bestn = nothing, false, typemax(Int)
            for o in (:key, :first, :last)
                d = drop_total_derivatives(t; onto=o)
                s = simplify(_sub(d))
                fired, n = s != simplify(d), length(summands(s))
                if best === nothing || (fired && !bestfired) || (fired == bestfired && n < bestn)
                    best, bestfired, bestn = s, fired, n
                end
            end
            acc = acc + best
        end
        nxt = simplify(acc)
        nxt == cur && break
        cur = nxt
    end
    return drop_total_derivatives(cur)
end
reduce_bp_with(t::Tensorial, rules; kw...) = reduce_bp_with(TensorExpression(t), rules; kw...)

##############################################################################################
# The general case: reduction modulo total derivatives by linear algebra.

# A total derivative preserves both the field content of a term and its derivative degree, 
# so the problem splits into finite sectors. Within one, the monomials reachable from those of the
# expression are enumerated, expanded, and are row-reduced.
##############################################################################################

"""
The canonical form of a monomial: the term `simplify` reduces it to, paired with the coefficient
that costs. `nothing` when the monomial vanishes or is not a single term.
"""
function _canonical_monomial(tt::TensorTerm)
    e = simplify(TensorExpression(tt))
    isZeroTens(e) && return nothing
    s = summands(e)
    length(s) == 1 || return nothing
    return (s[1][2], s[1][1])
end

"""Index a canonical monomial into `basis`, returning its position."""
function _basis_slot!(index::Dict{TensorTerm,Int}, basis::Vector{TensorTerm}, t::TensorTerm)
    haskey(index, t) && return index[t]
    push!(basis, t)
    return index[t] = length(basis)
end

"""
Expand an expression into a coefficient vector over `basis`, extending the basis with any
monomial not yet seen. Returns `nothing` if some summand is not a usable monomial.
"""
function _coefficient_vector!(index, basis, e::TensorExpression)
    acc = Dict{Int,Any}()
    for (c, tt) in summands(e)
        cm = _canonical_monomial(tt)
        cm === nothing && continue
        (t, k) = cm
        i = _basis_slot!(index, basis, t)
        acc[i] = get(acc, i, 0) + c * k
    end
    return acc
end

"""
Every total derivative of the sector, as coefficient vectors over a shared monomial basis.

Generators are obtained by lifting one derivative off a factor of a known monomial and differentiating it again. 
Each such `∂_μ(g)` is a total derivative by construction, and expanding it can reveal monomials not present in the original
expression, so the process repeats until done.
"""
function _total_derivative_rows(exp::TensorExpression; budget::Int=512)
    index = Dict{TensorTerm,Int}()
    basis = TensorTerm[]
    target = _coefficient_vector!(index, basis, exp)
    rows = Vector{Dict{Int,Any}}()
    seen = Set{TensorTerm}()
    frontier = copy(basis)
    while !isempty(frontier)
        length(basis) > budget && return (nothing, nothing, nothing)
        next = TensorTerm[]
        for m in frontier
            m in seen && continue
            push!(seen, m)
            parts = [strip_derivatives(f) for f in terms(m)]
            for (i, (base, ops)) in enumerate(parts)
                for j in eachindex(ops)
                    lowered = copy(parts)                       # one derivative lifted off factor i
                    lowered[i] = (base, ops[[k for k in eachindex(ops) if k != j]])
                    g = TensorTerm(TensorSuperType[rebuild_nested(o, b) for (b, o) in lowered])
                    before = length(basis)
                    row = _coefficient_vector!(index, basis,
                                               simplify(derive(_with_operand(ops[j], TensorExpression(g)))))
                    isempty(row) || push!(rows, row)
                    append!(next, basis[(before+1):end])
                end
            end
        end
        frontier = next
    end
    return (target, rows, basis)
end

"""Exact Gaussian elimination on sparse rows. gives the echelon rows keyed by pivot"""
function _echelon(rows::Vector{Dict{Int,Any}})
    pivots = Dict{Int,Dict{Int,Any}}()
    for row in rows
        r = copy(row)
        while true
            filter!(p -> !_is_zero_coefficient(p.second), r)
            isempty(r) && break
            p = minimum(keys(r))
            haskey(pivots, p) || (pivots[p] = _scale_row(r, inv(r[p])); break)
            _subtract_row!(r, pivots[p], r[p])
        end
    end
    return pivots
end

_is_zero_coefficient(c) = c isa Number && iszero(c)
_scale_row(r, s) = Dict{Int,Any}(k => v * s for (k, v) in r)
function _subtract_row!(r, piv, factor)
    for (k, v) in piv
        r[k] = get(r, k, 0) - factor * v
    end
    return r
end

"""
    drop_total_derivatives_general(exp; budget = 512) -> TensorExpression

Reduce `exp` modulo total derivatives at any degree, by constructing the total-derivative
subspace of each sector explicitly and reducing against it.

Returns `exp` unchanged if the sector grows past `budget` monomials, or if a coefficient cannot be divided exactly
"""
function drop_total_derivatives_general(exp::TensorExpression; budget::Int=512)
    e = simplify(derive(exp))
    (target, rows, basis) = _total_derivative_rows(e; budget=budget)
    target === nothing && return exp
    all(c -> c isa Number, values(target)) || return exp
    all(r -> all(c -> c isa Number, values(r)), rows) || return exp
    pivots = _echelon(rows)
    r = copy(target)
    while true
        filter!(p -> !_is_zero_coefficient(p.second), r)
        isempty(r) && break
        p = minimum(keys(r))
        haskey(pivots, p) || break
        _subtract_row!(r, pivots[p], r[p])
        filter!(q -> !_is_zero_coefficient(q.second), r)
        haskey(r, p) && break                       # elimination stalled
    end
    filter!(p -> !_is_zero_coefficient(p.second), r)
    isempty(r) && return zero()
    return simplify(TensorExpression(CoeffPair[CoeffPair(c, basis[i]) for (i, c) in r]))
end
drop_total_derivatives_general(t::Tensorial; kw...) =
    drop_total_derivatives_general(TensorExpression(t); kw...)
