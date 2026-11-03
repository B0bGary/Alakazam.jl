
"""Computes the Euler-Lagrange equations for a given lagrangian and tensor field."""
function EL(lag::Tensorial, A::TensorSuperType, DA::DerivativeSuperType)
    lag = TensorExpression(lag)
    #compute the variational derivative of the lagrangian with respect to the tensor field A
    v=VariationalDerivative(lag, A)
    DA = copy(DA)                    #the caller's derivative is theirs: it is given A on a copy
    set_operand!(DA, A)
    #apply the derivative operator to the variational derivative
    dv=apply_derivative(v)
    v∂=VariationalDerivative(lag, DA)
    dv∂=apply_derivative(v∂)
    dn=apply_derivative(simplify(PD(name(DA), simplify(dv∂), coord(DA), op_indices(DA))))
    #simplify the result
    return simplify(dv-dn)
end


FD(L::Tensorial, φ::TensorSuperType)=functional_derivative(L, φ)
"""
    functional_derivative(L, φ)

The functional derivative `δS/δφ` of the action `S = ∫L`, to any order in derivatives:

    δS/δφ = ∂L/∂φ - D_μ(∂L/∂(D_μφ)) + D_ν D_μ(∂L/∂(D_μ D_νφ)) - …

Boundary terms are dropped, which for a covariant derivative assumes the integration by parts holds for it. 
"""
function functional_derivative(L::Tensorial, φ::TensorSuperType)
    L = product_rule(copy(TensorExpression(L)))   #a copy: product_rule sets operands in place
    claimed = IndexTypes[]
    collect_indices!(claimed, L)
    collect_indices!(claimed, φ)
    #one variable per distinct chain: an occurrence that an earlier variable already takes is not new
    vars = TensorSuperType[]
    for (_, tt) in summands(L), t in terms(tt)
        t isa DerivativeSuperType || continue
        ops = get_operator_chain(t, DerivativeSuperType)
        cur = isTensor(operand(last(ops)))
        (cur isa TensorSuperType && name(cur) == name(φ)) || continue
        any(v -> _field_variation(t, v) !== nothing, vars) && continue
        #the same chain on φ itself, its derivative indices fresh
        inner = TensorExpression(φ)
        for op in reverse(ops)
            fresh = IndexPair[]
            for (i, p) in op_indices(op)
                n = isLightcone(i) ? i : _fresh_index(i.index_set, claimed)
                push!(claimed, n)
                push!(fresh, n => p)
            end
            new = copy(op)
            set_op_indices!(new, fresh)
            set_operand!(new, inner)
            inner = TensorExpression(new)
        end
        v = isTensor(inner)
        _field_variation(t, v) === nothing || push!(vars, v)
    end
    out = CoeffPair[]
    append!(out, summands(apply_derivative(VariationalDerivative(L, φ))))
    for v in vars
        P = apply_derivative(VariationalDerivative(L, v))
        isZeroTens(P) && continue
        #∫ P D₁D₂…δφ = (-1)ᵏ ∫ (…D₂D₁P) δφ: the outermost derivative comes off first
        chain = get_operator_chain(v, DerivativeSuperType)
        for d in chain
            P = apply_derivative(TensorExpression(_with_operand(d, TensorExpression(P))))
        end
        append!(out, summands(isodd(length(chain)) ? -P : P))
    end
    isempty(out) && return zero()
    return simplify(TensorExpression(out))
end

copy(V::VariationalSuperType) = typeof(V)(symcopy(V.data), symcopy(V.var))

"""
    VariationalDerivative(L, φ)
"""
mutable struct VariationalDerivative <: VariationalSuperType
    data::OperatorData
    var::TensorSuperType
    function VariationalDerivative(operand::Tensorial, var::TensorSuperType, args...)
        new(OperatorData("δ", TensorExpression(operand), args...), var)
    end
    function VariationalDerivative(data::OperatorData, var)
        new(data, var)
    end
end

"""Custom grading for variational derivatives to account for the w/r variable"""
function grading(D::VariationalSuperType)
    grade = (grading(operand(D))+grading(D.var))%2
    return grade
end
"""Custom grading for variational derivatives to account for the w/r variable"""
function op_grading(D::VariationalSuperType)
    grade =grading(D.var)
    return grade % 2
end

coord(V::VariationalSuperType)=get_operand_function_of(operand(V))

function apply_derivative(D::VariationalSuperType)
    var = D.var
    vnames = [q[1] for q in indices(var)]
    L = operand(D)
    any(q -> q[1] in vnames, free_indices_ordered(L)) &&
        throw(ArgumentError("the field $var shares an index name with a free index of $L, which the " *
                            "variation would then carry twice; name the field's indices apart"))
    #dummies of L that the field's indices would meet in a delta are renamed out of the way, on a copy:
    #L is the caller's, and product_rule too sets its operators' operands in place
    L = product_rule(rename_dummies_excluded(L, vnames))
    #the Leibniz rule, graded: each occurrence of the field is replaced by deltas joining its slots to the
    #field's, the sign being that of moving the variation past the factors before it
    gvar = grading(var)
    out = CoeffPair[]
    for (coef, tt) in summands(L)
        ts = terms(tt)
        grad = 0
        for j in eachindex(ts)
            d = _field_variation(ts[j], var)
            if d !== nothing
                f = TensorSuperType[ts[1:j-1]; d; ts[j+1:end]]
                append!(out, summands((isodd(grad) ? -coef : coef) * TensorTerm(f)))
            end
            grad += grading(ts[j]) * gvar
        end
    end
    isempty(out) && return zero()
    return _project_onto_field(merge_sum(TensorExpression(out)), var)
end

"""
The variation of the factor `t` by the field `var`: a delta joining each slot of `t` to the same slot
of `var`, or a metric where their levels differ. `nothing` when `t` is not that field, a derivative of
it included, which is a variable of its own; and at a different level over a set with no metric.
"""
function _field_variation(t::TensorSuperType, var::TensorSuperType)
    (typeof(t) == typeof(var) && name(t) == name(var)) || return nothing
    if t isa DerivativeSuperType
        #the same chain all the way down: ∂_μ∂_ν A is not ∂_a∂_b B
        a, b = isTensor(operand(t)), isTensor(operand(var))
        (a isa TensorSuperType && b isa TensorSuperType && _field_variation(a, b) !== nothing) || return nothing
    end
    ti, vi = indices(t), indices(var)
    length(ti) == length(vi) || return nothing
    d = TensorSuperType[]
    for ((i, p), (j, q)) in zip(ti, vi)
        if p != q
            g = generate_metric(i => p, j => p)     #nothing over a set that cannot raise or lower
            g === nothing && return nothing
            push!(d, g)
        else
            push!(d, KroneckerDelta(i => p, j => op_pos(q)))
        end
    end
    return d
end

"""
`R`, a variation by the field `var`, projected by the transpose of the field's Young projector, so that
it pairs correctly with every variation the field's symmetry allows. For a symmetric or antisymmetric
field that is the field's own symmetry; for a mixed one it is the transposed symmetry.
"""
function _project_onto_field(R::TensorExpression, var::TensorSuperType)
    s = var isa DerivativeSuperType ? isTensor(root_oprnd(var)) : var
    (s isa TensorSuperType && !isempty(tableaux(s))) || return R
    #the field's projector, tableau by tableau, read back as the index list each of its terms puts
    #in the field's slots
    perms = [(1 // 1, s)]
    for tab in tableaux(s)
        proj, norm = compute_young_projector(tab), 1 // compute_hook_product(tab)
        perms = [(c * norm * c2, terms(tt)[1]) for (c, t) in perms for (c2, tt) in summands(apply_young_projector(t, proj))]
    end
    inds = indices(s)
    out = CoeffPair[]
    for (c, t) in perms
        rep = Dict{IndexPair,IndexPair}()
        for (v, u) in zip(inds, indices(t)), p in (upper, lower)
            rep[u[1]=>p] = v[1] => p
        end
        Rk = copy(R)
        set_indices_deep!(Rk, rep)
        append!(out, summands(c * Rk))
    end
    isempty(out) && return zero()
    return merge_sum(TensorExpression(out))
end
