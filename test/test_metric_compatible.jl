###############################################################################
# Identities between metrics and nested derivatives, each with its answer worked out by hand, for
# ∂ and ∇, under both default positions, in every kind of metric: constant, varying in what the
# derivatives act on, varying in something else, and compatible with ∇.
#   N: a derivative's own level (opposite the default position), M: the other.
#   P: the metric passes through the derivative. C: the derivatives commute with each other (∂ does;
#   ∇ has no commutator declared, so must not be taken to).
###############################################################################
function _metric_derivative_matrix(conv, default, setting, P_of)
    N, M = op_pos(default), default
    st = SET_REGISTRY["mx_st_$conv"]
    a, b, c = (get_index(n * conv, st) for n in ("a", "b", "c"))
    x = COORD_REGISTRY["mx_x$conv"]
    u = Tensor("u", x)
    A(i, l) = Tensor("A", [i => l], x)
    g(l, i, j) = Metric([i => l, j => l])
    same(l, r) = isZeroTens(simplify(l - r))
    ∂(i, l, T) = PartialDerivative("∂", T, [x], [i => l])
    ∇(i, l, T) = Derivative("∇", T, [x], [i => l])
    @testset "$conv default $(default == upper ? "up" : "down"), $setting, $dn" for (dn, D, C) in (("∂", ∂, true), ("∇", ∇, false))
        P = P_of(dn)
        @test same(g(M, a, b) * D(a, N, D(b, N, u)), D(a, M, D(a, N, u)))                     #g^{ab}D_aD_b u = D^aD_a u
        @test same(g(M, a, b) * D(a, N, D(b, N, u)), D(a, N, D(a, M, u))) == P                #  = D_aD^a u
        @test same(D(a, N, D(a, M, u)), D(a, M, D(a, N, u))) == P
        @test same(D(a, N, D(b, N, u)), D(b, N, D(a, N, u))) == C
        @test same(D(a, M, D(b, M, u)), D(b, M, D(a, M, u))) == (P && C)
        @test same(D(a, N, g(M, a, b) * D(b, N, u)), D(a, N, D(a, M, u)))                     #D_a(g^{ab}D_b u) = D_aD^a u
        @test same(g(N, a, b) * D(a, M, D(b, M, u)), D(a, N, D(a, M, u)))                     #g_{ab}D^aD^b u = D_aD^a u
        @test same(D(a, N, A(a, M)), D(a, M, A(a, N))) == P
        @test same(D(a, N, D(b, N, A(b, M))), D(b, N, D(a, N, A(b, M)))) == C
        @test same(D(a, N, D(a, M, D(b, N, u))), D(b, N, D(a, N, D(a, M, u)))) == (P && C)
        #related by relabelling and the symmetry of g alone, so no commutation is needed
        @test same(D(a, N, D(b, N, D(a, M, D(b, M, u)))), D(a, M, D(b, M, D(a, N, D(b, N, u))))) == P
        @test isZeroTens(apply_derivative(TensorExpression(D(c, N, TensorExpression(g(N, a, b)))))) == P
        @test same(D(a, N, g(M, a, b) * A(b, N)), g(M, a, b) * D(a, N, A(b, N))) == P
        #a pair running into a product inside the operand: seen only where the metric passes the derivative
        @test same(D(a, N, u * D(a, M, u)), D(a, M, u * D(a, N, u))) == P
        @test same(D(a, N, D(a, M, u)), D(b, N, D(b, M, u)))
        #a metric outside contracts into the operand only if it passes the derivative
        @test same(g(M, a, b) * D(c, N, A(b, N)), D(c, N, A(a, M))) == P
        #already in canonical form: simplifying again changes nothing
        let s1 = simplify(D(a, N, D(b, N, D(a, M, D(b, M, u)))))
            @test isZeroTens(simplify(s1 - simplify(s1)))
        end
    end
    @test !same(∂(a, N, ∇(b, N, u)), ∇(b, N, ∂(a, N, u)))                                     #no commutator declared
end

for (conv, default) in (("U", upper), ("L", lower))
    st = IndexSet("mx_st_$conv", 4, default, Metric)
    Index(["a$conv", "b$conv", "c$conv"], st)
    x = Coordinate("mx_x$conv", st)
    y = Coordinate("mx_y$conv", IndexSet("mx_aux_$conv", 3, upper, Metric))
    metric_function!(st)
    _metric_derivative_matrix(conv, default, "constant", D -> true)
    metric_function!(st, x)
    _metric_derivative_matrix(conv, default, "g(x)", D -> false)
    metric_function!(st, y)
    _metric_derivative_matrix(conv, default, "g(y), derivatives in x", D -> true)
end

@compatible Derivative Metric

for (conv, default) in (("U", upper), ("L", lower))
    metric_function!(SET_REGISTRY["mx_st_$conv"], COORD_REGISTRY["mx_x$conv"])
    _metric_derivative_matrix(conv, default, "g(x), ∇ compatible", D -> D == "∇")
end

@testset "metric compatible derivatives" begin
    mst = IndexSet("mc_spacetime", 4, lower, Metric)
    mμ, mν, mρ, mσ = Index(["mμ", "mν", "mρ", "mσ"], mst)
    mx = Coordinate("mc_x", mst)
    metric_function!(mst, mx)
    mets(e) = [t for t in terms(summands(TensorExpression(e))[1][2]) if t isa Metric]

    A = Tensor("A", [mμ => upper], mx)
    g = Metric(mμ => lower, mν => lower)
    #at their own level, upper for this set (opposite its default), so they carry no metric of their own
    ∇A = Derivative("∇", A, [mx], [mν => upper])
    ∂A = PartialDerivative("∂", A, [mx], [mν => upper])

    #the macro defines the function, for the declared types and their subtypes only
    @test metric_compatible(Derivative, Metric) && !metric_compatible(PartialDerivative, Metric)
    @test metric_compatible(∇A, g) && Alakazam.commutes_with_metric(∇A, g)
    #the wrong way round is refused, and defines nothing
    @test_throws ArgumentError eval(:(@compatible Metric Derivative))
    @test !metric_compatible(Metric, Derivative)
    @test_throws LoadError eval(:(@compatible Derivative))              #refused at expansion

    #popping: ∇ annihilates g although g depends on x, so g comes out
    @test function_of(only(mets(pop_metric(∇A)))) == [mx]
    @test isempty(mets(operand(terms(summands(TensorExpression(pop_metric(∇A)))[1][2])[end])))
    #declared for ∇ only: ∂ still keeps it
    @test !metric_compatible(∂A, g) && isempty(mets(pop_metric(∂A)))
    #∇ carries it out of itself but not out of a ∂ inside it
    let p = isTensor(pop_metric(Derivative("∇", ∂A, [mx], [mρ => upper])))
        @test p isa Derivative && isempty(mets(operand(p)))
        @test !isempty(mets(operand(isTensor(operand(p)))))
    end

    #applying: ∇g = 0, so ∇(gA) = g∇A, while ∂g survives
    @test isZeroTens(apply_derivative(Derivative("∇", TensorExpression(g), [mx], [mρ => lower])))
    @test !isZeroTens(apply_derivative(PartialDerivative("∂", TensorExpression(g), [mx], [mρ => lower])))
    gA = Metric(mμ => lower, mσ => lower) * A
    @test apply_derivative(Derivative("∇", gA, [mx], [mν => lower])) ==
          TensorExpression(Metric(mμ => lower, mσ => lower) * Derivative("∇", A, [mx], [mν => lower]))
    @test length(summands(apply_derivative(PartialDerivative("∂", gA, [mx], [mν => lower])))) == 2

    #a derivative's own raised index hides a metric, ∇^a = g^{ab}∇_b, and ∇ being compatible lets it
    #through: g^{ab}∇_a∇_b u, ∇_a∇^a u and ∇^a∇_a u are one thing although g depends on x
    u = Tensor("u", mx)
    gab = Metric([mμ => upper, mν => upper])
    @test isZeroTens(simplify(gab * Derivative("∇", Derivative("∇", u, [mx], [mν => lower]), [mx], [mμ => lower]) -
                              Derivative("∇", Derivative("∇", u, [mx], [mμ => upper]), [mx], [mμ => lower])))
    @test isZeroTens(simplify(Derivative("∇", Derivative("∇", u, [mx], [mμ => upper]), [mx], [mμ => lower]) -
                              Derivative("∇", Derivative("∇", u, [mx], [mμ => lower]), [mx], [mμ => upper])))
    #while ∂ is not compatible, so the same with partials is not
    @test !isZeroTens(simplify(gab * PartialDerivative("∂", PartialDerivative("∂", u, [mx], [mν => lower]), [mx], [mμ => lower]) -
                               PartialDerivative("∂", PartialDerivative("∂", u, [mx], [mμ => upper]), [mx], [mμ => lower])))

    #and out of a ∂ it was not in: a flat set's metric leaves ∂ as before
    metric_function!(mst)
    @test !isempty(mets(pop_metric(PartialDerivative("∂", A, [mx], [mν => lower]))))
end

Base.delete_method(which(metric_compatible, (Type{Derivative}, Type{Metric})))
@test !metric_compatible(Derivative, Metric)

@testset "metrics into derivative operands" begin
    ost = IndexSet("op_st", 4, upper)
    ofl = IndexSet("op_fl", 2, EpsilonTensor, upper)
    oa, ob, oc, od, oe = Index(["oa", "ob", "oc", "od", "oe"], ost)
    oi, oj = Index(["oi", "oj"], ofl)
    ox = Coordinate("op_x", ost)
    T(n, ps...) = Tensor(n, collect(ps), ox)
    ∂c(X) = PartialDerivative("∂", X, [ox], [oc => lower])
    g = Metric([oa => lower, ob => lower])
    omets(e) = [t for t in terms(summands(TensorExpression(e))[1][2]) if t isa Alakazam.MetricSuperType]
    expands_to(l, r) = isZeroTens(simplify(apply_derivative(TensorExpression(l)) - apply_derivative(TensorExpression(r))))

    #a constant metric passes ∂ into a product, not only into a single tensor
    @test contract_metrics(g * ∂c(T("X", oa => upper) * T("Y", oc => upper))) == TensorExpression(∂c(T("X", ob => lower) * T("Y", oc => upper)))
    #ε's sign goes in with it, summand by summand when the summands take different ends
    @test simplify(EpsilonTensor([oi => lower, oj => lower]) * ∂c(T("Θ", oi => upper) * T("Y", oc => upper))) ==
          simplify(-∂c(T("Θ", oj => lower) * T("Y", oc => upper)))
    εsum = EpsilonTensor([oi => lower, oj => lower]) * ∂c(T("Θ", oi => upper) * T("Φ", oj => upper) * T("Y", oc => upper) +
                                                       T("Φ", oj => upper) * T("Θ", oi => upper) * T("Y", oc => upper))
    @test isempty(omets(simplify(εsum))) && expands_to(simplify(εsum), -2 * ∂c(T("Θ", oj => lower) * T("Φ", oj => upper) * T("Y", oc => upper)))
    #a summand whose own dummy has the name the metric brings in has it moved, onto one the term does not use
    clash = g * ∂c(T("X", oa => upper) * T("X", od => lower) * T("Z", od => upper) + T("X", oa => upper) * T("W", ob => lower) * T("V", ob => upper))
    @test isempty(omets(simplify(clash))) &&
          expands_to(simplify(clash), ∂c(T("X", ob => lower) * T("X", od => lower) * T("Z", od => upper) + T("X", ob => lower) * T("W", oe => lower) * T("V", oe => upper)))
    #any other operator is opaque: what it does need not commute with raising an index
    @test !isempty(omets(contract_metrics(g * Alakazam.Operator("O", T("X", oa => upper)))))
    #and a metric that ∂ does not pass stays outside, product or not
    metric_function!(ost, ox)
    @test !isempty(omets(contract_metrics(Metric([oa => lower, ob => lower]) * ∂c(T("X", oa => upper) * T("Y", oc => upper)))))
    metric_function!(ost)
end
