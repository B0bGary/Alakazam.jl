@testset "Trace" begin

    #setup
    if 1 == 1
        trace_st = IndexSet("trace_spacetime", 4, upper)
        @indices μ:τ ∈ trace_st
        trace_fl = IndexSet("trace_flavour", 2, upper, Alakazam.fermionic)
        @indices a:h ∈ trace_fl
        x = Coordinate("x", [trace_st])
        Symbolics.@syms cc hh
    end

    @matrix TrMx #for testing
    A = TrMx("A"); B = TrMx("B"); C = TrMx("C")
    Af = TrMx("A", [a => lower]); Bf = TrMx("B", [b => lower])   # odd-graded matrices
    V(i) = Tensor("V", [i => lower])
    F(i) = Tensor("F", [i => lower])                             # odd, since a:h are fermionic

    @test grading(A) == 0 && grading(Af) == 1 && grading(F(c)) == 1

    ###############################################
    ### linearity and pulling scalars out
    ###############################################
    # a trace of matrices alone has nothing to factor
    @test factor_trace(Trace(A * B)) == TensorExpression(Trace(A * B))
    # numeric and symbolic coefficients come out in front
    @test factor_trace(Trace(2 * A * B)) == 2 * Trace(A * B)
    @test factor_trace(Trace(cc * A * B)) == cc * Trace(A * B)
    # a non-matrix tensor comes out wherever it sat
    @test factor_trace(Trace(V(μ) * A * B)) == V(μ) * Trace(A * B)
    @test factor_trace(Trace(A * V(μ) * B)) == V(μ) * Trace(A * B)
    @test factor_trace(Trace(A * B * V(μ))) == V(μ) * Trace(A * B)
    # linearity across summands, each handled independently
    @test factor_trace(Trace(A * B + cc * B * C)) == Trace(A * B) + cc * Trace(B * C)
    # a summand with no matrices leaves an explicit identity behind rather than an empty term,
    # which would be indistinguishable from a literal zero
    @test factor_trace(Trace(V(μ))) == V(μ) * Trace(IdMatrix())
    @test IdMatrix() isa MatrixSuperType && grading(IdMatrix()) == 0
    @test isempty(indices(IdMatrix()))
    @test IdMatrix() != Mat("𝟙")                      # a distinct type, not just a name
    @test !isZeroTens(TensorExpression(IdMatrix()))    # and not confusable with zero
    # an empty trace stays zero
    @test isZeroTens(factor_trace(Trace(zero())))
    # the identity absorbs into an adjacent matrix
    @test contract_adjacent(TensorExpression(IdMatrix() * A * B), MatrixSuperType) ==
          TensorExpression(A * B)
    @test contract_adjacent(TensorExpression(A * IdMatrix() * B), MatrixSuperType) ==
          TensorExpression(A * B)
    @test contract_adjacent(TensorExpression(A * B * IdMatrix()), MatrixSuperType) ==
          TensorExpression(A * B)

    ###############################################
    ### supertrace
    ###############################################
    # F_c and F_d are odd and Af is odd, so pulling F_d left past Af costs a sign
    pulled = factor_trace(Trace(F(c) * Af * F(d)))
    @test Alakazam.simplify(pulled) ==
          Alakazam.simplify(F(c) * Trace(TensorExpression(Af)) * F(d))
    @test first_coef(pulled) == -1

    # with an even matrix in the middle nothing is crossed, so no sign appears
    @test first_coef(factor_trace(Trace(F(c) * A * F(d)))) == 1

    ###############################################
    ### (graded) cyclicity
    ###############################################
    # even matrices: every rotation is free
    @test Trace(A * B) == Trace(B * A)
    @test Trace(A * B * C) == Trace(B * C * A)
    @test Trace(A * B * C) == Trace(C * A * B)
    @test Alakazam.signed_cyclic_equals(Trace(A * B * C), Trace(C * A * B)) == (true, 1)
    # a permutation that is not a rotation is not matched
    @test Alakazam.signed_cyclic_equals(Trace(A * B * C), Trace(B * A * C)) == (false, 0)
    # different content is not matched either
    @test Trace(A * B) != Trace(A * C)

    # two odd matrices: str(XY) = -str(YX), so they are related but not equal
    @test Alakazam.signed_cyclic_equals(Trace(Af * Bf), Trace(Bf * Af)) == (true, -1)
    @test Trace(Af * Bf) != Trace(Bf * Af)
    # one odd and one even: the rotation is free again
    @test Trace(Af * B) == Trace(B * Af)
    @test Alakazam.signed_cyclic_equals(Trace(Af * B), Trace(B * Af)) == (true, 1)

    # rotations of a term are generated with their signs, and rotating all the way round returns the term itself
    rots = Alakazam.cyclic_rotations(Alakazam.isTensorTerm(TensorExpression(A * B * C)))
    @test length(rots) == 3
    @test [r[1] for r in rots] == [1, 1, 1]
    @test rots[1][2] == Alakazam.isTensorTerm(TensorExpression(A * B * C))

    ###############################################
    ### traces nested in larger expressions
    ###############################################
    # a trace inside a product is reached
    @test factor_trace(V(μ) * Trace(cc * A * B)) == cc * V(μ) * Trace(A * B)
    # and inside another operator's operand
    inner = factor_trace(PartialDerivative("∂", Trace(cc * A * B), x, μ => lower))
    @test inner == PartialDerivative("∂", cc * Trace(A * B), x, μ => lower)
    # something with no trace in it is returned untouched
    @test factor_trace(V(μ) * V(ν)) == V(μ) * V(ν)

    ###############################################
    ### a trace nested under a derivative must not be walked as part of the derivative chain
    ###############################################
    # only partial derivatives may be exchanged by operator index.
    Mv = TrMx("M", [ν => lower]); Nv = TrMx("N", [σ => lower])
    dtr = PartialDerivative("∂", Trace(Mv * Nv), x, μ => lower)

    @test length(Alakazam.get_operator_chain(dtr)) == 2# encoding walks both
    @test length(Alakazam.get_operator_chain(dtr, PartialSuperType)) == 1# reordering does not
    @test Alakazam.symmetrise_partials(dtr) == TensorExpression(dtr)
    @test Alakazam.simplify(dtr) == TensorExpression(dtr)
    @test Alakazam.hard_simplify(dtr) == TensorExpression(dtr)
    @test dtr ≅ dtr

    # nested partials are still exchanged
    Aq = Tensor("A", [ν => lower], x)
    dd = PartialDerivative("∂", PartialDerivative("∂", Aq, x, ρ => lower), x, μ => lower)
    @test length(Alakazam.get_operator_chain(dd, PartialSuperType)) == 2

    ###############################################
    ### trace relations, registered with @mrule and matched up to cyclicity
    ###############################################
    Symbolics.@syms kk
    @matrix TrMg
    @matrix TrM5
    G(i) = TrMg("G", [i => lower])
    @mrule Trace(TrMg("G", [μ => lower])*TrMg("G", [ν => lower])) = 4*Metric([μ => lower, ν => lower])

    # dispatched on type and matched by index, exactly as the reduce_matrix rules are
    @test reduce_trace(Trace(G(μ)*G(ν))) == 4*Metric([μ => lower, ν => lower])
    @test reduce_trace(Trace(G(ρ)*G(σ))) == 4*Metric([ρ => lower, σ => lower])
    # a product with no registered relation is left alone
    @test reduce_trace(Trace(C)) == TensorExpression(Trace(C))
    @test reduce_trace(Trace(A * B)) == TensorExpression(Trace(A * B))

    # rules apply wherever a trace sits
    @test reduce_trace(V(ρ) * Trace(G(μ)*G(ν))) ==
          4 * V(ρ) * Metric([μ => lower, ν => lower])
    @test reduce_trace(Trace(G(μ)*G(ν)) + Trace(G(ν)*G(μ))) ==
          8 * Metric([μ => lower, ν => lower])
    @test reduce_trace(PartialDerivative("∂", Trace(G(μ)*G(ν)), x, ρ => lower)) ==
          TensorExpression(PartialDerivative("∂", 4*Metric([μ => lower, ν => lower]), x, ρ => lower))

    # cyclicity when the factors are of different types, so that the registered
    # ordering is not the one dispatch sees.
    @mrule Trace(TrM5("P")*TrMg("G", [μ => lower])*TrMg("G", [ν => lower])) =
        kk*Metric([μ => lower, ν => lower])
    P = TrM5("P")
    want = kk*Metric([μ => lower, ν => lower])
    @test reduce_trace(Trace(P*G(μ)*G(ν))) == want
    @test reduce_trace(Trace(G(μ)*G(ν)*P)) == want
    @test reduce_trace(Trace(G(ν)*P*G(μ))) == want
    # and a product matching no rotation is untouched
    @test reduce_trace(Trace(A * B * C)) == TensorExpression(Trace(A * B * C))

    # The identity is transparent under a trace: it is stripped before the registered relations are
    # matched, so it cannot break up a run that dispatch would otherwise have to see contiguously.
    Id = IdMatrix(spinor_indexset(trace_st))
    tr4(x) = Alakazam.simplify(eliminate_trace(Trace(TensorExpression(x))))
    @test tr4(Id) == TensorExpression(4)
    @test tr4(Id * Id * Id) == TensorExpression(4)
    @test tr4(Id * G(μ) * G(ν)) == 4 * Metric([μ => lower, ν => lower])
    @test tr4(G(μ) * Id * G(ν)) == 4 * Metric([μ => lower, ν => lower])
    @test tr4(G(μ) * G(ν) * Id) == 4 * Metric([μ => lower, ν => lower])
    # an identity carrying no index set has no dimension to report, so the trace stands
    @test tr4(IdMatrix()) == TensorExpression(Trace(IdMatrix()))
end