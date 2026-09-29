@testset "Index generation" begin

    #names are generated when an operation needs a dummy and the set has none
    #spare, so these surface inside simplify rather than at construction

    ###############################################
    #own pool is still used first, so nothing changes when one alphabet is enough
    small = IndexSet("gen_small", 4, lower)
    Index(["a", "b"], small)
    nxt = Alakazam.generate_indices(small)
    @test nxt.name ∈ Alakazam.LATIN_LOWER
    @test nxt.name ∉ ["a", "b"]
    @test nxt.name == "c"

    ###############################################
    #an exhausted alphabet falls through to the next pool
    full = IndexSet("gen_full", 4, lower)
    Index(copy(Alakazam.LATIN_LOWER), full)
    @test length(full.indices) == length(Alakazam.LATIN_LOWER)

    beyond = Alakazam.generate_indices(full)
    @test beyond isa Alakazam.IndexSuperType
    @test beyond.name ∉ Alakazam.LATIN_LOWER
    @test beyond.name ∉ [i.name for i in full.indices[1:end-1]]

    #and it keeps going, so the supply is not one alphabet
    for _ in 1:30
        @test_nowarn Alakazam.generate_indices(full)
    end
    names = [i.name for i in full.indices]
    @test length(names) == length(unique(names))
    @test length(names) > length(Alakazam.LATIN_LOWER)

    ###############################################
    #a template name in no pool. the Index constructor accepts one, so
    #generation has to cope rather than raise
    odd = IndexSet("gen_odd", 4, lower)
    Index(["i1", "i2"], odd)
    @test all(n -> !any(p -> n in p, Alakazam.ORDERED_POOLS), ["i1", "i2"])

    fresh = Alakazam.generate_indices(odd)
    @test fresh isa Alakazam.IndexSuperType
    @test fresh.name ∉ ["i1", "i2"]

    ###############################################
    #every pool used up returns nothing, which is what the callers all test for
    #with `nid === nothing && throw(OutOfIndices(...))`
    every = IndexSet("gen_every", 4, lower)
    Index(unique(reduce(vcat, Alakazam.ORDERED_POOLS)), every)
    @test Alakazam.generate_indices(every) === nothing

    #so safe_pop! reaches the OutOfIndices it has always tested for
    drained = copy(every.indices)
    empty!(drained)
    @test_throws Alakazam.OutOfIndices safe_pop!(drained, every)

    ###############################################
    #safe_pop!, which is the path the simplifier takes
    popset = IndexSet("gen_pop", 4, lower)
    Index(copy(Alakazam.LATIN_LOWER), popset)
    pool = copy(popset.indices)
    popped = [safe_pop!(pool, popset) for _ in 1:(length(Alakazam.LATIN_LOWER) + 12)]
    @test length(popped) == length(Alakazam.LATIN_LOWER) + 12
    @test length(unique(p.name for p in popped)) == length(popped)

    ###############################################
    #an expression with more contracted pairs than one alphabet holds
    # `let` so the loop updates `expr` rather than shadowing it: at module top level (as when this
    # file is used as a precompile workload) a bare `expr = 1` is a global and the for-loop body
    # would introduce a new local instead.
    let
        wide = IndexSet("gen_wide", 4, lower)
        widx = Index(copy(Alakazam.LATIN_LOWER), wide)
        expr = 1
        for k in 1:12
            p, q = widx[2k-1], widx[2k]
            expr = expr * SymmetricTensor("A", [p => upper, q => upper]) *
                   SymmetricTensor("B", [p => lower, q => lower])
        end
        @test_nowarn Alakazam.simplify(expr)
    end


    ###############################################
    ### raise_all / lower_all / flip_all: free indices only, dummies preserved
    ###############################################
    pset = IndexSet("posn_test", 4, upper)
    pm, pn = Index("pm", pset), Index("pn", pset)
    A1 = Tensor("A", [pm => upper])
    B1 = Tensor("B", [pm => lower, pn => lower])      # pm is a dummy, pn is free
    trm = A1 * B1

    @test Set(Alakazam.name.(dummies(trm))) == Set(["pm"])
    # the free index moves, the dummy pair does not
    @test indices(terms(raise_all(trm))[2]) == [pm => lower, pn => upper]
    @test indices(terms(lower_all(trm))[2]) == [pm => lower, pn => lower]
    @test indices(terms(flip_all(trm))[2])  == [pm => lower, pn => upper]
    # op_pos on a whole expression is an alias for flip_all
    @test op_pos(trm) == flip_all(trm)
    # ... but remains the primitive on a bare position
    @test op_pos(upper) == lower && op_pos(lower) == upper
    # flip_all is its own inverse
    @test flip_all(flip_all(trm)) == trm

    # building M_{ij} from M^{ij} so the two contract
    Mup = Tensor("M", [pm => upper, pn => upper])
    @test indices(flip_all(Mup)) == [pm => lower, pn => lower]
    @test length(free_indices(Alakazam.simplify(flip_all(Mup) * Mup))) == 0

    # a contracted spinor pair is left alone: under an antisymmetric metric reversing it would
    # change the sign, and Alakazam does not relate the two forms
    sset = IndexSet("posn_spin", 2, upper, EpsilonTensor, Alakazam.fermionic)
    sa = Index("sa", sset)
    spair = Tensor("psi", [sa => upper]) * Tensor("chi", [sa => lower])
    @test flip_all(spair) == spair

    @testset "an index set with no metric" begin
        #An algebra with no quadratic invariant - SU(N≥3), whose ε is invariant but not a metric -
        #gives a set with no metric. Its indices cannot be raised or lowered, and whatever would need
        #to either declines or says why, rather than failing somewhere else.
        #not `nm`: the precompile workload replays these files at top level, where that name would
        #become a global that test_latex.jl's loop variable then collides with
        nomet = IndexSet("nometric_su3", upper, SLAlgebra(3))
        na, nb, nc = Index(["na", "nb", "nc"], nomet)
        @test nomet.metric_type === nothing

        #asked for directly, a metric is refused with a reason rather than handed back as nothing
        @test_throws ArgumentError Metric(na => upper, nb => upper)
        @test_throws ArgumentError Metric([na => upper, nb => upper])

        #there is nothing to pop an index with, so it stays where it is
        Tl = Tensor("T", [na => lower])
        @test pop_metric(Tl) == TensorExpression(Tl)
        #and T^a S_a and T_a S^a are different contractions, not one term written twice
        @test length(summands(simplify(Tensor("T", [na => upper]) * Tensor("S", [na => lower]) -
                                       Tensor("T", [na => lower]) * Tensor("S", [na => upper])))) == 2

        #the invariant ε still contracts, although no metric does; its sign used to test
        #nothing <: AntisymMetricSuperType, a TypeError
        @test !isZeroTens(eliminate_epsilon(TensorExpression(
            EpsilonTensor([na => lower, nb => lower, nc => lower]) *
            EpsilonTensor([na => upper, nb => upper, nc => upper]))))

        #flavour, not spacetime: beside a real vector index, :dim is taken from the vector one
        @test !Alakazam._is_spacetime_set(nomet)
        vset = IndexSet("nometric_vec", 7, lower)
        vμ = Index("vμ", vset)
        @test Alakazam._matched_dimension([Tensor("G", [na => upper, vμ => upper])]) == 7

        #σσ = δ + iεσ needs a metric, so with none the product is left alone rather than erroring
        σσ = Pauli(na => upper) * Pauli(nb => upper)
        @test contract_adjacent(σσ, Pauli) == σσ

        #A matrix rule that disagrees with the term on a level cannot bridge it, so it does not match -
        #and must not have relabelled the term on the way out, as the matcher writes in place
        P1, P2 = Mat("X", na => upper), Mat("Y", na => lower)
        A1, A2 = Mat("X", nb => lower), Mat("Y", nb => upper)
        before = string(A1) * string(A2)
        (hit, _) = Alakazam.match_and_apply_rule(Alakazam.TensorSuperType[P1, P2],
            Alakazam.TensorSuperType[A1, A2], TensorExpression(1), A1 * A2)
        @test !hit
        @test string(A1) * string(A2) == before

        #varying by the other level of a field is zero with no metric: they are different fields
        xm = Coordinate("x_nometric", vset)
        L = Tensor("φ", [na => lower], xm) * Tensor("J", [na => upper], xm)
        @test !isZeroTens(simplify(apply_derivative(VariationalDerivative(L, Tensor("φ", [nb => lower])))))
        @test isZeroTens(simplify(apply_derivative(VariationalDerivative(L, Tensor("φ", [nb => upper])))))

        #and so has nothing to depend on a coordinate
        @test_throws ErrorException metric_function!(nomet, xm)
    end

    @testset "a metric depending on coordinates" begin
        cst = IndexSet("mf_spacetime", 4, lower, Metric)
        cμ, cν, cρ = Index(["cμ", "cν", "cρ"], cst)
        cx = Coordinate("mf_x", cst)
        cfl = IndexSet("mf_flat", 3, lower, Metric)
        ca, cb = Index(["ca", "cb"], cfl)
        cy = Coordinate("mf_y", cfl)

        #flat by default
        @test isempty(metric_function(cst))
        η = Metric(cμ => lower, cν => lower)
        @test name(η) == "η" && isempty(function_of(η))

        @test metric_function!(cst, cx) === cst
        @test metric_function(cst) == [cx] && metric_function(cμ) == [cx]
        #every metric the set builds now varies, and is named apart from a constant one
        g = Metric(cμ => lower, cν => lower)
        @test name(g) == "g" && function_of(g) == [cx]
        @test g != η
        #the metrics in the first term of a term or expression
        mets(e) = [t for t in terms(summands(TensorExpression(e))[1][2]) if t isa Metric]
        #popping puts g(x) beside the tensor
        @test any(t -> function_of(t) == [cx], mets(pop_metric(Tensor("A", [cμ => upper], cx))))
        #and a derivative in x no longer moves it out, ∂(gA) ≠ g∂A, so it is popped inside, next to A.
        #The derivatives here sit at their own level, upper for this set (opposite its default), so
        #that they carry no metric of their own to pop
        cA = Tensor("A", [cμ => upper], cx)
        ∂A = PartialDerivative("∂", cA, [cx], [cν => upper])
        let p = isTensor(pop_metric(∂A))
            @test p isa PartialDerivative && function_of(mets(operand(p))[1]) == [cx]
            @test isempty(mets(pop_metric(∂A)))                     #nothing outside it
            @test Set(free_indices(p)) == Set(free_indices(∂A))
            @test all(t -> t isa Metric || in_default_pos(indices(t)[1]), terms(summands(operand(p))[1][2]))
        end
        @test isZeroTens(simplify(∂A - ∂A)) && isZeroTens(simplify(pop_metric(∂A) - ∂A))
        #nested, it stays in the innermost derivative
        let p = isTensor(pop_metric(PartialDerivative("∂", ∂A, [cx], [cρ => upper])))
            @test isempty(mets(operand(p))) && !isempty(mets(operand(isTensor(operand(p)))))
        end
        #with several summands, each is popped where it stands
        let p = isTensor(pop_metric(PartialDerivative("∂", cA + Tensor("C", [cμ => upper], cx), [cx], [cρ => upper])))
            @test p isa PartialDerivative && all(s -> any(t -> t isa Metric, terms(s[2])), summands(operand(p)))
        end
        #undeclared, a covariant derivative is not compatible, so it keeps the metric too
        #(declaring one is tested in test_metric_compatible.jl, which is not replayed at precompile)
        ∇A = Derivative("∇", cA, [cx], [cν => upper])
        @test !metric_compatible(∇A, g) && isTensor(pop_metric(∇A)) isa Derivative
        #at the default level a covariant derivative hides a metric of its own, ∇_ν = g_{νρ}∇^ρ, which
        #is exact with the metric outside, so that one does come out (a partial's chain accounts for it)
        let p = pop_metric(Derivative("∇", cA, [cx], [cν => lower]))
            @test length(mets(p)) == 1 && Set(free_indices(p)) == Set(free_indices(Derivative("∇", cA, [cx], [cν => lower])))
        end
        #copy on write: popping and contracting leave what they were given as it was
        let D = Derivative("∇", Tensor("A", [cμ => upper], cx) * Tensor("B", [cρ => upper], cx), [cx], [cν => lower])
            oi, ii, s = copy(op_indices(D)), copy(indices(D)), string(D)
            pop_metric(D)
            Alakazam.contract_metrics(TensorExpression(Metric(cν => upper, Index("cσ", cst) => upper) * D), Metric)
            @test op_indices(D) == oi && indices(D) == ii && string(D) == s
        end
        #nor treats it as a constant
        @test !isZeroTens(apply_derivative(PartialDerivative("∂", TensorExpression(g), [cx], [cρ => lower])))

        #another set is untouched, and its metric still pops out of its own derivative
        @test name(Metric(ca => lower, cb => lower)) == "η"
        dB = PartialDerivative("∂", Tensor("B", [ca => upper], cy), [cy], [cb => upper])
        @test !isempty(mets(pop_metric(dB)))

        #a partial's own level is opposite to the default position: cst is lower by default, so ∂^μ.
        #At the default position it hides a varying metric, ∂_μ = g_{μν}∂^ν, and does not swap
        cφ = Tensor("φ", cx)
        ∂mf(i, T) = PartialDerivative("∂", T, [cx], [i])
        @test !isZeroTens(simplify(∂mf(cμ => lower, ∂mf(cμ => upper, cφ)) - ∂mf(cμ => upper, ∂mf(cμ => lower, cφ))))
        @test !isZeroTens(simplify(∂mf(cμ => upper, ∂mf(cν => lower, cφ)) - ∂mf(cν => lower, ∂mf(cμ => upper, cφ))))
        @test !isZeroTens(simplify(∂mf(cμ => lower, ∂mf(cν => lower, cφ)) - ∂mf(cν => lower, ∂mf(cμ => lower, cφ))))
        #at their own level they always do
        @test isZeroTens(simplify(∂mf(cμ => upper, ∂mf(cν => upper, cφ)) - ∂mf(cν => upper, ∂mf(cμ => upper, cφ))))
        #and relabelling a dummy is not a swap
        @test isZeroTens(simplify(∂mf(cμ => lower, ∂mf(cμ => upper, cφ)) - ∂mf(cν => lower, ∂mf(cν => upper, cφ))))
        #the swap rule itself: one partial hiding g(x) blocks it, two at their own level do not
        @test !Alakazam._partials_swap(∂mf(cμ => lower, cφ), ∂mf(cν => upper, cφ))
        @test !Alakazam._partials_swap(∂mf(cμ => lower, cφ), ∂mf(cν => lower, cφ))
        @test Alakazam._partials_swap(∂mf(cμ => upper, cφ), ∂mf(cν => upper, cφ))
        #so a chain is reordered past only what it may swap with: ∂^μ∂^ρ may, ∂_ν between them may not
        @test isZeroTens(simplify(∂mf(cμ => upper, ∂mf(cρ => upper, ∂mf(cν => lower, cφ))) -
                                  ∂mf(cρ => upper, ∂mf(cμ => upper, ∂mf(cν => lower, cφ)))))
        @test !isZeroTens(simplify(∂mf(cμ => upper, ∂mf(cν => lower, ∂mf(cρ => upper, cφ))) -
                                   ∂mf(cρ => upper, ∂mf(cν => lower, ∂mf(cμ => upper, cφ)))))
        #a derivative in another coordinate leaves the metric alone, so the swap is allowed
        @test Alakazam._partials_swap(∂mf(cμ => lower, cφ), PartialDerivative("∂", cφ, [cy], [ca => lower]))

        #a derivative's own raised index is a metric too, ∇^a = g^{ab}∇_b. With g(x) and ∇ not
        #declared compatible, ∇_a∇^a u = ∇_a(g^{ab}∇_b u) differs from g^{ab}∇_a∇_b u by ∇_a g^{ab}
        ∇mf(i, T) = Derivative("∇", T, [cx], [i])
        g2 = Metric([cμ => upper, cν => upper])
        @test !isZeroTens(simplify(g2 * ∇mf(cμ => lower, ∇mf(cν => lower, cφ)) - ∇mf(cμ => lower, ∇mf(cμ => upper, cφ))))
        #but with the metric outside every derivative the two are the same
        @test isZeroTens(simplify(g2 * ∇mf(cμ => lower, ∇mf(cν => lower, cφ)) - ∇mf(cμ => upper, ∇mf(cμ => lower, cφ))))

        #declaring no coordinates makes it constant again
        metric_function!(cst)
        #and then the metric passes through, so all three are the same
        @test isZeroTens(simplify(g2 * ∇mf(cμ => lower, ∇mf(cν => lower, cφ)) - ∇mf(cμ => lower, ∇mf(cμ => upper, cφ))))
        @test isZeroTens(simplify(∇mf(cμ => lower, ∇mf(cμ => upper, cφ)) - ∇mf(cμ => upper, ∇mf(cμ => lower, cφ))))
        #and with a constant metric every level swaps
        @test Alakazam._partials_swap(∂mf(cμ => lower, cφ), ∂mf(cν => upper, cφ))
        @test isZeroTens(simplify(∂mf(cμ => lower, ∂mf(cμ => upper, cφ)) - ∂mf(cμ => upper, ∂mf(cμ => lower, cφ))))
        @test isempty(metric_function(cst)) && name(Metric(cμ => lower, cν => lower)) == "η"
    end

end
