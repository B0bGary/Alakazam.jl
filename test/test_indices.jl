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

end
