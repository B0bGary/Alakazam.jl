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
