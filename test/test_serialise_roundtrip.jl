# A type the serialiser has never seen
struct _RoundTripThing
    label::String
    n::Int
    where_it_lives::IndexSet
end

###############################################################################
# Self-contained save -> clear -> load cycles, asserting SEMANTICS rather than presence.
#
# The commutator key motivates several of them: it is now (name, concrete type, index signature),
# so a relation only resolves after loading if serialisation restores the operand's exact concrete
# type. Reconstructing a GammaM as a Mat would leave every relation silently unfindable, and
# nothing previously would have noticed.
###############################################################################

@testset "Serialisation round trip" begin

    _RT = joinpath(tempdir(), "alakazam_roundtrip_test")

    ###########################################################################
    # index sets and their properties
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    st  = IndexSet("spacetime", 4, lower, Metric)
    sp  = IndexSet("spinor", 4, lower, fermionic, SymplecticForm)
    su2 = IndexSet("su2", 2, upper, EpsilonTensor)
    μ, ν, ρ = Index(["μ", "ν", "ρ"], st)
    α, β    = Index(["α", "β"], sp)
    i, j    = Index(["i", "j"], su2)
    x = Coordinate("x", st)

    save(TensorExpression(Tensor("Q", [μ => lower], x)), _RT, true)
    clear_global!()
    _ = @load _RT

    @test haskey(SET_REGISTRY, "spacetime") && haskey(SET_REGISTRY, "spinor")
    let s = SET_REGISTRY["spacetime"], p = SET_REGISTRY["spinor"], u = SET_REGISTRY["su2"]
        @test dimension(s) == 4 && dimension(p) == 4 && dimension(u) == 2
        @test s.metric_type <: SymMetricSuperType        # a genuine metric, not an epsilon
        @test p.metric_type <: AntisymMetricSuperType    # symplectic form survived
        @test p.grading == fermionic.grading             # grading is not lost
        @test s.grading == bosonic.grading
    end

    ###########################################################################
    # concrete types, which the commutator key depends on
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    st2 = IndexSet("spacetime", 4, lower, Metric)
    sp2 = IndexSet("spinor", 4, lower, fermionic, SymplecticForm)
    a, b = Index(["a", "b"], st2)
    γi   = Index("γi", sp2)
    y    = Coordinate("y", st2)

    objs = Dict("mat"    => Mat("M", [a => upper]),
                "spinor" => Spinor("S", y),
                "gam1"   => Gamma([a => upper]),
                "gam2"   => Gamma([a => upper, b => upper]),
                "gamT"   => Gamma([a => upper, γi => lower, Index("δi", sp2) => lower]),
                "ket"    => Ket("K", [a => lower]),
                "bra"    => Bra("B", [a => upper]))
    for (nm, o) in objs
        save(TensorExpression(o), _RT * "_" * nm, true)
    end
    typenames = Dict(nm => nameof(typeof(o)) for (nm, o) in objs)

    clear_global!()
    for (nm, want) in typenames
        got = @load (_RT * "_" * nm)
        t = terms(first_summand(got))[1]
        @test nameof(typeof(t)) === want
    end

    ###########################################################################
    # number fields survive, including the complex default on matrices
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    st3 = IndexSet("spacetime", 4, lower, Metric)
    c = Index("c", st3)
    save(TensorExpression(Mat("Mc", [c => upper]) * Tensor("Tr", [c => lower])), _RT, true)
    clear_global!()
    let back = @load _RT
        fs = terms(first_summand(back))
        m = only(t for t in fs if name(t) == "Mc")
        r = only(t for t in fs if name(t) == "Tr")
        @test numberfield(m) isa ℂ       # matrices default complex
        @test numberfield(r) isa ℝ       # plain tensors do not
    end

    ###########################################################################
    # commutators: the relation must RESOLVE after loading, not merely be present
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    st4 = IndexSet("spacetime", 4, lower, Metric)
    p, q, r = Index(["p", "q", "r"], st4)
    A  = Mat("A", [p => upper]); B = Mat("B", [q => upper]); D = Mat("D", [r => upper])
    G1 = Gamma([p => upper])                      # rank 1
    G2 = Gamma([p => upper, q => upper])          # rank 2, SAME NAME, distinct key
    @commutator [A, B] = zero()
    @commutator [G1, D] = zero()
    @anticommutator {G2, D} = zero()
    keys_before = Set(keys(Alakazam.commutators))

    save(TensorExpression(A), _RT, true)
    clear_global!()
    _ = @load _RT

    let s = SET_REGISTRY["spacetime"]
        p2, q2, r2 = Index(["p", "q", "r"], s)
        A2 = Mat("A", [p2 => upper]); B2 = Mat("B", [q2 => upper]); D2 = Mat("D", [r2 => upper])
        G1b = Gamma([p2 => upper]); G2b = Gamma([p2 => upper, q2 => upper])

        @test has_commutation_rel(A2, B2)
        @test Alakazam.can_swap(A2, B2)
        # the two same-named gammas kept SEPARATE entries through the round trip
        @test has_commutation_rel(G1b, D2)
        @test has_commutation_rel(G2b, D2)
        @test get_commutator(G1b, D2)[1] ==  1     # stayed a commutator
        @test get_commutator(G2b, D2)[1] == -1     # stayed an anticommutator
        @test Alakazam.can_swap(G1b, D2)
        @test !Alakazam.can_swap(G2b, D2)
        # keys are restored identically, not just repopulated with something equivalent
        @test keys_before ⊆ Set(keys(Alakazam.commutators))
    end

    ###########################################################################
    # symmetries, coefficients and operators
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    st5 = IndexSet("spacetime", 4, lower, Metric)
    u, v, w = Index(["u", "v", "w"], st5)
    z = Coordinate("z", st5)
    @syms ħ
    # opposite positions so u,v form genuine dummy pairs, and a spare w on both summands so the
    # sum is well formed
    sym  = SymmetricTensor("Sy", [u => lower, v => lower], z)
    anti = AntisymmetricTensor("An", [u => upper, v => upper], z)
    der  = PD("∂", Tensor("f", z), z, [w => lower])
    expr = (2 * ħ) * TensorExpression(sym * anti * Tensor("g", [w => lower], z)) +
           TensorExpression(der)

    save(expr, _RT, true)
    clear_global!()
    let back = @load _RT
        @test length(summands(back)) == length(summands(expr))
        allt = [t for (_, tt) in summands(back) for t in terms(tt)]
        sy = only(t for t in allt if name(t) == "Sy")
        an = only(t for t in allt if name(t) == "An")
        @test isSymmetric(sy)
        @test isAntisymmetric(an)
        @test any(t isa DerivativeSuperType for t in allt)     # the operator came back as one
    end

    ###########################################################################
    # fields OUTSIDE (data, indices), which the tensor serialisation path drops
    ###########################################################################
    # `to_dict(::TensorSuperType)` stores only "data", and `reconstruct` calls T(data, indices).
    # A type carrying anything else loses it
    clear_global!(); Alakazam.clear_commutators!()
    st6 = IndexSet("spacetime", 4, lower, Metric)
    sp6 = IndexSet("spinor", 2, lower, fermionic, EpsilonTensor)
    save(TensorExpression(IdMatrix(st6)), _RT, true)
    clear_global!()
    let back = @load _RT
        idm = terms(first_summand(back))[1]
        @test idm isa IdMatrix
        @test idm.index_set !== nothing                      # the space must survive
        @test idm.index_set.name == "spacetime"
    end

    ###########################################################################
    # non-default signature
    ###########################################################################
    # `signature` is a Tuple{Int,Int}, which JSON has no representation for; its `from_field`
    # handler is commented out, so it round trips through the generic vector path.
    clear_global!(); Alakazam.clear_commutators!()
    st7 = IndexSet("st", 4, lower, Metric, (2, 2))
    @test st7.signature == (2, 2)
    save(TensorExpression(Tensor("W", [Index("w", st7) => lower])), _RT, true)
    clear_global!()
    _ = @load _RT
    @test SET_REGISTRY["st"].signature == (2, 2)

    ###########################################################################
    # conjugate index sets, which reference each other
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    ca = IndexSet("sp",  2, lower, fermionic, EpsilonTensor)
    cb = IndexSet("spb", 2, lower, fermionic, EpsilonTensor)
    conjugates!(ca, cb)
    @test ca.conjugate === cb && cb.conjugate === ca
    save(TensorExpression(Tensor("V", [Index("v1", ca) => lower])), _RT, true)
    clear_global!()
    _ = @load _RT
    let a = SET_REGISTRY["sp"], b = SET_REGISTRY["spb"]
        @test a.conjugate !== nothing && b.conjugate !== nothing
        @test a.conjugate.name == "spb" && b.conjugate.name == "sp"
        @test a.conjugate === b && b.conjugate === a     # the registry objects, not copies
    end

    ###########################################################################
    # a metric depending on coordinates, which are built after the sets
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    let cs = IndexSet("curvedst", 4, lower, Metric)
        cx = Coordinate("cx", cs)
        metric_function!(cs, cx)
        save(TensorExpression(Tensor("V", [Index("cv", cs) => lower])), _RT, true)
    end
    clear_global!()
    _ = @load _RT
    let back = SET_REGISTRY["curvedst"]
        @test [c.name for c in metric_function(back)] == ["cx"]
        @test metric_function(back)[1] === COORD_REGISTRY["cx"]   # the registry object, not a copy
    end
    # a file written before the field existed loads flat
    let d = Alakazam.to_dict(SET_REGISTRY["curvedst"])
        delete!(d, "metric_function")
        clear_global!()
        #through JSON, as a file would be: the raw dict still holds metric_type as a Type
        @test isempty(metric_function(Alakazam.from_dict(Alakazam.JSON.parse(Alakazam.JSON.json(d)))))
    end

    ###########################################################################
    # an index-set ALGEBRA, whose signature lives in TYPE PARAMETERS
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    sa = IndexSet("spin5", 4, lower, Spin(4, 1))
    @test sa.algebra !== nothing
    save(TensorExpression(Tensor("Sp", [Index("s1", sa) => lower])), _RT, true)
    clear_global!()
    _ = @load _RT
    let back = SET_REGISTRY["spin5"]
        @test back.algebra !== nothing
        @test typeof(back.algebra) === typeof(sa.algebra)      # parameters and all
        @test dimension(back) == dimension(sa)
    end

    # the literal-type-parameter case directly, independent of any algebra
    @test Alakazam._resolve_type("Val{:x}") === Val{:x}
    @test Alakazam._resolve_type("Array{Int64, 2}") === Array{Int64,2}

    ###########################################################################
    # a USER-DEFINED type carrying a new field, including a referential one
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    let us = IndexSet("userset", 3, lower, Metric)
        m = _RoundTripThing("tag", 7, us)
        d = Alakazam.to_dict(m)
        @test haskey(d, "where_it_lives")
        @test d["where_it_lives"]["__ref__"] == "IndexSet"     # stored by reference
        back = Alakazam.from_json_value(d, Main)
        @test back.where_it_lives === us                        # identity, not a copy
        @test back.label == "tag" && back.n == 7
    end

    ###########################################################################
    # sort_order: global state that silently changes RESULTS if lost
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    st8 = IndexSet("spacetime", 4, lower, Metric)
    o1, o2 = Index(["o1", "o2"], st8)
    custom_sort_order!(["ZZ", "AA"])          # deliberately NOT alphabetical
    before_order = copy(get_sort_order())
    ordered = Alakazam.simplify(TensorExpression(Tensor("AA", [o1 => lower], Coordinate("q", st8)) *
                                                 Tensor("ZZ", [o2 => lower])))
    save(ordered, _RT, true)
    clear_global!()
    let back = @load _RT
        @test get_sort_order() == before_order          # the order itself survived
        # and it still governs: the reloaded expression sorts the same way it did when saved
        @test Alakazam.simplify(back) == back
        @test name(terms(first_summand(back))[1]) == name(terms(first_summand(ordered))[1])
    end

    ###########################################################################
    # derivative_actions: how a derivative acts on a basis tensor
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    stt = IndexSet("spacetime", 4, lower, Metric)
    sptt = IndexSet("spin", 2, lower, fermionic, EpsilonTensor)
    aa, bb = Index(["aa", "bb"], sptt)
    θc = Coordinate("θ", sptt)
    ∂θ = PD("∂", [θc], [aa => lower])
    @action ∂θ(Tensor("θ", [bb => upper], θc)) = KroneckerDelta([aa => lower, bb => upper])
    @test !isempty(derivative_actions)

    save(TensorExpression(Tensor("θ", [bb => upper], θc)), _RT, true)
    clear_global!()
    _ = @load _RT
    # the registration must survive AND still fire: an action that reloads but does nothing is
    # indistinguishable from a lost one until a derivative silently fails to evaluate
    @test !isempty(derivative_actions)
    let sp2 = SET_REGISTRY["spin"], θ2 = COORD_REGISTRY["θ"]
        a2, b2 = Index(["aa", "bb"], sp2)
        applied = derive(PD("∂", Tensor("θ", [b2 => upper], θ2), θ2, [a2 => lower]))
        @test !isZeroTens(applied)                       # the action fired
        @test applied == TensorExpression(KroneckerDelta([a2 => lower, b2 => upper]))
    end

    ###########################################################################
    # derivative_commutators: the CovariantCommutator table
    ###########################################################################
    # Distinct from the matrix `commutators` dict, and stored under its own header.
    clear_global!(); Alakazam.clear_commutators!()
    stc = IndexSet("spacetime", 4, lower, Metric)
    spc2 = IndexSet("spin", 2, lower, fermionic, EpsilonTensor)
    mc = Index("mc", stc); ac, bc = Index(["ac", "bc"], spc2)
    xc = Coordinate("xc", stc); θd = Coordinate("θd", spc2)
    Dc = Derivative("Dc", [xc, θd], [ac => lower])
    Dd = Derivative("Dc", [xc, θd], [bc => lower])
    ∂c = PD("∇", [xc], [mc => lower])
    # the result must carry the same free indices as the product, so a sigma supplies the spinor
    # pair and contracts the vector index of the derivative 
    @anticommutator {Dc, Dd} =
        im * TensorExpression(Tensor("σ", [mc => upper, ac => lower, bc => lower]) * ∂c)
    @test !isempty(Alakazam.derivative_algebra)
    n_before = length(Alakazam.derivative_algebra)

    save(TensorExpression(Tensor("Zc", [mc => lower], xc)), _RT, true)
    clear_global!()
    _ = @load _RT
    @test !isempty(Alakazam.derivative_algebra)
    @test length(Alakazam.derivative_algebra) == n_before    # restored, not duplicated

    ###########################################################################
    # covariant_expansions: what `@covariant` registers
    ###########################################################################
    # A lost expansion surfaces as KeyError("No expansion registered for ...")
    clear_global!(); Alakazam.clear_commutators!()
    ste = IndexSet("spacetime", 4, lower, Metric)
    spe = IndexSet("spin", 2, lower, fermionic, EpsilonTensor)
    me = Index("me", ste); ae = Index("ae", spe)
    xe = Coordinate("xe", ste); θe = Coordinate("θe", spe)
    ∂e = PD("∂", [θe], [ae => lower])
    ∂xe = PD("∇", [xe], [me => lower])
    De = Derivative("De", [xe, θe], [ae => lower])
    # every summand of the expansion must carry the derivative's own free index, so the vector
    # derivative is dressed with a factor holding it
    @covariant De = ∂e + im * Tensor("Ψ", [ae => lower, me => upper], θe) * ∂xe
    @test !isempty(Alakazam.DERIVATIVE_EXPANSIONS)

    save(TensorExpression(Tensor("Ze", [me => lower], xe)), _RT, true)
    clear_global!()
    _ = @load _RT
    @test !isempty(Alakazam.DERIVATIVE_EXPANSIONS)
    @test haskey(Alakazam.DERIVATIVE_EXPANSIONS, "De")

    ###########################################################################
    # TensorData payload: weights, labels, function_of
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    stw = IndexSet("st", 4, lower, Metric)
    aw = Index("aw", stw)
    xw = Coordinate("xw", stw)
    weighted = Tensor("Wt", [aw => lower], xw, Dict{String,Number}("conformal" => 2, "ghost" => -1))
    # built directly rather than through a T-string: the macro resolves against a default
    # "spacetime" set, which the clear above removes
    labelled = Tensor("Lb", [aw => lower], xw, Alakazam.LabelPair["tg" => lower])
    @test !isempty(labels(labelled))

    save(TensorExpression(weighted), _RT * "_w", true)
    save(TensorExpression(labelled), _RT * "_l", true)
    clear_global!()
    let back = @load (_RT * "_w")
        t = terms(first_summand(back))[1]
        @test weights(t)["conformal"] == 2
        @test weights(t)["ghost"] == -1
        @test !isempty(function_of(t))                     # the coordinate dependence survived
        @test function_of(t)[1].name == "xw"
    end
    clear_global!()
    let back = @load (_RT * "_l")
        t = terms(first_summand(back))[1]
        @test labels(t) == labels(labelled)                # a genuine label round trips
    end
    # and the empty-label case, which is the shape the {} vs [] handling exists for
    clear_global!()
    let stq = IndexSet("st", 4, lower, Metric)
        plain = Tensor("Pl", [Index("aq", stq) => lower])
        @test isempty(labels(plain))
        save(TensorExpression(plain), _RT, true)
        clear_global!()
        t = terms(first_summand(@load _RT))[1]
        @test isempty(labels(t))                           # not {}-shaped garbage
    end

    ###########################################################################
    # coefficients: complex, and symbolic
    ###########################################################################
    # `_tree_from_dict` has a dedicated branch for Complex (rebuilt from "re"/"im"), and the whole
    # `_tree_to_dict` path for symbolic scalars. Matrices now default to ℂ, so an imaginary
    # coefficient is ordinary rather than exotic.
    clear_global!(); Alakazam.clear_commutators!()
    stk = IndexSet("st", 4, lower, Metric)
    ak = Index("ak", stk); xk = Coordinate("xk", stk)
    @syms κ
    cexpr = (3 + 2im) * TensorExpression(Tensor("Cx", [ak => lower], xk))
    sexpr = (2 * κ) * TensorExpression(Tensor("Sx", [ak => lower], xk))
    mixed = (im * κ) * TensorExpression(Tensor("Mx", [ak => lower], xk))

    for (nm, e) in (("c", cexpr), ("s", sexpr), ("m", mixed))
        save(e, _RT * "_" * nm, true)
    end
    clear_global!()
    @test (@load (_RT * "_c")) == cexpr        # complex coefficient exact
    clear_global!()
    @test (@load (_RT * "_s")) == sexpr        # symbolic coefficient exact
    clear_global!()
    @test (@load (_RT * "_m")) == mixed        # symbolic TIMES imaginary

    ###########################################################################
    # lightcone indices, which carry their own flag and may legitimately repeat
    ###########################################################################
    clear_global!(); Alakazam.clear_commutators!()
    let lc = IndexSet("lightcone", 2, lower, Metric)
        @test lc.islightcone
        p1 = Index("p1", lc)
        save(TensorExpression(Tensor("Lc", [p1 => lower])), _RT, true)
        clear_global!()
        _ = @load _RT
        @test SET_REGISTRY["lightcone"].islightcone         # the flag is not lost
    end

    ###########################################################################
    # the file-level API: header-only saves, headerless saves, skip_header,
    # idempotency, and repeated loads
    ###########################################################################
    # `clear_global!` alone here, WITHOUT `clear_commutators!`
    clear_global!()
    sti = IndexSet("st", 4, lower, Metric)
    ai = Index("ai", sti); xi = Coordinate("xi", sti)
    ei = 2 * TensorExpression(Tensor("Ei", [ai => lower], xi))

    save_header(_RT * "_hdr")
    @test isfile(_RT * "_hdr.json")

    save(ei, _RT * "_nohdr")
    @test load(_RT * "_nohdr") == ei

    # `skip_header` on a file that has one: a documented argument of `load`
    save(ei, _RT * "_hdr2", true)
    @test load(_RT * "_hdr2", true) == ei

    # Idempotency
    j1 = read(_RT * "_hdr2.json", String)
    clear_global!()
    reloaded = load(_RT * "_hdr2")
    @test reloaded == ei
    save(reloaded, _RT * "_hdr3", true)
    @test read(_RT * "_hdr3.json", String) == j1        # byte-identical round trip

    # loading twice must not duplicate global state
    clear_global!()
    first_load = load(_RT * "_hdr2");  n_sets = length(SET_REGISTRY)
    second_load = load(_RT * "_hdr2")
    @test length(SET_REGISTRY) == n_sets
    @test first_load == second_load

    clear_global!(); Alakazam.clear_commutators!()
end
