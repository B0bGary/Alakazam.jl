@testset "Complex Conjugation" begin

    #setup
    if 1 == 1
        clear_global!()
        # spinor and its conjugate are DIFFERENT index sets -> complex representation
        spinor = IndexSet("spinor", 2, lower, EpsilonTensor, fermionic)
        spinordot = IndexSet("spinordot", 2, lower, EpsilonTensor, fermionic)
        spinor.conjugate = spinordot
        spinordot.conjugate = spinor
        # SU(2) is self-conjugate with an ANTISYMMETRIC form -> quaternionic, the case that signs
        su2 = IndexSet("su2", 2, lower, EpsilonTensor)
        su2.conjugate = su2
        # spacetime is self-conjugate with a SYMMETRIC metric -> fixed
        spacetime = IndexSet("spacetime", 4, upper, Metric)
        # spacetime.conjugate = spacetime

        α = Index("α", spinor)
        αd = Index("̇α", spinordot)
        β = Index("β", spinor)
        i = Index("i", su2)
        j = Index("j", su2)
        a = Index("a", spacetime)

        θ(x...) = Tensor("θ", ℂ, collect(x))
        θbar(x...) = Tensor(Alakazam.togglebar("θ"), ℂ, collect(x))
    end
    # @show Alakazam.togglebar("a")
    # @test Alakazam.togglebar("θ")=="̄θ"

    # @show t=θ(α=>upper,i=>upper)
    # @show t'
    #  @show t=θ(α=>upper,i=>lower)
    # @show t'
    # @show t=θ(α=>lower,i=>upper)
    # @show t'
    # @show t=θ(α=>lower,i=>lower)
    # @show t'
    # @show t=θbar(α=>upper,i=>upper)
    # @show t'
    #  @show t=θbar(α=>upper,i=>lower)
    # @show t'
    # @show t=θbar(α=>lower,i=>upper)
    # @show t'
    # @show t=θbar(α=>lower,i=>lower)
    # @show t'

    @testset "index set classification" begin
        # read off `conjugate` and `metric_type`, no dedicated field
        @test Alakazam.conj_class(spinor) === :complex
        @test Alakazam.conj_class(spinordot) === :complex
        @test Alakazam.conj_class(su2) === :quaternionic
        @test Alakazam.conj_class(spacetime) === :fixed
        # self-conjugate with no metric cannot be told apart from quaternionic, and must not be
        # guessed.
        # The constructor always supplies a metric, so the field is cleared directly
        bare = IndexSet("bare", 2, lower)
        bare.conjugate = bare
        bare.metric_type = nothing
        @test_throws ErrorException Alakazam.conj_class(bare)
    end

    @testset "number field decides whether conjugation acts" begin
        F = Tensor("F", [a => upper])                    # real by default
        G = Tensor("G", iℝ, [a => upper])
        @test F' == F
        @test isZeroTens(G' + G)                          # imaginary: conj = -self
        @test isReal(F) && !isReal(G)
        @test isImaginary(G)
    end

    @testset "4D N=2 reality conditions" begin
        # the four conditions, which between them fix the whole sign structure. The minus signs
        # come from (eps^{ij})* = -eps_{ij}: an SU(2) index displaced from its field's reference.
        @test θ(α => upper, i => lower)'  == θbar(αd => upper, i => upper)
        @test θ(α => upper, i => upper)'  == -1 * θbar(αd => upper, i => lower)
        @test θbar(αd => lower, i => upper)' == θ(α => lower, i => lower)
        @test θbar(αd => lower, i => lower)' == -1 * θ(α => lower, i => upper)
    end

    @testset "involution" begin
        for T in (θ(α => upper, i => lower), θ(α => upper, i => upper),
                  θbar(αd => lower, i => upper), θbar(αd => lower, i => lower))
            @test (T')' == T
        end
        # the two SU(2) signs cancel, which is what makes conj o conj the identity rather than -1
        @test (θ(α => upper, i => upper)')' == θ(α => upper, i => upper)
    end

    @testset "positions and index sets" begin
        c = θ(α => upper, i => lower)'
        # spinor position PRESERVED, index set swapped for its partner
        @test indices(c)[1][2] == upper
        @test indices(c)[1][1].index_set.name == "spinordot"
        # SU(2) position FLIPPED, index set unchanged
        @test indices(c)[2][2] == upper
        @test indices(c)[2][1].index_set.name == "su2"
        # a real index set moves neither
        v = Tensor("V", ℂ, [a => upper])'
        @test indices(v)[1][2] == upper
        @test indices(v)[1][1].index_set.name == "spacetime"
    end

    @testset "barred names toggle" begin
        @test Alakazam.togglebar("psi") != "psi"
        @test Alakazam.togglebar(Alakazam.togglebar("psi")) == "psi"
        @test Alakazam.isbarred(Alakazam.togglebar("psi"))
        @test !Alakazam.isbarred("psi")
        @test name(θ(α => upper, i => lower)') == Alakazam.togglebar("θ")
    end

    @testset "invariant forms" begin
        # epsilon on a quaternionic set: both positions flip and the pair costs one sign
        eps_su2 = SymplecticForm("ε", [i => upper, j => upper])
        r = eps_su2'
        @test indices(r)[1][2] == lower && indices(r)[2][2] == lower
        # epsilon on spinors: partner set, positions kept, NO sign
        eps_sp = SymplecticForm("ε", [α => lower, β => lower])
        rs = eps_sp'
        @test indices(rs)[1][1].index_set.name == "spinordot"
        @test indices(rs)[1][2] == lower
        # Kronecker moves to the partner set and keeps its positions; no sign, it is the identity
        d = KroneckerDelta(α => lower, β => upper)
        dc = d'
        @test dc isa KroneckerSuperType
        @test indices(dc)[1][1].index_set.name == "spinordot"
        @test indices(dc)[1][2] == lower && indices(dc)[2][2] == upper
    end

    @testset "antimorphism reverses a product" begin
        P = Tensor("psi", ℂ, [a => upper])
        X = Tensor("chi", ℂ, [a => lower])
        # (AB)^dagger = B^dagger A^dagger: the factors come back in the opposite order
        pd = (P * X)'
        @test terms(pd)[1] == X'
        @test terms(pd)[2] == P'
        # a real product is unchanged apart from that reordering
        R1 = Tensor("R", [a => upper])
        R2 = Tensor("S", [a => lower])
        @test isZeroTens((R1 * R2)' - R2 * R1)
    end

    @testset "sums distribute with conjugated coefficients" begin
        P = Tensor("psi", ℂ, [a => upper])
        Q = Tensor("phi", ℂ, [a => upper])
        @test isZeroTens((P + Q)' - (P' + Q'))
        @test isZeroTens((2 * P)' - 2 * P')
    end
    @testset "partner indices follow the declared convention" begin
        # The dot is conventional, not mandatory. Partners are matched on the undecorated base, so
        # a user who marks the conjugate set with a bar or a hat gets those back, not a dot.
        CD, CM, CH = Alakazam.COMBINING_DOT, Alakazam.COMBINING_MACRON, Alakazam.COMBINING_HAT
        function cpair(nm)
            a = IndexSet(nm, 2, lower, EpsilonTensor, fermionic)
            b = IndexSet(nm * "c", 2, lower, EpsilonTensor, fermionic)
            a.conjugate = b; b.conjugate = a
            (a, b)
        end
        conj1(idx) = indices(Tensor("ψ", ℂ, [idx => lower])')[1][1]

        # a barred partner is found, not shadowed by a freshly dotted one
        s1, c1 = cpair("cv1")
        x, xb = Index("α", s1), Index(string(CM) * "α", c1)
        @test conj1(x) === xb
        @test (Tensor("ψ", ℂ, [x => lower])')' == Tensor("ψ", ℂ, [x => lower])

        # likewise a hat
        s2, c2 = cpair("cv2")
        y, yh = Index("β", s2), Index(string(CH) * "β", c2)
        @test conj1(y) === yh

        # generating a partner infers the mark from the set instead of assuming a dot
        s3, c3 = cpair("cv3")
        Index("α", s3); Index(string(CM) * "α", c3)
        @test conj1(Index("γ", s3)).name == string(CM) * "γ"

        # an empty fermionic partner set offers nothing to infer from: the dot is the default
        s4, _ = cpair("cv4")
        @test conj1(Index("δ", s4)).name == string(CD) * "δ"

        # and the reverse direction strips the mark rather than adding another
        s5, c5 = cpair("cv5")
        p, pd = Index("ε", s5), Index(string(CD) * "ε", c5)
        @test conj1(pd) === p

        @test Alakazam.indexmark(string(CD) * "α") == (CD, "α")
        @test Alakazam.indexmark("α") == (nothing, "α")
    end
    @testset "antimorphism is generic in its leaf and scalar actions" begin
        # The reversal is structural and belongs to AntimorphismSuperType as a whole; dagger only
        # supplies what it does to one factor. Both parameters must actually be honoured, or a
        # future transpose/reversion could not reuse this.
        A = Tensor("A", ℂ, [a => upper])
        B = Tensor("B", ℂ, [a => lower])
        tt = A * B

        # leaf action = identity: factors reverse, nothing is conjugated
        rev = Alakazam.apply_antimorphism(identity, tt)
        @test terms(rev)[1] == B && terms(rev)[2] == A
        # leaf action = dagger: reversed AND conjugated
        dag = Alakazam.apply_antimorphism(Alakazam._dagger, tt)
        @test terms(dag)[1] == B' && terms(dag)[2] == A'
        @test dag == tt'

        # the scalar action is independent of the leaf action
        e = (2 + 3im) * A
        @test first_coef(Alakazam.apply_antimorphism(identity, e, conj)) == 2 - 3im
        @test first_coef(Alakazam.apply_antimorphism(identity, e, identity)) == 2 + 3im
    end
end
