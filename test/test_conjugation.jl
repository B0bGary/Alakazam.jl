struct MajC <: AntisymMetricSuperType
    data::Alakazam.TensorData
    indices::Vector{Alakazam.IndexPair}
end
Alakazam.conj_convention(::Type{MajC}, S::IndexSet) = (upper, -1)

# A form following the convention of Gates Superspace, where C is imaginary:
# (C_{ab})* = -C_{adot bdot}, from which (psi_a)^dag = -psibar_{adot} follows. Declared here exactly
# as a user would declare it, so this doubles as a check that one method really is enough.
struct GatesC <: AntisymMetricSuperType
    data::Alakazam.TensorData
    indices::Vector{Alakazam.IndexPair}
end
Alakazam.conj_convention(::Type{GatesC}, S::IndexSet) = (upper, -1)

@testset "Complex Conjugation" begin

    #setup
    if 1 == 1
        clear_global!()
        # spinor and its conjugate are DIFFERENT index sets -> complex representation
        spinor = IndexSet("spinor", 2, lower, EpsilonTensor, fermionic)
        spinordot = IndexSet("spinordot", 2, lower, EpsilonTensor, fermionic)
        spinor.conjugate = spinordot
        spinordot.conjugate = spinor
        # SU(2) is self-conjugate and quaternionic. That is carried by the ALGEBRA: an
        # antisymmetric form alone does not imply it, so a set declaring only the form is real.
        su2 = IndexSet("su2", 2, lower, EpsilonTensor, SUAlgebra{ℂ}(2))
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
        # A self-conjugate set with nothing else declared is real. The invariant form is not
        # consulted at all: its symmetry does not determine the reality, so a set carrying an
        # antisymmetric form but no algebra is real too, and quaternionic must be declared.
        bare = IndexSet("bare", 2, lower)
        bare.conjugate = bare
        bare.metric_type = nothing
        @test Alakazam.conj_class(bare) === :fixed
        antisym_only = IndexSet("antisym_only", 4, lower, EpsilonTensor, fermionic)
        antisym_only.conjugate = antisym_only
        @test Alakazam.conj_class(antisym_only) === :fixed
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
            local su = IndexSet(nm, 2, lower, EpsilonTensor, fermionic)
            local sc = IndexSet(nm * "c", 2, lower, EpsilonTensor, fermionic)
            su.conjugate = sc; sc.conjugate = su
            (su, sc)
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
    @testset "conj_convention is the single override point" begin
        # the defaults, which every other testset above is implicitly relying on
        @test Alakazam.conj_convention(su2.metric_type, su2) == (su2.default_position, -1)
        @test Alakazam.conj_convention(spinor.metric_type, spinor) == (spinor.default_position, 1)
        @test Alakazam.conj_convention(spacetime.metric_type, spacetime)[2] == 1

        # the barred flip lives in the engine, not the convention: theta and thetabar share one
        # index set but carry their doublet index on opposite levels
        (ref_u, _) = Alakazam._conj_rule(su2, false)
        (ref_b, _) = Alakazam._conj_rule(su2, true)
        @test ref_u == su2.default_position
        @test ref_b == Alakazam._flip(su2.default_position)

        # form and fields read the same pair, so they cannot disagree. The epsilon takes the sign
        # outright; a field slot takes it only when displaced from the reference
        s = Alakazam.conj_convention(su2.metric_type, su2)[2]
        @test SymplecticForm("ε", [i => upper, j => upper])' ==
              s * SymplecticForm("ε", [i => lower, j => lower])
        atref = θ(α => upper, i => su2.default_position)'
        @test first_coef(TensorExpression(atref)) == 1
        away  = θ(α => upper, i => Alakazam._flip(su2.default_position))'
        @test first_coef(TensorExpression(away)) == s
    end

    @testset "a user convention reaches the fields" begin
        # Gates' rule, declared once at the top of this file as `conj_convention(::Type{GatesC},…)`.
        # Nothing else is stated: the signs on psi follow from the form's.
        g  = IndexSet("gates", 2, upper, GatesC, fermionic)
        gd = IndexSet("gatesdot", 2, upper, GatesC, fermionic)
        g.conjugate = gd; gd.conjugate = g
        a  = Index("α", g); Index("α", gd)      # the partner, matched on the undecorated base
        ψ(x...)    = Tensor("ψ", ℂ, collect(x))
        ψbar(x...) = Tensor(Alakazam.togglebar("ψ"), ℂ, collect(x))

        @test Alakazam.conj_class(g) === :complex
        @test Alakazam.conj_convention(g.metric_type, g) == (upper, -1)

        # reference is upper, so the raised index is unsigned and the lowered one is not:
        # (psi^a)^dag = psibar^adot,   (psi_a)^dag = -psibar_adot
        cu = ψ(a => upper)'
        cl = ψ(a => lower)'
        @test first_coef(TensorExpression(cu)) == 1
        @test first_coef(TensorExpression(cl)) == -1
        # the index still moves to the partner set, and keeps its level
        @test indices(first_summand(TensorExpression(cu)))[1][1].index_set.name == "gatesdot"
        @test indices(first_summand(TensorExpression(cu)))[1][2] == upper
        @test indices(first_summand(TensorExpression(cl)))[1][2] == lower

        # and it is still an involution, the two signs cancelling
        @test (ψ(a => lower)')' == ψ(a => lower)
    end
    @testset "class x field x position" begin
        # Every combination of conjugation class, number field and index level, checked for all
        # four things conjugation can do to a leaf: the overall sign, whether the name is barred,
        # which index set the slot lands in, and at which level.
        #
        # The field and the indices are independent, so the sign is their product: iℝ contributes
        # -1 on top of whatever the slot contributes, and only ℂ bars the name.
        st  = IndexSet("mx_st", 4, upper, Metric)                        # :fixed,  (upper, +1)
        # Spin(2,1) is what makes this real: an antisymmetric form alone would classify as
        # quaternionic, which is the ambiguity the algebra route exists to settle. The
        # explicit MajC still wins over the metric the algebra would supply.
        mj  = IndexSet("mx_maj", upper, MajC, fermionic, Spin(2, 1))    # :fixed,  (upper, -1)
        mj.conjugate = mj
        wl  = IndexSet("mx_wl",  2, lower, EpsilonTensor, fermionic)     # :complex
        wld = IndexSet("mx_wld", 2, lower, EpsilonTensor, fermionic)
        wl.conjugate = wld; wld.conjugate = wl
        su  = IndexSet("mx_su2", 2, lower, EpsilonTensor, SUAlgebra{ℂ}(2))  # :quaternionic
        su.conjugate = su
        for (S, nm) in ((st,"a"), (mj,"α"), (wl,"α"), (su,"i"))
            Index(nm, S)
        end
        Index("α", wld)                       # the partner, matched on the undecorated base

        @test Alakazam.conj_class(st) === :fixed
        @test Alakazam.conj_class(mj) === :fixed
        @test Alakazam.conj_class(wl) === :complex
        @test Alakazam.conj_class(su) === :quaternionic

        """(coefficient, name barred?, index set, level) after conjugating a one-index tensor."""
        function probe(S, pos, field)
            t = field === nothing ? Tensor("ψ", [first(S.indices) => pos]) :
                                    Tensor("ψ", field, [first(S.indices) => pos])
            e = TensorExpression(adjoint(t))
            tens = terms(first_summand(e))[1]
            (first_coef(e), Alakazam.isbarred(name(tens)),
             indices(tens)[1][1].index_set.name, indices(tens)[1][2])
        end

        # set, level, field  =>  sign, barred, resulting set, resulting level.
        # Fields are explicit here; the inferred default is checked separately below.
        cases = [
            # :fixed, default convention: no sign from the index at either level
            (st, upper, ℝ,   1, false, "mx_st", upper),
            (st, lower, ℝ,   1, false, "mx_st", lower),
            (st, upper, iℝ, -1, false, "mx_st", upper),
            (st, lower, iℝ, -1, false, "mx_st", lower),
            (st, upper, ℂ,   1, true,  "mx_st", upper),
            (st, lower, ℂ,   1, true,  "mx_st", lower),
            # Majorana: index and level untouched, but the lower slot is imaginary
            (mj, upper, ℝ,   1, false, "mx_maj", upper),
            (mj, lower, ℝ,  -1, false, "mx_maj", lower),
            (mj, upper, iℝ, -1, false, "mx_maj", upper),
            (mj, lower, iℝ,  1, false, "mx_maj", lower),
            (mj, upper, ℂ,   1, true,  "mx_maj", upper),
            (mj, lower, ℂ,  -1, true,  "mx_maj", lower),
            # complex: moves to the partner set, keeps its level. ℂ is the only legal field
            (wl, upper, ℂ,   1, true,  "mx_wld", upper),
            (wl, lower, ℂ,   1, true,  "mx_wld", lower),
            # quaternionic: keeps its set, flips level, signs away from the reference (lower).
            # ℝ is legal -- a real field on a quaternionic index is a symplectic Majorana condition
            (su, lower, ℝ,   1, false, "mx_su2", upper),
            (su, upper, ℝ,  -1, false, "mx_su2", lower),
            (su, lower, iℝ, -1, false, "mx_su2", upper),
            (su, upper, iℝ,  1, false, "mx_su2", lower),
            (su, lower, ℂ,   1, true,  "mx_su2", upper),
            (su, upper, ℂ,  -1, true,  "mx_su2", lower),
        ]
        for (S, pos, field, sgn, bar, set, npos) in cases
            @test probe(S, pos, field) == (sgn, bar, set, npos)
        end

        # A representation not isomorphic to its conjugate admits no reality condition, so a tensor
        # carrying one cannot be declared real or imaginary. Quaternionic is not caught: that is a
        # symplectic Majorana condition, which is legitimate.
        @test_throws ErrorException probe(wl, lower, ℝ)
        @test_throws ErrorException probe(wl, upper, iℝ)

        # with nothing declared, the field is read off the indices
        idx(S) = first(S.indices)
        @test numberfield(Tensor("t", [idx(st) => upper])) isa Alakazam.RealNumberSet
        @test numberfield(Tensor("t", [idx(mj) => upper])) isa Alakazam.RealNumberSet
        @test numberfield(Tensor("t", [idx(wl) => lower])) isa Alakazam.ComplexNumberSet
        @test numberfield(Tensor("t", [idx(su) => lower])) isa Alakazam.ComplexNumberSet
        @test numberfield(Tensor("t")) isa Alakazam.RealNumberSet             # a scalar
        # ...and an explicit field still wins, in both directions
        @test numberfield(Tensor("t", ℂ, [idx(st) => upper])) isa Alakazam.ComplexNumberSet
        @test numberfield(Tensor("t", ℝ, [idx(su) => lower])) isa Alakazam.RealNumberSet
    end

    @testset "conjugation class across the spin groups" begin
        # the class must follow the division algebra for every signature, not just the four the
        # exceptional isomorphisms pin
        for d in 2:11, p in (1, 2)
            q = d - p
            q >= 1 || continue
            alg = Spin(p, q)
            da  = division_algebra(alg)
            S = IndexSet("sp_$(p)_$(q)", lower, fermionic, alg, EpsilonTensor)
            if da === ℂ
                # a complex representation needs a partner; declaring it self-conjugate is an error
                S.conjugate = S
                @test_throws ErrorException Alakazam.conj_class(S)
            else
                S.conjugate = S
                @test Alakazam.conj_class(S) === (da === ℍ ? :quaternionic : :fixed)
            end
        end
    end

    @testset "several indices at once" begin
        # Everything above is rank 1. With more than one slot the signs multiply, classes can mix
        # inside one tensor, and the field is decided by all the slots together.
        st  = IndexSet("mi_st", 4, upper, Metric);            a = Index("a", st)
        su  = IndexSet("mi_su2", 2, lower, EpsilonTensor, SUAlgebra{ℂ}(2)); su.conjugate = su
        i, j = Index("i", su), Index("j", su)
        wl  = IndexSet("mi_wl",  2, lower, EpsilonTensor, fermionic)
        wld = IndexSet("mi_wld", 2, lower, EpsilonTensor, fermionic)
        wl.conjugate = wld; wld.conjugate = wl
        al = Index("α", wl); Index("α", wld)
        mj = IndexSet("mi_mj", upper, MajC, fermionic, Spin(2, 1)); mj.conjugate = mj
        m = Index("μ", mj)

        """(coefficient, name, [(set, level)...]) after conjugating."""
        function shape(t)
            e = TensorExpression(adjoint(t))
            tn = terms(first_summand(e))[1]
            (first_coef(e), name(tn),
             [(x[1].index_set.name, x[2]) for x in indices(tn)])
        end

        # two quaternionic slots: one sign per slot away from the reference, so they multiply
        @test shape(Tensor("T", ℂ, [i=>lower, j=>lower]))[1] ==  1
        @test shape(Tensor("T", ℂ, [i=>upper, j=>lower]))[1] == -1
        @test shape(Tensor("T", ℂ, [i=>upper, j=>upper]))[1] ==  1   # (-1)^2
        # ...and both levels flip
        @test shape(Tensor("T", ℂ, [i=>upper, j=>upper]))[3] ==
              [("mi_su2", lower), ("mi_su2", lower)]

        # two Majorana slots, both at the reference: no sign, nothing moves, no bar
        @test shape(Tensor("X", [m=>upper, m=>upper])) ==
              (1, "X", [("mi_mj", upper), ("mi_mj", upper)])

        # mixed classes in one tensor: each slot follows its own rule
        @test shape(Tensor("ψ", ℂ, [al=>lower, i=>lower])) ==
              (1, Alakazam.togglebar("ψ"), [("mi_wld", lower), ("mi_su2", upper)])
        @test shape(Tensor("ψ", ℂ, [al=>lower, a=>upper])) ==
              (1, Alakazam.togglebar("ψ"), [("mi_wld", lower), ("mi_st", upper)])
        # Majorana sign survives alongside a slot that contributes nothing, and stays unbarred
        @test shape(Tensor("X", [m=>lower, a=>upper])) ==
              (-1, "X", [("mi_mj", lower), ("mi_st", upper)])

        # the inferred field looks at every slot: one complex index is enough
        @test numberfield(Tensor("q", [al=>lower, a=>upper])) isa Alakazam.ComplexNumberSet
        @test numberfield(Tensor("q", [m=>lower,  a=>upper])) isa Alakazam.RealNumberSet
        @test numberfield(Tensor("q", [a=>upper,  a=>lower])) isa Alakazam.RealNumberSet

        # and so does the contradiction check, wherever the complex index sits
        @test_throws ErrorException adjoint(Tensor("q", ℝ, [a=>upper, al=>lower]))
        @test_throws ErrorException adjoint(Tensor("q", ℝ, [al=>lower, a=>upper]))

        # involution with several slots: the signs cancel pairwise across the two applications
        for t in (Tensor("T", ℂ, [i=>upper, j=>lower]),
                  Tensor("T", ℂ, [i=>upper, j=>upper]),
                  Tensor("ψ", ℂ, [al=>lower, i=>lower]),
                  Tensor("X", [m=>lower, a=>upper]),
                  Tensor("X", [m=>lower, m=>upper]))
            @test adjoint(adjoint(t)) == t
        end
    end

    @testset "@conjugates declares the pairing" begin
        # two sets become each other's conjugate; one becomes its own. Both forms forward to
        # conjugates!, so the function is equally usable where a macro would be awkward.
        clear_global!()
        w  = IndexSet("cw", upper, fermionic, Spin(4))
        wb = IndexSet("cwbar", upper, fermionic, Spin(4))
        @conjugates w wb
        @test w.conjugate === wb
        @test wb.conjugate === w
        @test Alakazam.conj_class(w) === :complex

        # a comma between them reads the same
        a = IndexSet("ca", upper, fermionic, Spin(4))
        b = IndexSet("cb", upper, fermionic, Spin(4))
        @conjugates a, b
        @test a.conjugate === b && b.conjugate === a

        maj = IndexSet("cmaj", upper, fermionic, Spin(3))
        @conjugates maj
        @test maj.conjugate === maj
        @test Alakazam.conj_class(maj) === :fixed

        su2 = IndexSet("csu2", lower, SUAlgebra{ℂ}(2))
        @conjugates su2
        @test Alakazam.conj_class(su2) === :quaternionic

        # naming one set twice is the self-conjugate case, not an error
        m2 = IndexSet("cm2", upper, fermionic, Spin(3))
        @conjugates m2 m2
        @test m2.conjugate === m2

        p = IndexSet("cp", upper, fermionic, Spin(4))
        q = IndexSet("cq", upper, fermionic, Spin(4))
        conjugates!(p, q)
        @test p.conjugate === q
    end

    @testset "@conjugates refuses a pairing that cannot hold" begin
        clear_global!()
        # a representation and its conjugate share a dimension and their statistics
        r4 = IndexSet("cr4", 4, lower, EpsilonTensor)
        r2 = IndexSet("cr2", 2, lower, EpsilonTensor)
        @test_throws ErrorException conjugates!(r4, r2)
        f2 = IndexSet("cf2", 2, lower, EpsilonTensor, fermionic)
        b2 = IndexSet("cb2", 2, lower, EpsilonTensor, bosonic)
        @test_throws ErrorException conjugates!(f2, b2)

        # a complex representation cannot be its own conjugate, and is caught at the declaration
        cx = IndexSet("ccx", upper, fermionic, Spin(4))
        @test_throws ErrorException conjugates!(cx)
        @test cx.conjugate === nothing        # and the set is left as it was
    end

end
