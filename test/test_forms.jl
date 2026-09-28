@testset "Differential Forms" begin

    #setup
    if 1==1
        euclidean = IndexSet("spacetime", 4, (0,4), lower)
        eu = Index(["a","b","c","d","e","f","g","h","i","j","k","l"], euclidean)
        ea, eb = eu[1], eu[2]

        lorentz = IndexSet("lorentzian", 4, lower)
        lo = Index(["p","q","r","s","t","u","v","w","x","y","z","ø"], lorentz)
        la, lb = lo[1], lo[2]

        spinor = IndexSet("spinor_form", 4, upper, fermionic)
        α = Index("α", spinor)
        β = Index("β", spinor)

        flavour = IndexSet("flavour_form", 3, lower)
        fi = Index("i", flavour)
        fj = Index("j", flavour)
    end

    #a one-form wedged with itself vanishes, being bosonic
    @test isZeroTens(simplify(Tensor("A",[ea=>lower]) ∧ Tensor("A",[eb=>lower])))

    #and two distinct one-forms give the antisymmetrisation
    @test isZeroTens(simplify(
        (Tensor("A",[ea=>lower]) ∧ Tensor("B",[eb=>lower]))
        - (Tensor("A",[ea=>lower])*Tensor("B",[eb=>lower])
           - Tensor("A",[eb=>lower])*Tensor("B",[ea=>lower]))))

    #the wedge lives on the cotangent space and never sees the metric, so the same
    #antisymmetrisation must come back in Lorentzian signature as in Euclidean. The
    #epsilon pair used to build it carries a det(η), which is divided back out.
    @test isZeroTens(simplify(
        (Tensor("A",[la=>lower]) ∧ Tensor("B",[lb=>lower]))
        - (Tensor("A",[la=>lower])*Tensor("B",[lb=>lower])
           - Tensor("A",[lb=>lower])*Tensor("B",[la=>lower]))))

    #an index from another set rides through the wedge untouched
    if 1==1
        ψ1 = Tensor("ψ",[ea=>lower, α=>upper])
        ψ2 = Tensor("ψ",[eb=>lower, β=>upper])
        w = simplify(ψ1 ∧ ψ2)
        @test any(ind[1] == α for t in terms(summands(w)[1][2]) for ind in indices(t))
        @test any(ind[1] == β for t in terms(summands(w)[1][2]) for ind in indices(t))
    end

    #so a spinor-valued one-form wedged with itself does NOT vanish, the two
    #sign flips cancelling. every gravitino bilinear rests on that
    @test !isZeroTens(simplify(Tensor("ψ",[ea=>lower, α=>upper])
                               ∧ Tensor("ψ",[eb=>lower, β=>upper])))

    #and it is the antisymmetrisation carrying the graded sign
    @test isZeroTens(simplify(
        (Tensor("ψ",[ea=>lower,α=>upper]) ∧ Tensor("ψ",[eb=>lower,β=>upper]))
        - (Tensor("ψ",[ea=>lower,α=>upper])*Tensor("ψ",[eb=>lower,β=>upper])
           - Tensor("ψ",[eb=>lower,α=>upper])*Tensor("ψ",[ea=>lower,β=>upper]))))

    #a flavour index survives the same way
    if 1==1
        wf = simplify(Tensor("A",[ea=>lower, fi=>lower]) ∧ Tensor("B",[eb=>lower, fj=>lower]))
        @test any(ind[1] == fi for t in terms(summands(wf)[1][2]) for ind in indices(t))
        @test any(ind[1] == fj for t in terms(summands(wf)[1][2]) for ind in indices(t))
    end

    #the Hodge star already kept them, and still does
    if 1==1
        s = ⋆Tensor("ψ",[ea=>lower, α=>upper])
        @test any(ind[1] == α for t in terms(summands(s)[1][2]) for ind in indices(t))
    end

    #the form operations generate the dummies they need rather than indexing
    #into whatever the caller happened to declare
    if 1==1
        thin = IndexSet("spacetime", 4, (0,4), lower)
        tp, tq = Index(["p1","q1"], thin)
        @test !isZeroTens(simplify(Tensor("A",[tp=>lower]) ∧ Tensor("B",[tq=>lower])))
    end

    #a form of degree above the dimension vanishes
    @test isZeroTens(
        AntisymmetricTensor("F",[eu[1]=>lower,eu[2]=>lower,eu[3]=>lower])
        ∧ AntisymmetricTensor("G",[eu[4]=>lower,eu[5]=>lower,eu[6]=>lower]))

    #A 0-form is a scalar, and wedging with one is plain multiplication. The index set has
    #to be read from both factors for this to work at all: a scalar carries no slot to
    #name one with, so `f ∧ ω` has to find the set on the other side.
    if 1==1
        xs = Coordinate("x_scal", euclidean)
        f0 = Tensor("f", xs)
        g0 = Tensor("g", xs)
        A1s = Tensor("A", eu[1] => lower, xs)
        F2s = AntisymmetricTensor("F", [eu[1] => lower, eu[2] => lower], xs)
        @test isZeroTens(simplify((f0 ∧ g0) - f0 * g0))
        @test isZeroTens(simplify((f0 ∧ A1s) - f0 * A1s))
        @test isZeroTens(simplify((A1s ∧ f0) - A1s * f0))
        @test isZeroTens(simplify((f0 ∧ F2s) - f0 * F2s))

        #A scalar carrying dummies of its own has them renamed clear of the form's slot,
        #which is why the product cannot simply be taken first: S_{ab}S^{ab}A_{a} puts `a`
        #on three slots and is rejected outright. The wedge renames, so the slot survives
        #and the answer is the one the same scalar written on free names would give.
        sdirty = AntisymmetricTensor("S", [eu[1] => lower, eu[2] => lower], xs) *
                 AntisymmetricTensor("S", [eu[1] => upper, eu[2] => upper], xs)
        sclean = AntisymmetricTensor("S", [eu[3] => lower, eu[4] => lower], xs) *
                 AntisymmetricTensor("S", [eu[3] => upper, eu[4] => upper], xs)
        @test Set(pr[1] for pr in Alakazam.free_indices_ordered(sdirty ∧ A1s)) ==
              Set([eu[1]])
        @test isZeroTens(simplify((sdirty ∧ A1s) - sclean * A1s))
    end

    ###########################################################################
    ### the Hodge star
    ###########################################################################

    #⋆ hands back the argument's own slots. 
    if 1==1
        Fcd = AntisymmetricTensor("F",[eu[3]=>lower, eu[4]=>lower])
        free = Set(pr[1] for pr in Alakazam.free_indices_ordered(⋆(Fcd)))
        @test eu[3] in free
        @test eu[4] in free
    end

    #⋆ tops up a declared pool 
    if 1==1
        short = IndexSet("spacetime", 4, (0,4), lower)
        sm, sn = Index(["m1","m2"], short)
        @test !isZeroTens(⋆(AntisymmetricTensor("F",[sm=>lower, sn=>lower])))
    end

    #an epsilon carries all its indices at one level
    @test_throws Alakazam.InvalidDifferentialForm ⋆(
        AntisymmetricTensor("F",[ea=>upper, eb=>lower]))

    #⋆⋆ω = (-1)^{p(n-p)} sgn(det η) ω
    if 1==1
        #Euclidean, p=2, n=4:  (+1)^{} · (+1) = +1
        F2e = AntisymmetricTensor("F",[ea=>lower, eb=>lower])
        @test isZeroTens(simplify(⋆(⋆(F2e)) - F2e))

        #Lorentzian, p=2, n=4:  (-1)^{2·2} · (-1) = -1
        F2l = AntisymmetricTensor("F",[la=>lower, lb=>lower])
        @test isZeroTens(simplify(⋆(⋆(F2l)) + F2l))

        #Lorentzian, p=1, n=4:  (-1)^{1·3} · (-1) = +1
        V1l = Tensor("V",[la=>lower])
        @test isZeroTens(simplify(⋆(⋆(V1l)) - V1l))
    end

    #⋆ keeps the names it was given even when they are not the ones that sort first.

    if 1==1
        late4 = IndexSet("st_late4", 4, (0,4), lower)
        q4 = Index(["a_l4","b_l4","c_l4","zz_l4"], late4)      #zz_l4 sorts last
        Vl4 = Tensor("V", [q4[4] => lower])
        @test q4[4] in Set(pr[1] for pr in Alakazam.free_indices_ordered(⋆(Vl4)))
        @test isZeroTens(simplify(⋆(⋆(Vl4)) + Vl4))            #(-1)^{1·3} = -1

        #an odd dimension, where reused and fresh slots divide up differently again
        late5 = IndexSet("st_late5", 5, (0,5), lower)
        q5 = Index(["a_l5","b_l5","c_l5","d_l5","zz_l5"], late5)
        @test isZeroTens(simplify(⋆(⋆(Tensor("V", [q5[5] => lower])))
                                  - Tensor("V", [q5[5] => lower])))         #(-1)^{1·4} = +1
        Fl5 = AntisymmetricTensor("F", [q5[4] => lower, q5[5] => lower])
        @test isZeroTens(simplify(⋆(⋆(Fl5)) - Fl5))                          #(-1)^{2·3} = +1
    end

    ###########################################################################
    ### the wedge against the star
    ###########################################################################

    #F ∧ ⋆F is expressible
    if 1==1
        F1 = AntisymmetricTensor("F",[eu[1]=>lower, eu[2]=>lower])
        F2 = AntisymmetricTensor("F",[eu[3]=>lower, eu[4]=>lower])
        @test !isZeroTens(simplify(F1 ∧ ⋆(F2)))

        # gives the Maxwell term. Contracting the top form with an epsilon gives
        # (1/2)(F·F) ϵ_{abcd} ϵ^{abcd} = 12 sgn(det η) F·F, and this set is Euclidean.
        scal = simplify(EpsilonTensor([eu[1]=>upper, eu[2]=>upper,
                                       eu[3]=>upper, eu[4]=>upper]) * (F1 ∧ ⋆(F2)))
        FdotF = AntisymmetricTensor("F",[eu[1]=>lower, eu[2]=>lower]) *
                AntisymmetricTensor("F",[eu[1]=>upper, eu[2]=>upper])
        @test isZeroTens(simplify(scal - 12 * FdotF))
    end

    #`expand` contracts the wedge's own pair of epsilons eagerly
    if 1==1
        A1 = Tensor("A",[ea=>lower])
        B1 = Tensor("B",[eb=>lower])
        @test isZeroTens(simplify(∧(A1, B1; expand=true) - (A1 ∧ B1)))
    end

    ###########################################################################
    ### the exterior derivative against the wedge
    ###########################################################################

    #d composes with ∧.
    if 1==1
        xf = Coordinate("x_forms", euclidean)
        A1 = Tensor("A", eu[1]=>lower, xf)
        B1 = Tensor("B", eu[2]=>lower, xf)
        @test !isZeroTens(simplify(apply_derivative(exterior_derivative(A1 ∧ B1, xf))))
    end

    #The wedge's slot convention, checked against the textbook formula rather than
    #against another of our own results. For a one-form α and a two-form δ,
    #    (α∧δ)_{μνρ} = α_μ δ_{νρ} + α_ν δ_{ρμ} + α_ρ δ_{μν}.
    #A form's slots are held in canonical (name-sorted) order, so α ∧ dβ with α on "c"
    #and dβ on "d","z_conv" is that sum on (c,d,z_conv), whichever factor each came from.
    if 1==1
        xc = Coordinate("x_conv", euclidean)
        ca, cb = eu[3], eu[4]                       #"c" and "d"
        Ac(i) = Tensor("A", i => lower, xc)
        Bc(i) = Tensor("B", i => lower, xc)
        ∂c(i, T) = PartialDerivative("∂", xc, i => lower)(T)
        cz = Index("z_conv", euclidean)             #sorts after both
        dBc(i, j) = ∂c(i, Bc(j)) - ∂c(j, Bc(i))     #(dβ)_{ij}, written out by hand
        ours = simplify(apply_derivative(
                   Ac(ca) ∧ exterior_derivative(Bc(cb), xc, cz)))
        hand = simplify(apply_derivative(
                   Ac(ca) * dBc(cb, cz) + Ac(cb) * dBc(cz, ca) + Ac(cz) * dBc(ca, cb)))
        @test isZeroTens(simplify(ours - hand))
        #and the opposite sign is genuinely wrong, not merely a different arrangement
        @test !isZeroTens(simplify(ours + hand))
    end

    #Graded Leibniz for two one-forms, d(α∧β) = dα∧β + (-1)^{deg α} α∧dβ, which for a
    #one-form α is a minus - written here exactly as a textbook has it. Both the wedge
    #and d hold their slots in canonical order, so all three terms are built on the same
    #one and no transposition is hidden in the signs.
    if 1==1
        xl = Coordinate("x_leib", euclidean)
        A1 = Tensor("A", eu[1]=>lower, xl)
        B1 = Tensor("B", eu[2]=>lower, xl)
        w  = A1 ∧ B1                        #form the wedge first, so its dummies are drawn
        zz = Index("z_leib", euclidean)     #then take a slot neither side is using
        lhs = simplify(apply_derivative(exterior_derivative(w, xl, zz)))
        rhs = simplify(apply_derivative(exterior_derivative(A1, xl, zz) ∧ B1
                                        - A1 ∧ exterior_derivative(B1, xl, zz)))
        @test isZeroTens(simplify(lhs - rhs))
        #the other sign really is wrong, so the minus above is being tested
        @test !isZeroTens(simplify(lhs -
            simplify(apply_derivative(exterior_derivative(A1, xl, zz) ∧ B1
                                      + A1 ∧ exterior_derivative(B1, xl, zz)))))
    end

    #The (-1)^{deg α} is carried explicitly, so it has to be right at each degree.
    #Checked at p=1, q=2, where the sign is a minus, and at p=2, q=1, where it is a plus.
    if 1==1
        xg = Coordinate("x_grade", euclidean)

        #p = 1, q = 2: (-1)^p = -1
        a1 = Tensor("A", eu[5] => lower, xg)
        b2 = AntisymmetricTensor("B", [eu[6] => lower, eu[7] => lower], xg)
        w1 = a1 ∧ b2
        z1 = Index("z_grade1", euclidean)
        @test isZeroTens(simplify(
            simplify(apply_derivative(exterior_derivative(w1, xg, z1)))
            - simplify(apply_derivative(exterior_derivative(a1, xg, z1) ∧ b2
                                        - a1 ∧ exterior_derivative(b2, xg, z1)))))

        #p = 2, q = 1: (-1)^p = +1
        c2 = AntisymmetricTensor("C", [eu[8] => lower, eu[9] => lower], xg)
        d1 = Tensor("D", eu[10] => lower, xg)
        w2 = c2 ∧ d1
        z2 = Index("z_grade2", euclidean)
        @test isZeroTens(simplify(
            simplify(apply_derivative(exterior_derivative(w2, xg, z2)))
            - simplify(apply_derivative(exterior_derivative(c2, xg, z2) ∧ d1
                                        + c2 ∧ exterior_derivative(d1, xg, z2)))))
    end

    ###########################################################################
    ### grading, across degrees and across repeated wedges
    ###########################################################################

    if 1==1
        six = IndexSet("sixdim", 6, (0,6), lower)
        s6 = Index(["s1","s2","s3","s4","s5","s6","s7","s8","s9"], six)
        five = IndexSet("fivedim", 5, (0,5), lower)
        s5 = Index(["v1","v2","v3","v4","v5","v6","v7"], five)
        #a p-form on the given slots, coordinate-dependent only where d needs it
        pform(nm, is) = length(is) == 1 ? Tensor(nm, [is[1] => lower]) :
                        AntisymmetricTensor(nm, [k => lower for k in is])
        pformx(nm, is, x) = length(is) == 1 ? Tensor(nm, is[1] => lower, x) :
                            AntisymmetricTensor(nm, [k => lower for k in is], x)

        #α∧β = (-1)^{pq} β∧α. Both orderings come back on the same canonical slots, so
        #the grading sign has nowhere to hide in a relabelling.
        for (α, β, pq) in [
                (pform("A", eu[1:1]), pform("B", eu[2:2]), 1),
                (pform("A", eu[1:1]), pform("B", eu[2:3]), 2),
                (pform("A", eu[1:2]), pform("B", eu[3:4]), 4),
                (pform("A", s6[1:1]), pform("B", s6[2:4]), 3),
                (pform("A", s6[1:2]), pform("B", s6[3:5]), 6),
                (pform("A", s6[1:3]), pform("B", s6[4:6]), 9)]
            @test isZeroTens(simplify((α ∧ β) - (-1)^pq * (β ∧ α)))
        end

        #(α∧β)∧γ = α∧(β∧γ). This is what catches a wedge renaming a slot it was given:
        #α∧β draws its dummies from the pool and they can land on the index γ carries,
        #and the two bracketings then come back on different slots.
        #Higher dimensions and a top form are in reach because `reduce_epsilon` expands
        #only the free slots; over the full n! determinant an n=6 triple did not finish.
        for (α, β, γ) in [
                (pform("A", eu[1:1]), pform("B", eu[2:2]), pform("C", eu[3:3])),
                (pform("A", eu[1:1]), pform("B", eu[2:2]), pform("C", eu[3:4])),
                (pform("A", eu[1:2]), pform("B", eu[3:3]), pform("C", eu[4:4])),
                (pform("A", s5[1:1]), pform("B", s5[2:2]), pform("C", s5[3:3])),
                (pform("A", s5[1:1]), pform("B", s5[2:3]), pform("C", s5[4:5])),  #r = n
                (pform("A", s6[1:1]), pform("B", s6[2:2]), pform("C", s6[3:3])),
                (pform("A", s6[1:1]), pform("B", s6[2:2]), pform("C", s6[3:4]))]
            @test isZeroTens(simplify(((α ∧ β) ∧ γ) - (α ∧ (β ∧ γ))))
        end
        #the slots really are the ones handed in, not fresh names picked to dodge a dummy
        @test Set(pr[1] for pr in Alakazam.free_indices_ordered(
                      (pform("A", eu[1:1]) ∧ pform("B", eu[2:2])) ∧ pform("C", eu[3:3]))) ==
              Set(eu[1:3])

        #α∧α vanishes for odd degree by the grading, and need not for even - F∧F is the
        #Pontryagin density, so a wedge that forced it to zero would be wrong
        @test isZeroTens(simplify(pform("A", eu[1:1]) ∧ pform("A", eu[2:2])))
        @test isZeroTens(simplify(pform("A", s6[1:3]) ∧ pform("A", s6[4:6])))
        @test !isZeroTens(simplify(pform("A", eu[1:2]) ∧ pform("A", eu[3:4])))

        #no form of degree above the dimension
        @test isZeroTens(simplify(pform("A", eu[1:2]) ∧ pform("B", eu[3:5])))
        @test isZeroTens(simplify(pform("A", eu[1:3]) ∧ pform("B", eu[4:6])))

        #d² = 0 at each degree, the property that fixes d's signs
        x4 = Coordinate("x_grade4", euclidean)
        x6 = Coordinate("x_grade6", six)
        for (ω, xx, s1, s2) in [
                (Tensor("f", x4),              x4, eu[11], eu[12]),
                (pformx("A", eu[1:1], x4),     x4, eu[11], eu[12]),
                (pformx("A", eu[1:2], x4),     x4, eu[11], eu[12]),
                (pformx("A", s6[1:3], x6),     x6, s6[7],  s6[8])]
            @test isZeroTens(simplify(apply_derivative(
                exterior_derivative(exterior_derivative(ω, xx, s1), xx, s2))))
        end

        #graded Leibniz at each degree pair, with the (-1)^p written out
        for (nm, α, β, p, xx, st) in [
                ("g1", pformx("A", eu[1:1], x4), pformx("B", eu[2:2], x4), 1, x4, euclidean),
                ("g2", pformx("A", eu[1:1], x4), pformx("B", eu[2:3], x4), 1, x4, euclidean),
                ("g3", pformx("A", eu[1:2], x4), pformx("B", eu[3:3], x4), 2, x4, euclidean),
                ("g4", pformx("A", s6[1:2], x6), pformx("B", s6[3:4], x6), 2, x6, six),
                ("g5", pformx("A", s6[1:1], x6), pformx("B", s6[2:4], x6), 1, x6, six),
                ("g6", pformx("A", s6[1:3], x6), pformx("B", s6[4:5], x6), 3, x6, six)]
            #names of their own: `w` and `zz` are already taken above, and assigning to
            #those from inside a loop is a soft-scope ambiguity Julia warns about
            wg = α ∧ β                      #wedge first, so its dummies are drawn
            zg = Index("z_" * nm, st)       #then a slot neither side is using
            @test isZeroTens(simplify(
                simplify(apply_derivative(exterior_derivative(wg, xx, zg)))
                - simplify(apply_derivative(
                      exterior_derivative(α, xx, zg) ∧ β
                      + (-1)^p * (α ∧ exterior_derivative(β, xx, zg))))))
        end

        #and over three factors, where the second sign compounds:
        #d(α∧β∧γ) = dα∧β∧γ + (-1)^p α∧dβ∧γ + (-1)^{p+q} α∧β∧dγ
        if 1==1
            α = pformx("A", eu[1:1], x4)
            β = pformx("B", eu[2:2], x4)
            γ = pformx("C", eu[3:3], x4)
            w = (α ∧ β) ∧ γ
            zt = Index("z_trip", euclidean)
            @test isZeroTens(simplify(
                simplify(apply_derivative(exterior_derivative(w, x4, zt)))
                - simplify(apply_derivative(
                      (exterior_derivative(α, x4, zt) ∧ β) ∧ γ
                      - (α ∧ exterior_derivative(β, x4, zt)) ∧ γ
                      + (α ∧ β) ∧ exterior_derivative(γ, x4, zt)))))
        end
    end

    #The Hodge star of a summand with no indices to infer the form index set from.
    #`⋆` on a lone tensor says so plainly; on a PRODUCT of index-free tensors the
    #index set came back as `nothing` and reached `_star_term`, which has no
    #method for it, so the user saw a MethodError naming a private function
    #instead of the message the same situation gives one factor along.
    if 1==1
        xstar = Coordinate("x_star", euclidean)
        sq = Tensor("q_star", xstar)
        sh = Tensor("h_star", xstar)

        #a lone index-free tensor already asked for the index set
        @test_throws Alakazam.InvalidDifferentialForm ⋆(sq)
        #and a product of them must ask the same way, not throw a MethodError
        @test_throws Alakazam.InvalidDifferentialForm ⋆(sq * sh)
        @test_throws Alakazam.InvalidDifferentialForm ⋆(2 * sq * sh)

        #given the index set, both take the star and agree with each other
        @test !isZeroTens(simplify(⋆(sq, euclidean)))
        @test isZeroTens(simplify(⋆(sq * sh, euclidean) - ⋆(sh * sq, euclidean)))
        #and starring a 0-form twice recovers it, which is the round trip that
        #ran into the MethodError: the middle result is a top form, and starring
        #that leaves an index-free summand for the third star to look at
        @test isZeroTens(simplify(⋆(⋆(sq, euclidean)) - sq))
        @test_throws Alakazam.InvalidDifferentialForm ⋆(simplify(⋆(⋆(sq, euclidean))) * sh)
    end

end
