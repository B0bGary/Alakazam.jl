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

    #the overall sign is the signature: ϵϵ carries a minus in Lorentzian
    @test isZeroTens(simplify(
        (Tensor("A",[la=>lower]) ∧ Tensor("B",[lb=>lower]))
        + (Tensor("A",[la=>lower])*Tensor("B",[lb=>lower])
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

end
