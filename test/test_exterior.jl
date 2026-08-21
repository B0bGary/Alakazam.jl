@testset "Exterior Derivative" begin

    #setup
    if 1==1
        ext = IndexSet("spacetime", 4, (0,4), lower)
        ei = Index(["xa","xb","xc","xd","xe","xf","xg","xh","xi","xj","xk","xl"], ext)
        xa, xb, xc, xd, xe = ei[1], ei[2], ei[3], ei[4], ei[5]
        xx = Coordinate("xext", [ext])

        ff = Tensor("f", xx)                                  #a zero-form
        AA = Tensor("A", [xb=>lower], xx)                     #a one-form
        FF = AntisymmetricTensor("F", [xc=>lower, xd=>lower], xx)  #a two-form
    end

    #d of a zero-form is the gradient
    @test exterior_derivative(ff, xx, xa) == TensorExpression(PartialDerivative("∂", xx, xa=>lower)(ff))

    #d of a one-form is the curl, two terms with opposite signs
    @test length(summands(exterior_derivative(AA, xx, xa))) == 2

    #d of a two-form is three terms, alternating
    @test length(summands(exterior_derivative(FF, xx, xa))) == 3

    if 1==1
        dab = exterior_derivative(Tensor("A", [xb=>lower], xx), xx, xa)
        dba = exterior_derivative(Tensor("A", [xa=>lower], xx), xx, xb)
        @test isZeroTens(simplify(apply_derivative(dab - dba)))

        #and that component is the curl, which is manifestly antisymmetric in xa, xb
        hand = PartialDerivative("∂", xx, xa=>lower)(Tensor("A", [xb=>lower], xx)) -
               PartialDerivative("∂", xx, xb=>lower)(Tensor("A", [xa=>lower], xx))
        @test isZeroTens(simplify(apply_derivative(dab - hand)))
        @test !isZeroTens(simplify(apply_derivative(dab + hand)))
    end

    #d is nilpotent, the property that fixes the signs. the derivative has to
    #be applied for the terms to meet
    @test isZeroTens(simplify(apply_derivative(
        exterior_derivative(exterior_derivative(ff, xx, xb), xx, xa))))
    @test isZeroTens(simplify(apply_derivative(
        exterior_derivative(exterior_derivative(AA, xx, xa), xx, xe))))
    @test isZeroTens(simplify(apply_derivative(
        exterior_derivative(exterior_derivative(FF, xx, xb), xx, xa))))

    #so F = dA satisfies dF = 0, the Bianchi identity of a gauge field
    if 1==1
        Fda = exterior_derivative(Tensor("A", [xb=>lower], xx), xx, xa)
        @test isZeroTens(simplify(apply_derivative(exterior_derivative(Fda, xx, xe))))
    end

    #a form of degree above the dimension has nowhere to go
    if 1==1
        top = AntisymmetricTensor("T", [ei[1]=>lower, ei[2]=>lower,
                                        ei[3]=>lower, ei[4]=>lower], xx)
        @test isZeroTens(exterior_derivative(top, xx, ei[5]))
    end

    #the index for the new slot is generated when none is given
    @test length(summands(exterior_derivative(AA, xx))) == 2

    #the unicode alias is the same function
    @test 𝖽(AA, xx, xa) == exterior_derivative(AA, xx, xa)

    #a derivative may be named, so several can sit in one expression
    @test name(terms(summands(exterior_derivative(ff, xx, xa; dname="D"))[1][2])[1]) == "D"

    #a Grassmann-odd form: the form slots are spacetime, a spinor index rides along as a spectator.
    #The exterior derivative's signs are indexed by form degree, and Grassmann parity is carried
    #through untouched, so d^2=0 must still hold and the parity must not change.
    if 1==1
        spin = IndexSet("spinor", 2, lower, EpsilonTensor, fermionic)
        sal  = Index(["α","β"], spin)[1]

        #declared with the SPACETIME slot first, and with the SPINOR slot first: the same form
        psi_tf = Tensor("ψ", [xb=>lower, sal=>lower], xx)
        psi_sf = Tensor("ψ", [sal=>lower, xb=>lower], xx)

        @test grading(psi_tf) == 1 && grading(psi_sf) == 1

        #the form set is chosen from the slots that are directions of the coordinate, so the spinor
        #set is never nominated however the tensor was declared.
        @test Alakazam._form_index_set(psi_sf, xx, nothing) == ext
        @test Alakazam._form_index_set(psi_tf, xx, nothing) == ext

        #d is parity preserving, and leaves the spinor slot where it was
        dpsi = exterior_derivative(psi_sf, xx, xa)
        @test grading(dpsi) == 1
        @test length(summands(dpsi)) == 2
        @test all(pr[1].index_set == spin for pr in indices(terms(summands(dpsi)[1][2])[1])
                  if pr[1].index_set == spin)

        #the tensor and expression paths agree, neither keying on declaration order
        @test exterior_derivative(TensorExpression(psi_sf), xx, xa) == dpsi

        #nilpotency survives the grading, for degree 1 and degree 2
        @test isZeroTens(simplify(apply_derivative(
            exterior_derivative(exterior_derivative(psi_sf, xx, xa), xx, xc))))

        #an odd 2-form: antisymmetric on its spacetime slots only, spinor slot a spectator
        Psi = Tensor("Ψ", [sal=>lower, xa=>lower, xb=>lower], xx, YoungTableau([1,1], [2,3], 3))
        @test isForm(Psi, ext) && grading(Psi) == 1
        @test length(summands(exterior_derivative(Psi, xx, xc))) == 3
        @test isZeroTens(simplify(apply_derivative(
            exterior_derivative(exterior_derivative(Psi, xx, xc), xx, xd))))
    end

    #d with the coordinate left implicit, so that 𝖽(α∧β) can be written directly. The set
    #comes off the slots and the coordinate from what the form is a function of, falling
    #back to the registry when it is a function of nothing.
    if 1==1
        Ai = Tensor("Ai", [xb=>lower], xx)
        Bi = Tensor("Bi", [xc=>lower], xx)

        #inferring the coordinate gives the same answer as naming it
        @test isZeroTens(simplify(apply_derivative(𝖽(Ai, xa)
                                                   - 𝖽(Ai, xx, xa))))
        @test isZeroTens(simplify(apply_derivative(exterior_derivative(Ai, xa)
                                                   - exterior_derivative(Ai, xx, xa))))

        #the reported case: a wedge, with no coordinate to hand
        @test !isZeroTens(simplify(apply_derivative(𝖽(Ai ∧ Bi))))
        @test isZeroTens(simplify(apply_derivative(𝖽(𝖽(Ai)))))

        #a form that is a function of nothing differentiates to zero rather than refusing:
        #there is no dependence to read, and every coordinate gives the same answer
        @test isZeroTens(simplify(apply_derivative(𝖽(Tensor("Ci", [xb=>lower])))))

        #but a set carrying no coordinate at all has nothing to differentiate along
        bare = IndexSet("bare_ext", 4, (0,4), lower)
        @test_throws Alakazam.InvalidDifferentialForm 𝖽(
            Tensor("Di", [Index("bi", bare)=>lower]))
    end

end
