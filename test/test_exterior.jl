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

    #the result is totally antisymmetric, so swapping two slots flips the sign
    if 1==1
        dab = exterior_derivative(Tensor("A", [xb=>lower], xx), xx, xa)
        dba = exterior_derivative(Tensor("A", [xa=>lower], xx), xx, xb)
        @test isZeroTens(simplify(apply_derivative(dab + dba)))
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

end
