@testset "Tensor Equalities" begin


    ####setup

    spacetime_indices = IndexSet("spacetimelc", 2, lower)
    flavour_indices = IndexSet("flavour_indices", 2, upper,fermionic)
    lightcone_indices = IndexSet("lightcone", 2, lower)
    spinor_indices=IndexSet("spinor_indices", 2, lower)
    spinor_conj_indices=IndexSet("spinor_conj_indices", 2, lower)
    @indices a:z ∈ flavour_indices
    x = Coordinate("x", [spacetime_indices, spacetime_indices])
    θ = Coordinate("θ", [spacetime_indices])
    θbar = Coordinate("Θ", [spacetime_indices])

    #############

    ####
    #simple tensors
    @test Tensor("A", [i=>upper, j=>upper])≅Tensor("A", [i=>upper, j=>upper]) #should be exactly equal
    @test ! (Tensor("A", [i=>upper, j=>upper])≅Tensor("A", [j=>upper, i=>upper])) #not equal since no syms
    @test SymmetricTensor("A", [i=>upper, j=>upper])≅SymmetricTensor("A", [j=>upper, i=>upper]) # equal since syms
    @test !(SymmetricTensor("A", [i=>upper, j=>upper])≅SymmetricTensor("B", [j=>upper, i=>upper])) # not equal since different names
    @test !(AntisymmetricTensor("A", [i=>upper, j=>upper]) ≅ AntisymmetricTensor("A", [j=>upper, i=>upper])) # not equal since antisym
    @test !(AntisymmetricTensor("A", [i=>upper, i=>lower]) ≅ AntisymmetricTensor("A", [i=>lower, i=>upper])) # not equal since antisym
    #over an ε the see-saw itself costs the sign
    let es = IndexSet("tterm_eps", 2, upper, EpsilonTensor), e1 = Index("e1", es)
        @test !(AntisymmetricTensor("A", [e1=>upper, e1=>lower]) ≅ AntisymmetricTensor("A", [e1=>lower, e1=>upper]))
        @test Alakazam.same_tensor_mod_dummies(AntisymmetricTensor("A", [e1=>upper, e1=>lower]), AntisymmetricTensor("A", [e1=>lower, e1=>upper])) == -1
    end

    ####
    #simple products and contractions
    #should be not equal by antisym
    @test !((AntisymmetricTensor("A", [i=>upper, j=>upper], x)*SymmetricTensor("B", [i=>lower, k=>lower], x)) ≅ (AntisymmetricTensor("A", [j=>upper, i=>upper], x)*SymmetricTensor("B", [i=>lower, k=>lower], x)))
    #should be equal since contractions imply equal to  zero(), even tho they look unequal by antisymmetry.
    @test (AntisymmetricTensor("A", [i=>upper, j=>upper], x)*SymmetricTensor("B", [i=>lower, j=>lower], x))≅(AntisymmetricTensor("A", [j=>upper, i=>upper], x)*SymmetricTensor("B", [i=>lower, j=>lower], x))
    #equal by double antisym
    @test (AntisymmetricTensor("A", [i=>upper, j=>upper], x)*AntisymmetricTensor("B", [i=>lower, j=>lower], x)) ≅ (AntisymmetricTensor("A", [j=>upper, i=>upper], x)*AntisymmetricTensor("B", [j=>lower, i=>lower], x))


    ####
    #testing the grading and swaps 
    #should be equal by grading (even*odd=even)
    @test (Tensor("A", [i=>upper], x)*Tensor("B", [j=>lower, k=>lower], x))≅(Tensor("B", [j=>lower, k=>lower], x) * Tensor("A", [i=>upper], x))
    #should be equal by grading (even*even=even)
    @test (Tensor("A", [i=>upper, t=>upper], x)*Tensor("B", [j=>lower, k=>lower], x))≅(Tensor("B", [j=>lower, k=>lower], x) * Tensor("A", [i=>upper, t=>upper], x))
    #should be not equal by grading (odd*odd=odd) 
    @test ! ((Tensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [j=>lower, k=>lower, l=>lower], x))≅(Tensor("B", [j=>lower, k=>lower, l=>lower], x) * Tensor("A", [i=>upper, t=>upper, s=>upper], x)))

    #should be equal by grading plus antisymm (odd*odd=odd) 
    @test (AntisymmetricTensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [j=>lower, k=>lower, l=>lower], x))≅(Tensor("B", [j=>lower, k=>lower, l=>lower], x) * AntisymmetricTensor("A", [i=>upper, s=>upper, t=>upper], x))
    #should be equal by two graded swaps
    @test (Tensor("A", [i=>upper], x)*Tensor("B", [j=>lower], x)*Tensor("C", [k=>lower], x))≅(Tensor("C", [k=>lower], x)*Tensor("A", [i=>upper], x)*Tensor("B", [j=>lower], x))
    #should be not equal
    @test ! ((Tensor("A", [i=>upper], x)*Tensor("B", [j=>lower], x)*Tensor("C", [k=>lower], x))≅(Tensor("C", [k=>lower], x)*Tensor("B", [j=>lower], x)*Tensor("A", [i=>upper], x)))

    #should be equal by grading plus antisymm thanks to the contraction relabelling
    @test (AntisymmetricTensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [t=>lower], x)*Tensor("B", [s=>lower], x))≅(AntisymmetricTensor("A", [i=>upper, s=>upper, t=>upper], x)*Tensor("B", [s=>lower], x)*Tensor("B", [t=>lower], x))

    #should be equal since equal to  zero() via contraction
    #sort of have to use sort here... i think this is a bad case, not sure how to get around  zero() both sides are really equal to  zero()
    @test (sort(SymmetricTensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [t=>lower], x)*Tensor("B", [s=>lower], x)))≅(SymmetricTensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [s=>lower], x)*Tensor("B", [t=>lower], x))

    #should be equal since the tableaux encodes totally symmetric
    @test (Tensor("A", [i=>upper, j => upper, k=>upper, l=>upper], x, YoungTableau([4], [1, 2, 3, 4]))≅Tensor("A", [j=>upper, i => upper, l=>upper, k=>upper], x, YoungTableau([4], [1, 2, 3, 4])))

    @test (1*Tensor("A", [i=>upper, j => upper, k=>upper, l=>upper], x, YoungTableau([1, 1, 1, 1], [1, 2, 3, 4]))≅-Tensor("A", [j=>upper, i => upper, k=>upper, l=>upper], x, YoungTableau([1, 1, 1, 1], [1, 2, 3, 4])))


    ###############################################
    ### type is a discriminator on every "same tensor" path, not only ==
    ###############################################
    # A Boson, a Fermion, a Gamma and a plain Tensor of the same name and indices are different by type
    let ty = IndexSet("ttype", 4, upper), a = Index("ttya", ty), y = Coordinate("ttyy", ty),
        Tn = Tensor("Q", [a => lower], y), Bo = Boson("Q", [a => lower], y),
        Fe = Fermion("Q", [a => lower], y)

        @test typeof(Tn) != typeof(Bo) != typeof(Fe)
        @test Tn != Bo && Tn != Fe && Bo != Fe
        @test hash(Tn) != hash(Bo)
        @test Alakazam.same_tensor_mod_syms(Tn, Bo) == (false, 0)
        @test Alakazam.same_tensor_mod_syms(Tn, Tn)[1]
        @test !Alakazam.same_tensorterm_mod_syms(TensorTerm(Tn), TensorTerm(Bo))
        @test Alakazam._graph_vertex_name(Tn) != Alakazam._graph_vertex_name(Bo)
        # the merge paths: a difference of two different types is not zero
        @test length(summands(Alakazam.simplify(TensorExpression(Tn) + TensorExpression(Bo)))) == 2
        @test !isZeroTens(Alakazam.simplify(TensorExpression(Tn) - TensorExpression(Bo)))
        @test !isZeroTens(Alakazam.hard_simplify(TensorExpression(Tn) - TensorExpression(Bo)))
        # while the same type still merges, so nothing has been made over-cautious
        @test Alakazam.simplify(TensorExpression(Tn) + TensorExpression(Tn)) == 2 * Tn
        @test isZeroTens(Alakazam.simplify(TensorExpression(Bo) - TensorExpression(Bo)))
        # a Gamma is not a tensor that happens to be called γ
        let g = Gamma(a => lower), tg = Tensor("γ", [a => lower])
            @test !isZeroTens(Alakazam.hard_simplify(TensorExpression(g) - TensorExpression(tg)))
        end
    end

    #a covariant derivative is not a partial, and root_oprnd used to have a
    #method for the partial alone, so it raised on anything carrying one
    if 1==1
        cdx = Coordinate("xcd", [spacetime_indices])
        uu = Tensor("ucd", cdx)
        cd(ii, pp, arg) = Derivative("∇cd", cdx, ii=>pp)(arg)
        #the tensor at the bottom of a covariant nest
        @test Alakazam.root_oprnd(cd(i, lower, cd(i, upper, uu))) == TensorExpression(uu)
        #a nest of two different derivative kinds resolves to its innermost operand
        @test Alakazam.root_oprnd(cd(i, lower, PartialDerivative("∂cd", cdx, j=>lower)(uu))) ==
              TensorExpression(uu)
        #≅ answers rather than raising
        @test cd(i, lower, cd(i, upper, uu)) ≅ cd(i, lower, cd(i, upper, uu))
    end

    #renaming dummies never moves a slot: X_i^{si} and X_i^{is} hold the free s in different slots
    let su = IndexSet("tterm_su2", upper, SLAlgebra(2)), (iq, jq, sq, mq, nq, oq) = Index(["iq", "jq", "sq", "mq", "nq", "oq"], su)
        xq = Coordinate("tterm_x", [spacetime_indices])
        μd = Index("μd", spacetime_indices)
        X(p...) = Tensor("X", Pair{IndexSuperType,IndexPosition}[p...], xq)
        ∂q(T) = PartialDerivative("∂", T, [xq], [μd => lower])
        Xa, Xb = X(iq => lower, sq => upper, iq => upper), X(iq => lower, iq => upper, sq => upper)
        @test Alakazam.same_tensor_mod_dummies(Xa, Xb) == 0
        @test Alakazam.same_tensor_mod_dummies(Xa, X(jq => lower, sq => upper, jq => upper)) == 1
        @test Alakazam.same_tensor_mod_dummies(Xa, X(iq => upper, sq => upper, iq => lower)) == -1   #ε see-saw
        #which used to cancel inside a derivative, and in a pair of ε's contracted into one
        @test !isZeroTens(Alakazam.simplify(∂q(Xa) - ∂q(Xb)))
        @test isZeroTens(Alakazam.simplify(∂q(Xa) - ∂q(X(jq => lower, sq => upper, jq => upper))))
        @test !isZeroTens(Alakazam.simplify(EpsilonTensor([sq => upper, mq => upper]) * EpsilonTensor([nq => lower, oq => lower]) *
                                            ∂q(X(mq => lower, oq => upper, nq => upper))))
    end

    #a relabelling the simplify fast path sees only through a slot symmetry
    let st = IndexSet("tterm_st", 4, upper, Metric), (μq, νq) = Index(["μq", "νq"], st)
        A(i) = Tensor("A", [i => upper]); B(i) = Tensor("B", [i => upper])
        S = SymmetricTensor("S", [μq => lower, νq => lower]); F = AntisymmetricTensor("F", [μq => lower, νq => lower])
        @test isZeroTens(Alakazam.simplify(S * A(μq) * B(νq) - S * A(νq) * B(μq)))
        @test isZeroTens(Alakazam.simplify(F * A(μq) * B(νq) + F * A(νq) * B(μq)))
        @test Alakazam.simplify(F * A(μq) * B(νq) - F * A(νq) * B(μq)) == Alakazam.simplify(2 * F * A(μq) * B(νq))
        T = Tensor("T", [μq => lower, νq => lower])          #with no symmetry nothing merges
        @test length(summands(Alakazam.simplify(T * A(μq) * B(νq) - T * A(νq) * B(μq)))) == 2
    end

    #a trace is over one set: T^i_i over an ε is not T^a_a over a δ, which the signature once merged
    let es = IndexSet("tterm_tre", upper, SLAlgebra(2)), ds = IndexSet("tterm_trd", 3, upper, Metric, (0, 3)),
        ie = Index("ie", es), ad = Index("ad", ds)
        T(p...) = TensorExpression(Tensor("T", Pair{IndexSuperType,IndexPosition}[p...]))
        @test length(summands(Alakazam.simplify(T(ie => upper, ie => lower) - T(ad => upper, ad => lower)))) == 2
        #and the term graph, which named no dummy, let hard_simplify and ≅ merge them too
        @test !isZeroTens(Alakazam.hard_simplify(T(ie => upper, ie => lower) - T(ad => upper, ad => lower)))
        @test !(Tensor("T", [ie => upper, ie => lower]) ≅ Tensor("T", [ad => upper, ad => lower]))
        @test isZeroTens(Alakazam.simplify(T(ie => upper, ie => lower) + T(ie => lower, ie => upper)))
        #a pair repeated at one level is not a contraction at all
        @test_throws SameLevelDummyError Alakazam.same_tensor_mod_dummies(Tensor("T", [ie => upper, ie => upper]),
                                                                          Tensor("T", [ie => lower, ie => lower]))
    end

end