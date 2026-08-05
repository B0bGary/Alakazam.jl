@testset "Manipulations and Simplifications" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4,lower)
        spinor_indices=IndexSet("spinor_indices", 2,lower,fermionic)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2,lower,fermionic)
        @indices α:γ ∈ spinor_indices
        χ = Index("χ", spinor_indices)
        αd = Index("̇α", spinor_conj_indices)
        βd = Index("̇β", spinor_conj_indices)
        γd = Index("̇γ", spinor_conj_indices)
        χd = Index("̇χ", spinor_conj_indices)
        μ = Index("μ", spacetime_indices)
        ν = Index("ν", spacetime_indices)
        ρ = Index("ρ", spacetime_indices)
        τ = Index("τ", spacetime_indices)
        σ = Index("σ", spacetime_indices)
        κ = Index("κ", spacetime_indices)
        λ = Index("λ", spacetime_indices)
        ζ = Index("ζ", spacetime_indices)
        ϕ = Index("ϕ", spacetime_indices)

        x = Coordinate("x", [spacetime_indices])
        θ = Coordinate("θ", [spinor_indices])
        θbar = Coordinate("Θ", [spinor_conj_indices])
        custom_sort_order!(["θ","Θ","ϵ"])
    end
        
    #needs to commute partials, and relabel dummies inside partials, and commute Bs
        A=Tensor("A", [μ=>upper],x)
        B=Tensor("B", [ν=>upper,σ=>upper],x)
        B2=Tensor("B", [ρ=>upper,κ=>upper],x)
        B3=Tensor("B", [ρ=>upper,σ=>upper],x)
        B4=Tensor("B", [ν=>upper,κ=>upper],x)
        e1=B2*B*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        e2=B3*B4*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        es= Alakazam.simplify(e1+e2)
        ex=2*Tensor("B", [ρ=>upper,σ=>upper],x)*Tensor("B", [ν=>upper,κ=>upper],x)*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        # @test Alakazam.hard_simplify(es-ex)==Alakazam. zero()
        @test es≅ex

        #same, but now Bs are SUSY, so the swap should make everything  zero().
        A=Tensor("A", [μ=>upper],x)
        B=Tensor("B", [ν=>upper,α=>upper],x)
        B2=Tensor("B", [ρ=>upper,β=>upper],x)
        B3=Tensor("B", [ρ=>upper,α=>upper],x)
        B4=Tensor("B", [ν=>upper,β=>upper],x)
        e1=B2*B*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        e2=B3*B4*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        es= Alakazam.simplify(e1+e2)
        @test es≅Alakazam. zero()


        #commute an ordinary PD with SUSY PD
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        e1=PartialDerivative("∂",PartialDerivative("∇",A, [x], [ν => lower]), [θ], [α => lower])
        e2=PartialDerivative("∇",PartialDerivative("∂",A, [θ], [α => lower]), [x], [ν => lower])
        @test e1≅e2

        #commute a SUSY with SUSY PD
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        e1=PartialDerivative("∂",PartialDerivative("∂",A, [θ], [β => lower]), [θ], [α => lower])
        e2=-PartialDerivative("∂",PartialDerivative("∂",A, [θ], [α => lower]), [θ], [β => lower])
        @test e1≅e2


        #commute a SUSY with different SUSY PD
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        e1=PartialDerivative("𝛛",PartialDerivative("∂",A, [θ], [β => lower]), [θbar], [αd => lower])
        e2=-PartialDerivative("∂",PartialDerivative("𝛛",A, [θbar], [αd => lower]), [θ], [β => lower])
        @test e1≅e2

        #commute a tensor past susy diff op
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [ϕ=>upper],[x,θ,θbar])
        e1=B*PartialDerivative("𝛛",PartialDerivative("∂",A, [θ], [β => lower]), [θbar], [αd => lower])
        e2=-PartialDerivative("∂",PartialDerivative("𝛛",A, [θbar], [αd => lower]), [θ], [β => lower])*B
        @test e1≅e2

        #commute a susy tensor past susy diff op
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [γ=>upper],[x,θ,θbar])
        e1=B*PartialDerivative("𝛛",PartialDerivative("∂",A, [θ], [β => lower]), [θbar], [αd => lower])
        e2=-PartialDerivative("∂",PartialDerivative("𝛛",A, [θbar], [αd => lower]), [θ], [β => lower])*B
        @test e1≅e2

        #commute a susy tensor past differnt susy diff op
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [γ=>upper],[x,θ,θbar])
        e1=B*PartialDerivative("∂",PartialDerivative("∇",A, [x], [ν => lower]), [θ], [α => lower])
        e2=-PartialDerivative("∇",PartialDerivative("∂",A, [θ], [α => lower]), [x], [ν => lower])*B
        @test e1≅e2


        
        #products of diff ops
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [γ=>upper],[x,θ,θbar])
        e1=PartialDerivative("∇",A, [x], [ν => lower])*PartialDerivative("∂",B, [θ], [α => lower])
        e2=PartialDerivative("∂",B, [θ], [α => lower])*PartialDerivative("∇",A, [x], [ν => lower])
        @test e1≅e2

        #products of nested diff ops, with different levels
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [γ=>upper],[x,θ,θbar])
        e1=PartialDerivative("∇",A, [x], [ν => lower])*PartialDerivative("∇",PartialDerivative("∂",B, [θ], [α => lower]), [x], [μ => lower])
        e2=PartialDerivative("∇",PartialDerivative("∂",B, [θ], [α => lower]), [x], [μ => lower])*PartialDerivative("∇",A, [x], [ν => lower])
        @test e1≅e2

        ###############################################
        ### isomorphisms must not swap matrices that carry a commutation relation
        ###############################################
        @matrix CmMx
        Symbolics.@syms hcm
        Pc = CmMx("P"); Qc = CmMx("Q"); Rc = CmMx("R"); Sc = CmMx("S")
        @commutator [Pc, Qc] = hcm*Sc

        @test !Alakazam.can_swap(Pc, Qc)
        # R has no relation registered with P, so it is NOT assumed to commute with it
        @test !Alakazam.can_swap(Pc, Rc)
        @test !Alakazam.has_commutation_rel(Pc, Rc)
        # declaring the commutator to vanish is what makes a pair freely exchangeable
        Tc = CmMx("T")
        @commutator [Pc, Tc] = zero()
        @test Alakazam.can_swap(Pc, Tc)
        # the graph records names, indices and tableau symmetries but not commutation relations, so
        # an isomorphism could otherwise identify PQ with QP and cancel a nonzero commutator
        # PQ - QP is the commutator, and sorting emits it rather than cancelling the two terms
        @test Alakazam.simplify(Pc*Qc - Qc*Pc) == expand_commutators(Alakazam.Lie(Pc, Qc, 0))
        @test Alakazam.hard_simplify(Pc*Qc - Qc*Pc) == hcm*Sc
        # an undeclared pair does not cancel: PR and RP are different products
        @test !isZeroTens(Alakazam.hard_simplify(Pc*Rc - Rc*Pc))
        # one declared to vanish does
        @test isZeroTens(Alakazam.hard_simplify(Pc*Tc - Tc*Pc))
        @test isZeroTens(Alakazam.hard_simplify(Pc*Qc - Pc*Qc))
        @test Alakazam.hard_simplify(Pc*Qc + Pc*Qc) == 2*Pc*Qc

        ###############################################
        ### sorting orders the matrices and keeps the remainder
        ###############################################
        # QP = PQ - [P,Q]. This used to depend on whether a custom sort order happened to be set.
        @test sort(TensorExpression(Qc*Pc)) == Pc*Qc - hcm*Sc
        @test sort(TensorExpression(Pc*Qc)) == TensorExpression(Pc*Qc)   # already ordered
        # a pair declared to commute is reordered with no remainder
        @test sort(TensorExpression(Tc*Pc)) == TensorExpression(Pc*Tc)
        # an undeclared pair is left exactly as written
        @test sort(TensorExpression(Rc*Pc)) == TensorExpression(Rc*Pc)
        # QPR = PQR - [P,Q]R, and the remainder keeps S where the commutator was formed: S sits
        # to the left of R because nothing licenses moving it past
        @test Alakazam.simplify(Qc*Pc*Rc) == Alakazam.simplify(Pc*Qc*Rc - hcm*Sc*Rc)

end
