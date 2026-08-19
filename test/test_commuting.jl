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

###############################################################################
# Edge cases for the commutation-relation key.
###############################################################################

@testset "Commutator keys" begin

    clear_global!()
    clear_commutators!()

    st  = IndexSet("spacetime", 4, lower, Metric)
    st2 = IndexSet("spacetime_other", 4, lower, Metric)      # same kind, different name
    sp  = IndexSet("spinor", 4, lower, fermionic, SymplecticForm)
    spc = IndexSet("spinor_conj", 4, lower, fermionic, SymplecticForm)

    μ, ν, ρ, σ = Index(["μ", "ν", "ρ", "σ"], st)
    p, q       = Index(["p", "q"], st2)
    α, β       = Index(["α", "β"], sp)
    # conjugate-set indices, accented as users spell them
    αdot = Index("α" * Alakazam.COMBINING_DOT, spc)
    αbar = Index("α" * Alakazam.COMBINING_MACRON, spc)

    key = Alakazam._ckey

    ###########################################################################
    # what must NOT be distinguished, or a declared relation is lost
    ###########################################################################

    # index names are irrelevant: a relation for X^μ must hold for X^ν
    @test key(Mat("X", [μ => upper])) == key(Mat("X", [ν => upper]))

    # index POSITION is irrelevant: raising or lowering does not change what an object is
    @test key(Mat("X", [μ => upper])) == key(Mat("X", [μ => lower]))

    # index-set NAME is irrelevant when the kind agrees. This is the one that broke the built-in
    # Clifford identities: they are registered against one set, and a gamma built on a set the
    # user named themselves must still find them.
    @test key(Gamma([μ => upper])) == key(Gamma([p => upper]))

    # writing order is irrelevant, since the signature is a sorted multiset
    @test key(Gamma([μ => upper, α => lower, β => lower])) ==
          key(Gamma([α => lower, β => lower, μ => upper]))

    ###########################################################################
    # what MUST be distinguished, or one declaration overwrites another
    ###########################################################################

    # rank: γ^a and γ^{ab} are different objects sharing a name
    @test key(Gamma([μ => upper])) != key(Gamma([μ => upper, ν => upper]))
    @test key(Gamma([μ => upper, ν => upper])) !=
          key(Gamma([μ => upper, ν => upper, ρ => upper]))

    # kind: a spacetime index is not a spinor index
    @test key(Mat("X", [μ => upper])) != key(Mat("X", [α => upper]))

    # matrix form vs component form of one symbol
    @test key(Gamma([μ => upper])) != key(Gamma([μ => upper, α => lower, β => lower]))

    # ACCENTS: σ_{αα̇} carries one index from each bundle and σ_{αβ} two from one. These support
    # genuinely different relations. σ^{μν} is self-dual on undotted indices while σ̄^{μν} is
    # anti-self-dual on dotted ones.
    @test key(Pauli([μ => upper, α => lower, αdot => lower])) !=
          key(Pauli([μ => upper, α => lower, β => lower]))

    # and different accents stay apart, for a user who writes both α̇ and ᾱ with distinct meanings
    @test key(Mat("Y", [αdot => lower])) != key(Mat("Y", [αbar => lower]))
    @test key(Mat("Y", [αdot => lower]) ) == key(Mat("Y", [αdot => lower]))

    ###########################################################################
    # end to end: declare, look up, and confirm the lookups agree with each other
    ###########################################################################

    A1 = Mat("A", [μ => upper])
    B1 = Mat("B", [ν => upper])
    @commutator [A1, B1] = zero()

    @test has_commutation_rel(A1, B1)
    @test Alakazam.can_swap(A1, B1)
    @test get_commutator(A1, B1)[1] == 1

    # the same relation is found through differently-named indices of the same set
    @test has_commutation_rel(Mat("A", [ρ => upper]), Mat("B", [σ => upper]))
    @test Alakazam.can_swap(Mat("A", [ρ => lower]), Mat("B", [σ => lower]))

    ###########################################################################
    # two relations for one name coexist rather than overwrite
    ###########################################################################

    G1 = Gamma([μ => upper])                    # rank 1
    G2 = Gamma([μ => upper, ν => upper])        # rank 2, same name
    D1 = Mat("D", [ρ => upper])     # an index of its own: γ^{μ}D^{μ} is not a valid dummy pair
    # the two differ in grading, so they are distinguishable after the fact even though both
    # right-hand sides vanish.
    @commutator [G1, D1] = zero()          # commutator, so G1 and D1 freely exchange
    @anticommutator {G2, D1} = zero()      # ANTIcommutator: same names, different entry

    @test get_commutator(G1, D1)[1] == 1                # recovered as a commutator
    @test get_commutator(G2, D1)[1] == -1               # rank-2 entry survived independently
    @test Alakazam.can_swap(G1, D1)
    @test !Alakazam.can_swap(G2, D1)                    # anticommuting: not freely swappable

    # removing one leaves the other in place
    remove_commutator!(G2, D1)
    @test !has_commutation_rel(G2, D1)
    @test has_commutation_rel(G1, D1)

    ###########################################################################
    # same-name pairs still resolve, as `generate_conformal_algebra` relies on for [P,P] = 0
    ###########################################################################

    P1 = Mat("P", [μ => upper])
    P2 = Mat("P", [ν => upper])
    @commutator [P1, P2] = zero()
    @test has_commutation_rel(P1, P2)
    @test Alakazam.can_swap(Mat("P", [ρ => upper]), Mat("P", [σ => upper]))

    clear_commutators!()
    clear_global!()
end
