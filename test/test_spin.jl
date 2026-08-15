@testset "Spin algebras" begin

    @testset "signature, dimension and the classical realisations" begin
        @test spacetime_dimension(Spin(2, 1)) == 3
        @test spacetime_dimension(Spin(1, 9)) == 10
        @test spacetime_dimension(Spin(1, 10)) == 11
        # Spin(d) is mostly-plus Lorentzian, matching the IndexSet default
        @test Alakazam.signature(Spin(4)) == (1, 3)
        @test spacetime_dimension(Spin(11)) == 11

        # the exceptional isomorphisms, in both orderings since Spin(p,q) = Spin(q,p)
        @test Alakazam.classical(Spin(2, 1)) isa SLAlgebra{Alakazam.RealNumberSet}
        @test Alakazam.classical(Spin(1, 3)) isa SLAlgebra{Alakazam.ComplexNumberSet}
        @test Alakazam.classical(Spin(3, 1)) isa SLAlgebra{Alakazam.ComplexNumberSet}
        @test Alakazam.classical(Spin(5, 1)) isa SLAlgebra{Alakazam.QuaternionNumberSet}
        @test Alakazam.classical(Spin(9, 1)) isa SLAlgebra{Alakazam.OctonionNumberSet}
        # and they run out: no classical algebra for the 32-component spinor of 11D
        @test Alakazam.classical(Spin(1, 10)) === nothing
    end

    @testset "reality of the minimal spinor" begin
        # These four are fixed independently by the exceptional isomorphisms SL(2,𝔽): the division
        # algebra IS the endomorphism algebra, so they pin the mod-8 formula rather than relying on
        # it. d = 2 + dim_R(𝔽) for 𝔽 = R, C, H, O of real dimension 1, 2, 4, 8.
        @test division_algebra(Spin(1, 2)) === ℝ    # d=3,  Majorana
        @test division_algebra(Spin(1, 3)) === ℂ    # d=4,  Weyl
        @test division_algebra(Spin(1, 5)) === ℍ    # d=6,  symplectic Majorana-Weyl
        @test division_algebra(Spin(1, 9)) === ℝ    # d=10, Majorana-Weyl

        # the classical realisation must agree with the mod-8 answer wherever both exist
        for (p, q) in ((1, 2), (1, 3), (1, 5), (1, 9))
            @test division_algebra(Spin(p, q)) === division_algebra(Alakazam.classical(Spin(p, q)))
        end

        # Spin(p,q) = Spin(q,p), so the reality cannot depend on which way the signature is written
        for (p, q) in ((1, 2), (1, 3), (1, 5), (1, 9), (1, 10), (2, 5))
            @test division_algebra(Spin(p, q)) === division_algebra(Spin(q, p))
        end

        # periodicity: the type repeats with d in steps of 8 at fixed timelike count
        for q in 2:6
            @test division_algebra(Spin(1, q)) === division_algebra(Spin(1, q + 8))
        end
    end

    @testset "minimal spinor dimension" begin
        # Dirac is 2^(d/2), halved by the Weyl projection in even d
        @test Alakazam.dimension(Spin(1, 2)) == 2     # d=3,  Dirac
        @test Alakazam.dimension(Spin(1, 3)) == 2     # d=4,  Weyl
        @test Alakazam.dimension(Spin(1, 5)) == 4     # d=6,  Weyl
        @test Alakazam.dimension(Spin(1, 9)) == 16    # d=10, Weyl
        @test Alakazam.dimension(Spin(1, 10)) == 32   # d=11, Dirac
    end

    @testset "conjugation class follows from the algebra" begin
        clear_global!()
        # 2+1D Majorana: self-conjugate and real, the case the metric alone gets wrong. An
        # antisymmetric form makes the metric route say quaternionic; the algebra says real.
        maj = IndexSet("majorana",lower, fermionic, Spin(2, 1))
        maj.conjugate = maj
        @test Alakazam.conj_class(maj) === :fixed
        @test maj.range==2
        @test maj.metric_type==EpsilonTensor

        # SU(2) doublet: also self-conjugate with an antisymmetric form, but genuinely quaternionic
        su2 = IndexSet("su2", lower, SUAlgebra{ℂ}(2))
        su2.conjugate = su2
        @test Alakazam.conj_class(su2) === :quaternionic
        @test su2.range==2
        @test su2.metric_type==EpsilonTensor

        # a complex representation cannot be its own conjugate
        w = IndexSet("weyl", 2, lower, EpsilonTensor, fermionic, Spin(3, 1))
        w.conjugate = w
        @test_throws ErrorException Alakazam.conj_class(w)

        # declared as a partner pair, as a complex representation must be
        u = IndexSet("u", 2, lower, EpsilonTensor, fermionic, Spin(3, 1))
        v = IndexSet("v", 2, lower, EpsilonTensor, fermionic, Spin(3, 1))
        u.conjugate = v; v.conjugate = u
        @test Alakazam.conj_class(u) === :complex
    end
    @testset "the full Lorentzian table, d = 2..11" begin
        # The type of the minimal chiral spinor: Weyl where chirality exists, Dirac in odd
        # dimensions. Checked against the standard table rather than against the four cases the
        # exceptional isomorphisms pin, so the mod-8 formula is asserted across its range.
        #
        # In odd d this is the Dirac spinor, so it reads off the Majorana condition directly:
        # R at 3, 9, 11 for Majorana, H at 5, 7 for symplectic Majorana. In even d it is the Weyl
        # spinor, and the type predicts Majorana-Weyl exactly: R at 2 and 10 where MW exists,
        # H at 6 for the symplectic case, C at 4 and 8 where MW does not exist.
        expected = Dict(2 => ℝ, 3 => ℝ, 4 => ℂ, 5 => ℍ, 6 => ℍ,
                        7 => ℍ, 8 => ℂ, 9 => ℝ, 10 => ℝ, 11 => ℝ)
        for (d, want) in expected
            @test division_algebra(Spin(1, d - 1)) === want
        end

        # Majorana-Weyl exists exactly where the chiral spinor is real, which is what makes d = 8
        # complex rather than real: 8D has Majorana spinors and Weyl spinors but not both at once.
        mw_dims = (2, 10)
        for d in 2:11
            iseven(d) || continue
            @test (division_algebra(Spin(1, d - 1)) === ℝ) == (d in mw_dims)
        end

        # the signature convention cannot matter, since Spin(p,q) and Spin(q,p) are one group
        for d in 2:11
            @test division_algebra(Spin(1, d - 1)) === division_algebra(Spin(d - 1, 1))
        end
    end

end
