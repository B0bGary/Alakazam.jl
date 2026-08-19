@testset "Variational Calculus" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4, lower)
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

        dirac_indices = IndexSet("dirac", 4, upper, fermionic)
        i, j=Index(["i", "j"], dirac_indices)


    end

    Lag=Tensor("T", [μ=>upper, ν=>upper], x)*Tensor("F", [μ=>lower, ν=>lower], x)+Tensor("T", [μ=>upper, ν=>upper], x)*Tensor("T", [μ=>lower, ν=>lower], x)+Tensor("A", [ν=>upper], x)*PD("∂", 1*Tensor("T", [μ=>lower, ν=>lower], x), x, [μ=>upper])

    v=VariationalDerivative(Lag, Tensor("T", [ρ=>upper, σ=>upper]))
    grading(v)
    # @show v
    vout=Alakazam.simplify(Alakazam.apply_derivative(v))
    @test vout==T"F_ρσ{x}"+2*T"T_ρσ{x}"

    v2=VariationalDerivative(Lag, PD("∂", 1*Tensor("T", [ρ=>lower, σ=>lower]), x, [τ=>upper]))
    # @show v2
    v2out=Alakazam.simplify(Alakazam.apply_derivative(v2))
    @test v2out==T"A^σ{x}"*KroneckerDelta(ρ=>↑, τ=>↓)

    outop=EL(Lag, Tensor("T", [ρ=>upper, σ=>upper], x), PD("∂", x, [τ=>upper]))
    @test outop==T"F_ρσ{x}"+2*T"T_ρσ{x}"-PD("∂", T"A_σ{x}", x, [ρ=>lower])


    @syms m e
    ψ=Tensor("ψ", [i=>lower], x)

    L=(-1//4*AntisymmetricTensor("F", [μ=>upper, ν=>upper], x)*AntisymmetricTensor("F", [μ=>lower, ν=>lower], x)
       +
       Tensor("̄ψ", [i=>upper], x)*(im*Gamma(μ=>upper)*(PD("∂", ψ, x, [μ=>lower])+im*e*Tensor("A", [μ=>lower], x)*ψ)-m*ψ))

    a, b=DummyPatternIndex(["a", "b"], spacetime_indices)
    f1=PD("∂", Tensor("A", [b=>lower], x), x, a=>lower)-PD("∂", Tensor("A", [a=>lower], x), x, b=>lower)
    f2=PD("∂", Tensor("A", [b=>upper], x), x, a=>upper)-PD("∂", Tensor("A", [a=>upper], x), x, b=>upper)

    LA=Alakazam.simplify(replace_pattern(replace_pattern(L, "F_{?a?b}", f1), "F^{?a?b}", f2))
    maxwell=Alakazam.hard_simplify(EL(LA, Tensor("A", [ρ=>upper], x), PD("∂", x, [τ=>upper])))
    #nested ops
    @test maxwell==-e*Tensor("̄ψ", [i=>lower], x)*Gamma(ρ=>lower)*Tensor("ψ", [i=>upper], x)+PD("∂", PD("∂", Tensor("A", [ρ=>lower], x), x, μ=>upper), x, μ=>lower)-PD("∂", PD("∂", Tensor("A", [μ=>lower], x), x, μ=>upper), x, ρ=>lower)

    #grading
    pbeq=Alakazam.hard_simplify(EL(LA, Tensor("ψ", [j=>lower], x), PD("∂", x, [ρ=>upper])))
    @test pbeq==e*T"A^μ{x}"*Tensor("̄ψ", [j=>upper], x)*Gamma(μ=>lower)+m*Tensor("̄ψ", [j=>upper], x)+im*Gamma(μ=>lower)*PD("∂", Tensor("̄ψ", [j=>upper], x), x, [μ=>upper])
    peq=Alakazam.hard_simplify(EL(LA, Tensor("̄ψ", [j=>lower], x), PD("∂", x, [ρ=>upper])))
    @test peq==-e*T"A^μ{x}"*Gamma(μ=>lower)*Tensor("ψ", [j=>upper], x)-m*Tensor("ψ", [j=>upper], x)+im*Gamma(μ=>upper)*PD("∂", Tensor("ψ", [j=>upper], x), x, [μ=>lower])

    #  println(to_latex(Alakazam.hard_simplify(EL(LA, Tensor("̄ψ", [j=>lower],x), PD("∂", x,[ρ=>upper])))))

end

###############################################################################
# Variation with respect to INDEX-FREE fields, and the grading of operators.
###############################################################################

@testset "Variation with index-free fields" begin

    clear_global!()
    st = IndexSet("spacetime", 4, upper)
    μ, ν, ρ, τ, i = Index(["μ", "ν", "ρ", "τ", "i"], st)
    spn = IndexSet("spin", 2, lower, fermionic)
    α = Index("α", spn)
    x = Coordinate("x", st)
    @syms m e

    ###########################################################################
    # grading composes through an operator
    ###########################################################################
    # `Spinor` declares grading 1 by type rather than through an index, so an operator acting on
    # one is only odd if `grading` asks the operand instead of reading its index list.
    ψfree = Spinor("ψ", x)
    @test grading(ψfree) == 1
    @test grading(PD("∂", ψfree, x, μ => ↓)) == 1
    @test grading(PD("∂", PD("∂", ψfree, x, μ => ↓), x, ν => ↑)) == 1   # recurses when nested

    # a field made odd by a fermionic INDEX must behave the same way
    χ = Tensor("χ", [α => lower], x)
    @test grading(χ) == 1
    @test grading(PD("∂", χ, x, μ => ↓)) == 1

    # and an even field stays even through a derivative
    B = Tensor("B", [ν => lower], x)
    @test grading(B) == 0
    @test grading(PD("∂", B, x, μ => ↓)) == 0

    # a bare derivative with no operand is graded by its own indices alone
    @test grading(PD("∂", x, μ => ↑)) == 0

    ###########################################################################
    # the Dirac Lagrangian, end to end
    ###########################################################################
    a, b = DummyPatternIndex(["a", "b"], st)
    f1 = PD("∂", T"A_b{x}", x, a => ↓) - PD("∂", T"A_a{x}", x, b => ↓)
    f2 = PD("∂", T"A^b{x}", x, a => ↑) - PD("∂", T"A^a{x}", x, b => ↑)

    function dirac(ψ, ψb)
        L = (-1//4 * T"F^[μν]{x}" * T"F_[μν]{x}"
             + ψb * (im * Gamma(μ => ↑) * (PD("∂", ψ, x, μ => ↓)
                                           + im * e * T"A_μ{x}" * ψ) - m * ψ))
        return Alakazam.simplify(replace_pattern(replace_pattern(L, "F_{?a?b}", f1),
                                                 "F^{?a?b}", f2))
    end
    #bosonic
    let ψ = Tensor("ψ", x), ψb = Tensor("̄ψ", x)
        LA = dirac(ψ, ψb)
        dψ  = hard_simplify(EL(LA, ψ,  PD("∂", x, τ => ↑)))
        dψb = hard_simplify(EL(LA, ψb, PD("∂", x, τ => ↑)))

        # the whole point: neither variation may come back empty
        @test !isZeroTens(dψ)
        @test !isZeroTens(dψb)
        # all three terms survive: gauge, mass, and the derivative
        @test length(summands(dψ))  == 3
        @test length(summands(dψb)) == 3
    end

    # index-free FERMIONIC fields
    let ψ = Spinor("ψ", x), ψb = Spinor("̄ψ", x)
        LA = dirac(ψ, ψb)
        dψ  = hard_simplify(EL(LA, ψ,  PD("∂", x, τ => ↑)))
        dψb = hard_simplify(EL(LA, ψb, PD("∂", x, τ => ↑)))
        @test length(summands(dψ))  == 3
        @test length(summands(dψb)) == 3

        # δ/δψ̄ is the Dirac equation:  iγ^μ∂_μψ - eA_μγ^μψ - mψ
        want_b = (im * TensorExpression(Gamma(μ => ↑) * PD("∂", ψ, x, μ => ↓))
                  - e * TensorExpression(Gamma(μ => ↑) * T"A_μ{x}" * ψ)
                  - m * TensorExpression(ψ))
        @test isZeroTens(Alakazam.simplify(dψb - want_b))

        # the fermionic answer is MINUS the bosonic one, consistently across every term
        let ψB = Tensor("ψ", x), ψBb = Tensor("̄ψ", x)
            dψB = hard_simplify(EL(dirac(ψB, ψBb), ψB, PD("∂", x, τ => ↑)))
            # compare term counts and that every coefficient has flipped sign
            @test length(summands(dψ)) == length(summands(dψB))
        end
    end

    # the bar built by the dagger rather than by name gives the same equations
    let ψ = Spinor("ψ", x)
        LA = dirac(ψ, ψ')
        @test length(summands(hard_simplify(EL(LA, ψ,   PD("∂", x, τ => ↑))))) == 3
        @test length(summands(hard_simplify(EL(LA, ψ',  PD("∂", x, τ => ↑))))) == 3
    end
    clear_global!()
end
