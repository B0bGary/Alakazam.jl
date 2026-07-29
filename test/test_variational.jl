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