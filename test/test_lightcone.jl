@testset "Lightcone" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4, lower)
        spinor_indices=IndexSet("spinor_indices", 2, lower, Alakazam.fermionic, EpsilonTensor)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2, Alakazam.fermionic, lower, EpsilonTensor)
        α, β, γ, χ = Index(["α", "β", "γ", "χ"], spinor_indices)
        αd, βd, γd, χd = Index(["̇α", "̇β", "̇γ", "̇χ"], spinor_conj_indices)
        μ, ν, ρ, τ = Index(["μ", "ν", "ρ", "τ"], spacetime_indices)
        flavour_indices = IndexSet("flavour_indices", 2, upper, Alakazam.fermionic, EpsilonTensor)
        i, j, k, l, m, n, p, q, r = Index(["i", "j", "k", "l", "m", "n", "p", "q", "r"], flavour_indices)
        x = Coordinate("x", spacetime_indices)
        θ = Coordinate("θ", spinor_indices)
        θbar = Coordinate("Θ", spinor_conj_indices)
        custom_sort_order!(["θ", "Θ", "ϵ"])

    end
    lp=ℓ⁺
    lm=ℓ⁻
    @test (T=Tensor("T", [ℓ⁺=>upper, ℓ⁺=>upper]);true)
    @test (S=Tensor("S", [ℓ⁻=>upper, ℓ⁻=>upper]);true)
    @test (S=Tensor("S", [ℓ⁺⁻ => upper, ℓ⁺⁻ => upper])*Tensor("S", [ℓ⁺⁻ => upper, ℓ⁻⁺ => upper]);true)
    @test expand_lc(S)==Tensor("S", [lp =>upper, lp=>upper])*Tensor("S", [lp =>upper, lm=>upper])+Tensor("S", [lm =>upper, lm=>upper])*Tensor("S", [lm =>upper, lp=>upper])
    #check overridden appension behaviour
    @test T×T==T*T

    sf=generate_superfield_naive("X", [i=>upper], [Tensor("θ", [ℓ⁺=>lower, i=>upper], θ), Tensor("θ", [ℓ⁻ => lower, i=>upper], θ)], x)
    @test length(sf)==9

    S=generate_superfield_naive("s", [ℓ⁺=>lower, ℓ⁺=>lower], [Tensor("θ", [ℓ⁺=>upper], θ), Tensor("Θ", [ℓ⁺ => upper], θbar)], x)
    T=generate_superfield_naive("t", [ℓ⁻ => lower, ℓ⁻=>lower, ℓ⁻=>lower, ℓ⁻=>lower], [Tensor("θ", [ℓ⁺=>upper], θ), Tensor("Θ", [ℓ⁺ => upper], θbar)], x)
    W=generate_superfield_naive("w", [ℓ⁻=>lower], [Tensor("θ", [ℓ⁺=>upper], θ), Tensor("Θ", [ℓ⁺ => upper], θbar)], x)
    Wbar=generate_superfield_naive("̄w", [ℓ⁻=>lower], [Tensor("θ", [ℓ⁺=>upper], θ), Tensor("Θ", [ℓ⁺ => upper], θbar)], x)


    ∇=PD("∇", [x], [ℓ⁺=>lower,ℓ⁺=>lower])
    ∂=PD("∂", [θ], [ℓ⁺=>lower])
    𝛛=PD("𝛛", [θbar], [ℓ⁺=>lower])
    @action ∂(Tensor("θ", [ℓ⁺=>upper], θ)) = 1
    @action ∂(Tensor("θ", [ℓ⁺=>lower], θ)) = 1
    @action 𝛛(Tensor("Θ", [ℓ⁺=>upper], θbar)) = 1
    @action 𝛛(Tensor("Θ", [ℓ⁺=>lower], θbar)) = 1
    covar_bar=Derivative("̄D", [x, θbar], [ℓ⁺=>lower])
    covar=Derivative("D", [x, θ], [ℓ⁺=>lower])
    @covariant covar=∂-im//2*Tensor("Θ", [ℓ⁺=>upper], θbar)*∇
    @covariant covar_bar=-𝛛+im//2*Tensor("θ", [ℓ⁺=>upper], θ)*∇


    con1=derive(covar(Wbar))
    projcon=project_to_component(con1,[covar,covar_bar],1,1)
    inner=isTensor(operand(isTensor(first_summand(projcon))))
    ref=Tensor("̄w",[lm=>lower,lp=>lower],"(1,0)"=>lower,x) #make sure function_of matches
    @test inner==ref
    @test name(inner)==name(ref)
    @test indices(inner)==indices(ref)
    @test tableaux(inner)==tableaux(ref)

    @commutator {covar,covar_bar}=im*∇
    @commutator [covar,∇]=0
    @commutator [covar_bar,∇]=0
    T=Tensor("T",x)
    nest=covar(covar_bar(T))
    @test commute(nest)==-covar_bar(covar(T))+im*∇(T)

end
