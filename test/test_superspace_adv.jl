@testset "Superspace Advanced" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4, lower)
        spinor_indices=IndexSet("spinor_indices", lower, Alakazam.fermionic, SLAlgebra(2))
        spinor_conj_indices=IndexSet("spinor_conj_indices", Alakazam.fermionic, lower, SLAlgebra(2))
        α, β, γ, χ = Index(["α", "β", "γ", "χ"], spinor_indices)
        αd, βd, γd, χd = Index(["̇α", "̇β", "̇γ", "̇χ"], spinor_conj_indices)
        μ, ν, ρ, τ = Index(["μ", "ν", "ρ", "τ"], spacetime_indices)
        flavour_indices = IndexSet("flavour_indices", upper, SpAlgebra(4))
        i, j, k, l, m, n, p, q, r = Index(["i", "j", "k", "l", "m", "n", "p", "q", "r"], flavour_indices)
        x = Coordinate("x", spacetime_indices)
        θ = Coordinate("θ", spinor_indices)
        θbar = Coordinate("Θ", spinor_conj_indices)
        custom_sort_order!(["θ", "Θ", "ϵ"])

    end

    ∇=PartialDerivative("∇", [x], [μ => lower])
    ∂=PartialDerivative("∂", [θ], [α => lower])
    𝛛=PartialDerivative("𝛛", [θbar], [αd => lower])
    @action ∂(Tensor("θ", [β=>lower], θ)) = EpsilonTensor([β => lower, α=>lower])
    @action ∂(Tensor("θ", [β=>upper], θ)) = KroneckerDelta([α => lower, β=>upper])
    @action 𝛛(Tensor("Θ", [βd=>lower], θbar)) = EpsilonTensor([βd => lower, αd=>lower])
    @action 𝛛(Tensor("Θ", [βd=>upper], θbar)) = KroneckerDelta([αd => lower, βd=>upper])
    covar_bar_expansion=𝛛+im//2*Tensor("θ", [α=>upper], θ)*Tensor("σ", [μ=>upper, α=>lower, αd=>lower])*∇#+Am
    covar_expansion=∂+im//2*Tensor("σ", [μ=>upper, α=>lower, αd=>lower])*Tensor("Θ", [αd=>upper], θbar)*∇
    covar_bar=Derivative("̄D", [x, θbar], [αd => lower])
    covar=Derivative("D", [x, θ], [α => lower])
    @covariant covar=covar_expansion
    @covariant covar_bar=covar_bar_expansion

    #test generating products
    base_grassmann=Tensor("θ", [α=>lower, i=>upper], θ)
    @test (X=generate_superfield("X", [i=>upper], [Tensor("θ", [α=>lower, i=>upper], θ)], x); true)
    # tex=load("superfield")
    # @show tex

    # @show measure(tex,[covar,covar_bar])


    #easy tests
    flavour_indices = IndexSet("flavour_indices", upper, SLAlgebra(2))
    i, j, k, l, m, n, p, q, r = Index(["i", "j", "k", "l", "m", "n", "p", "q", "r"], flavour_indices)
    X=generate_superfield("X", [Tensor("θ", [α=>lower, i=>upper], θ)], x)
    @test length(X)==6

    X=generate_superfield("X", [i=>upper], [Tensor("θ", [α=>lower], θ)], x)
    covar_expansion=∂
    covar=Derivative("D", [x, θ], [α => lower])
    @covariant covar=covar_expansion
    @test integrate(X, [covar], false)==-2*Tensor("X", i=>upper, "(2)"=>lower, x)




    # flavour_indices = IndexSet("flavour_indices",upper, SpAlgebra(4))
    # i, j, k, l, m, n, p, q, r = Index(["i", "j", "k", "l", "m", "n", "p", "q", "r"], flavour_indices)
    # X=generate_superfield("X", [Tensor("θ", [α=>lower, i=>upper], θ)], x)
    # println(to_latex(X))
    # @show length(X)



    ∇=PartialDerivative("∇", [x], [μ => lower])
    ∂=PartialDerivative("∂", [θ], [α => lower])
    𝛛=PartialDerivative("𝛛", [θbar], [αd => lower])
    @action ∂(Tensor("θ", [β=>lower], θ)) = EpsilonTensor([β => lower, α=>lower])
    @action ∂(Tensor("θ", [β=>upper], θ)) = KroneckerDelta([α => lower, β=>upper])
    @action 𝛛(Tensor("Θ", [βd=>lower], θbar)) = EpsilonTensor([βd => lower, αd=>lower])
    @action 𝛛(Tensor("Θ", [βd=>upper], θbar)) = KroneckerDelta([αd => lower, βd=>upper])
    covar_bar_expansion=𝛛+im//2*Tensor("θ", [α=>upper], θ)*Tensor("σ", [μ=>upper, α=>lower, αd=>lower])*∇#+Am
    covar_expansion=∂+im//2*Tensor("σ", [μ=>upper, α=>lower, αd=>lower])*Tensor("Θ", [αd=>upper], θbar)*∇
    Dbar=Derivative("̄D", [x, θbar], [αd => lower])
    D=Derivative("D", [x, θ], [α => lower])
    @covariant D=covar_expansion
    @covariant Dbar=covar_bar_expansion

    grassmann=[Tensor("θ", [α=>lower], θ), Tensor("Θ", [αd=>lower], θbar)]
    sspace=Superspace(grassmann, [D, Dbar])
    X=generate_superfield("X", [i=>upper], sspace, x)
    xout=hard_simplify(integrate(X, sspace, false))
    # @show xout
    # println(to_latex(xout))

end
