
@testset "Naive Superspace" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4, lower)
        spinor_indices=IndexSet("spinor_indices", 2, lower, Alakazam.fermionic, EpsilonTensor)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2, Alakazam.fermionic, lower, EpsilonTensor)
        α, β, γ, χ = Index(["α", "β", "γ", "χ"], spinor_indices)
        αd, βd, γd, χd = Index(["̇α", "̇β", "̇γ", "̇χ"], spinor_conj_indices)
        μ, ν, ρ, τ = Index(["μ", "ν", "ρ", "τ"], spacetime_indices)
        flavour_indices = IndexSet("flavour_indices", 2, upper, EpsilonTensor)
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
    Dbar=Derivative("̄D", [x, θbar], [αd => lower])
    D=Derivative("D", [x, θ], [α => lower])
    @covariant D=covar_expansion
    @covariant Dbar=covar_bar_expansion
    #############################
    #test algebra of derivatives
    #############################
    X=Tensor("X", [x, θ, θbar])
    a1=Derivative("D", Derivative("̄D", X, [x, θbar], [γd=>lower]), [x, θ], [γ=>lower])
    b1=Derivative("̄D", Derivative("D",X, [x, θ], [γ=>lower]), [x, θbar], [γd=>lower])
    @test out2=derive(a1+b1)== im*Tensor("σ", [μ=>upper, γ=>lower, γd=>lower])∇(X)

    #test generating products
    base_grassmann=Tensor("θ", [α=>lower, i=>upper], θ)
    sf=generate_superfield_naive("X", [Tensor("θ", [α=>lower, i=>upper], θ)], x)
    # @show sf
    # show_function_of(true)
    # println(to_latex(sf))
    @test length(sf)==5

    grassmann=[Tensor("θ", [α=>lower], θ),Tensor("Θ", [αd=>upper], θ)]
    sspace=Superspace(grassmann,[D,Dbar])

    # @show sf2=generate_superfield_naive("Z",sspace,x)
    # println(to_latex(sf2))

    @test (generate_superfield_naive("X", [i=>upper], [Tensor("θ", [α=>lower, i=>upper], θ)], x);true)
    @test length(generate_superfield_naive("X", [i=>upper], [Tensor("θ", [α=>lower, i=>upper], θ),Tensor("Θ", [αd=>upper, i=>lower], θ)], x))==25

end
