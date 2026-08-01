@testset "Serialisation Export" begin
    clear_global!()
    #test
    flavour_indices = IndexSet("flavour_indices", 2, upper)
    a, b, c, e = Index(["a", "b", "c", "e"], flavour_indices)
    x=Coordinate("x", flavour_indices)
    spacetime_indices = IndexSet("spacetime", 4, lower)
    spinor_indices=IndexSet("spinor_indices", 2, lower, Alakazam.fermionic, EpsilonTensor)
    spinor_conj_indices=IndexSet("spinor_conj_indices", 2, Alakazam.fermionic, lower, EpsilonTensor)
    α, β, γ, χ = Index(["α", "β", "γ", "χ"], spinor_indices)
    αd, βd, γd, χd = Index(["̇α", "̇β", "̇γ", "̇χ"], spinor_conj_indices)
    μ, ν, ρ, τ = Index(["μ", "ν", "ρ", "τ"], spacetime_indices)
    θ = Coordinate("θ", spinor_indices)
    θbar = Coordinate("Θ", spinor_conj_indices)

    ∇=PartialDerivative("∇", [x], [μ => lower])
    ∂=PartialDerivative("∂", [θ], [α => lower])
    𝛛=PartialDerivative("𝛛", [θbar], [αd => lower])
    @action ∂(Tensor("θ", [β=>lower], θ)) = EpsilonTensor([β => lower, α=>lower])
    @action ∂(Tensor("θ", [β=>upper], θ)) = KroneckerDelta([α => lower, β=>upper])
    @action 𝛛(Tensor("Θ", [βd=>lower], θbar)) = EpsilonTensor([βd => lower, αd=>lower])
    @action 𝛛(Tensor("Θ", [βd=>upper], θbar)) = KroneckerDelta([αd => lower, βd=>upper])
    covar_bar_expansion=𝛛+im/2*Tensor("θ", [α=>upper], θ)*Tensor("σ", [μ=>upper, α=>lower, αd=>lower])*∇#+Am
    covar_expansion=∂+im/2*Tensor("σ", [μ=>upper, α=>lower, αd=>lower])*Tensor("Θ", [αd=>upper], θbar)*∇
    covar_bar=Derivative("̄D", [x, θbar], [αd => lower])
    covar=Derivative("D", [x, θ], [α => lower])

    @covariant covar=covar_expansion
    @covariant covar_bar=covar_bar_expansion
    custom_sort_order!(["θ", "Θ", "ϵ"])



    #####testing

    d=SymmetricTensor("T", [a=>upper, b=>upper])
    d2=SymmetricTensor("T", [a=>lower, b=>lower])

    @syms ħ
    expr=(2*ħ)*d*d2
    _TESTFILE1 = joinpath(tempdir(), "test_expression")
    save(expr, _TESTFILE1)

    # println(to_JSON(expr))
    #version with headers


    @commutator {covar, covar_bar}=im*Tensor("σ", [μ=>upper, α=>lower, αd=>lower])*∇
    @commutator [covar, ∇]=0
    @commutator [covar_bar, ∇]=0

    X=Alakazam.Matrix("x", [μ=>upper])
    P=Alakazam.Matrix("p", [ν=>upper])
    @syms ħ
    result=im*ħ*Metric([μ=>upper, ν=>upper])

    add_commutator!(1, X, P, result)

    _TESTFILE2 = joinpath(tempdir(), "test_expression_headers")

    save(expr, _TESTFILE2, true)

end