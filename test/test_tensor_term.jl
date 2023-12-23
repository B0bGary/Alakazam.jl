@testset "Tensor Term" begin
    @testset "Tensor Construction" begin
        spinor_indices = IndexSet("spinor", 2)
        spacetime_indices = IndexSet("spacetime", 4,lower)

        μ = Index("̇μ", spinor_indices, 1)
        ν = Index("̇ν", spinor_indices, 1)
        ρ = Index("̇ρ", spinor_indices, 1)
        τ = Index("̇τ", spinor_indices, 1)

        a = Index("a", spacetime_indices, 0)
        b = Index("b", spacetime_indices)

        # Shorthand
        ↑ = upper
        ↓ = lower

        # Underscore expr::Symbol cases
        @test Tensor("θ", [a => lower, b => lower]) == @Tensor θ_ab
        @test Tensor("T", [μ => lower, ν => lower]) == @Tensor T_μν
        @test Tensor("G") == @Tensor G

        # Caret expr::Expr cases
        @test Tensor("G", [μ => ↑, ν => ↑, μ => ↓, ν => ↓, ρ => ↑, τ => ↑]) == 
            @Tensor G^μν_μν^ρτ
        @test Tensor("T", [τ => ↓, ρ => ↑, μ => ↓, μ => ↑, τ => ↑]) ==
            @Tensor T_τ^ρ_μ^μ^τ
        
        # Unicode Tensor name which complicates string manipulation, see
        # https://docs.julialang.org/en/v1/manual/strings/
        @test Tensor("θ", [μ => ↑]) == @Tensor θ^μ
        @test Tensor("θ", [μ => ↑, ν => ↑]) == @Tensor θ^μν
        @test Tensor("θ", [μ => ↓, μ => ↑, τ => ↑]) == @Tensor θ_μ^μτ
        
        # This should probably throw an error of some kind
        # @test_throws IndexSummationError Tensor("θ", [a => lower, a => upper]) 
    end

    # @testset "TensorExpression Construction"
    #     #...
    # end
end