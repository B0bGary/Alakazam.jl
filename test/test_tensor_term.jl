@testset "Tensor Term" begin
    @testset "Tensor Construction" begin
        spinor_indices = IndexSet("spinor", 2)
        spacetime_indices = IndexSet("spacetime", 4,lower)

        μ = Index("̇μ", spinor_indices, 1)
        ν = Index("̇ν", spinor_indices, 1)

        a = Index("a", spacetime_indices, 0)
        b = Index("b", spacetime_indices)

        @test Tensor("θ", [a => lower, b => lower]) == @Tensor θ_ab
        @test Tensor("T", [μ => lower, ν => lower]) == @Tensor T_μν
        
        # This should probably throw an error of some kind
        # @test_throws IndexSummationError Tensor("θ", [a => lower, a => upper]) 
    end

    # @testset "TensorExpression Construction"
    #     #...
    # end
end