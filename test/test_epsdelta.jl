@testset "Epsilon and Delta Contraction Identities" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 2, SymMetricSuperType,(1,1),lower)
        spinor_indices=IndexSet("spinor_indices", 2,lower,fermionic)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2,lower,fermionic)
        α = Index("α", spinor_indices)
        β = Index("β", spinor_indices)
        γ = Index("γ", spinor_indices)
        χ = Index("χ", spinor_indices)
        αd = Index("̇α", spinor_conj_indices)
        βd = Index("̇β", spinor_conj_indices)
        γd = Index("̇γ", spinor_conj_indices)
        χd = Index("̇χ", spinor_indices)
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
        

    #basic 2 index contractions in 2D
        e1=EpsilonTensor([μ=>upper,ν=>upper])
        e2=EpsilonTensor([μ=>lower,ρ=>lower])
        delta=KroneckerDelta([ν=>upper,ρ=>lower])
        @test eliminate_epsilon(e1*e2)≅-delta

        e1=EpsilonTensor([μ=>upper,ν=>upper])
        e2=EpsilonTensor([μ=>lower,ν=>lower])
        # for summand in summands(eliminate_epsilon(e1*e2))
        #     @show summand
        #     @show length(summand[2])
        # end
        # @show length(summands(TensorExpression(-2))[1][2])
        @test eliminate_epsilon(e1*e2)==-2

        e1=EpsilonTensor([ν=>upper,μ=>upper])
        e2=EpsilonTensor([μ=>lower,ν=>lower])
        delta=KroneckerDelta([ρ=>lower,ν=>upper])
        @test eliminate_epsilon(e1*e2)==2

        e1=EpsilonTensor([ν=>upper,μ=>upper])
        e2=EpsilonTensor([μ=>lower,ρ=>lower])
        delta=KroneckerDelta([ν=>upper,ρ=>lower])
        @test eliminate_epsilon(e1*e2)≅delta
      
        e1=EpsilonTensor([μ=>upper,ν=>upper])
        e2=EpsilonTensor([ρ=>lower,μ=>lower])
        delta=KroneckerDelta([ν=>upper,ρ=>lower])
        @test eliminate_epsilon(e1*e2)≅delta

         e1=EpsilonTensor([ν=>upper,μ=>upper])
        e2=EpsilonTensor([ρ=>lower,μ=>lower])
        delta=KroneckerDelta([ν=>upper,ρ=>lower])
        @test eliminate_epsilon(e1*e2)==-delta

        #basic 3 index
        spacetime_indices.range=3 #need bigger range.
         e2=EpsilonTensor([ν=>upper,μ=>upper,ρ=>upper])
        e1=EpsilonTensor([ν=>lower,μ=>lower,ρ=>lower])
    @test eliminate_epsilon(e1*e2)==-6

     spacetime_indices.range=4#need bigger range.
         e2=EpsilonTensor([ν=>upper,μ=>upper,ρ=>upper,τ=>upper])
        e1=EpsilonTensor([ν=>lower,μ=>lower,ρ=>lower,σ=>lower])
        d=-6*KroneckerDelta([τ=>upper,σ=>lower])
    @test Alakazam.simplify(eliminate_epsilon(e1*e2))==d

end