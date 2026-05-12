@testset "Superspace" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetimelc", 4,lower)
        # flavour_indices = IndexSet("flavour_indices", 2,upper)
        spinor_indices=IndexSet("spinor_indices", 2,lower)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2,lower)


        α = Index("α", spinor_indices, 1)
        β = Index("β", spinor_indices, 1)
        γ = Index("γ", spinor_indices, 1)
        χ = Index("χ", spinor_indices, 1)
        αd = Index("̇α", spinor_conj_indices, 1)
        βd = Index("̇β", spinor_conj_indices, 1)
        γd = Index("̇γ", spinor_conj_indices, 1)
        χd = Index("̇χ", spinor_indices, 1)
        

        μ = Index("μ", spacetime_indices)
        ν = Index("ν", spacetime_indices)
        ρ = Index("ρ", spacetime_indices)
        τ = Index("τ", spacetime_indices)
        σ = Index("σ", spacetime_indices)
        κ = Index("κ", spacetime_indices)

        x = Coordinate("x", [spacetime_indices])
        custom_sort_order(["θ","Θ","ϵ"])

        A=Tensor("A", [μ=>upper],x)
        B=Tensor("B", [ν=>upper,σ=>upper],x)
        B2=Tensor("B", [ρ=>upper,κ=>upper],x)
        B3=Tensor("B", [ρ=>upper,σ=>upper],x)
        B4=Tensor("B", [ν=>upper,κ=>upper],x)
        e1=B2*B*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        e2=B3*B4*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])

        @show e1
        @show e2

        # Alakazam.symmetrise_partials(e1)
        
        @show Alakazam.hard_simplify(e1+e2)

    end



    
end