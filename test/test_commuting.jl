@testset "Manipulations and Simplifications" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4,lower)
        spinor_indices=IndexSet("spinor_indices", 2,lower,fermionic)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2,lower,fermionic)
        @indices α:γ ∈ spinor_indices
        χ = Index("χ", spinor_indices)
        αd = Index("̇α", spinor_conj_indices)
        βd = Index("̇β", spinor_conj_indices)
        γd = Index("̇γ", spinor_conj_indices)
        χd = Index("̇χ", spinor_conj_indices)
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
        
    #needs to commute partials, and relabel dummies inside partials, and commute Bs
        A=Tensor("A", [μ=>upper],x)
        B=Tensor("B", [ν=>upper,σ=>upper],x)
        B2=Tensor("B", [ρ=>upper,κ=>upper],x)
        B3=Tensor("B", [ρ=>upper,σ=>upper],x)
        B4=Tensor("B", [ν=>upper,κ=>upper],x)
        e1=B2*B*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        e2=B3*B4*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        es= Alakazam.simplify(e1+e2)
        ex=2*Tensor("B", [ρ=>upper,σ=>upper],x)*Tensor("B", [ν=>upper,κ=>upper],x)*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        # @test Alakazam.hard_simplify(es-ex)==Alakazam. zero()
        @test es≅ex

        #same, but now Bs are SUSY, so the swap should make everything  zero().
        A=Tensor("A", [μ=>upper],x)
        B=Tensor("B", [ν=>upper,α=>upper],x)
        B2=Tensor("B", [ρ=>upper,β=>upper],x)
        B3=Tensor("B", [ρ=>upper,α=>upper],x)
        B4=Tensor("B", [ν=>upper,β=>upper],x)
        e1=B2*B*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        e2=B3*B4*PartialDerivative("∇",PartialDerivative("∇",A, [x], [ν => lower]), [x], [ρ => lower])
        es= Alakazam.simplify(e1+e2)
        @test es≅Alakazam. zero()


        #commute an ordinary PD with SUSY PD
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        e1=PartialDerivative("∂",PartialDerivative("∇",A, [x], [ν => lower]), [θ], [α => lower])
        e2=PartialDerivative("∇",PartialDerivative("∂",A, [θ], [α => lower]), [x], [ν => lower])
        @test e1≅e2

        #commute a SUSY with SUSY PD
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        e1=PartialDerivative("∂",PartialDerivative("∂",A, [θ], [β => lower]), [θ], [α => lower])
        e2=-PartialDerivative("∂",PartialDerivative("∂",A, [θ], [α => lower]), [θ], [β => lower])
        @test e1≅e2


        #commute a SUSY with different SUSY PD
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        e1=PartialDerivative("𝛛",PartialDerivative("∂",A, [θ], [β => lower]), [θbar], [αd => lower])
        e2=-PartialDerivative("∂",PartialDerivative("𝛛",A, [θbar], [αd => lower]), [θ], [β => lower])
        @test e1≅e2

        #commute a tensor past susy diff op
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [ϕ=>upper],[x,θ,θbar])
        e1=B*PartialDerivative("𝛛",PartialDerivative("∂",A, [θ], [β => lower]), [θbar], [αd => lower])
        e2=-PartialDerivative("∂",PartialDerivative("𝛛",A, [θbar], [αd => lower]), [θ], [β => lower])*B
        @test e1≅e2

        #commute a susy tensor past susy diff op
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [γ=>upper],[x,θ,θbar])
        e1=B*PartialDerivative("𝛛",PartialDerivative("∂",A, [θ], [β => lower]), [θbar], [αd => lower])
        e2=-PartialDerivative("∂",PartialDerivative("𝛛",A, [θbar], [αd => lower]), [θ], [β => lower])*B
        @test e1≅e2

        #commute a susy tensor past differnt susy diff op
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [γ=>upper],[x,θ,θbar])
        e1=B*PartialDerivative("∂",PartialDerivative("∇",A, [x], [ν => lower]), [θ], [α => lower])
        e2=-PartialDerivative("∇",PartialDerivative("∂",A, [θ], [α => lower]), [x], [ν => lower])*B
        @test e1≅e2


        
        #products of diff ops
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [γ=>upper],[x,θ,θbar])
        e1=PartialDerivative("∇",A, [x], [ν => lower])*PartialDerivative("∂",B, [θ], [α => lower])
        e2=PartialDerivative("∂",B, [θ], [α => lower])*PartialDerivative("∇",A, [x], [ν => lower])
        @test e1≅e2

        #products of nested diff ops, with different levels
        A=Tensor("A", [μ=>upper],[x,θ,θbar])
        B=Tensor("B", [γ=>upper],[x,θ,θbar])
        e1=PartialDerivative("∇",A, [x], [ν => lower])*PartialDerivative("∇",PartialDerivative("∂",B, [θ], [α => lower]), [x], [μ => lower])
        e2=PartialDerivative("∇",PartialDerivative("∂",B, [θ], [α => lower]), [x], [μ => lower])*PartialDerivative("∇",A, [x], [ν => lower])
        @test e1≅e2
    
end