@testset "Lie Brackets" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4, lower)
        spinor_indices=IndexSet("spinor_indices", 2, lower,fermionic)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2, lower,fermionic)
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
        custom_sort_order!(["θ", "Θ", "ϵ"])
    end
    #basic proof of concept com rel
    X=Alakazam.Matrix("x", [μ=>upper])
    P=Alakazam.Matrix("p", [ν=>upper])
    Q=Alakazam.Matrix("q", [ρ=>upper])
    clear_commutators!()
    br = @Lie [X, P]
    br2= @Lie {X,P}
    
    @test swap(br) ≅ br
    @test swap(swap(br)) ≅ br
    @test swap(br2) ≅ br2
    @test swap(swap(br2)) ≅ br2

    # @show br
    br2 = @Lie [X, {P, Q}]
    # @show product_rule(br2)
    br3 = @Lie [X, P*Q]
    # @show product_rule(br3)
    # @show expand_commutators(br2) #graded comms

    remove_commutator!(X, P)
    remove_commutator!(br)
    jac=@Lie[X, [P, Q]] + @Lie[P, [Q, X]] + @Lie[Q, [X, P]]
    @test expand_commutators(jac)==Alakazam. zero()
    
    
    expr= @Lie[X, P]
    to_latex(expr)
    @syms ħ
    result=im*ħ*Metric([μ=>upper, ν=>upper])
    add_commutator!(1, X, P, result)
    @test expand_commutators(@Lie[P,X]) == -result
    @test Alakazam.simplify(expand_commutators(expr))== result
    # @show expand_commutators(@Lie[X, P*Q])

    add_commutator!(1,X,P,0)
    # @show expr
    # println(to_latex(expr))
    # println(to_latex(@Lie{X, P}))
end