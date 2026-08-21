@testset "Lie Brackets and Matrix Reduction" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4, lower, Metric)
        spinor_indices=IndexSet("spinor_indices", 2, lower, Alakazam.fermionic)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2, lower, Alakazam.fermionic)
        α, β, γ, χ = Index(["α", "β", "γ", "χ"], spinor_indices)
        αd, βd, γd, χd = Index(["̇α", "̇β", "̇γ", "̇χ"], spinor_conj_indices)
        μ, ν, ρ, τ, σ, κ, λ, ζ, ϕ = Index("μ", "ν", "ρ", "τ", "σ", "κ", "λ", "ζ", "ϕ", spacetime_indices)

        x = Coordinate("x", [spacetime_indices])
        θ = Coordinate("θ", [spinor_indices])
        θbar = Coordinate("Θ", [spinor_conj_indices])
        custom_sort_order!(["θ", "Θ", "ϵ"])
    end
    #basic proof of concept com rel
    spacetime_indices.range=4
    spacetime_indices.signature=(1, 3)
    g1=Gamma(μ=>upper)
    g2=Gamma(ν=>upper)
    g5=Gamma5()
    # @show μ.index_set.metric_type
    @anticommutator {g1, g2}=(2*Metric([μ=>upper, ν=>upper]))
    @anticommutator {g5, g1}=Alakazam.zero()
    custom_sort_order!(["γ5", "γ"])
    exp=g2*g5
    @test sort(exp)==-g5*g2
    @test sort(g2*g1)==-g1*g2+2*Metric([μ=>upper, ν=>upper])


    X=Alakazam.Matrix("x", [μ=>upper])
    P=Alakazam.Matrix("p", [ν=>upper])
    @syms ħ
    result=im*ħ*Metric([μ=>upper, ν=>upper])
    @commutator [X, P]=result
    @test get_commutator(X, P)[2]==result

    X2=Alakazam.Matrix("x", [ρ=>upper])
    P2=Alakazam.Matrix("p", [σ=>upper])
    custom_sort_order!(["x", "p"])
    @test sort(P2*X2)==X2*P2-im*ħ*Metric([ρ=>upper, σ=>upper])




    su3_indices = IndexSet("su3_flavour", (0, 8), lower)
    @indices a:f∈su3_indices
    @matrix GMMatrix
    # mutable struct GMMatrix <: MatrixSuperType
    #     data::TensorData
    #     GMMatrix(args...) = new(TensorData(args...))
    # end
    λa=GMMatrix("λ", a=>lower)
    λb=GMMatrix("λ", b=>lower)
    rhs=(2//3*Metric(a=>lower, b=>lower)
         +(T"D_(abc)"+im*T"F_[abc]")*GMMatrix("λ", c=>upper))
    # @show rhs

    @mrule λa*λb = rhs

    λc=GMMatrix("λ", c=>lower)
    λd=GMMatrix("λ", d=>lower)
    λe=GMMatrix("λ", e=>lower)
    λf=GMMatrix("λ", f=>upper)

    @test contract_adjacent(λc*λd*λe, GMMatrix)==(2//3*Metric(c=>lower, d=>lower)*λe
                                                  + 2//3*T"D_(cde)"
                                                  + T"D_(cd^b)"*T"D_(bef)"*λf
                                                  + im*T"D_(cd^b)"*T"F_[bef]"*λf
                                                  + 2//3*im*T"F_[cde]" + im*T"F_[cd^b]"*T"D_(bef)"*λf
                                                  -
                                                  T"F_[cd^b]"*T"F_[bef]"*λf)



    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4, lower)
        spinor_indices=IndexSet("spinor_indices", 2, lower, fermionic)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2, lower, fermionic)
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
    @syms ħ
    result=im*ħ*Metric([μ=>upper, ν=>upper])
    add_commutator!(1, X, P, result)

    X2=Alakazam.Matrix("x", [ρ=>upper])
    P2=Alakazam.Matrix("p", [σ=>upper])
    result2=im*ħ*Metric([σ=>upper, ρ=>upper])
    @test get_commutator(X2, P2)==(1, result2) #com rel
    @test get_commutator(P2, X2)==(1, -result2)#com rel

    #canonical ordering
    custom_sort_order!(["x", "p", "γ5", "γ"])
    t=Tensor("a", [ν=>upper])#some other tens
    #test on a tensor expression
    expresp=3*t*(X2*P2-result2)
    @test sort(3*P2*X2*t)==expresp

    #now, lets work in 3D euclidean space for pauli matrices.
    spacetime_indices.signature=(0, 3)
    spacetime_indices.range=3
    #pauli matrices.
    Pσ=Pauli(μ=>upper)
    Pσ2=Pauli(ν=>upper)
    @commutator [Pσ, Pσ2]=(2*im*EpsilonTensor([μ=>upper, ν=>upper, ρ=>upper])*Pauli(ρ=>lower))
    # @show Pσ
    @anticommutator {Pσ, Pσ2}=(2*im*Metric([μ=>upper, ν=>upper]))
    @commutator {Pσ, Pσ2}=2*im*Metric([μ=>upper, ν=>upper]) #these do the same thing
    t=Tensor("a")
    Pσ3=Pauli(ρ=>upper)
    @test Alakazam.simplify(contract_sigma(Pσ*Pσ2*Pσ3))==Metric([μ=>upper, ν=>upper])*Pσ3+im*EpsilonTensor([μ=>upper, ν=>upper, ρ=>upper])+Metric([ν=>upper, ρ=>upper])*Pσ-Metric([μ=>upper, ρ=>upper])*Pσ2
    #basic to lorentzian for gamma matrices
    spacetime_indices.range=4
    spacetime_indices.signature=(1, 3)
    g1=Gamma(μ=>upper)
    g2=Gamma(ν=>upper)

    g5=Gamma5()
    @anticommutator {g1, g2}=(2*Metric([μ=>upper, ν=>upper]))
    @anticommutator {g5, g1}=Alakazam.zero()
    exp=g2*g5
    @test sort(exp)==-g5*g2 #anticomrels


    #test an artificial example of a commutator in the way
    spacetime_indices.signature=(0, 3)
    spacetime_indices.range=3
    t=Alakazam.Matrix("T", σ=>upper)
    custom_sort_order!(["σ", "T"])
    @commutator [Pσ, t]=2*im*Tensor("A", [μ=>upper, σ=>upper]) #these do the same thing
    @test Alakazam.simplify(contract_sigma(Pσ*t*Pσ2)) == t*Metric([μ=>upper, ν=>upper])+im*Pauli(ρ=>lower)*t*EpsilonTensor([μ=>upper, ν=>upper, ρ=>upper])-2*im*Tensor("A", [ν=>upper, σ=>upper])*Pσ
    # @show Alakazam.simplify(contract_sigma(Pσ*t*Pσ2*Pσ3))

    # @show Alakazam.simplify(Alakazam.simplify(contract_adjacent(Pσ*Pσ2*Pσ3,PauliSuperType)))
    #looks good!

    T1=Alakazam.Matrix("T", σ=>upper)
    T2=Alakazam.Matrix("T", μ=>upper)
    @mrule T1*T2 = 2*Metric(σ=>upper, μ=>upper)
    T3=Alakazam.Matrix("T", ν=>lower)
    exp=T1*T3
    ## check generic matrix product reduction rule, where we have to use a metric contraction and relabeling
    @test contract_adjacent(exp, Alakazam.Matrix)==2*KroneckerDelta(ν=>lower, σ=>upper)


    spacetime_up = IndexSet("spacetime_up", (1, 3), upper, Metric)
    μ, ν, ρ, σ = Index("μ", "ν", "ρ", "σ", spacetime_up)
    g1=Gamma(μ=>upper)
    g2=Gamma(ν=>upper)
    g3=Gamma(μ=>lower)
    # @show g1*g2*g3
    @test sort(g1*g2*g3)==g1*g2*g3
    # @test sort(g1*g3*g2)==-g1*g2*g3+2*g1

    #:dim, not spacetime_up.range: naming the set captures its dimension at definition, and
    #these rules stay registered for the rest of the session, so gammas from any other set
    #would then be reduced with this set's dimension
    @mrule g1*g2*g3 = (2 - :dim)*g2
    @mrule g1*g3 = :dim


    g4=Gamma(ρ=>upper)
    g5=Gamma(σ=>upper)


    @mrule g1*g2*g4*g3 = 2*g4*g2+(:dim - 2)*g2*g4 # =4Metric in 4D
    @test contract_adjacent(g1*g2*g3, Alakazam.Gamma)==-2*g2
    @test contract_adjacent(g1*g3, Alakazam.Gamma)==4
    @test sort(contract_adjacent(g1*g2*g4*g3, Alakazam.Gamma))==4*Metric([ν=>upper, ρ=>upper])#great!!!
    @test sort(contract_adjacent(g1*g5*g4*g2*g3, Alakazam.Gamma))==-2*g2*g4*g5

    @mrule Gamma5()*Gamma5() = 1

    gfive=Gamma5()
    @test contract_adjacent(gfive*gfive, Alakazam.Gamma5)==1
    @test contract_adjacent(gfive*g1*gfive, Gamma5)==-g1


end