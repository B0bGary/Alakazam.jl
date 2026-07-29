
@testset "Advanced Tensor Arithmetic" begin

    #setup
    if 1 == 1
        flavour_indices = IndexSet("flavour_indices", 2, upper,Alakazam.fermionic)
        @indices a:z ∈ flavour_indices
        x = Coordinate("x", [flavour_indices])
        θ = Coordinate("θ", [flavour_indices])
        θbar = Coordinate("Θ", [flavour_indices])
    end
    ###############################################
    ### begin test cases
    ###############################################



    ###############################################
    ### ADVANCED ARITHMETIC
    ###############################################

    #TENSORS APPENSION

    #Tensor appension basics
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    S1 = Tensor("S", [a => lower, b => lower, c => lower], x)
    S2 = Tensor("S", [d => lower, e => lower, f => lower], x)
    @test T1 × S1 == T1 * S2 #appension should automatically get the next indices in the index set

    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    S1 = Tensor("S", [a => lower, b => lower, b => upper], x)
    S2 = Tensor("S", [d => lower, e => lower, e => upper], x)
    @test T1 × S1 == T1 * S2 #appension should also work with dummy indices

    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    V1 = Tensor("V", [t => lower, u => lower, v => lower], x)
    S1 = Tensor("S", [a => lower, b => lower, b => upper], x)
    S2 = Tensor("S", [d => lower, e => lower, e => upper], x)
    @test (T1 * V1) × S1 == (T1 * V1) * S2 #appension of tensor to a tensorterm

    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    V1 = Tensor("V", [t => lower, u => lower, v => lower], x)
    S1 = Tensor("S", [a => lower, b => lower, b => upper], x)
    W1 = Tensor("W", [t => lower, a => upper, c => upper], x)

    S2 = Tensor("S", [d => lower, e => lower, e => upper], x)
    W2 = Tensor("W", [f => lower, d => upper, g => upper], x)
    @test (T1 * V1) × (S1 * W1) == (T1 * V1) * (S2 * W2) #appension of tensorterm to a tensorterm


    #hodge star
    spacetime_indices = IndexSet("spacetime", 4, lower)
    μ = Index("μ", spacetime_indices)
    ν = Index("ν", spacetime_indices)
    ρ = Index("ρ", spacetime_indices)
    σ = Index("σ", spacetime_indices)
    τ = Index("τ", spacetime_indices)
    λ = Index("λ", spacetime_indices)

    #Basic hodge dual of a 2-form in 4D
    F=AntisymmetricTensor("F",[μ=>lower,ν=>lower])
    # @show tableaux(EpsilonTensor([μ=>upper,ν=>upper,τ=>upper,λ=>upper]))
    Fd=1//2*EpsilonTensor([μ=>upper,ν=>upper,τ=>upper,λ=>upper])*AntisymmetricTensor("F",[τ=>lower,λ=>lower])
    # @show  tableaux(terms(summands(Fd)[1][2])[1])
    sf=star(F)
    @test star(F)==Fd
    @test ⋆(F)==Fd

    A=AntisymmetricTensor("A",[μ=>lower])
    Ad=EpsilonTensor([μ=>upper,ν=>upper,ρ=>upper,τ=>upper])*Tensor("A",[τ=>lower])
    # @test ⋆(A)==Ad #failed for some reason... need to fix equality logic
     @test ⋆(A)-Ad==Alakazam. zero()
    f=Tensor("f")
    fd=EpsilonTensor([μ=>upper,ν=>upper,ρ=>upper,σ=>upper])*Tensor("f")
    # @test ⋆(f,spacetime_indices)==fd #failed also...
    @test ⋆(f,spacetime_indices)-fd==Alakazam. zero()

    #######
    #now, suppose colour indices for example
    #######
    colour_indices = IndexSet("colour", 3, upper)
    a = Index("a", colour_indices)
    b = Index("b", colour_indices)
    c = Index("c", colour_indices)
    d = Index("d", colour_indices)
    #Basic hodge dual of a 2-form in 4D
    F=Tensor("F",[a=>upper,μ=>lower,ν=>lower],YoungTableau([1,1],[2,3]))
    Fd=1//2*EpsilonTensor([μ=>upper,ν=>upper,τ=>upper,λ=>upper])*Tensor("F",[a=>upper,τ=>lower,λ=>lower],YoungTableau([1,1],[2,3]))
    sf=⋆(F)
    @test sf==Fd
    # @show ⋆(F)
    A=Tensor("A",[a=>upper,b=>upper,μ=>lower])
    Ad=EpsilonTensor([μ=>upper,ν=>upper,ρ=>upper,τ=>upper])*Tensor("A",[a=>upper,b=>upper,τ=>lower])
    # @test ⋆(A)==Ad #failed for some reason... need to fix equality logic
     @test ⋆(A)-Ad==Alakazam. zero()


    #now, in 6D for example
    spacetime_indices = IndexSet("spacetime", 6, lower)
    μ = Index("μ", spacetime_indices)
    ν = Index("ν", spacetime_indices)
    ρ = Index("ρ", spacetime_indices)
    σ = Index("σ", spacetime_indices)
    τ = Index("τ", spacetime_indices)
    λ = Index("λ", spacetime_indices)
    α = Index("α", spacetime_indices)
    β = Index("β", spacetime_indices)
    γ = Index("γ", spacetime_indices)
    δ = Index("δ", spacetime_indices)
    χ = Index("χ", spacetime_indices)
    ϕ = Index("ϕ", spacetime_indices)
    κ = Index("κ", spacetime_indices)
    ο = Index("ο", spacetime_indices)
    ω = Index("ω", spacetime_indices)
    θ = Index("θ", spacetime_indices)
    υ = Index("υ", spacetime_indices)
    # 2-form in 6D
    F=Tensor("F",[a=>upper,μ=>lower,ν=>lower],YoungTableau([1,1],[2,3]))#only antisym in ST inds
    Fd=1//2*EpsilonTensor([μ=>upper,ν=>upper,ρ=>upper,σ=>upper,α=>upper,β=>upper])*Tensor("F",[a=>upper,α=>lower,β=>lower],YoungTableau([1,1],[2,3]))
    @test ⋆(F)-Fd==Alakazam. zero()
    #3-form H in 6D
    H=Tensor("H",[a=>upper,μ=>lower,ν=>lower,σ=>lower],YoungTableau([1,1,1],[2,3,4]))
    Hd=1//factorial(3)*EpsilonTensor([μ=>upper,ν=>upper,ρ=>upper,α=>upper,β=>upper,γ=>upper])*Tensor("H",[a=>upper,α=>lower,β=>lower,γ=>lower],YoungTableau([1,1,1],[2,3,4]))
    @test ⋆(H)==Hd

   
    V=Tensor("v",[μ=>lower])
    ω=AntisymmetricTensor("ω",[τ=>lower,ν=>lower,ρ=>lower])
    contr=Tensor("v",[μ=>lower])*AntisymmetricTensor("ω",[μ=>upper,ν=>lower,ρ=>lower])
    @test V⨼ω==contr
    @test ι(V,ω)==contr
    @test contraction(V,ω)==contr

    F=AntisymmetricTensor("F",[μ=>lower,ν=>lower])
    H=AntisymmetricTensor("H",[λ=>lower,ρ=>lower,σ=>lower])
    #antisym of wedge
    @test Alakazam.hard_simplify(F∧H-H∧F)==Alakazam. zero()#want this to simplify....

    @test (1×1)==1
    @test F∧H==wedge(F,H)
end