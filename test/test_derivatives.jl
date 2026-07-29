@testset "Derivatives and Calculus" begin

    #setup
    if 1 == 1
        flavour_indices = IndexSet("flavour", 2, upper,fermionic)
        spacetime_indices = IndexSet("spacetime", 4, upper)
       @indices λ:τ ∈ spacetime_indices
       @indices i:z ∈ flavour_indices
        x = Coordinate("x", [spacetime_indices])
        θ = Coordinate("θ", [flavour_indices])
    end
    ###############################################
    ### begin test cases
    ###############################################

    ###############################################
    #streamlined deriv rego
    @test_nowarn D=PartialDerivative("∂", x,μ=>lower)
    X=Tensor("X",x,ν=>upper)
    res=KroneckerDelta(μ=>lower,ν=>upper)
    @action D(X)=res
    #basic check registration worked
    @test apply_derivative(D(X))==res
    @test res==apply_derivative(D(X))

    #shorthand Partial
    D=Partial("∂", x,ρ=>lower)
    X=Tensor("X",x,σ=>upper)
    res=Delta(ρ=>lower,σ=>upper)
    #test that deriv works symbolically on the basis
    @test apply_derivative(D(X))==res


    ###############################################
    #basic product rule on a multivar polynomial
    eq=(SymmetricTensor("A", [μ => lower, ν => lower])*
    Tensor("X", [μ => upper], x)*Tensor("X", [ν => upper], x)
    +Tensor("B", [μ => lower])*Tensor("X", [μ => upper], x)
    +Tensor("C"))

    derivedeq=(2*SymmetricTensor("A", [μ => lower, ρ => lower])*
    Tensor("X", [μ => upper], x)
    +Tensor("B", [ρ=> lower]))
    @test Alakazam.simplify(apply_derivative(D(eq)))==derivedeq

    @test (simplify(D(X));true)


end