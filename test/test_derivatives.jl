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

    ###############################################
    # apply_derivative on a single PRODUCT, not a sum. `TensorTerm` matched no
    # method, so an identity written as one product raised a MethodError while
    # the same identity written as a difference went through.
    p = Tensor("p", x)
    q = Tensor("q", x)
    Dm = Partial("∂", x, μ => lower)
    @test apply_derivative(p * Dm(q)) == p * Dm(q)
    Dn = Partial("∂", x, ν => lower)
    @test apply_derivative(Dm(p) * Dn(q)) == Dm(p) * Dn(q)
    # a metric or a delta alongside the derivative is the case that surfaced it
    @test (apply_derivative(Metric([μ => upper, ν => upper]) * Dm(q)); true)
    # a term with no derivative in it at all is returned unchanged
    @test apply_derivative(p * q) == p * q
    @test apply_derivative(p) == p
    # and the sum form keeps working, so the new methods shadow nothing
    @test apply_derivative(p * Dm(q) - p * Dm(q)) == zero()

end