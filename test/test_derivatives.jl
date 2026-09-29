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
    #nested partials compare up to any allowed reordering of the chain, not one swap at a time,
    #whatever the chain acts on
    ∂dv(i, T) = PartialDerivative("∂", T, [x], [i])
    ∂θv(i, T) = PartialDerivative("∂", T, [θ], [i])
    Xs, Ys, Zs = Tensor("X", x), Tensor("Y", x), Tensor("Z", x)
    #deep enough that single swaps no longer connect every ordering
    @test isZeroTens(simplify(∂dv(λ => lower, ∂dv(μ => lower, ∂dv(ν => lower, ∂dv(ρ => lower, Xs)))) -
                              ∂dv(ρ => lower, ∂dv(ν => lower, ∂dv(μ => lower, ∂dv(λ => lower, Xs))))))
    box2 = ∂dv(μ => lower, ∂dv(ν => lower, ∂dv(μ => upper, ∂dv(ν => upper, Xs))))
    @test isZeroTens(simplify(box2 - ∂dv(ν => upper, ∂dv(μ => upper, ∂dv(ν => lower, ∂dv(μ => lower, Xs))))))
    #a single term stays a single term, rather than an average of its orderings
    let s = summands(simplify(box2))
        @test length(s) == 1 && isequal(s[1][1], 1)
    end
    #the root is an arbitrary expression
    @test isZeroTens(simplify(∂dv(μ => lower, ∂dv(ν => lower, Xs * Ys + Zs)) - ∂dv(ν => lower, ∂dv(μ => lower, Ys * Xs + Zs))))
    #and may hold chains of its own
    @test isZeroTens(simplify(∂dv(μ => lower, Ys * ∂dv(ν => lower, ∂dv(λ => lower, Xs))) -
                              ∂dv(μ => lower, Ys * ∂dv(λ => lower, ∂dv(ν => lower, Xs)))))
    #fermionic partials anticommute
    ψs = Tensor("ψ", [x, θ])
    @test isZeroTens(simplify(∂θv(i => lower, ∂θv(j => lower, ψs)) + ∂θv(j => lower, ∂θv(i => lower, ψs))))
    @test !isZeroTens(simplify(∂θv(i => lower, ∂θv(j => lower, ψs)) - ∂θv(j => lower, ∂θv(i => lower, ψs))))
    #and a bosonic one commutes past them
    @test isZeroTens(simplify(∂dv(μ => lower, ∂θv(i => lower, ψs)) - ∂θv(i => lower, ∂dv(μ => lower, ψs))))
    #reorderings that are equal hash alike, as merging relies on
    let p = ∂dv(λ => lower, ∂dv(μ => lower, ∂dv(ν => lower, Xs))), q = ∂dv(ν => lower, ∂dv(λ => lower, ∂dv(μ => lower, Xs)))
        @test p == q && hash(p) == hash(q)
    end
    #and symmetrise_partials puts both orders of a chain in the same canonical one
    @test symmetrise_partials(∂dv(μ => lower, ∂dv(ρ => lower, Tensor("A", [ν => lower], x)))) ==
          symmetrise_partials(∂dv(ρ => lower, ∂dv(μ => lower, Tensor("A", [ν => lower], x))))

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

    ###############################################
    # a metric raising a derivative's own index must refresh THAT derivative's cached indices. The
    # refresh used to land on the factor after it, so ∂^b kept the cache of ∂_a: the term then had no
    # dummies, was never relabelled, and two equal terms never merged
    let s2 = IndexSet("c2_contract", (0, 2), lower)
        ka, kb, kc, kd = Index(["ka", "kb", "kc", "kd"], s2)
        kx = Coordinate("kx", [s2])
        kp = Tensor("kp", kx)
        Dk(i) = PartialDerivative("∂", kp, [kx], [i => lower])
        gu(i, j) = Metric([i => upper, j => upper])
        gl(i, j) = Metric([i => lower, j => lower])
        w = gu(kc, kd) * Dk(kc) * Dk(kd)
        @test isZeroTens(simplify(gu(ka, kb) * gl(ka, kb) * w - 2 * w))
        @test isZeroTens(simplify(gu(ka, kb) * (-2 * Dk(ka) * Dk(kb)) + gu(ka, kb) * gl(ka, kb) * w))
        let tt = summands(Alakazam.contract_metrics(TensorExpression(gu(ka, kb) * Dk(ka) * Dk(kb)), Metric))[1][2]
            @test all(t -> indices(t) == op_indices(t), terms(tt))
            @test length(dummies(tt)) == 1
        end
    end
end