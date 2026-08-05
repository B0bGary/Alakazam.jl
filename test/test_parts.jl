#Integration by parts: reduction modulo total derivatives.
@testset "Integration by parts" begin

    x = Coordinate("ibp_x", lightcone_indices)
    ∂pp = PD("∇", [x], [ℓ⁺ => lower, ℓ⁺ => lower])
    ∂mm = PD("∇", [x], [ℓ⁻ => lower, ℓ⁻ => lower])

    A = Tensor("A", [ℓ⁺ => lower, ℓ⁺ => lower], x)                              # even
    B = Tensor("B", [ℓ⁻ => lower, ℓ⁻ => lower, ℓ⁻ => lower, ℓ⁻ => lower], x)    # even
    ψ = Tensor("ψ", [ℓ⁻ => lower], x)                                           # odd
    χ = Tensor("χ", [ℓ⁻ => lower], x)                                           # odd
    @test grading(A) == 0 && grading(B) == 0
    @test grading(ψ) == 1 && grading(χ) == 1

    dtd(e) = drop_total_derivatives(e)

    ###############################################
    ### peeling and rebuilding
    ###############################################
    (b, ops) = strip_derivatives(∂pp(∂mm(A)))
    @test b == A
    @test length(ops) == 2
    @test rebuild_nested(ops, b) == ∂pp(∂mm(A))
    # outermost first
    @test op_indices(ops[1]) == op_indices(∂pp)
    # a bare tensor peels to itself
    @test strip_derivatives(A) == (A, Alakazam.PartialSuperType[])
    # a derivative over a product is not of the form ∂ᵏφ, so it is left whole
    @test isempty(strip_derivatives(∂pp(A * B))[2])

    ###############################################
    ### total derivatives vanish
    ###############################################
    # the defining case: ∂(AB) = (∂A)B + A(∂B)
    @test isZeroTens(dtd(A * ∂pp(B) + ∂pp(A) * B))
    @test isZeroTens(dtd(∂pp(A * B)))
    # with a coefficient carried along
    @test isZeroTens(dtd(im//2 * A * ∂pp(B) + im//2 * B * ∂pp(A)))
    # two odd fields: the graded reordering has to be picked up for this to cancel
    @test isZeroTens(dtd(im//2 * ψ * ∂pp(χ) - im//2 * χ * ∂pp(ψ)))
    @test isZeroTens(dtd(∂pp(ψ * χ)))
    # second derivatives, sign (-1)² = +1
    @test isZeroTens(dtd(A * ∂mm(∂mm(B)) - ∂mm(∂mm(A)) * B))
    # third, sign (-1)³ = -1
    @test isZeroTens(dtd(A * ∂mm(∂mm(∂mm(B))) + ∂mm(∂mm(∂mm(A))) * B))
    # derivatives of two different kinds
    @test isZeroTens(dtd(A * ∂pp(∂mm(B)) - ∂mm(∂pp(A)) * B))
    # a lone total derivative, of one factor or of a product
    @test isZeroTens(dtd(TensorExpression(∂pp(A))))
    @test isZeroTens(dtd(TensorExpression(∂pp(∂mm(A)))))

    ###############################################
    ### a field against itself
    ###############################################
    # A∂A = ½∂(A²) for even A, so it drops
    @test isZeroTens(dtd(TensorExpression(A * ∂pp(A))))
    @test isZeroTens(dtd(TensorExpression(B * ∂mm(∂mm(∂mm(B))))))
    # but ψ∂ψ with ψ odd is not a total derivative: ∂(ψψ) = (∂ψ)ψ + ψ∂ψ = 0 identically
    @test !isZeroTens(dtd(TensorExpression(ψ * ∂pp(ψ))))
    # even A against an even number of derivatives is likewise not one
    @test !isZeroTens(dtd(TensorExpression(A * ∂pp(∂pp(A)))))

    ###############################################
    ### what must survive
    ###############################################
    @test dtd(TensorExpression(A * ∂pp(B))) == TensorExpression(A * ∂pp(B))
    @test isZeroTens(dtd(TensorExpression(A * B)) - TensorExpression(A * B))   # no derivatives at all
    @test dtd(TensorExpression(3)) == TensorExpression(3)                      # a bare constant
    # a genuine remainder is kept with the right coefficient: the first two cancel
    @test dtd(im * A * ∂pp(B) + im * B * ∂pp(A) - 2im * A * ∂pp(B)) ==
          simplify(-2im * A * ∂pp(B))
    # idempotent, as a normal form must be
    r = dtd(im * A * ∂pp(B) + im * B * ∂pp(A) - 2im * A * ∂pp(B))
    @test dtd(r) == r

    ###############################################
    ### linearity
    ###############################################
    e1 = A * ∂pp(B) + ∂pp(A) * B
    e2 = 5 * A * ∂pp(B)
    @test dtd(e1 + e2) == dtd(e1) + dtd(e2)

    ###############################################
    ### graded derivatives
    ###############################################
    # A partial is not always even-graded so integrating by parts
    # picks up a graded sign as well as the minus: 
    # (∂X)Y ≡ -(-1)^{|∂||X|} X(∂Y). Which grading  enters depends on direction,
    # moving right the derivative passes what remains of the left factor, moving left it passes it all.
    let θc = Coordinate("ibp_θ", lightcone_indices),
        ∂θ = PD("∂", [θc], [ℓ⁺ => lower]),
        Aθ = Tensor("A", [ℓ⁺ => lower, ℓ⁺ => lower], [x, θc]),      # even, depends on θ
        Bθ = Tensor("B", [ℓ⁻ => lower, ℓ⁻ => lower], [x, θc]),      # even
        ψθ = Tensor("ψ", [ℓ⁻ => lower], [x, θc]),                   # odd
        χθ = Tensor("χ", [ℓ⁻ => lower], [x, θc])                    # odd

        @test Alakazam.op_grading(∂θ) == 1        # the derivative itself is odd
        @test Alakazam.op_grading(∂pp) == 0       # while ∂₊₊ carries two indices and is even

        # ∂(XY) = (∂X)Y + (-1)^{|∂||X|}X(∂Y), so this combination is a total derivative
        for (X, Y) in ((Aθ, Bθ), (ψθ, Bθ), (Bθ, ψθ), (ψθ, χθ))
            s = iseven(grading(X)) ? 1 : -1
            @test isZeroTens(dtd(∂θ(X) * Y + s * (X * ∂θ(Y))))
            # and the orthogonal combination is not one
            @test !isZeroTens(dtd(∂θ(X) * Y - s * (X * ∂θ(Y))))
        end
        # an even derivative is unaffected by any of this
        @test isZeroTens(dtd(∂pp(ψθ) * χθ + ψθ * ∂pp(χθ)))
    end

    ###############################################
    ### beyond quadratic: handled by the linear-algebra route
    ###############################################
    # The quadratic normal form does not reach three factors, so these fall through to
    # `drop_total_derivatives_general`, which builds the total-derivative subspace explicitly.
    # A(∂A)B ≡ -½A²∂B, since ∂(A²B) = 2A(∂A)B + A²(∂B).
    @test dtd(TensorExpression(A * B * ∂pp(A))) ==
          simplify(-1//2 * A * A * ∂pp(B))
    # a genuine cubic total derivative goes to zero
    @test isZeroTens(dtd(∂pp(A)*B*B + A*∂pp(B)*B + A*B*∂pp(B)))
    @test isZeroTens(drop_total_derivatives_general(∂pp(A)*B*B + A*∂pp(B)*B + A*B*∂pp(B)))
    # one that is not stays put
    @test !isZeroTens(dtd(TensorExpression(∂pp(A)*B*B)))
    # quadratic terms are unaffected by the presence of cubic ones: different sectors cannot mix
    @test isZeroTens(dtd(∂pp(A)*B + A*∂pp(B) + ∂pp(A)*B*B + A*∂pp(B)*B + A*B*∂pp(B)))

end
