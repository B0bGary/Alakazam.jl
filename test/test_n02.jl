#The N=(0,2) supercurrent-squared calculation reproduced in the paper.
@testset "N=(0,2) supercurrent squared" begin

    x  = Coordinate("n02x", lightcone_indices)
    θ  = Coordinate("n02θ", lightcone_indices)
    θb = Coordinate("n02Θ", lightcone_indices)
    custom_sort_order!(["n02θ", "n02Θ", "ϵ"])
    θp, θbp = Tensor("n02θ", [ℓ⁺ => upper], θ), Tensor("n02Θ", [ℓ⁺ => upper], θb)
    ∂pp = PD("∇", [x], [ℓ⁺ => lower, ℓ⁺ => lower])
    ∂mm = PD("∇", [x], [ℓ⁻ => lower, ℓ⁻ => lower])
    ∂θ, ∂θb = PD("∂", [θ], [ℓ⁺ => lower]), PD("𝛛", [θb], [ℓ⁺ => lower])
    @action ∂θ(θp) = 1
    @action ∂θb(θbp) = 1
    D    = Derivative("nD",  [x, θ],  [ℓ⁺ => lower])
    Dbar = Derivative("̄nD", [x, θb], [ℓ⁺ => lower])
    @covariant D    =  ∂θ  - im//2 * θbp * ∂pp
    @covariant Dbar = -∂θb + im//2 * θp  * ∂pp
    @commutator {D, Dbar} = im * ∂pp

    sf(n, inds) = generate_superfield_naive(n, inds, [θp, θbp], x)
    S, T = sf("nS", [ℓ⁺ => lower, ℓ⁺ => lower]), sf("nT", [ℓ⁻ => lower for _ in 1:4])
    W, Wb = sf("nW", [ℓ⁻ => lower]), sf("̄nW", [ℓ⁻ => lower])

    c(F, a, b) = isTensor(grassmann_coefficient(F, [θp, θbp], [a, b]))
    sub(e, rs) = foldl((a, r) -> replace_pattern(a, TensorTerm(r[1]), r[2]), rs; init=TensorExpression(e))

    ###############################################
    ### grassmann_order and drop_nilpotent
    ###############################################
    # a single lightcone component is nilpotent at order 1; a ± spans two
    @test grassmann_order(θp) == 1 && grassmann_order(θbp) == 1
    @test grassmann_order(Tensor("n02pm", [ℓ⁺⁻ => upper], θ)) == 2
    let A = Tensor("nA", x)
        @test drop_nilpotent(TensorExpression(θp * θp * A) + TensorExpression(θp * A), [θp, θbp]) ==
              TensorExpression(θp * A)
    end
    # reading a component off the expansion agrees with hunting for its label
    @test c(W, 1, 0) == only(t for (_, tt) in summands(W) for t in terms(tt)
                               if any(l -> l[1] == "(1,0)", labels(t)))

    ###############################################
    ### solving the constraints
    ###############################################
    # each constraint is an equation  for the superfield on its left, so that is the field whose
    # components it determines
    cons = [(determines = W,  vanishes = (F) -> Dbar(F[3])),
            (determines = Wb, vanishes = (F) -> D(F[4])),
            (determines = T,  vanishes = (F) -> D(F[2])    - 1//2 * ∂mm(F[4])),
            (determines = T,  vanishes = (F) -> Dbar(F[2]) - 1//2 * ∂mm(F[3])),
            (determines = S,  vanishes = (F) -> ∂mm(F[1])  - D(F[3]) + Dbar(F[4]))]

    sols, diffs = Pair[], Pair[]
    # `derive` with the Grassmann basis drops the nilpotent tails the covariant derivatives make
    expand(i) = derive(cons[i].vanishes([sub(F, vcat(sols, diffs)) for F in (S, T, W, Wb)]), [θp, θbp])

    function unknown_in(eqn, F)
        target = Alakazam.name(c(F, 0, 0))
        for (_, tt) in summands(eqn)
            length(terms(tt)) == 1 || continue
            factor = terms(tt)[1]
            Alakazam.name(first(strip_derivatives(factor))) == target && return factor
        end
        error("no lone factor of $target appears")
    end

    for (i, con) in enumerate(cons), (a, b) in ((0, 0), (1, 0), (0, 1), (1, 1))
        equation = grassmann_coefficient(expand(i), [θp, θbp], [a, b])
        isZeroTens(equation) && continue                 # nothing new at this order
        (unknown, value) = solve_linear(equation, unknown_in(equation, con.determines))
        push!(unknown isa PartialSuperType ? diffs : sols, unknown => value)
    end

    ###############################################
    ### the twelve relations that come out
    ###############################################
    # seven fix a component outright, five relate only a derivative of one
    @test length(sols) == 7 && length(diffs) == 5
    @test isZeroTens(Dict(sols)[c(W, 0, 1)])                                     # D̄𝒲 = 0 at θ⁰
    @test Dict(sols)[c(W, 1, 1)] == -im//2 * ∂pp(c(W, 0, 0))
    @test Dict(sols)[c(T, 1, 0)] == 1//2 * ∂mm(c(Wb, 0, 0))                      # B₊₋₋₋₋ = ½∂₋₋w̄₋
    @test Dict(diffs)[∂mm(c(S, 1, 0))] == -im * ∂pp(c(Wb, 0, 0))                 # ∂₋₋S₊₊₊ = −i∂₊₊w̄₋
    @test Dict(diffs)[∂mm(c(S, 0, 0))] == c(W, 1, 0) + c(Wb, 0, 1)               # ∂₋₋j₊₊ = F + F̄
    # the split is exactly by whether the solved-for object carries a derivative
    @test all(r[1] isa PartialSuperType for r in diffs)
    @test !any(r[1] isa PartialSuperType for r in sols)

    ###############################################
    ### every constraint then holds identically
    ###############################################
    for i in 1:5
        @test isZeroTens(Alakazam.simplify(symmetrise_partials(
            Alakazam.simplify(sub(expand(i), vcat(sols, diffs))))))
    end

    ###############################################
    ### the operator itself
    ###############################################
    raw = integrate(TensorExpression(T * S - Wb * W), [D, Dbar], false)
    @test length(summands(raw)) == 12
    # the four ∂₊₊ terms the covariant measure contributes are a total derivative
    @test length(summands(drop_total_derivatives(raw))) == 8

    O = reduce_bp_with(raw, vcat(sols, diffs))
    # O_{T²} = −S₍₁,₁₎T₍₀,₀₎ + ¼(W₍₁,₀₎ − W̄₍₀,₁₎)², i.e. T₊₊₊₊T₋₋₋₋ − (Im F₊₋)², since the
    # multiplet is written 𝒮₊₊ = j₊₊ + θ⁺S₊₊₊ − θ̄⁺S̄₊₊₊ − θ⁺θ̄⁺T₊₊₊₊ so that S₍₁,₁₎ = −T₊₊₊₊
    @test isZeroTens(Alakazam.simplify(O -
        (-c(S, 1, 1) * c(T, 0, 0) + 1//4 * (c(W, 1, 0) - c(Wb, 0, 1)) * (c(W, 1, 0) - c(Wb, 0, 1)))))
    # the fermions cancel outright
    @test all(grading(f) == 0 for (_, tt) in summands(O) for f in terms(tt))

    ###############################################
    ### in its usual form
    ###############################################
    # ∂₋₋j₊₊ = F₊₋ + F̄₊₋ fixes Re F₊₋ = 1/2∂₋₋j₊₊, so the independent piece of 𝒲₋ is Im F₊₋. Naming
    # it Θ₊₋ leaves O_{T²} = T₊₊₊₊T₋₋₋₋ − Θ₊₋², Zamolodchikov's TT̄ − Θ².
    Θ  = Tensor("nΘ", indices(c(W, 1, 0)), x)
    Tp = Tensor("n𝐓", indices(c(S, 1, 1)), x)
    trace_form = Alakazam.simplify(sub(O, [c(W, 1, 0)  => 1//2 * ∂mm(c(S, 0, 0)) + im * Θ,
                                           c(Wb, 0, 1) => 1//2 * ∂mm(c(S, 0, 0)) - im * Θ,
                                           c(S, 1, 1)  => -Tp]))
    @test length(summands(trace_form)) == 2
    @test isZeroTens(Alakazam.simplify(trace_form - (Tp * c(T, 0, 0) - Θ * Θ)))

    clear_sort_order!()
end
