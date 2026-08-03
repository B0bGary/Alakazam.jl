@testset "Tensor Replacement Patterns" begin

    #setup
    if 1 == 1
        spacetime_indices = IndexSet("spacetime", 4, upper)
        @indices μ:τ ∈ spacetime_indices

        flavour_indices = IndexSet("flavour_indices", 2, upper, Alakazam.fermionic)
        @indices a:z ∈ flavour_indices

        x = Coordinate("x", [flavour_indices])
        θ = Coordinate("θ", [flavour_indices])
        θbar = Coordinate("Θ", [flavour_indices])
    end
    ###############################################
    ### begin test cases
    ###############################################

    #TENSORS

    #Tensor addition basics
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    S1 = Tensor("S", [a => lower, b => lower, c => lower], x)
    #tensor exp
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    R1 = Tensor("R", [a => lower, b => lower, c => lower], x)
    T4 = Tensor("T", [s => lower, t => upper, u => lower], x)
    T5 = Tensor("T", [d => lower, e => upper, t => lower], x)

    #check pattern matching parsing and behaviour with grading and wildcards
    term=T4*T5*T1
    @test replace_pattern(term, "?S_{?a?e?f}", R1) == TensorExpression(R1*T4*T5)
    term=T4*T1
    @test replace_pattern(term, "?S_{?a?e?f}", R1) == -R1*T4

    f1 = AntisymmetricTensor("F", [μ => lower, ν=>lower], x)
    f2 = AntisymmetricTensor("F", [μ => upper, ν=>upper], x)
    L=-1//4*f1*f2

    a=DummyPatternIndex("a", spacetime_indices)
    b=DummyPatternIndex("b", spacetime_indices)
    da=PartialDerivative("∂", Tensor("A", [b=>lower], x), x, a=>lower)
    db=PartialDerivative("∂", Tensor("A", [a=>lower], x), x, b=>lower)
    dau=PartialDerivative("∂", Tensor("A", [b=>upper], x), x, a=>upper)
    dbu=PartialDerivative("∂", Tensor("A", [a=>upper], x), x, b=>upper)

    F2inA=-1//4*(2*PartialDerivative("∂", Tensor("A", [ν=>lower], x), x, μ=>lower)*PartialDerivative("∂", Tensor("A", [ν=>upper], x), x, μ=>upper)-2*PartialDerivative("∂", Tensor("A", [μ=>upper], x), x, ν=>upper)*PartialDerivative("∂", Tensor("A", [ν=>lower], x), x, μ=>lower))

    #nontrivial replacement pattern with products of operators in the replacement
    #wildcard
    @test Alakazam.simplify(replace_pattern(replace_pattern(L, "?T_{?a?b}", da-db), "?T^{?a?b}", dau-dbu)) == F2inA
    #exact
    @test Alakazam.simplify(replace_pattern(replace_pattern(L, "F_{?a?b}", da-db), "F^{?a?b}", dau-dbu)) == F2inA
    #mixed exact and wildcard indices
    da2=PartialDerivative("∂", DummyPatternTensor("A", [b=>lower], x), x, μ=>lower)
    db2=PartialDerivative("∂", DummyPatternTensor("A", [μ=>lower], x), x, b=>lower)
    @test Alakazam.simplify(replace_pattern(replace_pattern(L, "F_{μ?b}", da2-db2), "F^{?a?b}", dau-dbu)) == F2inA
    da3=PartialDerivative("∂", DummyPatternTensor("A", [ν=>lower], x), x, μ=>lower)
    db3=PartialDerivative("∂", DummyPatternTensor("A", [μ=>lower], x), x, ν=>lower)
    @test Alakazam.simplify(replace_pattern(replace_pattern(L, "F_{μν}", da3-db3), "F^{?a?b}", dau-dbu)) == F2inA
    @test Alakazam.simplify(replace_pattern(replace_pattern(L, "?G_{μν}", da3-db3), "F^{?a?b}", dau-dbu)) == F2inA

    #replace in a single go by popping off the metrics.
    upperF2=pop_metric(L)
    @test Alakazam.hard_simplify(contract_metrics(replace_pattern(upperF2, "?T^{?a?b}", dau-dbu))) == F2inA
    @test Alakazam.hard_simplify(contract_metrics(replace_pattern(upperF2, "F^{?a?b}", dau-dbu))) == F2inA


    ### test bianchi...
    R1 = RiemannTensor([μ => lower, ν => lower, σ => lower, τ => lower], x)
    D1 = Partial("∇", R1, x, ρ => lower)   # +abe
    D2 = Partial("∇", RiemannTensor([ν => lower, ρ => lower, σ => lower, τ => lower], x), x, μ => lower)   # +bea
    D3 = Partial("∇", RiemannTensor([ρ => lower, μ => lower, σ => lower, τ => lower], x), x, ν => lower)   # +eab
    c=DummyPatternIndex("c", spacetime_indices)
    d=DummyPatternIndex("d", spacetime_indices)
    ddg=1//2*(PD("∇", (PD("∇", SymmetricTensor("g", [a=>lower, d=>lower], x), x, c=>lower)), x, b=>lower)
              +
              PD("∇", (PD("∇", SymmetricTensor("g", [b=>lower, c=>lower], x), x, d=>lower)), x, a=>lower)
              -
              PD("∇", (PD("∇", SymmetricTensor("g", [a=>lower, c=>lower], x), x, d=>lower)), x, b=>lower)
              -
              PD("∇", (PD("∇", SymmetricTensor("g", [b=>lower, d=>lower], x), x, c=>lower)), x, a=>lower))
    bi = replace_pattern(D1+D2+D3, "R_{?a?b?c?d}", ddg)
    @test Alakazam.hard_simplify(apply_derivative(bi))==Alakazam. zero()

    ###############################################
    ### μlti-tensor patterns
    ### A pattern string parses into a TensorTerm, so it may name several tensors at once.
    ###############################################
    Zs = Tensor("Z")
    Am(p) = Tensor("A", [μ => p]); Bm(p) = Tensor("B", [μ => p])
    Cn = Tensor("C", [ν => upper]); Km(p) = Tensor("K", [μ => p])

    # a contracted pair collapses to a scalar, by explicit index name or by wildcard
    @test replace_pattern(Am(lower) * Bm(upper), "A_{μ}B^{μ}", Zs) == TensorExpression(Zs)
    @test replace_pattern(Am(lower) * Bm(upper), "A_{?x}B^{?x}", Zs) == TensorExpression(Zs)
    # matching is order-independent: the pattern need not appear contiguously or in order
    @test replace_pattern(Bm(upper) * Am(lower), "A_{?x}B^{?x}", Zs) == TensorExpression(Zs)
    # factors outside the match are preserved
    @test replace_pattern(Am(lower) * Bm(upper) * Cn, "A_{?x}B^{?x}", Zs) == Zs * Cn
    # a self-contraction can be sent to zero, e.g. imposing an on-shell condition k_a k^a = 0
    @test replace_pattern(Km(lower) * Km(upper) + Am(lower) * Am(upper), "K_{?x}K^{?x}", zero()) ==
          Am(lower) * Am(upper)
    # three tensors in one pattern (two was the most that used to work)
    @test replace_pattern(Am(lower) * Bm(upper) * Cn, "A_{?x}B^{?x}C^{?y}", Zs) == TensorExpression(Zs)
    # wildcard tensor names inside a μlti-tensor pattern
    @test replace_pattern(Am(lower) * Bm(upper), "?P_{?x}?Q^{?x}", Zs) == TensorExpression(Zs)
    # the replacement may itself be a sum
    @test replace_pattern(Am(lower) * Bm(upper), "A_{?x}B^{?x}", Zs + 2 * Tensor("W")) ==
          Zs + 2 * Tensor("W")
    # applied to every summand
    @test replace_pattern(Am(lower) * Bm(upper) * Cn + Km(lower) * Km(upper) * Cn, "A_{?x}B^{?x}", Zs) ==
          Zs * Cn + Km(lower) * Km(upper) * Cn
    # a pattern that matches nothing leaves the expression untouched
    @test replace_pattern(Am(lower) * Cn, "A_{?x}B^{?x}", Zs) == Am(lower) * Cn

    # grading: matched factors are dragged together and pick up the graded sign
    Ff(p) = Tensor("F", [e => p]); Gf(p) = Tensor("G", [e => p]); Hf = Tensor("H", [f => upper])
    @test grading(Ff(lower)) == 1 && grading(Hf) == 1
    @test replace_pattern(Ff(lower) * Gf(upper) * Hf, "F_{?x}G^{?x}", Zs) == Zs * Hf        # adjacent
    @test replace_pattern(Ff(lower) * Hf * Gf(upper), "F_{?x}G^{?x}", Zs) == -(Zs * Hf)     # G past H


    ###############################################
    ### labels, and the DT"..." pattern macro
    ###############################################
    # labels parse from the same string syntax as indices: a digit, or anything not naming a
    # variable in scope, is a label; 'quotes' force a label even when that name IS in scope.
    @test labels(T"k_1^μ") == [("1" => lower)] && length(indices(T"k_1^μ")) == 1
    @test labels(T"A_'μ'ν") == [("μ" => lower)] && length(indices(T"A_'μ'ν")) == 1
    @test labels(T"A_'bare'μ") == [("bare" => lower)]
    @test T"k_1^μ" != T"k_2^μ"

    # DT"..." shares that parser but emits pattern objects
    @test DT"?A_?x?y" isa DummyPatternTensor
    @test all(p[1] isa DummyPatternIndex for p in indices(DT"?A_?x?y"))
    @test DT"A_?x" isa DummyPatternTensor && Alakazam.name(DT"A_?x") == "A"

    # ... and can be used directly as a pattern, singly or as a product
    Pk(p) = Tensor("P", [μ => p]); Qk(p) = Tensor("Q", [μ => p]); Rk = Tensor("R", [ν => upper])
    Zd = Tensor("Z")
    @test replace_pattern(TensorExpression(Pk(lower) * Qk(upper) * Rk),
                          TensorTerm(DT"P_?x"), TensorExpression(Zd)) == Zd * Qk(upper) * Rk
    @test replace_pattern(TensorExpression(Pk(lower) * Qk(upper) * Rk),
                          DT"P_?x" * DT"Q^?x", TensorExpression(Zd)) == Zd * Rk

    # a pattern naming no label is a wildcard over labels; one that names a label must match it
    k1e = TensorExpression(T"k_1^μ" * Tensor("m", [μ => lower]))
    k2e = TensorExpression(T"k_2^μ" * Tensor("m", [μ => lower]))
    @test replace_pattern(k1e, TensorTerm(DT"k^?x"), TensorExpression(Zd)) != k1e   # unlabelled
    @test replace_pattern(k2e, TensorTerm(DT"k^?x"), TensorExpression(Zd)) != k2e   # matches both
    @test replace_pattern(k1e, TensorTerm(DT"k_1^?x"), TensorExpression(Zd)) != k1e # labelled
    @test replace_pattern(k2e, TensorTerm(DT"k_1^?x"), TensorExpression(Zd)) == k2e # k_2 untouched


    ###############################################
    ### wildcard indices must bind consistently across the whole pattern
    ### (a repeated ?x demands a genuine contraction, not two unrelated indices)
    ###############################################
    Wq(x, p) = Tensor("W", [x => p])
    Vv(x, p) = Tensor("V", [x => p])
    Zw = Tensor("Z")
    @test replace_pattern(Wq(μ, lower) * Wq(μ, upper), "W_{?x}W^{?x}", Zw) == TensorExpression(Zw)
    @test replace_pattern(Vv(μ, lower) * Wq(μ, upper), "V_{?x}W^{?x}", Zw) == TensorExpression(Zw)
    # different indices: ?x cannot bind to both, so these must NOT match
    @test replace_pattern(Wq(μ, lower) * Wq(ν, upper), "W_{?x}W^{?x}", Zw) ==
          Wq(μ, lower) * Wq(ν, upper)
    @test replace_pattern(Vv(μ, lower) * Wq(ν, upper), "V_{?x}W^{?x}", Zw) ==
          Vv(μ, lower) * Wq(ν, upper)
    # distinct wildcards may of course bind to distinct indices
    @test replace_pattern(Vv(μ, lower) * Wq(ν, upper), "V_{?x}W^{?y}", Zw) == TensorExpression(Zw)


    ###############################################
    ### the replacement is spliced in where the match began
    ### (for non-commuting factors the position is part of the answer)
    ###############################################
    @matrix RepMx
    Ax = RepMx("A"); Bx = RepMx("B"); Cx = RepMx("C"); Dx = RepMx("D"); Xx = RepMx("X")

    @test replace_pattern(TensorExpression(Ax*Bx*Cx),    "B",  Xx) == Ax*Xx*Cx
    @test replace_pattern(TensorExpression(Ax*Bx*Cx*Dx), "C",  Xx) == Ax*Bx*Xx*Dx
    @test replace_pattern(TensorExpression(Ax*Bx*Cx),    "A",  Xx) == Xx*Bx*Cx
    @test replace_pattern(TensorExpression(Ax*Bx*Cx),    "C",  Xx) == Ax*Bx*Xx
    # a contiguous multi-tensor match collapses to the anchor, in place
    @test replace_pattern(TensorExpression(Ax*Bx*Cx*Dx), "BC", Xx) == Ax*Xx*Dx
    # a split match: both factors are consumed, the one between them stays put
    @test replace_pattern(TensorExpression(Ax*Bx*Cx*Dx), "BD", Xx) == Ax*Xx*Cx

    # the shape a spin sum takes: u(p1) ubar(p1) -> projector, inside a Dirac chain.
    # the projector must land between the gamma matrices, not at the front of the term.
    Ub = RepMx("U", "2" => lower); Uu = RepMx("u", "2" => lower)
    u1 = RepMx("u", "1" => lower); U1 = RepMx("U", "1" => lower)
    ga = Gamma(μ => upper); gb = Gamma(ν => upper)
    chain = TensorExpression(Ub * ga * u1 * U1 * gb * Uu)
    @test replace_pattern(chain, DT"u_1" * DT"U_1", TensorExpression(Xx)) == Ub * ga * Xx * gb * Uu


    ###############################################
    ### non-commuting factors between the match are crossed with the registered
    ### commutation relation, not merely a sign
    ###############################################
    @matrix DragMx
    Ad = DragMx("A"); Bd = DragMx("B"); Cd = DragMx("C"); Xd = DragMx("X"); Ed = DragMx("E")
    Symbolics.@syms hh
    @commutator [Cd, Bd] = hh * Ed          # so B cannot be moved past C for free

    @test !Alakazam.can_swap(Bd, Cd)
    # adjacent: nothing moves, no remainder
    @test replace_pattern(TensorExpression(Ad * Bd * Cd), "AB", Xd) == Xd * Cd
    # separated by C:  A C B = A(BC + [C,B]) = ABC + h A E,  and only the first part matches
    @test replace_pattern(TensorExpression(Ad * Cd * Bd), "AB", Xd) == Xd * Cd + hh * Ad * Ed
    # a factor that DOES commute is still crossed for free, with no remainder
    Fd = DragMx("F")
    @test Alakazam.can_swap(Bd, Fd)
    @test replace_pattern(TensorExpression(Ad * Fd * Bd), "AB", Xd) == Xd * Fd

end
