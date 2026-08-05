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
    da2=PartialDerivative("∂", Tensor("A", [b=>lower], x), x, μ=>lower)
    db2=PartialDerivative("∂", Tensor("A", [μ=>lower], x), x, b=>lower)
    @test Alakazam.simplify(replace_pattern(replace_pattern(L, "F_{μ?b}", da2-db2), "F^{?a?b}", dau-dbu)) == F2inA
    da3=PartialDerivative("∂", Tensor("A", [ν=>lower], x), x, μ=>lower)
    db3=PartialDerivative("∂", Tensor("A", [μ=>lower], x), x, ν=>lower)
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
    ### multi-tensor patterns
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
    # three tensors in one pattern
    @test replace_pattern(Am(lower) * Bm(upper) * Cn, "A_{?x}B^{?x}C^{?y}", Zs) == TensorExpression(Zs)
    # wildcard tensor names inside a multi-tensor pattern
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
    ### labels and wildcards through the T"..." macro
    ###############################################
    # labels parse from the same string syntax as indices: a digit, or anything not naming a
    # variable in scope, is a label; 'quotes' force a label even when that name IS in scope.
    @test labels(T"k_1^μ") == [("1" => lower)] && length(indices(T"k_1^μ")) == 1
    @test labels(T"A_'μ'ν") == [("μ" => lower)] && length(indices(T"A_'μ'ν")) == 1
    @test labels(T"A_'bare'μ") == [("bare" => lower)]
    @test T"k_1^μ" != T"k_2^μ"

    # patterns and tensors come from the SAME macro: a wildcard is just a leading ?
    @test Alakazam.name(T"?A_?x?y") == "?A"
    @test all(p[1] isa DummyPatternIndex for p in indices(T"?A_?x?y"))
    @test T"A_?x" isa Tensor && Alakazam.name(T"A_?x") == "A"

    # ... and can be used directly as a pattern, singly or as a product
    Pk(p) = Tensor("P", [μ => p]); Qk(p) = Tensor("Q", [μ => p]); Rk = Tensor("R", [ν => upper])
    Zd = Tensor("Z")
    @test replace_pattern(TensorExpression(Pk(lower) * Qk(upper) * Rk),
                          TensorTerm(T"P_?x"), TensorExpression(Zd)) == Zd * Qk(upper) * Rk
    @test replace_pattern(TensorExpression(Pk(lower) * Qk(upper) * Rk),
                          T"P_?x" * T"Q^?x", TensorExpression(Zd)) == Zd * Rk

    # a pattern naming no label is a wildcard over labels; one that names a label must match it
    k1e = TensorExpression(T"k_1^μ" * Tensor("m", [μ => lower]))
    k2e = TensorExpression(T"k_2^μ" * Tensor("m", [μ => lower]))
    mlow = Tensor("m", [μ => lower])
    @test replace_pattern(k1e, TensorTerm(T"k^?x"), TensorExpression(Zd)) == Zd * mlow   # unlabelled
    @test replace_pattern(k2e, TensorTerm(T"k^?x"), TensorExpression(Zd)) == Zd * mlow   # matches both
    @test replace_pattern(k1e, TensorTerm(T"k_1^?x"), TensorExpression(Zd)) == Zd * mlow # labelled
    @test replace_pattern(k2e, TensorTerm(T"k_1^?x"), TensorExpression(Zd)) == k2e       # k_2 untouched


    ###############################################
    ### wildcard names in a replacement, and concrete ones left alone
    ###############################################
    # NB a wildcard is ? plus a SINGLE character; "?bw" parses as index ?b carrying the label w.
    # built without coordinate dependence, since T"..." patterns carry none
    gw = DummyPatternIndex("g", spacetime_indices)
    ep = Tensor("P", [μ => lower]) * Tensor("C", [ν => upper])
    epB = Tensor("B", [μ => lower]) * Tensor("C", [ν => upper])

    # a wildcard-named replacement is bound to whatever the pattern matched
    @test replace_pattern(TensorExpression(ep), TensorTerm(T"?A_?g"), TensorExpression(T"?A_?g")) == ep

    # a replacement whose name is concrete is NOT looked up in the bindings. This used to throw a
    # KeyError, because a pattern-typed replacement was taken to imply a wildcard name.
    @test replace_pattern(TensorExpression(ep), TensorTerm(T"?A_?g"), TensorExpression(T"B_?g")) == epB

    # and the same written as a plain Tensor, which is now the only representation
    @test replace_pattern(TensorExpression(ep), TensorTerm(T"?A_?g"),
                          TensorExpression(Tensor("B", [gw => lower]))) == epB

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
    # a split match cannot be gathered across C, which does not commute with B or D, so the
    # pattern does not apply and the term stands
    @test replace_pattern(TensorExpression(Ax*Bx*Cx*Dx), "BD", Xx) == Ax*Bx*Cx*Dx

    # the shape a spin sum takes: u(p1) ubar(p1) -> projector, inside a Dirac chain.
    # the projector must land between the gamma matrices, not at the front of the term.
    Ub = RepMx("U", "2" => lower); Uu = RepMx("u", "2" => lower)
    u1 = RepMx("u", "1" => lower); U1 = RepMx("U", "1" => lower)
    ga = Gamma(μ => upper); gb = Gamma(ν => upper)
    chain = TensorExpression(Ub * ga * u1 * U1 * gb * Uu)
    @test replace_pattern(chain, T"u_1" * T"U_1", TensorExpression(Xx)) == Ub * ga * Xx * gb * Uu


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
    # a factor with no relation registered is not assumed to commute, so it cannot be crossed
    # at all and the match is never gathered
    Fd = DragMx("F")
    @test !Alakazam.can_swap(Bd, Fd)
    @test replace_pattern(TensorExpression(Ad * Fd * Bd), "AB", Xd) == Ad * Fd * Bd
    # declaring it to commute restores the free crossing
    Gd = DragMx("G")
    @commutator [Ad, Gd] = zero()
    @commutator [Bd, Gd] = zero()
    @test Alakazam.can_swap(Bd, Gd)
    @test replace_pattern(TensorExpression(Ad * Gd * Bd), "AB", Xd) == Xd * Gd


    ###############################################
    ###############################################
    ### EDGE CASES
    ### Silent non-matching is the dominant failure mode here, so most of these pin down that
    ### something specific does NOT happen.
    ###############################################
    ###############################################

    Ze = Tensor("Z"); We = Tensor("W")
    Ae(p) = Tensor("A", [μ => p]); Be(p) = Tensor("B", [μ => p]); Ce = Tensor("C", [ν => upper])
    Aq(i, p) = Tensor("A", [i => p])

    ### --- parser limits: these fail SILENTLY, which is why they are pinned ---------------------
    # a tensor name in a pattern string is a single character
    @test replace_pattern(Ae(lower) * Ce, "AB_{?x}", Ze) == Ae(lower) * Ce
    # a wildcard is ? plus a single character; "?abc" is not a three-letter wildcard
    @test replace_pattern(Ae(lower) * Ce, "?abc", Ze) == Ae(lower) * Ce
    # a pattern with no indices matches only a tensor with no indices
    @test replace_pattern(Ae(lower) * Ce, "C", Ze) == Ae(lower) * Ce
    @test replace_pattern(Ae(lower) * Tensor("C"), "C", Ze) == Ae(lower) * Ze
    # a pattern longer than the term cannot match
    @test replace_pattern(Ae(lower), "A_{?x}B^{?x}C^{?y}", Ze) == TensorExpression(Ae(lower))
    # whitespace is insignificant, and a lone index needs no braces
    @test replace_pattern(Ae(lower) * Ce, " A_{?x} ", Ze) == Ze * Ce
    @test replace_pattern(Ae(lower) * Ce, "A_μ", Ze) == Ze * Ce
    # degenerate patterns are rejected rather than silently matching everything
    @test_throws Exception replace_pattern(Ae(lower) * Ce, "?", Ze)
    @test_throws Exception replace_pattern(Ae(lower) * Ce, "", Ze)

    ### --- mixed upper/lower slots need the macro, not a pattern string -------------------------
    # the string parser restarts its slot counter at each ^/_ group, so "A_{?x}^{?y}" describes a
    # ONE-index tensor and cannot match A_{μ}{}^{ν}.
    Mx2 = Tensor("A", [μ => lower, ν => upper])
    @test length(indices(Alakazam.parse_pattern2tensor("A_{?x}^{?y}"))) == 1
    @test replace_pattern(TensorExpression(Mx2), "A_{?x}^{?y}", Ze) == TensorExpression(Mx2)
    # T"..." goes through the T"..." parser, which tracks positions properly
    @test length(indices(T"A_?x^?y")) == 2
    @test replace_pattern(TensorExpression(Mx2), TensorTerm(T"A_?x^?y"), TensorExpression(Ze)) ==
          TensorExpression(Ze)

    ### --- wildcard binding ---------------------------------------------------------------------
    # a repeated wildcard NAME demands the same tensor throughout
    @test replace_pattern(Ae(lower) * Be(upper), "?A_{?x}?A^{?x}", Ze) == Ae(lower) * Be(upper)
    @test replace_pattern(Aq(μ, lower) * Aq(μ, upper), "?A_{?x}?A^{?x}", Ze) == TensorExpression(Ze)
    # distinct wildcards are NOT required to bind distinct tensors
    @test replace_pattern(Aq(μ, lower) * Aq(μ, upper), "?P_{?x}?Q^{?x}", Ze) == TensorExpression(Ze)

    ### --- names and special tensors ------------------------------------------------------------
    # dotted and barred names survive the parser
    @test Alakazam.name(T"Ȧ_?x") == "Ȧ" && Alakazam.name(T"Ā_?x") == "Ā"
    @test replace_pattern(Tensor("Ȧ", [μ => lower]) * Ce, "Ȧ_{?g}", Ze) == Ze * Ce
    # built-in tensors are matched by their printed name, and by wildcards
    @test Alakazam.name(Metric([μ => lower, ν => lower])) == "η"
    @test replace_pattern(Metric([μ => lower, ν => lower]) * Tensor("C", [ρ => upper]), "η_{?g?h}", Ze) ==
          Ze * Tensor("C", [ρ => upper])
    @test replace_pattern(Metric([μ => lower, ν => lower]) * Tensor("C", [ρ => upper]), "?P_{?g?h}", Ze) ==
          Ze * Tensor("C", [ρ => upper])
    # a Kronecker has mixed slots, so again it needs the macro
    @test replace_pattern(TensorExpression(KroneckerDelta(μ => lower, ν => upper) * Tensor("C", [ρ => upper])),
                          TensorTerm(T"δ_?g^?h"), TensorExpression(Ze)) == Ze * Tensor("C", [ρ => upper])
    # matching is by name, so a Fermion is matched by a plain pattern even though == now checks type
    @test replace_pattern(Fermion("A", [μ => lower]) * Ce, "A_{?g}", Ze) == Ze * Ce

    ### --- coefficients are carried through -----------------------------------------------------
    @test replace_pattern((3//4) * Ae(lower), "A_{?g}", Ze) == (3//4) * Ze
    @test replace_pattern(Ae(lower), "A_{?g}", 2 * Ze) == 2 * Ze
    @test replace_pattern(Ae(lower) * Ce, "A_{?g}", Ze + We) == Ze * Ce + We * Ce

    ### --- zero, and total annihilation ---------------------------------------------------------
    @test isZeroTens(replace_pattern(Ae(lower) * Ce, "A_{?x}", zero()))
    @test isZeroTens(replace_pattern(Ae(lower), "A_{?x}", zero()))

    ### --- graded signs accumulate one per crossing -----------------------------------------------
    Fe(p) = Tensor("F", [e => p]); Ge(p) = Tensor("G", [e => p])
    Hp = Tensor("H", [f => upper]); Jp = Tensor("J", [g => upper])
    @test grading(Fe(lower)) == 1 && grading(Hp) == 1 && grading(Jp) == 1
    @test replace_pattern(Fe(lower) * Ge(upper) * Hp * Jp, "F_{?x}G^{?x}", Ze) == Ze * Hp * Jp    # 0
    @test replace_pattern(Fe(lower) * Hp * Ge(upper) * Jp, "F_{?x}G^{?x}", Ze) == -(Ze * Hp * Jp) # 1
    @test replace_pattern(Fe(lower) * Hp * Jp * Ge(upper), "F_{?x}G^{?x}", Ze) == Ze * Hp * Jp    # 2

    ### --- commutator remainders, one per factor crossed ------------------------------------------
    # NB commutators are keyed by tensor NAME, so these must not reuse A/B/C from the block above.
    @matrix EdgeMx
    Pm = EdgeMx("P"); Qm = EdgeMx("Q"); Rm = EdgeMx("R"); Sm = EdgeMx("S"); Xm = EdgeMx("X")
    Symbolics.@syms hs
    @commutator [Rm, Qm] = hs * Sm
    @test !Alakazam.can_swap(Qm, Rm)
    @test replace_pattern(TensorExpression(Pm * Qm * Rm), "PQ", Xm) == Xm * Rm          # adjacent
    @test replace_pattern(TensorExpression(Pm * Rm * Qm), "PQ", Xm) == Xm * Rm + hs * Pm * Sm
    # two crossings give two remainders. They do not merge: S lands where each commutator was
    # formed, giving P S R from the outer crossing and P R S from the inner, and nothing licenses
    # reordering S past R
    @test replace_pattern(TensorExpression(Pm * Rm * Rm * Qm), "PQ", Xm) ==
          Xm * Rm * Rm + hs * Pm * Sm * Rm + hs * Pm * Rm * Sm
    # factors after the match are untouched
    @test replace_pattern(TensorExpression(Pm * Rm * Qm * EdgeMx("T")), "PQ", Xm) ==
          Xm * Rm * EdgeMx("T") + hs * Pm * Sm * EdgeMx("T")
    # a matrix with no registered relation cannot be crossed, so the match is not gathered
    @test replace_pattern(TensorExpression(Pm * EdgeMx("U") * Qm), "PQ", Xm) ==
          Pm * EdgeMx("U") * Qm

    ### --- operators --------------------------------------------------------------------------------
    gs = DummyPatternIndex("g", spacetime_indices)
    # a pattern reaches into an operator's operand
    @test replace_pattern(TensorExpression(PartialDerivative("∂", Tensor("A", [ν => lower], x), x, μ => lower)),
                          "A_{?g}", Tensor("B", [gs => lower], x)) ==
          TensorExpression(PartialDerivative("∂", Tensor("B", [ν => lower], x), x, μ => lower))
    # and a replacement may itself contain one
    @test replace_pattern(Tensor("Q", [ν => lower], x) * Tensor("C", [ρ => upper]), "Q_{?g}",
                          PartialDerivative("∂", Tensor("A", [gs => lower], x), x, μ => lower)) ==
          PartialDerivative("∂", Tensor("A", [ν => lower], x), x, μ => lower) * Tensor("C", [ρ => upper])

    ### --- repeat, and termination on self-referential replacements ---------------------------------
    @test replace_pattern(Aq(μ, lower) * Aq(ν, lower) * Aq(ρ, lower), "A_{?x}", Ze) == Ze * Ze * Ze
    @test replace_pattern(Aq(μ, lower) * Aq(ν, lower), "A_{?x}", Ze, false) == Ze * Aq(ν, lower)
    # replacing a tensor by itself, or by itself times something, must not loop forever
    @test replace_pattern(Ae(lower) * Ce, "A_{?x}", Ae(lower)) == Ae(lower) * Ce
    @test replace_pattern(Ae(lower) * Ce, "A_{?x}", Ae(lower) * We) == Ae(lower) * We * Ce

    ### --- hazard: a replacement may introduce an index already free in the term ---------------------
    # ν is free in C^{ν}, so this creates a contraction the user may not have intended. Pinned as
    # current behaviour: the replacement's indices are used as given, without renaming.
    @test replace_pattern(Ae(lower) * Ce, "A_{?g}", Tensor("B", [ν => lower])) ==
          Tensor("B", [ν => lower]) * Ce

    ### --- a decoy occurrence of the same name must not hide a match -------------------------------
    # a decoy of a DIFFERENT name is separated by the histogram
    @test replace_pattern(Tensor("D", [ν => lower]) * Aq(μ, lower) * Aq(μ, upper), "A_{?x}A^{?x}", Ze) ==
      Tensor("D", [ν => lower]) * Ze
    # a same-named decoy used to defeat the match: ?x bound to ν at position 1 and was never
    # reconsidered, though the pattern sits at positions 2,3. The scan now restarts from each
    # occurrence, so the later match is found.
    @test replace_pattern(Aq(ν, lower) * Aq(μ, lower) * Aq(μ, upper), "A_{?x}A^{?x}", Ze) ==
      Aq(ν, lower) * Ze
    @test replace_pattern(Aq(ν, lower) * Aq(ρ, lower) * Aq(μ, lower) * Aq(μ, upper), "A_{?x}A^{?x}", Ze) ==
      Aq(ν, lower) * Aq(ρ, lower) * Ze
    # a decoy in a slot the pattern cannot use never binds, so it was always harmless
    @test replace_pattern(Aq(ν, upper) * Aq(μ, lower) * Aq(μ, upper), "A_{?x}A^{?x}", Ze) ==
      Aq(ν, upper) * Ze
    # the matched factors need not be adjacent: here the pattern sits at positions 1 and 3
    @test replace_pattern(Aq(ν, lower) * Aq(μ, lower) * Aq(ν, upper), "A_{?x}A^{?x}", Ze) ==
      Ze * Aq(μ, lower)
    # two interleaved contractions, both replaced
    @test replace_pattern(Aq(μ, lower) * Aq(ν, lower) * Aq(ν, upper) * Aq(μ, upper), "A_{?x}A^{?x}", Ze) ==
      Ze * Ze
    # and a genuine non-match is still rejected
    @test replace_pattern(Aq(μ, lower) * Aq(ν, lower), "A_{?x}A^{?x}", Ze) ==
      Aq(μ, lower) * Aq(ν, lower)


    ###############################################
    ### ?A takes the matched NAME; ??A takes the matched OBJECT
    ###############################################
    firstt(v) = first(terms(first_summand(TensorExpression(v))))

    # a source carrying every kind of data a tensor can have
    Psrc = SymmetricTensor("P", [μ => lower, ν => lower], x, "lab" => lower)
    @test length(tableaux(Psrc)) == 1 && length(labels(Psrc)) == 1 && !isempty(Alakazam.function_of(Psrc))

    # ?A keeps only the name: coordinates, tableaux and labels come from the replacement as written
    viaName = firstt(replace_pattern(TensorExpression(Psrc), TensorTerm(T"?A_?g?h"), TensorExpression(T"?A_?g?h")))
    @test Alakazam.name(viaName) == "P"
    @test isempty(tableaux(viaName)) && isempty(labels(viaName)) && isempty(Alakazam.function_of(viaName))
    @test viaName != Psrc

    # ??A adopts the object, so replacing a tensor by itself is the identity
    viaObj = firstt(replace_pattern(TensorExpression(Psrc), TensorTerm(T"?A_?g?h"), TensorExpression(T"??A_?g?h")))
    @test viaObj == Psrc
    @test length(tableaux(viaObj)) == 1 && length(labels(viaObj)) == 1

    # the type comes with the object, so grading is preserved. Under ?A a Fermion would come back
    # as a Tensor and silently stop anticommuting.
    Fsrc = Fermion("P", [μ => lower, ν => lower], x)
    @test typeof(firstt(replace_pattern(TensorExpression(Fsrc), TensorTerm(T"?A_?g?h"),
                                        TensorExpression(T"?A_?g?h")))) == Tensor
    fObj = firstt(replace_pattern(TensorExpression(Fsrc), TensorTerm(T"?A_?g?h"), TensorExpression(T"??A_?g?h")))
    @test typeof(fObj) == Fermion && grading(fObj) == 1 && fObj == Fsrc

    # both forms name the same binding, and both parsers accept the doubled marker
    @test Alakazam.name(first(terms(Alakazam.parse_pattern2tensor("??A_{?g}")))) == "??A"
    @test Alakazam.name(first(terms(Alakazam.parse_pattern2tensor("?A_{?g}")))) == "?A"
    @test firstt(replace_pattern(TensorExpression(Psrc), "?A_{?g?h}", "??A_{?g?h}")) == Psrc


    ###############################################
    ### further edge cases
    ###############################################
    Pw(i, p) = Tensor("P", [i => p]); Qw(i, p) = Tensor("Q", [i => p])

    # a wildcard index will not bind across two different index sets
    @test replace_pattern(Pw(μ, lower) * Qw(e, upper), "P_{?g}Q^{?g}", Ze) == Pw(μ, lower) * Qw(e, upper)
    @test replace_pattern(Pw(μ, lower) * Qw(μ, upper), "P_{?g}Q^{?g}", Ze) == TensorExpression(Ze)
    # distinct wildcards may of course span different sets
    @test replace_pattern(Pw(μ, lower) * Qw(e, upper), "P_{?g}Q^{?h}", Ze) == TensorExpression(Ze)

    # one expression factor cannot fill two pattern slots
    @test replace_pattern(Aq(μ, lower), "A_{?x}A^{?y}", Ze) == TensorExpression(Aq(μ, lower))
    @test replace_pattern(Aq(μ, lower) * Aq(μ, upper), "A_{?x}A^{?x}A_{?y}", Ze) ==
          Aq(μ, lower) * Aq(μ, upper)
    # two disjoint copies are both replaced when repeating, one when not
    @test replace_pattern(Aq(μ, lower) * Aq(μ, upper) * Aq(ν, lower) * Aq(ν, upper), "A_{?x}A^{?x}", Ze) ==
          Ze * Ze
    @test replace_pattern(Aq(μ, lower) * Aq(μ, upper) * Aq(ν, lower) * Aq(ν, upper), "A_{?x}A^{?x}", Ze, false) ==
          Ze * Aq(ν, lower) * Aq(ν, upper)
    # the same wildcard twice on ONE tensor is a self-contraction
    @test replace_pattern(Tensor("A", [μ => lower, μ => lower]), "A_{?x?x}", Ze) == TensorExpression(Ze)

    # a pattern carrying no coordinates still matches a coordinate-dependent tensor
    @test replace_pattern(Tensor("P", [μ => lower], x) * Ce, "P_{?g}", Ze) == Ze * Ce

    # tableaux survive a replacement, and remain functionally live
    Asym = AntisymmetricTensor("A", [DummyPatternIndex("g", spacetime_indices) => lower,
                                     DummyPatternIndex("h", spacetime_indices) => lower])
    repl = firstt(replace_pattern(TensorExpression(Tensor("P", [μ => lower, ν => lower])), "P_{?g?h}", Asym))
    @test length(tableaux(repl)) == 1
    @test isZeroTens(Alakazam.hard_simplify(repl * SymmetricTensor("B", [μ => upper, ν => upper])))

    # a pattern reaches through two levels of operator nesting
    inner2 = PartialDerivative("∂", Tensor("A", [ν => lower], x), x, μ => lower)
    @test replace_pattern(TensorExpression(PartialDerivative("∂", inner2, x, ρ => lower)),
                          "A_{?g}", Tensor("B", [gs => lower], x)) ==
          TensorExpression(PartialDerivative("∂",
              PartialDerivative("∂", Tensor("B", [ν => lower], x), x, μ => lower), x, ρ => lower))

    # degenerate inputs are handled rather than crashing
    @test isZeroTens(replace_pattern(zero(), "A_{?g}", Ze))
    @test replace_pattern(Aq(μ, lower), "A_{?g}", Ze) == TensorExpression(Ze)
    # a replacement naming an index the pattern never bound is reported, not silently dropped
    @test_throws InvalidReplacementPattern replace_pattern(
        Aq(μ, lower), "A_{?g}", Tensor("B", [DummyPatternIndex("zz", spacetime_indices) => lower]))

    # replacement is idempotent once the pattern is gone
    once = replace_pattern(Aq(μ, lower) * Ce, "A_{?g}", Ze)
    @test replace_pattern(once, "A_{?g}", Ze) == once

    # a fermion sitting between two matched matrices is crossed, and the commutator still fires
    ferm = Tensor("F", [e => lower])
    @test replace_pattern(TensorExpression(Pm * ferm * Rm * Qm), "PQ", Xm) ==
          Xm * ferm * Rm + hs * Pm * Sm * ferm


    ###############################################
    ### wildcard names are resolved wherever they sit, and reported when unbound
    ###############################################
    Pco = Tensor("P", [μ => lower], x); Cup = Tensor("C", [ν => upper])

    # a wildcard nested in an operator operand used to leak the literal "?A"/"??A" into the result
    @test replace_pattern(TensorExpression(Pco * Cup), "?A_{?g}",
                          PartialDerivative("∂", Tensor("??A", [gs => lower]), x, ρ => lower)) ==
          PartialDerivative("∂", Pco, x, ρ => lower) * Cup
    # ?A takes only the name, so the result carries no coordinate dependence
    @test replace_pattern(TensorExpression(Pco * Cup), "?A_{?g}",
                          PartialDerivative("∂", Tensor("?A", [gs => lower]), x, ρ => lower)) ==
          PartialDerivative("∂", Tensor("P", [μ => lower]), x, ρ => lower) * Cup
    # ... at any depth
    @test replace_pattern(TensorExpression(Pco * Cup), "?A_{?g}",
              PartialDerivative("∂", PartialDerivative("∂", Tensor("??A", [gs => lower]), x, σ => lower),
                                x, ρ => lower)) ==
          PartialDerivative("∂", PartialDerivative("∂", Pco, x, σ => lower), x, ρ => lower) * Cup

    # an unbound wildcard name is reported rather than printed into the answer
    @test_throws InvalidReplacementPattern replace_pattern(
        TensorExpression(Pco), "?B_{?g}", Alakazam.parse_pattern2tensor("?A_{?g}"))
    @test_throws InvalidReplacementPattern replace_pattern(
        TensorExpression(Pco), "?B_{?g}", Alakazam.parse_pattern2tensor("??A_{?g}"))
    @test_throws InvalidReplacementPattern replace_pattern(
        TensorExpression(Pco), "P_{?g}", Alakazam.parse_pattern2tensor("?A_{?g}"))

    ###############################################
    ### ??A resolution details
    ###############################################
    # the adopted data wins over anything written on the ??A itself
    @test labels(firstt(replace_pattern(TensorExpression(Psrc), "?A_{?g?h}",
                                        TensorExpression(T"??A_?g?h_'own'")))) == ["lab" => lower]
    # whereas ?A keeps what was written, since it takes only the name
    @test labels(firstt(replace_pattern(TensorExpression(Psrc), "?A_{?g?h}",
                                        TensorExpression(T"?A_?g?h_'own'")))) == ["own" => lower]
    # built-in tensor types are adopted with their type intact
    @test typeof(firstt(replace_pattern(TensorExpression(Metric([μ => lower, ν => lower])), "?A_{?g?h}",
                        Alakazam.parse_pattern2tensor("??A_{?g?h}")))) == Metric
    # a replacement may use the same wildcard more than once
    @test replace_pattern(TensorExpression(Pco * Cup), "?A_{?g}",
                          Alakazam.parse_pattern2tensor("?A_{?g}?A^{?g}")) ==
          Tensor("P", [μ => lower]) * Tensor("P", [μ => upper]) * Cup
    # ??A -> ??A is the identity and must not loop under repeat
    @test firstt(replace_pattern(TensorExpression(Psrc), "?A_{?g?h}", TensorExpression(T"??A_?g?h"))) == Psrc


    ###############################################
    ### every combination of argument forms is accepted
    ### Expression, pattern and replacement may each be a bare tensor, a term, an expression, or
    ### (pattern/replacement) a pattern string. Object patterns used to work in one shape only.
    ###############################################
    gmatch = DummyPatternIndex("g", spacetime_indices)
    let exps = (Tensor("A", [μ => lower]),
                TensorTerm(Tensor("A", [μ => lower])),
                TensorExpression(Tensor("A", [μ => lower]))),
        pats = ("A_{?g}", TensorTerm(Tensor("A", [gmatch => lower])), Tensor("A", [gmatch => lower])),
        reps = (Tensor("Z"), TensorTerm(Tensor("Z")), TensorExpression(Tensor("Z")), "Z")

        @test all(replace_pattern(e, p, r) == TensorExpression(Tensor("Z"))
                  for e in exps, p in pats, r in reps)
    end

    ###############################################
    ### the call leaves its arguments alone
    ###############################################
    src_e = Tensor("A", [μ => lower]) * Tensor("C", [ν => upper])
    src_before = string(src_e)
    reusable = Tensor("B", [gmatch => lower])
    reusable_before = string(reusable)
    @test replace_pattern(src_e, "A_{?g}", reusable) == Tensor("B", [μ => lower]) * Tensor("C", [ν => upper])
    @test string(src_e) == src_before                      # the expression is not mutated
    @test string(reusable) == reusable_before              # nor is the replacement
    # so the same replacement object may be reused, as a rule table would
    @test replace_pattern(Tensor("A", [ν => lower]), "A_{?g}", reusable) ==
          TensorExpression(Tensor("B", [ν => lower]))

    ###############################################
    ### repeat terminates and reaches every occurrence
    ###############################################
    let idxs = [Index("rp$(k)", spacetime_indices) for k in 1:40]
        many = prod(Tensor("A", [i => lower]) for i in idxs)
        done = replace_pattern(many, "A_{?g}", Tensor("Z"))
        @test length(terms(first_summand(done))) == 40
        @test all(Alakazam.name(t) == "Z" for t in terms(first_summand(done)))
    end

    ###############################################
    ### solve_linear
    ###############################################
    let y = Coordinate("slv_y", spacetime_indices),
        p = Index("slva", spacetime_indices), q = Index("slvb", spacetime_indices),
        A = Tensor("A", [p => lower], y), B = Tensor("B", [p => lower], y),
        C = Tensor("C", [p => lower, q => lower], y), E = Tensor("E", [q => lower], y),
        ∂ = PD("∂", [y], [q => lower])

        @test solve_linear(2 * A + B, A) == (A, -1//2 * B)
        @test solve_linear(TensorExpression(A - B), A) == (A, TensorExpression(B))
        # coefficients of the target are gathered from wherever they occur
        @test solve_linear(A + A + B, A) == (A, -1//2 * B)
        # a complex coefficient stays exact rather than falling to floating point
        @test solve_linear(im//2 * A + B, A) == (A, Alakazam.simplify(2im * B))
        # an operator is as good a target as a tensor, and the two are distinct
        @test solve_linear(∂(A) - C, ∂(A)) == (∂(A), TensorExpression(C))

        # failure is reported, not raised, so a caller can sweep a list of candidates. The
        # equation comes back untouched alongside the `nothing`.
        @test solve_linear(∂(A) - C, A)[1] === nothing              # only ∂⟦A⟧ occurs, not A
        let (t, back) = solve_linear(TensorExpression(B), A)                        # absent
            @test t === nothing && back == TensorExpression(B)
        end
        @test solve_linear(TensorExpression(A * E), A)[1] === nothing # nonlinear
        @test solve_linear(TensorExpression(A * E), A)[2] == TensorExpression(A * E)

        # which is what makes a search over candidates work: the first that comes out wins
        eqn = 2 * A + B
        @test first(r for r in (solve_linear(eqn, t) for t in (E, A, B)) if r[1] !== nothing) ==
              (A, -1//2 * B)
    end

    ###############################################
    ### rewriting inside an operator's operand
    ###############################################
    # A rule has to reach into an operand wherever the operator sits.
    let w = Coordinate("opd_w", spacetime_indices),
        p = Index("opda", spacetime_indices), q = Index("opdb", spacetime_indices),
        ∂ = PD("∂", [w], [q => lower]),
        A = Tensor("A", [p => lower], w), Z = Tensor("Z", [p => lower], w),
        Q = Tensor("Q", w)

        @test replace_pattern(TensorExpression(∂(A)), TensorTerm(A), Z) == TensorExpression(∂(Z))
        # an operator sharing its term with another factor
        @test replace_pattern(TensorExpression(Q * ∂(A)), TensorTerm(A), Z) == Q * ∂(Z)
        #  nested chains
        @test replace_pattern(TensorExpression(Q * ∂(∂(A))), TensorTerm(A), Z) == Q * ∂(∂(Z))
        # a replacement introducing its own operators must not be descended into forever
        @test replace_pattern(TensorExpression(Q * A), TensorTerm(A), ∂(Z)) == Q * ∂(Z)
        # labels discriminate inside an operand
        let L1 = Tensor("A", [p => lower], w, ["1" => lower]),
            L2 = Tensor("A", [p => lower], w, ["2" => lower])

            @test replace_pattern(TensorExpression(Q * ∂(L1)), TensorTerm(L1), Z) == Q * ∂(Z)
            @test replace_pattern(TensorExpression(Q * ∂(L2)), TensorTerm(L1), Z) == Q * ∂(L2)
            # pattern naming no label is still a wildcard over them
            @test replace_pattern(TensorExpression(Q * ∂(L2)), TensorTerm(A), Z) == Q * ∂(Z)
        end
    end

    ###############################################
    ### labels are a discriminator everywhere, not only at the top level
    ###############################################
    # two components of one superfield share a name, a coordinate and an index structure, and are told apart by their label alone. 
    let y = Coordinate("lbl_y", spacetime_indices),
        ∂ = PD("∂", [y], [Index("lblm", spacetime_indices) => lower]),
        L10 = Tensor("Q", [Index("lbla", spacetime_indices) => lower], y, ["(1,0)" => lower]),
        L01 = Tensor("Q", [Index("lbla", spacetime_indices) => lower], y, ["(0,1)" => lower]),
        Z = Tensor("Z", [Index("lbla", spacetime_indices) => lower], y)
        @test labels(L10) != labels(L01)
        @test L10 != L01
        @test TensorTerm(L10) != TensorTerm(L01)
        # expression equality runs through `same_tensor_mod_syms`, which must also see the label
        @test TensorExpression(L10) != TensorExpression(L01)
        # `==` is only usable as a Dict key if the hash agrees with it
        @test hash(L10) != hash(L01)
        # sum of two components must not collapse to one
        @test length(summands(Alakazam.simplify(TensorExpression(L10) + TensorExpression(L01)))) == 2
        # an operator carries none of its operand's labels, so equality has to look through it
        @test ∂(L10) != ∂(L01)
        # same for pattern matching, or a rule for one component fires on the other
        @test replace_pattern(TensorExpression(∂(L01)), TensorTerm(∂(L10)), Z) ==
              TensorExpression(∂(L01))
        @test replace_pattern(TensorExpression(∂(L10)), TensorTerm(∂(L10)), Z) ==
              TensorExpression(Z)
        # a pattern naming no labels stays a wildcard over them, as it was before labels existed
        let U = Tensor("Q", [Index("lbla", spacetime_indices) => lower], y)
            @test replace_pattern(TensorExpression(∂(L01)), TensorTerm(∂(U)), Z) ==
                  TensorExpression(Z)
        end

        # operand's name needs the same treatment. The histogram classes factors by `name`,
        # which for ∂_μ⟦R⟧ is the operator's, so without a further check every derivative sits in
        # one class and a rule written for one field fires on any other of the same shape.
        let R = Tensor("R", [Index("lbla", spacetime_indices) => lower], y),
            r = Tensor("r", [Index("lbla", spacetime_indices) => lower], y)

            @test replace_pattern(TensorExpression(∂(r)), TensorTerm(∂(R)), Z) ==
                  TensorExpression(∂(r))
            @test replace_pattern(TensorExpression(∂(R)), TensorTerm(∂(R)), Z) == TensorExpression(Z)
            # nested chains too
            @test replace_pattern(TensorExpression(∂(∂(r))), TensorTerm(∂(∂(R))), Z) ==
                  TensorExpression(∂(∂(r)))
            # but a wildcard operand stays a wildcard: ∂⟦?A⟧ still takes any derivative.
            let wild = Tensor("?A", [Index("lbla", spacetime_indices) => lower], y)
                @test replace_pattern(TensorExpression(∂(r)), TensorTerm(∂(wild)), Z) ==
                      TensorExpression(Z)
                @test replace_pattern(TensorExpression(∂(R)), TensorTerm(∂(wild)), Z) ==
                      TensorExpression(Z)
            end
        end
    end

end
