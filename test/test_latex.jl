@testset "LT strings: reading TeX back into tensors" begin

    clear_global!()
    st  = IndexSet("spacetime", 4, (0, 4), lower)
    sp  = IndexSet("spinor", 2, lower, EpsilonTensor, fermionic)
    spb = IndexSet("spinorbar", 2, lower, EpsilonTensor, fermionic)
    sp.conjugate = spb
    spb.conjugate = sp
    a, b, m, n = Index(["a", "b", "μ", "ν"], st)
    al    = Index("α", sp)
    aldot = Index(string(Alakazam.COMBINING_DOT) * "α", spb)
    albar = Index(string(Alakazam.COMBINING_MACRON) * "α", spb)
    x = Coordinate("x", [st])

    @testset "the printed form round-trips" begin
        # parse_latex resolves against the registry rather than a caller's scope, which is what
        # makes this possible at all: nothing here has to be in scope under the same variable name.
        for t in [Tensor("A", m => lower),
                  Tensor("A", m => upper, n => lower),
                  Tensor("ψ", al => lower),
                  Tensor("ψ", aldot => upper),
                  Tensor("ψ", albar => lower),
                  Tensor("k", "1" => lower, m => upper),
                  Tensor("F", m => lower, n => lower),
                  Tensor("A", a => lower, b => upper),
                  Tensor("A", m => lower) * Tensor("B", m => upper)]
            @test parse_latex(to_latex(t)) == t
        end

        # the decoration on a NAME round-trips too, not only on an index
        psibar = Tensor(string(Alakazam.COMBINING_MACRON) * "ψ", aldot => lower)
        @test parse_latex(to_latex(psibar)) == psibar

        # and a whole expression, coefficient included
        E = TensorExpression(Tensor("A", m => lower)) + TensorExpression(Tensor("B", m => lower))
        @test parse_latex(to_latex(E)) == E
        H = 1 // 2 * TensorExpression(Tensor("A", m => lower))
        @test parse_latex(to_latex(H)) == H
    end

    @testset "what a person would actually type" begin
        # one script carrying several indices, which the printed form never produces
        @test LT"F_{\mu\nu}" == Tensor("F", m => lower, n => lower)
        @test LT"F_{ab}" == Tensor("F", a => lower, b => lower)
        # an unbraced script
        @test LT"A^\mu" == Tensor("A", m => upper)
        @test LT"A^{\mu}_{\nu}" == Tensor("A", m => upper, n => lower)
        # juxtaposition is a product, and the repeated index contracts
        @test LT"A^\mu B_\mu" == Tensor("A", m => upper) * Tensor("B", m => lower)
        @test length(Alakazam.free_indices_ordered(LT"A_\mu B^\mu")) == 0
    end

    @testset "accents, braced or not, name one index" begin
        @test LT"\psi_{\dot\alpha}" == Tensor("ψ", aldot => lower)
        @test LT"\psi_{\dot{\alpha}}" == Tensor("ψ", aldot => lower)
        @test LT"\psi_{\bar\alpha}" == Tensor("ψ", albar => lower)
        # \overline and \bar are the same mark; the accent may sit on the name as well
        psibar = Tensor(string(Alakazam.COMBINING_MACRON) * "ψ", aldot => lower)
        @test LT"\bar\psi_{\dot\alpha}" == psibar
        @test LT"\overline{\psi}_{\dot\alpha}" == psibar
    end

    @testset "labels, wildcards and coordinates" begin
        @test LT"k_1^\mu" == Tensor("k", "1" => lower, m => upper)
        @test LT"A_{?x}" == T"A_?x"
        @test LT"A_\mu(x)" == Tensor("A", m => lower, x)
        # a name that is neither in scope nor registered is a label, the rule T"..." follows
        @test labels(LT"A_q") == labels(T"A_q")
    end

    @testset "decoration that carries no meaning is dropped" begin
        @test LT"A^\mu \, B_\mu" == Tensor("A", m => upper) * Tensor("B", m => lower)
        @test LT"\left( A^\mu \right)" == Tensor("A", m => upper)
        @test LT"$A^\mu$" == Tensor("A", m => upper)
        @test LT"A^{\mu} \quad B_{\mu}" == Tensor("A", m => upper) * Tensor("B", m => lower)
    end

    @testset "sums and coefficients" begin
        A = TensorExpression(Tensor("A", m => lower))
        B = TensorExpression(Tensor("B", m => lower))
        @test LT"A_\mu + B_\mu" == A + B
        @test LT"A_\mu - B_\mu" == A - B
        @test LT"\frac{1}{2} A_\mu" == 1 // 2 * A
        @test LT"2 A_\mu" == 2 * A
    end

    @testset "the symbol table really is the inverse of the printer" begin
        # It is built by reading Latexify's own table backwards rather than being maintained here,
        # so the property to assert is that uni_to_TeX and the reader undo one another, over the
        # whole range rather than for a handful of letters.
        @test length(Alakazam._TEX_TO_UNI) > 500

        chars = ['α', 'β', 'γ', 'δ', 'ϵ', 'ζ', 'η', 'θ', 'ι', 'κ', 'λ', 'μ', 'ν', 'ξ',
                 'π', 'ρ', 'σ', 'τ', 'υ', 'ϕ', 'χ', 'ψ', 'ω',
                 'Γ', 'Δ', 'Θ', 'Λ', 'Ξ', 'Π', 'Σ', 'Υ', 'Φ', 'Ψ',
                 '∂', '∇', 'ℓ', '∞', 'a', 'Z']
        for ch in chars
            @test Alakazam._tex_name(Alakazam.uni_to_TeX(string(ch))) == string(ch)
        end

        # the same through whole tensors, decorated and not. These go in a set of their own: a name
        # is unique only within a set, and several of the letters above are already indices of the
        # spinor sets declared at the top of this file.
        wide = IndexSet("widecharset", 4, (0, 4), lower)
        for ch in chars
            for mark in ("", string(Alakazam.COMBINING_DOT), string(Alakazam.COMBINING_MACRON))
                nm = mark * string(ch)
                # skip anything already claimed elsewhere, which the markup cannot disambiguate
                Alakazam._find_registered_index(nm) === nothing || continue
                t = Tensor("A", Index(nm, wide) => lower)
                @test parse_latex(to_latex(t)) == t
            end
        end
    end

    @testset "a name living in two sets cannot be resolved from the markup" begin
        # to_latex records the name and not the set it belongs to, so this is a real limit rather
        # than an oversight. The answer is at least deterministic: sets are searched by name.
        other = IndexSet("aaa_first_by_name", 4, (0, 4), lower)
        shared = Index("μ", other)
        @test Alakazam._find_registered_index("μ") === shared      # "aaa..." sorts before "spacetime"
        @test Alakazam._find_registered_index("μ") ===
              Alakazam._find_registered_index("μ")                 # and stable across calls
        # with the index bound in scope, LT prefers the binding and the ambiguity does not arise
        @test LT"A_\mu" == Tensor("A", m => lower)
        # drop just this set again; clearing the whole registry would strand the later testsets,
        # which resolve μ and ν through it rather than through a binding
        delete!(Alakazam.SET_REGISTRY, "aaa_first_by_name")
        @test Alakazam._find_registered_index("μ") === m
    end

    @testset "symmetric and antisymmetric notation" begin
        # the brackets of the literature build the same tableaux as the tensor-string syntax
        @test tableaux(LT"A_{[ab]}") == tableaux(T"A_[ab]")
        @test tableaux(LT"A_{(ab)}") == tableaux(T"A_(ab)")
        # braces are optional, so the tensor-string spelling is accepted as well
        @test tableaux(LT"A_[ab]") == tableaux(T"A_[ab]")
        @test tableaux(LT"A_(ab)") == tableaux(T"A_(ab)")
        @test tableaux(LT"F_{[\mu\nu]}") == tableaux(T"F_[ab]")

        # and they are real symmetries, not decoration
        @test isZeroTens(simplify(LT"A_{[ab]}" + LT"A_{[ba]}"))
        @test isZeroTens(simplify(LT"A_{(ab)}" - LT"A_{(ba)}"))
        @test !isZeroTens(simplify(LT"A_{ab}" + LT"A_{ba}"))

        # bars hold slots out of the group, as in R_{a[bc|d|e]}
        @test tableaux(LT"A_{[a|b|c]}") == tableaux(T"A_[a|b|c]")
        @test tableaux(LT"A_{(a|b|c)}") == tableaux(T"A_(a|b|c)")
        @test isZeroTens(simplify(LT"A_{[a|b|c]}" + LT"A_{[c|b|a]}"))

        # a group may cover part of the slots, or there may be several
        @test tableaux(LT"R_{a[bc]}") == tableaux(T"R_a[bc]")
        @test tableaux(LT"R_{[ab][cd]}") == tableaux(T"R_[ab][cd]")
        # slots are numbered across the whole tensor, not within one script
        @test tableaux(LT"A^{[ab]}") == tableaux(T"A^[ab]")
        @test tableaux(LT"A^{a}_{[bc]}") == tableaux(T"A^a_[bc]")

        # malformed grouping is reported rather than guessed at
        @test_throws ErrorException parse_latex("A_{[ab}")
        @test_throws ErrorException parse_latex("A_{ab]}")
        @test_throws ErrorException parse_latex("A_{[ab)}")
        # antisymmetrisation across two factors is not a symmetry of either, and says so
        @test_throws ErrorException parse_latex("\\partial_{[a}A_{b]}")

        # NOTE: to_latex does not write the brackets, so symmetry is the one thing the round trip
        # cannot recover. This pins that, so it is noticed if to_latex ever starts emitting them.
        F = AntisymmetricTensor("F", a => lower, b => lower)
        @test to_latex(F) == to_latex(Tensor("F", a => lower, b => lower))
        @test tableaux(parse_latex(to_latex(F))) != tableaux(F)
    end

    @testset "a symbolic coefficient is indistinguishable from an indexless tensor" begin
        # to_latex prints c * A and Tensor("c") * A the same way up to a space, so the printed form
        # simply does not carry the difference. LT resolves every named factor as a tensor, which
        # is a choice rather than a recovery, and this pins it.
        @syms c
        A = Tensor("A", m => lower)
        @test parse_latex(to_latex(c * TensorExpression(A))) == Tensor("c") * A
        # the two only ever differed by whitespace
        @test replace(to_latex(c * TensorExpression(A)), " " => "") ==
              replace(to_latex(Tensor("c") * A), " " => "")
        # comparing a symbolically weighted expression against a term answers false, not raises
        @test (c * TensorExpression(A) == TensorTerm(A)) == false
    end

    @testset "run-time strings" begin
        # the same reader, without a caller's scope to consult
        @test parse_latex("F_{\\mu\\nu}") == Tensor("F", m => lower, n => lower)
        @test_throws ErrorException parse_latex("A_\\mu(nosuchcoord)")
        # a + with no factor after it contributes no term, and is dropped rather than complained of
        @test parse_latex("A_\\mu +") == Tensor("A", m => lower)
        # but something that cannot begin a name is reported rather than skipped
        @test_throws ErrorException parse_latex("A_\\mu ]")
    end

end
