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
                #a name of its own: a global `t` exists, and assigning to it from inside
                #this loop is a soft-scope ambiguity Julia warns about on every build
                tens = Tensor("A", Index(nm, wide) => lower)
                @test parse_latex(to_latex(tens)) == tens
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

    @testset "writing: coefficients, special objects, operators" begin
        T = Tensor("T")
        tex(x) = strip(to_latex(x))
        @syms q         # one letter: latexify sets a longer symbol name in \mathtt

        # A coefficient is written from its value. The display string attaches factors by
        # juxtaposition, `2//3i`, and read back as Julia that binds as 2//(3i).
        @test tex(2//3 * im * T) == "\\frac{2i}{3} T"
        @test tex(-2//3 * im * T) == "-\\frac{2i}{3} T"
        @test tex(2//3 * T) == "\\frac{2}{3} T"
        @test tex(im * T) == "i T"
        # the same flaw for a symbol, which would have come out as 1/(2m)
        @test tex(1//2 * q * T) == "\\frac{q}{2} T"
        @test tex((1//2 + 2//3 * im) * T) == "\\left(\\frac{1}{2} + \\frac{2i}{3}\\right) T"
        # and the display string itself is left as it was
        @test string(2//3 * im * T) == "2//3iT"

        # γ5 carries its 5 as a label, so only the letter is its name
        @test tex(Gamma5()) == "\\gamma { }^{5}"
        # a word, not the product T·r
        @test startswith(tex(Trace(Gamma(m => lower) * Gamma(n => lower))), "\\mathrm{Tr}")

        # An operator keeps the order of its raised and lowered slots, as a tensor does. They used
        # to be collected into one group each, and these two different operators printed the same.
        d1 = Derivative("D", Tensor("B", x), [x], [m => upper, n => lower])
        d2 = Derivative("D", Tensor("B", x), [x], [n => lower, m => upper])
        @test startswith(tex(d1), "D{ }^{\\mu}{ }_{\\nu}")
        @test startswith(tex(d2), "D{ }_{\\nu}{ }^{\\mu}")
        @test tex(d1) != tex(d2)

        # a unit inside a symbolic sum is a term, not a dropped coefficient: 1+q printed as ( + q)
        @test string((1 + q) * T) == "(1 + q)T"
        @test string((q - 1) * T) == "(-1 + q)T"
        @test tex((1 + q) * T) == "\\left(1 + q\\right) T"
        @test isequal(coefficients(3 * Tensor("A", m => lower) + (1 + q) * Tensor("B", m => lower)), [3, 1 + q])
        @test coefficients(Tensor("A", m => lower) * Tensor("B", m => upper)) == [1]
    end

    @testset "coefficients, console and TeX, case by case" begin
        @variables q r      # Num, which is itself a Number, unlike the bare symbols of @syms
        @syms s
        E = TensorExpression(Tensor("T"))
        #coefficient => (console, TeX) on one tensor; `nothing` where the console form is not pinned
        cases = Pair{Any,Tuple{Union{String,Nothing},String}}[
            # real numbers
            -1 => ("-", "- "), 4//2 => ("2", "2 "), -3//4 => ("-3//4", "-\\frac{3}{4} "),
            -2.5 => ("-2.5", "-2.5 "),
            # complex numbers: i throughout, a negative real part's sign outside its fraction
            1 + im => ("(1 + i)", "\\left(1 + i\\right) "), -2 - 3im => ("(-2 - 3i)", "\\left(-2 - 3i\\right) "),
            -im // 3 => ("-1//3i", "-\\frac{i}{3} "), 1.5 + 2.5im => ("(1.5 + 2.5i)", "\\left(1.5 + 2.5i\\right) "),
            -1//2 - im => ("(-1//2 - i)", "\\left(-\\frac{1}{2} - i\\right) "),
            1//2 - 2//3 * im => ("(1//2 - 2//3i)", "\\left(\\frac{1}{2} - \\frac{2i}{3}\\right) "),
            # a number times a symbol
            -(2//3) * q => ("-2//3q", "-\\frac{2q}{3} "), q^2 => (nothing, "q^{2} "),
            # symbolic complex parts, which Symbolics holds as complex(re, im)
            im * q => ("iq", "i q "), -im * q / 2 => ("-1//2iq", "-\\frac{i q}{2} "),
            (1 + im) * q => ("(q + iq)", "\\left(q + i q\\right) "), q - im => ("(q - i)", "\\left(q - i\\right) "),
            -1//2 + im * q => ("(-1//2 + iq)", "\\left(-\\frac{1}{2} + i q\\right) "),
            2q - 3im * r => ("(2q - 3ir)", "\\left(2q - 3i r\\right) "),
            im * (1 + q) => ("i(1 + q)", "i \\left(1 + q\\right) "),
            im * (q + r) / 2 => ("1//2i(q + r)", "\\frac{i \\left(q + r\\right)}{2} "),
            # a number times a sum is bracketed, but a sum alone over a denominator needs none
            2 * (q + r) => ("2(q + r)", "2\\left(q + r\\right) "), (1//2) * (q + r) => ("1//2(q + r)", "\\frac{q + r}{2} "),
            -(q + r) / 3 => ("-1//3(q + r)", "-\\frac{q + r}{3} "),
            # quotients: bracketed whole on the console, or `1 / qT` reads as 1/(qT); the sign in front
            1 / q => ("(1 / q)", "\\frac{1}{q} "), -q / r => ("-(q / r)", "-\\frac{q}{r} "),
            1 / (2q) => ("(1 / (2q))", "\\frac{1}{2q} "), (1 + q) / r => ("((1 + q) / r)", "\\frac{1 + q}{r} "),
            q / (1 + q) => ("(q / (1 + q))", "\\frac{q}{1 + q} "),
            # sums, a unit or a fraction among their terms
            -1 - q => ("(-1 - q)", "\\left(-1 - q\\right) "), q - 1//2 => ("(-1//2 + q)", "\\left(-\\frac{1}{2} + q\\right) "),
            # bare symbols of @syms take the same paths
            im * s => ("is", "i s "), (1 + im) * s => ("(1 + i)s", "\\left(1 + i\\right) s "), s - 1 => ("(-1 + s)", "\\left(-1 + s\\right) "),
        ]
        for (c, (con, tex)) in cases
            con === nothing || @test string(c * E) == con * "T"
            @test strip(to_latex(c * E)) == tex * "T"
        end

        # a later summand's sign becomes the operator before it
        A, B = TensorExpression(Tensor("A", m => lower)), TensorExpression(Tensor("B", m => lower))
        nospace(x) = replace(to_latex(x), " " => "")
        @test string(A - (q / r) * B) == "A_{μ}-(q / r)B_{μ}"
        @test nospace(A - (q / r) * B) == "A{}_{\\mu}-\\frac{q}{r}B{}_{\\mu}"
        @test string(A - im * q * B) == "A_{μ}-iqB_{μ}"
        @test nospace(A - im * q * B) == "A{}_{\\mu}-iqB{}_{\\mu}"
        @test string(A + (1 - im) * B) == "A_{μ}+(1 - i)B_{μ}"
        @test string(A - (1 + q) * B) == "A_{μ}+(-1 - q)B_{μ}"

        # an expression with no tensor in it writes even a unit coefficient
        one_ = TensorExpression(1)
        @test string(-1 * one_) == "-1" && strip(to_latex(-1 * one_)) == "-1"
        @test string((1//2) * one_) == "1//2" && strip(to_latex((1//2) * one_)) == "\\frac{1}{2}"
        @test string(im * q * one_) == "iq" && strip(to_latex(im * q * one_)) == "i q"
    end

    @testset "every decoration, next to another" begin
        # stored mark-first; printed with each mark behind its own letter, and as its TeX accent
        hs = IndexSet("latex_hatset", 4, (0, 4))
        ahat, btilde = Index(string(Alakazam.COMBINING_HAT) * "a", hs), Index(string(Alakazam.COMBINING_TILDE) * "b", hs)
        abar = Index(string(Alakazam.COMBINING_MACRON) * "a", spb)
        T2 = Tensor("T", abar => upper, ahat => upper, btilde => lower)
        @test string(T2) == "T^{āâ∘}_{∘∘b̃}"
        @test strip(to_latex(T2)) == "T { }^{\\overline{a}}{ }^{\\hat{a}}{ }_{\\tilde{b}}"
        @test parse_latex(to_latex(T2)) == T2
        # stacked marks nest, the first-written outermost
        @test strip(to_latex(Tensor(string(Alakazam.COMBINING_HAT, Alakazam.COMBINING_MACRON, "A")))) == "\\hat{\\overline{A}}"
    end

end
