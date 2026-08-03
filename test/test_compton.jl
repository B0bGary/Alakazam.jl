#=
Compton scattering in massless scalar QED -- the Alakazam equivalent of the introductory example
of the Redberry paper (arXiv:1302.1219, Sec. 1.2).   γ(k₁)φ(k₃) → γ(k₂)φ(k₄)
=#

@testset "Compton scattering, massless scalar QED (Redberry Sec. 1.2)" begin
    st = IndexSet("compton_spacetime", 4, upper); @indices i:j ∈ st
    Symbolics.@syms e s12 s13 s23

    k(n, x, p) = Tensor(n, [x => p])
    V(a, b) = -im * e * (a + b)                       # V_i[p,q] = -ie(p+q)_i

    # (k₃+k₁)² = 2·s13 and (k₃-k₂)² = -2·s23 on shell
    function M(x, y, p)                               # x,y: the two photon indices
        m(n, z) = k(n, z, p)
        V(m("r",x), m("r",x)+m("p",x)) * (-im/(2*s13))  * V(-m("t",y), -m("r",y)-m("p",y)) +
        V(m("r",y), m("r",y)-m("q",y)) * (-im/(-2*s23)) * V(-m("r",x)+m("q",x), -m("t",x)) +
        2*im*e^2*Metric([x => p, y => p])              # seagull V^ij = 2ie²g^ij
    end

    w = DummyPatternIndex("w", st)                     # momentum conservation k₄ = k₁ + k₃ - k₂
    con(X, p) = replace_pattern(X, p === upper ? "t^{?w}" : "t_{?w}",
                                k("p",w,p) + k("r",w,p) - k("q",w,p))
    Mu, Ml = con(M(i,j,upper), upper), con(M(i,j,lower), lower)

    @test length(summands(M(i,j,upper))) == 13         # 6 + 6 + seagull
    @test Mu != M(i,j,upper)                           # the substitution fired

    # k_i·k_i = 0, and each remaining scalar product becomes its invariant
    RULES = vcat([("$(n)_{?x}$(n)^{?x}", zero()) for n in ("p","q","r")],
                 [("$(n)^{?x}$(n)_{?x}", zero()) for n in ("p","q","r")],
                 [(pat, TensorExpression(s))
                  for (a,b,s) in (("p","q",s12), ("p","r",s13), ("q","r",s23))
                  for pat in ("$(a)_{?x}$(b)^{?x}", "$(a)^{?x}$(b)_{?x}")])

    function to_invariants(x)
        for _ in 1:8
            prev = x
            for (pat, rep) in RULES
                x = replace_pattern(x, pat, rep)
            end
            x == prev && break
        end
        isZeroTens(x) && return 0
        is_const(x) || return nothing                  # something failed to reduce
        return first_coef(x)
    end

    # s + t + u = 0 for four massless legs  =>  k₁·k₂ = k₁·k₃ - k₂·k₃
    onshell(v) = Symbolics.substitute(v, Dict(s12 => s13 - s23))
    iszerosym(v) = Alakazam.coeff_equal(Symbolics.simplify(Symbolics.expand(v)), 0)

    squared = to_invariants(Alakazam.simplify(-(Ml * Mu)))
    @test squared !== nothing                          # fully reduced to a scalar
    @test iszerosym(onshell(squared) - 8 * e^4)        # |M|² = 8e⁴

    # gauge invariance: k₁_i M^ij = 0 and M^ij k₂_j = 0. The sharpest check on the amplitude --
    # each diagram violates it alone, and only the seagull's normalisation restores it.
    ward = true
    for n in ("p", "q", "r")                           # p,q,r span the momenta appearing in M
        w1 = to_invariants(Alakazam.simplify(k("p",i,lower) * Mu * k(n,j,lower)))
        w2 = to_invariants(Alakazam.simplify(k(n,i,lower) * Mu * k("q",j,lower)))
        ward &= w1 !== nothing && w2 !== nothing &&
                iszerosym(onshell(w1)) && iszerosym(onshell(w2))
    end
    @test ward

    # the expression Redberry prints, in its own variables:
    #   A = (k₂k₃) = s23,  B = (k₂k₄) = k₁·k₃ = s13,  C = (k₃k₄) = k₁·k₂ = s12
    redberry = -(e^4/2) * (-18*s23*s13 + s12^2 + s23^2 + 2*(s23-s13)*s12 + s13^2) / (s23*s13)
    @test iszerosym(onshell(redberry) - 8 * e^4)
    @test iszerosym(onshell(squared) - onshell(redberry))
end
