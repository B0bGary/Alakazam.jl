#=
Compton scattering in massless scalar QED: the Alakazam equivalent of the introductory example
of the Redberry paper (arXiv:1302.1219, Sec. 1.2).   γ(k₁)φ(k₃) → γ(k₂)φ(k₄)
=#

@testset "Compton scattering, massless scalar QED (Redberry Sec. 1.2)" begin
    st = IndexSet("compton_spacetime", 4, upper); @indices i:j ∈ st
    Symbolics.@syms e s12 s13 s23                  # sᵢⱼ = kᵢ·kⱼ, scalar products

    # kₙ is one tensor k carrying n as a label, so the code reads as the physics does
    k(n, x, p=upper) = Tensor("k", [x => p], "$n" => lower)
    V(a, b) = -im * e * (a + b)                    # V^i(p,p') = -ie(p+p')^i

    # (k₃+k₁)² = 2·s13 and (k₃-k₂)² = -2·s23 on shell
    M = V(T"k_3^i", T"k_3^i"+T"k_1^i") * (-im/( 2*s13)) * V(-T"k_4^j", -T"k_3^j"-T"k_1^j") +
        V(T"k_3^j", T"k_3^j"-T"k_2^j") * (-im/(-2*s23)) * V(-T"k_3^i"+T"k_2^i", -T"k_4^i") +
        2*im*e^2*Metric([i => upper, j => upper])  # seagull V^ij = 2ie²g^ij

    Mu = replace_pattern(M, T"k_4^?w", T"k_1^?w" + T"k_3^?w" - T"k_2^?w")   # k₄ = k₁+k₃-k₂
    Ml = lower_all(Mu)                             # M has no dummies, so this is M with i,j down

    @test length(summands(M)) == 13                # 6 + 6 + seagull
    @test Mu != M                                  # the substitution fired
    @test T"k_1^i" == k(1, i) && T"k_1^i" != T"k_2^i"     # the macro and the constructor agree

    # kᵢ·kᵢ = 0, and each remaining scalar product becomes its invariant. Patterns are built as
    # objects rather than strings, so the momenta keep their physical names: a pattern string
    # parses a tensor name as a single character, but "k1" is fine as an object pattern.
    kp(n, p) = k(n, DummyPatternIndex("x", st), p)
    prods = ((1,1,0), (2,2,0), (3,3,0), (1,2,s12), (1,3,s13), (2,3,s23))   # kᵢ² = 0,  kᵢ·kⱼ = sᵢⱼ
    RULES = [(kp(a,p)*kp(b,op_pos(p)), r) for (a,b,r) in prods for p in (lower,upper)]
    # both index orders are needed: a pattern fixes which slot is raised, so k₁_x k₂^x and
    # k₁^x k₂_x are different patterns even though they are the same scalar product

    function to_invariants(x)
        for _ in 1:8
            prev = x
            for (pat, rep) in RULES
                x = replace_pattern(x, pat, rep)
            end
            x == prev && break
        end
        isZeroTens(x) && return 0
        is_const(x) || return nothing               # something failed to reduce
        return first_coef(x)
    end

    # k₄ is eliminated by momentum conservation before the on-shell rules run, so its own condition
    # k₄² = 0 survives only as the scalar relation s12 = s13 - s23, imposed here at the end.
    onshell(v) = Symbolics.substitute(v, Dict(s12 => s13 - s23))
    iszerosym(v) = Alakazam.coeff_equal(Symbolics.simplify(Symbolics.expand(v)), 0)

    # every term of M^ij is purely imaginary, so M* = -M and |M|² = -M_ij M^ij
    squared = to_invariants(Alakazam.simplify(-(Ml * Mu)))
    @test squared !== nothing                       # fully reduced to a scalar
    @test iszerosym(onshell(squared) - 8 * e^4)     # |M|² = 8e⁴, a constant

    # gauge invariance: k₁_i M^ij = 0 and M^ij k₂_j = 0. The sharpest check on the amplitude --
    # each diagram violates it alone, and only the seagull's normalisation restores it.
    ward = true
    for n in 1:3                                    # k₁,k₂,k₃ span the momenta appearing in M
        w1 = to_invariants(Alakazam.simplify(k(1,i,lower) * Mu * k(n,j,lower)))
        w2 = to_invariants(Alakazam.simplify(k(n,i,lower) * Mu * k(2,j,lower)))
        ward &= w1 !== nothing && w2 !== nothing &&
                iszerosym(onshell(w1)) && iszerosym(onshell(w2))
    end
    @test ward

    # the expression Redberry prints, in its own variables:
    #   A = (k₂k₃) = s23,  B = (k₂k₄) = k₁·k₃ = s13,  C = (k₃k₄) = k₁·k₂ = s12
    # It differs from ours term by term because k₄² = 0 has not been imposed there; the two agree
    # once s12 = s13 - s23 is substituted, i.e. they are representatives modulo that relation.
    redberry = -(e^4/2) * (-18*s23*s13 + s12^2 + s23^2 + 2*(s23-s13)*s12 + s13^2) / (s23*s13)
    @test iszerosym(onshell(redberry) - 8 * e^4)
    @test iszerosym(onshell(squared) - onshell(redberry))
end

#=
Compton scattering in spinor QED. The same process with the scalar replaced by an electron.
γ(k₁)e⁻(k₃) → γ(k₂)e⁻(k₄), summed over spins and photon polarisations, checked against Peskin (5.87).
=#
@testset "Compton scattering, spinor QED" begin
    st = IndexSet("compton_spinor_st", 4, upper)
    @indices μ:ν ∈ st
    @indices a:j ∈ st                          # one dummy per slashed momentum
    Symbolics.@syms m s u                      # s = k₁·k₃, u = k₂·k₃
    iszerosym(v) = Alakazam.coeff_equal(Symbolics.simplify(Symbolics.expand(v)), 0)

    γ(i, p=upper) = Gamma(i => p)
    slash(n, i) = Tensor("k", [i => upper], "$n" => lower) * Gamma(i => lower)
    Id = IdMatrix(spinor_indexset(st))   # the space the gammas act on, not st itself
    tr(x) = Alakazam.simplify(eliminate_trace(Trace(TensorExpression(x))))

    ###############################################
    ### the trace of a gamma run, and how it differs from the same run outside a trace
    ###############################################
    @test spinor_dimension(st) == 4 && spinor_indexset(st).range == 4
    st10 = IndexSet("compton_d10", 10, upper)
    @test spinor_dimension(st10) == 32 && spinor_indexset(st10).range == 32
    @test spinor_indexset(st) === spinor_indexset(st)           # built once, reused
    @test tr(Id) == TensorExpression(4)
    @test tr(IdMatrix(spinor_indexset(st10))) == TensorExpression(32)
    # a space that has nothing to do with spinors reports its own range
    @test tr(IdMatrix(IndexSet("compton_colour", 3, upper))) == TensorExpression(3)
    # identities over different spaces stay distinct, even at equal size
    @test IdMatrix(IndexSet("compton_cA", 3, upper)) != IdMatrix(IndexSet("compton_cB", 3, upper))
    @test IdMatrix(spinor_indexset(st)) != IdMatrix()
    @test tr(γ(μ)) == zero()                                   # odd runs vanish
    @test tr(γ(μ) * γ(ν) * γ(ρ)) == zero()
    @test tr(γ(μ, lower) * γ(ν, lower)) == 4 * Metric([μ => lower, ν => lower])
    # the (n-1)!! term counts of the recursion
    @test length(summands(tr(γ(a)*γ(b)*γ(c)*γ(d)*γ(e)*γ(f)))) == 15
    @test length(summands(tr(γ(a)*γ(b)*γ(c)*γ(d)*γ(e)*γ(f)*γ(g)*γ(h)))) == 105

    # A contracted run reduces to the dimension only inside a trace.
    @test eliminate_matrices(TensorExpression(γ(μ) * γ(μ, lower))) == TensorExpression(4)
    @test tr(γ(μ) * γ(μ, lower)) == TensorExpression(16)                       # 4d
    @test tr(γ(ν) * γ(μ) * γ(μ, lower) * γ(ν, lower)) == TensorExpression(64)  # 4d²
    # and an identity between the gammas is transparent, not an obstruction to matching
    @test tr(Id * γ(μ) * Id * γ(μ, lower)) == TensorExpression(16)

    ###############################################
    ### explicit matrix element Σ|M|² = e⁴ tr[(k̸₄+m) Γ^{μν} (k̸₃+m) Γ̄_{μν}],  Γ̄ = γ⁰Γ†γ⁰ reversing the gammas
    ###############################################
    # denominators (k₃+k₁)²-m² = 2s and (k₃-k₂)²-m² = -2u
    Γ    = γ(ν) * (slash(3, a) + slash(1, b) + m * Id) * γ(μ) * (1//2) / s -
           γ(μ) * (slash(3, c) - slash(2, d) + m * Id) * γ(ν) * (1//2) / u
    Γbar = γ(μ, lower) * (slash(3, e) + slash(1, f) + m * Id) * γ(ν, lower) * (1//2) / s -
           γ(ν, lower) * (slash(3, g) - slash(2, h) + m * Id) * γ(μ, lower) * (1//2) / u

    kernel = TensorExpression((slash(4, i) + m * Id) * Γ * (slash(3, j) + m * Id) * Γbar)
    traced = Alakazam.simplify(eliminate_trace(Trace(kernel)))

    # k₄ = k₁ + k₃ - k₂, in both index positions since the trace lowers some of them
    onshell = replace_pattern(traced,  T"k_4^?w", T"k_1^?w" + T"k_3^?w" - T"k_2^?w")
    onshell = replace_pattern(onshell, T"k_4_?w", T"k_1_?w" + T"k_3_?w" - T"k_2_?w")

    kx(n, p) = Tensor("k", [DummyPatternIndex("x", st) => p], "$n" => lower)
    prods = ((1, 1, 0), (2, 2, 0), (3, 3, m^2),      # k₁² = k₂² = 0, k₃² = m²
             (1, 3, s), (2, 3, u), (1, 2, s - u))    # k₁·k₃ = s, k₂·k₃ = u, k₁·k₂ = s - u
    rules = [(kx(x, p) * kx(y, op_pos(p)), r) for (x, y, r) in prods for p in (lower, upper)]
    for _ in 1:4, (pat, rep) in rules
        onshell = replace_pattern(onshell, pat, rep)
    end
    onshell = Alakazam.simplify(onshell)

    @test length(summands(onshell)) == 1          # every index contracted away
    M2 = first_coef(onshell)
    # Peskin (5.87) quotes the spin-averaged square, so it carries a 1/4 the sum does not
    peskin = 2 * (s/u + u/s + 2m^2 * (1/s - 1/u) + m^4 * (1/s - 1/u)^2)
    @test iszerosym(M2 / 4 - peskin)
    # the m → 0 limit is the familiar s/u + u/s
    @test iszerosym(Symbolics.substitute(M2 / 4, Dict(m => 0)) - 2 * (s/u + u/s))
end
