####################################################################################
#  Alakazam — a sample of the package
#
#      julia --project=. demo.jl
#
#  This follows the shape of the manual: fundamentals, special tensors, operations,
#  operators, pattern matching, variational calculus, differential geometry,
#  conjugation, superspace, serialisation. Each part is a taste rather than a
#  treatment. For complete calculations see the demo/ folder:
#      demo/compton_spinor.jl   Compton scattering, with trace reduction rules
#      demo/ttbar2DN02.jl       an N=(0,2) superspace constraint analysis
####################################################################################

using Alakazam, Symbolics

section(t) = println("\n", "─"^76, "\n  ", t, "\n", "─"^76)
out(l, v) = println("  ", rpad(l, 28), v)

####################################################################################
section("Fundamentals — index sets, indices, coordinates, tensors")
####################################################################################

# An index set fixes a dimension, a signature and where indices sit by default.
st = IndexSet("spacetime", 4, lower)
μ, ν, ρ, σ, τ, λ, κ, χ = Index(["μ", "ν", "ρ", "σ", "τ", "λ", "κ", "χ"], st)
x = Coordinate("x", st)                          # tensors may be fields of x

out("a tensor", Tensor("T", μ => upper, ν => lower, τ => upper))
out("a tensor field", Tensor("A", μ => lower, x))
# a repeated index is summed, and the name of a dummy is not observable
out("A^μ B_μ ≅ A^ν B_ν", Tensor("A", μ => upper) * Tensor("B", μ => lower) ≅
                          Tensor("A", ν => upper) * Tensor("B", ν => lower))

####################################################################################
section("Building tensors — macros, T-strings, LaTeX, and back out again")
####################################################################################

out("@Tensor", @Tensor T_μν^τ)
out("T-string", T"T_μν^τ")
out("from LaTeX", LT"T_{\mu\nu}^{\tau}")
out("to LaTeX", to_latex(T"T_μν^τ"))

####################################################################################
section("Special tensors")
####################################################################################

out("metric", Metric(μ => lower, ν => lower))
out("Kronecker", KroneckerDelta(μ => upper, ν => lower))
out("epsilon", EpsilonTensor(μ => lower, ν => lower, ρ => lower, σ => lower))
out("raising an index", contract_metrics(Metric(ν => upper, ρ => upper) *
                                         Tensor("T", μ => upper, ν => lower)))

####################################################################################
section("Symmetries — Young tableaux, and the Bianchi identity")
####################################################################################

# Rows symmetrise, columns antisymmetrise. This is the Riemann window.
riemann_tableau = YoungTableau([2, 2], [1, 3, 2, 4])
R(i, j, k, l) = Tensor("R", i => lower, j => lower, k => lower, l => lower,
                       riemann_tableau)
out("R_μνρσ", R(μ, ν, ρ, σ))
out("R_μνρσ = -R_νμρσ", isZeroTens(simplify(R(μ, ν, ρ, σ) + R(ν, μ, ρ, σ))))
# the first Bianchi identity is a multi-term symmetry: no single term is zero,
# only the cyclic sum, which hard_simplify finds
out("first Bianchi identity", hard_simplify(R(μ, ν, σ, τ) + R(μ, σ, τ, ν) +
                                            R(μ, τ, ν, σ); multiterm = true))

####################################################################################
section("Operations on expressions")
####################################################################################

W = AntisymmetricTensor("W", [μ => lower, ν => lower])
S = SymmetricTensor("S", [μ => upper, ν => upper])
out("symmetric × antisymmetric", hard_simplify(S * W))
out("raise every index", raise_all(Tensor("T", μ => lower, ν => lower)))
out("canonical dummy names", rename_dummies(Tensor("A", λ => upper) *
                                            Tensor("B", λ => lower)))

####################################################################################
section("Matrices — gamma, Pauli, traces, bras and kets")
####################################################################################

out("γ^μ γ_μ", contract_adjacent(Gamma(μ => upper) * Gamma(μ => lower), Gamma))
out("γ^μ γ^ν γ_μ", contract_adjacent(Gamma(μ => upper) * Gamma(ν => upper) *
                                     Gamma(μ => lower), Gamma))
out("tr(γ_μ γ_ν)", simplify(eliminate_trace(Trace(Gamma(μ => lower) *
                                                  Gamma(ν => lower)))))

st2 = IndexSet("spacetime10", 10, lower)
M,N,P,Q = Index(["M", "N", "P", "Q"], st2)
println("dimension aware:")
out("γ^M γ_M", contract_adjacent(Gamma(M => upper) * Gamma(M => lower), Gamma))
out("γ^M γ^N γ_M", contract_adjacent(Gamma(M => upper) * Gamma(N => upper) *
                                     Gamma(M => lower), Gamma))
out("tr(γ_M γ_N)", simplify(eliminate_trace(Trace(Gamma(M => lower) *
                                                  Gamma(N => lower)))))                                     

out("a Pauli matrix", Pauli(μ => upper))
out("an inner product", Bra("ψ") * Ket("ϕ"))

# Matrices do not commute until told. Declared relations feed back into simplify.
X, P = Mat("X", μ => upper), Mat("P", ν => upper)
@commutator [X, P] = im * Metric(μ => upper, ν => upper)
out("XP - PX", simplify(X * P - P * X))

####################################################################################
section("Lie brackets and algebra generators")
####################################################################################

out("[X,P]", @Lie [X, P])
out("{X,P}", @Lie {X, P})

# the conformal algebra in 4D, as generators that act on fields
println("Representations of algebras - conformal weights:")
P4, M4, D4, K4 = generate_conformal_algebra(4)
@syms Δ
Φ = Ket("Φ", x)
@action D4(Φ) = im * Δ * Φ
@action M4(Φ) = 0
out("|Φ⟩ is an eigenstate of D", derive(D4(Φ)))
out("[D, P]|Φ⟩  raises weight", derive(commute(D4(P4(Φ)))))

####################################################################################
section("Derivatives, the product rule, and integration by parts")
####################################################################################

∂ = PD("∂", [x], [μ => lower])
A, B = Tensor("A", [ν => lower], x), Tensor("B", [ρ => lower], x)
out("∂_μ(AB)", apply_derivative(∂(A * B)))
out("∂_μ of a constant", apply_derivative(PD("∂", Tensor("k", [ν => lower]), x,
                                             μ => lower)))
# modulo total derivatives, a derivative may be moved from one factor to the other
out("(∂A)B + A(∂B) ≃ 0", drop_total_derivatives(∂(A) * B + A * ∂(B)))
out("(∂A)B ≃ -A(∂B)", drop_total_derivatives(∂(A) * B))

####################################################################################
section("Pattern matching with wildcard indices")
####################################################################################

# ?α is a wildcard slot: the rule fires wherever F appears, whatever it is indexed on
α̂, β̂ = DummyPatternIndex(["α", "β"], st)
Fdef = (PD("∂", Tensor("A", [β̂ => lower], x), x, α̂ => lower)
        - PD("∂", Tensor("A", [α̂ => lower], x), x, β̂ => lower))
out("F_μν → ∂A - ∂A", replace_pattern(T"F_[μν]{x}", "F_{?α?β}", Fdef))

####################################################################################
section("Variational calculus — a Lagrangian gives back its field equations")
####################################################################################

dirac = IndexSet("dirac", 4, upper, fermionic)
i, j = Index(["i", "j"], dirac)
@syms m e
ψ = Tensor("ψ", [i => lower], x)
L = (-1//4 * AntisymmetricTensor("F", [μ => upper, ν => upper], x)
             * AntisymmetricTensor("F", [μ => lower, ν => lower], x)
     + Tensor("̄ψ", [i => upper], x) *
       (im * Gamma(μ => upper) * (PD("∂", ψ, x, [μ => lower])
                                  + im * e * Tensor("A", [μ => lower], x) * ψ) - m * ψ))
println("  L = ", L)
LA = simplify(replace_pattern(replace_pattern(L, "F_{?α?β}", Fdef), "F^{?α?β}",
    PD("∂", Tensor("A", [β̂ => upper], x), x, α̂ => upper)
    - PD("∂", Tensor("A", [α̂ => upper], x), x, β̂ => upper)))
println("  δL/δA_ρ  = ", Alakazam.hard_simplify(EL(LA, Tensor("A", [ρ => upper], x),
                                                   PD("∂", x, [τ => upper]))))
println("  δL/δ̄ψ    = ", Alakazam.hard_simplify(EL(LA, Tensor("̄ψ", [j => lower], x),
                                                   PD("∂", x, [ρ => upper]))))

####################################################################################
section("Differential geometry — ⋆, ⨼, ∧, 𝖽")
####################################################################################

V = Tensor("V", [μ => lower])
ω = AntisymmetricTensor("ω", [μ => lower, ν => lower])
out("⋆ω", ⋆(ω))
out("V ⨼ ω", V ⨼ ω)
out("A ∧ B", Tensor("A", [μ => lower]) ∧ Tensor("B", [ν => lower]))
out("𝖽A", 𝖽(Tensor("A", μ => lower, x)))
out("𝖽𝖽A = 0", isZeroTens(simplify(apply_derivative(𝖽(𝖽(Tensor("A", μ => lower, x)))))))

####################################################################################
section("Complex conjugation — reality follows from the indices")
####################################################################################

# the number field decides how a spacetime field conjugates
out("real field", Tensor("A", μ => lower, ℝ)')
out("imaginary field", Tensor("B", μ => lower, iℝ)')
out("complex field", Tensor("c", μ => lower, ℂ)')

# and for internal indices, the algebra decides
su2 = IndexSet("su2", 2, lower, SUAlgebra(2))
@conjugates su2
a, b = Index(["a", "b"], su2)               # internal indices stay Latin
out("SU(2) doublet", Tensor("λ", a => lower)')
out("...with index raised", Tensor("λ", a => upper)')

####################################################################################
section("Superspace — expanding a superfield in Grassmann coordinates")
####################################################################################

spinor = IndexSet("spinor", 2, lower, EpsilonTensor, fermionic)
α, β = Index(["α", "β"], spinor)            # spinor indices Greek
flavour = IndexSet("flavour", 2, upper)
k, l = Index(["k", "l"], flavour)           # internal indices Latin
θc = Coordinate("θ", spinor)

θα = Tensor("θ", [α => lower], θc)
out("grading(θ)", grading(θα))
# the expansion terminates on its own: θ is Grassmann, so it closes at θ²
out("superfield in θ_α", generate_superfield_naive("X", [θα], x))
# internal indices ride along beside the spinor ones, untouched
out("spinor × flavour", Tensor("ψ", [α => lower, k => upper], x))

####################################################################################
section("Serialisation — save the state, clear it, load it back")
####################################################################################

expr = Tensor("Q", [μ => lower], x) * Tensor("R", [μ => upper], x)
path = joinpath(tempdir(), "alakazam_demo_expr")
save(TensorExpression(expr), path, true)          # true: include the global header
out("saved", expr)
out("loaded back", load(path))
out("round trip agrees", load(path) == TensorExpression(expr))

println("\n", "─"^76)
println("  See demo/ for complete calculations, and the paper for the full manual.")
println("─"^76)
