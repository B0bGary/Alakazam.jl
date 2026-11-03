#A lightcone superspace whose Grassmanns carry an SU(2) index, θ^{+i} and Θ^{+}_i. The names are this
#file's own, so its actions and covariant derivatives cannot meet those registered elsewhere.
@testset "SU(2) indexed superspace" begin
    su2q = IndexSet("su2q", upper, SLAlgebra(2))
    IDXq = Index(["i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "a", "b", "c", "d"], su2q)
    iq, jq, kq = IDXq[1:3]
    aq, bq, cq, dq = IDXq[13:16]
    xq = Coordinate("su2q_x", [lightcone_indices, lightcone_indices])
    θq, Θq = Coordinate("su2q_θ", lightcone_indices), Coordinate("su2q_Θ", lightcone_indices)
    ∂pq, ∂mq = PD("∇q", xq, ℓ⁺ => lower, ℓ⁺ => lower), PD("∇q", xq, ℓ⁻ => lower, ℓ⁻ => lower)
    @commutator [∂pq, ∂mq] = 0
    ∂θ, ∂Θ = PD("∂q", θq, ℓ⁺ => lower, iq => lower), PD("𝛛q", Θq, ℓ⁺ => lower, iq => upper)
    @action ∂θ(Tensor("θ", [ℓ⁺ => upper, jq => upper], θq)) = KroneckerDelta([iq => lower, jq => upper])
    @action ∂θ(Tensor("θ", [ℓ⁺ => upper, jq => lower], θq)) = EpsilonTensor([jq => lower, iq => lower])
    @action ∂Θ(Tensor("Θ", [ℓ⁺ => upper, jq => lower], Θq)) = KroneckerDelta([iq => upper, jq => lower])
    @action ∂Θ(Tensor("Θ", [ℓ⁺ => upper, jq => upper], Θq)) = EpsilonTensor([jq => upper, iq => upper])
    Dc = Derivative("Dq", [xq, θq], ℓ⁺ => lower, iq => lower)
    D̄c = Derivative("D̄q", [xq, Θq], ℓ⁺ => lower, iq => upper)
    @covariant Dc = ∂θ - im//2 * Tensor("Θ", [ℓ⁺ => upper, iq => lower], Θq) * ∂pq
    @covariant D̄c = -∂Θ + im//2 * Tensor("θ", [ℓ⁺ => upper, iq => upper], θq) * ∂pq
    ss = Superspace([Tensor("θ", [ℓ⁺ => upper, iq => upper], θq), Tensor("Θ", [ℓ⁺ => upper, iq => lower], Θq)], [Dc, D̄c])
    Dq(ix) = F -> derive(Derivative("Dq", [xq, θq], ℓ⁺ => lower, ix => lower)(F), ss)
    D̄q(ix) = F -> derive(Derivative("D̄q", [xq, Θq], ℓ⁺ => lower, ix => upper)(F), ss)
    sfq(nm, inds) = generate_superfield(nm, Pair{IndexSuperType,IndexPosition}[inds...], ss, xq)
    J = sfq("J", [])

    #a rule written with distinct indices also acts where the derivative's index is contracted with
    #its target's: ∂_{+j}θ^{+j} = δ_j^j = 2, and D_cD̄^c𝒥 is δ_d^c D_cD̄^d𝒥
    ∂j = PD("∂q", θq, ℓ⁺ => lower, jq => lower)
    @test Alakazam.simplify(apply_derivative(TensorExpression(∂j(Tensor("θ", [ℓ⁺ => upper, jq => upper], θq))))) ==
          TensorExpression(2)
    @test isZeroTens(Alakazam.simplify(derive(Dq(cq)(D̄q(cq)(J))) -
                                       derive(KroneckerDelta([dq => lower, cq => upper]) * Dq(cq)(D̄q(dq)(J)))))

    #dummies private to a derivative's operand: ∂₋₋X shows only its first summand's indices, which
    #carry none, so × once kept Y's dummies and they met X's when the derivative was expanded. And a
    #dummy × had renamed to k met the derivative's own k
    X, Y = sfq("X", []), sfq("Y", [])
    @test D̄q(bq)(∂mq(X) × Y) isa TensorExpression
    @test D̄q(bq)(∂mq(∂mq(X)) × Y) isa TensorExpression
    @test Dq(kq)(X × Y) isa TensorExpression

    #a projection keeps the index of a single copy of D free, and pairs two copies with ε/2; and a
    #derivative of a whole superfield is applied before it is split by weight
    for (p, q, nlow, nup) in ((1, 0, 1, 0), (0, 1, 0, 1), (1, 1, 1, 1), (2, 0, 0, 0), (2, 1, 0, 1), (2, 2, 0, 0))
        c = project_to_component(TensorExpression(J), ss.covariant_basis, p, q)
        free = [f for f in keys(free_indices(c)) if !Alakazam.isLightcone(f)]
        @test count(f -> f[2] == lower, free) == nlow && count(f -> f[2] == upper, free) == nup
    end
    W = sfq("W", [ℓ⁻ => lower, aq => upper])
    @test project_to_component(TensorExpression(∂mq(W)), ss.covariant_basis, 0, 0) ==
          derive(∂mq(project_to_component(TensorExpression(W), ss.covariant_basis, 0, 0)))

    #and drop_nilpotent takes the superspace itself, which was ambiguous between two methods
    @test drop_nilpotent(TensorExpression(J), ss) == drop_nilpotent(TensorExpression(J), ss.grassmann_basis)
end
