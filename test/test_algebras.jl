@testset "Algebras" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetime", 4, lower)
        spinor_indices=IndexSet("spinor_indices", 2, lower, Alakazam.fermionic, EpsilonTensor)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2, Alakazam.fermionic, lower, EpsilonTensor)
        α, β, γ, χ = Index(["α", "β", "γ", "χ"], spinor_indices)
        αd, βd, γd, χd = Index(["̇α", "̇β", "̇γ", "̇χ"], spinor_conj_indices)
        μ, ν, ρ, τ = Index(["μ", "ν", "ρ", "τ"], spacetime_indices)
        flavour_indices = IndexSet("flavour_indices", 2, upper, Alakazam.fermionic, EpsilonTensor)
        i, j, k, l, m, n, p, q, r = Index(["i", "j", "k", "l", "m", "n", "p", "q", "r"], flavour_indices)
        x = Coordinate("x", spacetime_indices)
        θ = Coordinate("θ", spinor_indices)
        θbar = Coordinate("Θ", spinor_conj_indices)
        custom_sort_order!(["θ", "Θ", "ϵ"])
    end

    P, M, D, K=generate_conformal_algebra(4)

    @syms Δ
    X=Tensor("X", x)
    Y=Tensor("Y", x)
    # @action K(X) = 0
    @action D(X) = im*Δ*X
    @action M(X) = 0
    desc=P(X)
    deri=derive(commute(D(desc)))
    # @show deri
    @test deri==(im+im*Δ)*desc #increase conf wght
    Pp=set_op_indices!(copy(P), [ρ=>lower])
    Pp2=set_op_indices!(copy(P), [τ=>lower])

    # Mop=set_op_indices!(copy(M), op_pos(indices(M)))
    # @show Mop
    # @action Mop(X)=0
    #  @show spn=derive((commute_right(M(Pp2(Pp(X))-Pp(Pp2(X))),"M")))
    #  @show spn2=product_rule(M(Pp2(Pp(X))-Pp(Pp2(X))))
    #  @show out=Alakazam.simplify(Metric([ρ=>↑,τ=>↑])*derive((commute_right(spn2,"M"))))
    desc3=set_op_indices!(copy(P), [ρ=>lower])(set_op_indices!(copy(P), [ν=>lower])(P(2*X)))
    want = D(desc3)
    # @show want
    # @show  (commute_right(want,"D"))
#     a = derive(commute_right(want, "D"))
#     b = (3im + im*Δ) * desc3
#     sa = summands(a)[1];
#     sb = summands(b)[1]
#     ca = sa.first;
#     cb = sb.first
#     @show typeof(ca), typeof(cb)
#     @show ca isa Number, cb isa Number
#     @show isequal(ca, cb)
#     @show Alakazam._sym_equal(ca, cb)
#     @show which(Alakazam.coeff_equal, (typeof(ca), typeof(cb)))
#     @show Alakazam.coeff_equal(sa.first, sb.first)
#     @show sa.second == sb.second
#     using Symbolics,SymbolicUtils
#     d = Symbolics.simplify(Symbolics.expand(ca - cb))
# @show d
# @show typeof(d)
# @show d === 0
# @show SymbolicUtils.isconst(d)
# @show SymbolicUtils.isconst(d) ? SymbolicUtils.unwrap_const(d) : nothing
#     # @show Alakazam.same_tensorterm_mod_syms(sa.second, sb.second)
#     @show length(terms(sa.second)), length(terms(sb.second))
#     for (x, y) in zip(terms(sa.second), terms(sb.second))
#         @show x == y, name(x), name(y), indices(x), indices(y)
#     end
    @test derive(commute_right(want, "D"))==(3*im+im*Δ)*(desc3)
end
