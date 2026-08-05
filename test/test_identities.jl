#The identities the package ships with. Makes sure theyre registered
@testset "Default identities" begin

    st = IndexSet("ident_spacetime", (1, 3), upper, Metric)
    @indices μ:τ ∈ st
    g(i, p=upper) = Gamma(i => p)
    ca(x) = contract_adjacent(x, Gamma)

    @test st.range == 4

    ###############################################
    ### Clifford algebra
    ###############################################
    # gammas must not report as freely commuting: if they did, expressions would reduce towrong answers rather than simply failing to reduce
    @test !Alakazam.can_swap(g(μ), g(ν))
    @test sort(TensorExpression(g(ν) * g(μ))) == -g(μ) * g(ν) + 2 * Metric([μ => upper, ν => upper])
    @test sort(TensorExpression(Gamma5() * g(ν))) == -(g(ν) * Gamma5())
    # γ5 γ5 = 1, in any dimension and any representation
    @test contract_adjacent(TensorExpression(Gamma5() * Gamma5()), Gamma5) == TensorExpression(1)

    ###############################################
    ### contraction identities: dimension read at reduction time
    ###############################################
    @test ca(g(μ) * g(μ, lower)) == TensorExpression(4)
    @test ca(g(μ) * g(ν) * g(μ, lower)) == -2 * g(ν)
    @test Alakazam.simplify(ca(g(μ) * g(ν) * g(ρ) * g(μ, lower))) ==
          4 * Metric([ν => upper, ρ => upper])
    # the recursion continues to longer runs: γ^μγ^νγ^ργ^σγ_μ = -2γ^σγ^ργ^ν in four dimensions
    @test Alakazam.simplify(ca(g(μ) * g(ν) * g(ρ) * g(σ) * g(μ, lower))) ==
          Alakazam.simplify(-2 * g(σ) * g(ρ) * g(ν))

    ###############################################
    ### traces, with nothing loaded
    ###############################################
    # these need no setup at all, should work from fresh session
    @test isZeroTens(reduce_trace(Trace(g(μ))))
    @test isZeroTens(reduce_trace(Trace(g(μ) * g(ν) * g(ρ))))
    @test isZeroTens(reduce_trace(Trace(Gamma5())))
    @test isZeroTens(reduce_trace(Trace(Gamma5() * g(μ) * g(ν))))
    # matched up to cyclic rotation, so γ5 need not be written first
    @test isZeroTens(reduce_trace(Trace(g(μ) * g(ν) * Gamma5())))
    # the one carrying tr(𝟙), which takes the Dirac value 2^{d/2} from the index set
    @test reduce_trace(Trace(g(μ) * g(ν))) == 4 * Metric([μ => upper, ν => upper])
    # exactly 4, not 4.0: a Float here would silently destroy exact arithmetic downstream
    @test first_coef(reduce_trace(Trace(g(μ) * g(ν)))) isa Integer
    @test first_coef(reduce_trace(Trace(g(μ) * g(ν)))) === 4
    # a metric wants both indices at one level, so one raised against one lowered is a delta
    @test reduce_trace(Trace(g(μ) * g(ν, lower))) == 4 * KroneckerDelta([μ => upper, ν => lower])

    ###############################################
    ### two index sets at once, with no registration and no interference
    ###############################################
    big = IndexSet("bigspacetime", (1, 9), upper, Metric)
    @indices a:d ∈ big
    @test big.range == 10
    @test ca(Gamma(a => upper) * Gamma(a => lower)) == TensorExpression(10)
    @test ca(Gamma(a => upper) * Gamma(b => upper) * Gamma(a => lower)) == -8 * Gamma(b => upper)
    # the 4d set is unaffected
    @test ca(g(μ) * g(μ, lower)) == TensorExpression(4)
    # and the Clifford relation adapts to the other set too
    @test sort(TensorExpression(Gamma(b => upper) * Gamma(a => upper))) ==
          -Gamma(a => upper) * Gamma(b => upper) + 2 * Metric([a => upper, b => upper])
    # traces likewise: 2^{10/2} = 32 here and 2^{4/2} = 4 there, in the same session. 
    @test reduce_trace(Trace(Gamma(a => upper) * Gamma(b => upper))) ==
          32 * Metric([a => upper, b => upper])
    @test reduce_trace(Trace(g(μ) * g(ν))) == 4 * Metric([μ => upper, ν => upper])

    ###############################################
    ### clear_global! empties the registries but puts the defaults back
    ###############################################
    clear_global!()
    # the Clifford relations survive, because `clear_global!` reloads them on the way out
    @test !Alakazam.can_swap(Gamma(μ => upper), Gamma(ν => upper))
    # the contraction and trace identities are methods, so nothing had to be reloaded for them
    @test ca(Gamma(μ => upper) * Gamma(μ => lower)) == TensorExpression(4)
    @test reduce_trace(Trace(g(μ) * g(ν))) == 4 * Metric([μ => upper, ν => upper])
    # loading again by hand harmless
    load_gamma_defaults!()
    @test !Alakazam.can_swap(Gamma(μ => upper), Gamma(ν => upper))

end
