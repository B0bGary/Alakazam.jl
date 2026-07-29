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
end