# Not in the precompile workload: declaring compatibility adds a method to Alakazam, which must not
# outlive the test. Declared at top level, as a user would, and removed at the end.
@compatible Derivative Metric

@testset "metric compatible derivatives" begin
    mst = IndexSet("mc_spacetime", 4, lower, Metric)
    mμ, mν, mρ, mσ = Index(["mμ", "mν", "mρ", "mσ"], mst)
    mx = Coordinate("mc_x", mst)
    metric_function!(mst, mx)
    mets(e) = [t for t in terms(summands(TensorExpression(e))[1][2]) if t isa Metric]

    A = Tensor("A", [mμ => upper], mx)
    g = Metric(mμ => lower, mν => lower)
    ∇A = Derivative("∇", A, [mx], [mν => lower])
    ∂A = PartialDerivative("∂", A, [mx], [mν => lower])

    #the macro defines the function, for the declared types and their subtypes only
    @test metric_compatible(Derivative, Metric) && !metric_compatible(PartialDerivative, Metric)
    @test metric_compatible(∇A, g) && Alakazam.commutes_with_metric(∇A, g)
    #the wrong way round is refused, and defines nothing
    @test_throws ArgumentError eval(:(@compatible Metric Derivative))
    @test !metric_compatible(Metric, Derivative)
    @test_throws LoadError eval(:(@compatible Derivative))              #refused at expansion

    #popping: ∇ annihilates g although g depends on x, so g comes out
    @test function_of(only(mets(pop_metric(∇A)))) == [mx]
    @test isempty(mets(operand(terms(summands(TensorExpression(pop_metric(∇A)))[1][2])[end])))
    #declared for ∇ only: ∂ still keeps it
    @test !metric_compatible(∂A, g) && isempty(mets(pop_metric(∂A)))
    #∇ carries it out of itself but not out of a ∂ inside it
    let p = isTensor(pop_metric(Derivative("∇", ∂A, [mx], [mρ => lower])))
        @test p isa Derivative && isempty(mets(operand(p)))
        @test !isempty(mets(operand(isTensor(operand(p)))))
    end

    #applying: ∇g = 0, so ∇(gA) = g∇A, while ∂g survives
    @test isZeroTens(apply_derivative(Derivative("∇", TensorExpression(g), [mx], [mρ => lower])))
    @test !isZeroTens(apply_derivative(PartialDerivative("∂", TensorExpression(g), [mx], [mρ => lower])))
    gA = Metric(mμ => lower, mσ => lower) * A
    @test apply_derivative(Derivative("∇", gA, [mx], [mν => lower])) ==
          TensorExpression(Metric(mμ => lower, mσ => lower) * Derivative("∇", A, [mx], [mν => lower]))
    @test length(summands(apply_derivative(PartialDerivative("∂", gA, [mx], [mν => lower])))) == 2

    #and out of a ∂ it was not in: a flat set's metric leaves ∂ as before
    metric_function!(mst)
    @test !isempty(mets(pop_metric(PartialDerivative("∂", A, [mx], [mν => lower]))))
end

Base.delete_method(which(metric_compatible, (Type{Derivative}, Type{Metric})))
@test !metric_compatible(Derivative, Metric)
