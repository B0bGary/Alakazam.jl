@testset "Tensor Macros" begin
    ####setup
    flavour_indices = IndexSet("flavour_indices", 2, upper)
    @indices a:z∈flavour_indices
    αdot=Index("̇α", flavour_indices)
    α=Index("A", flavour_indices)
    x = Coordinate("x", [flavour_indices])
    θ = Coordinate("θ", [flavour_indices])
    θbar = Coordinate("Θ", [flavour_indices])


    # Underscore expr::Symbol cases
    @test Tensor("A", [a => lower, b => lower]) == @Tensor A_ab
    @test Tensor("A", [a => upper, b => upper]) == @Tensor A^ab
    @test Tensor("A", [a => upper, b => upper, c=>lower, d=>upper]) == @Tensor A^ab_c^d
    @test Tensor("A", [a => lower, b => upper, c=>lower, d=>upper]) == @Tensor A_a^b_c^d

    @test Tensor("A", [a => lower, b => upper, c=>lower, d=>upper], x) == @Tensor A_a^b_c^d{x}
    @test Tensor("A", [a => lower, b => upper, c=>lower, d=>upper], [x, θ]) == @Tensor A_a^b_c^d{x,θ}


    #string parsing
    @test Tensor("A", [a => lower, b => lower]) ≅ @Tensor "A_ab"
    @test Tensor("A", [a => lower, b => upper]) ≅ @Tensor "A_a^b"
    A=@Tensor "A_[ab]c(d|e|f)"
    #check that tableau is parsed faithfully
    A2=Tensor("A", [a => lower, b => lower, c=>lower, d => lower, e => lower, f => lower], [YoungTableau([1, 1], [1, 2]), YoungTableau([2], [4, 6])])
    @test tableaux(A)==tableaux(A2)
    A=T"A_[ab]c(d|e|f)"
    @test tableaux(A)==tableaux(A2)


    #multiple coords
    B=@Tensor "B_ab{x,θ,θbar}"
    @test function_of(B)==[x, θbar, θ]
    B2=Tensor("B", [a => lower, b => lower], [θbar])
    @test function_of(B2) == [θbar]
    @test !(Tensor("B", [a => lower, b => lower], x) ≅ B)

    @test (A=T"A_a,b,c,d,e,f"; true)
    @test (A=T"A_a,^b,_c,d,e,f"; true)
    @test (A=T"A_a,^b,_c,^d,e,f"; true)
    @test (A=T"A_a,^b,_c,^d,e,f{x}"; true)
    @test (A=T"A_a,^b,_c,^d,e,f{x,θ}"; true)
    @test (A=T"A_a,^b,_αdot,^d,e,f{x,θ}"; true)
    @test (A=T"A_a^b_αdot^def{x,θ}"; true)
    @test (A=T"A_a^b_αdotα^def{x,θ}"; true)

    @test (A=T"A_(a^b)_[c^d]ef{x,θ}"; true)
    #powers
    @test (T"F_ab")^2==T"F_ab"*T"F^ab"
    #higher powers
    @test (T"F_ab")^3==T"F_ab"*T"F_cd"*T"F^cd"
    K=T"K_ab"
    @test (K;true) #user scope


    ###############################################
    ### suffixed index names (dot / bar / hat) must keep resolving as INDICES, not labels
    ###############################################
    abar = Index("ā", flavour_indices)
    ahat = Index("â", flavour_indices)
    adot = Index("̇a", flavour_indices)

    @test indices(T"A_adot") == [adot => lower] && isempty(labels(T"A_adot"))
    @test indices(T"A_abar") == [abar => lower] && isempty(labels(T"A_abar"))
    @test indices(T"A_ahat") == [ahat => lower] && isempty(labels(T"A_ahat"))
    # a suffixed name followed by a plain one: the suffix must be consumed with its stem
    @test indices(T"A_adotb") == [adot => lower, b => lower]
    @test indices(T"A^adot_b") == [adot => upper, b => lower]
    # the pre-existing dotted index from the setup above
    @test indices(T"A^αdot") == [αdot => upper]
    # suffixes survive in pattern strings too
    # a wildcard has no scope variable to resolve, so the parser folds the suffix into its
    # combining mark itself: ?xdot denotes the dotted index, not one literally called "xdot"
    @test Alakazam.name(indices(T"?A_?xdot")[1][1]) == string('̇') * "x"
    @test Alakazam.name(indices(T"?A_?xbar")[1][1]) == string('̄') * "x"
    @test Alakazam.name(indices(T"?A_?xhat")[1][1]) == string('̂') * "x"
    @test Alakazam.name(indices(T"?A_?x")[1][1]) == "x"          # unsuffixed is untouched
    @test indices(T"?A_?xdot")[1][1] isa DummyPatternIndex
    # the converted name is what the binding is keyed on, so a dotted wildcard round-trips
    # between pattern and replacement. (It does not restrict matching: a wildcard matches any
    # index, and its own name is only a binding label.)
    let dt = Index(string('̇') * "q", flavour_indices)
        @test replace_pattern(TensorExpression(Tensor("V", [dt => lower])),
                              TensorTerm(T"V_?qdot"), TensorExpression(T"W_?qdot")) ==
              TensorExpression(Tensor("W", [dt => lower]))
        @test replace_pattern(TensorExpression(Tensor("V", [b => lower])),
                              TensorTerm(T"V_?qdot"), TensorExpression(T"W_?qdot")) ==
              TensorExpression(Tensor("W", [b => lower]))
    end

    ###############################################
    ### {...} is the coordinate group, not an index group
    ###############################################
    @test length(function_of(T"A_a{x}")) == 1              # valid: x is a Coordinate
    @test_throws ErrorException T"A_{a}"                   # a is an Index, not a Coordinate
    @test_throws ErrorException T"A_{qqqq}"                # no such variable at all

end
