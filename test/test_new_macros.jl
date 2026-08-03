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
    @test Alakazam.name(indices(DT"?A_?xdot")[1][1]) == "xdot"
    @test indices(DT"?A_?xdot")[1][1] isa DummyPatternIndex

    ###############################################
    ### {...} is the coordinate group, not an index group
    ###############################################
    @test length(function_of(T"A_a{x}")) == 1              # valid: x is a Coordinate
    @test_throws ErrorException T"A_{a}"                   # a is an Index, not a Coordinate
    @test_throws ErrorException T"A_{qqqq}"                # no such variable at all

end
