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

    # with no ^, A_b{x} is a bare curly rather than a call, and used to run off the end of it
    let t = @Tensor A_b{x}
        @test t == Tensor("A", [b => lower], x) && function_of(t) == [x]
    end
    let t = @Tensor A_bc{x,θ}
        @test t == Tensor("A", [b => lower, c => lower], [x, θ]) && Set(function_of(t)) == Set([x, θ])
    end
    # a bare name binds too, like every other form
    let t = @Tensor Z
        @test Z === t && t == Tensor("Z")
    end
    let t = @Tensor A{x}
        @test isempty(indices(t)) && function_of(t) == [x]
    end
    @test_throws LoadError eval(:(@Tensor A_{b}{x}))   # index groups are not braced


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

    ###############################################
    ### | bars hold slots out of the group being accumulated
    ###############################################
    # The bracket spans three slots but only two join the group, so the tableau must be built
    # from what was accumulated rather than from the span. Both bracket kinds count alike.
    @test length(collect(tableaux(T"A_[a|b|c]"))[1].index_assignment) == 2
    @test collect(tableaux(T"A_[a|b|c]"))[1].index_assignment ==
          collect(tableaux(T"A_(a|b|c)"))[1].index_assignment
    # antisymmetric over slots 1 and 3, with slot 2 held out
    @test isZeroTens(simplify(T"A_[a|b|c]" + T"A_[c|b|a]"))
    # and an unpaused group is unaffected
    @test isZeroTens(simplify(T"A_[ab]" + T"A_[ba]"))
    @test length(collect(tableaux(T"A_[abc]"))[1].index_assignment) == 3

    ###############################################
    ### a NAME may carry the same ASCII decorations an index may
    ###############################################
    # The mark leads, which is the form conjugation and to_latex already use, so a name written
    # with a suffix and a name produced by conjugating are the same object.
    let BAR = string(Alakazam.COMBINING_MACRON), DOT = string(Alakazam.COMBINING_DOT),
        HAT = string(Alakazam.COMBINING_HAT), TIL = string(Alakazam.COMBINING_TILDE)

        @test name(T"Abar_a") == BAR * "A"
        @test name(T"Adot_a") == DOT * "A"
        @test name(T"Ahat_a") == HAT * "A"
        # tilde was in INDEX_SUFFIXES but missing from the tensor-string table, so it used to be
        # read as part of a label rather than as a decoration
        @test name(T"Atilde_a") == TIL * "A"

        # a decorated name and a decorated index together, which is the spinor case
        let t = T"Abar_?xdot"
            @test name(t) == BAR * "A"
            @test first(indices(t))[1].name == DOT * "x"
        end

        # a decorated name with nothing after it is still a tensor
        @test name(T"Abar") == BAR * "A"
        @test isempty(indices(T"Abar"))

        # the wildcard marker stays outside the mark
        @test name(T"?Abar_a") == "?" * BAR * "A"

        # and an undecorated name is untouched
        @test name(T"A_a") == "A"
        @test first(indices(T"Abar^a"))[2] == upper
    end

    ###############################################
    ### @Tensor A_bc takes the suffixes too
    ###############################################
    # This form BINDS a variable named after the tensor, and a decorated name is not a valid
    # identifier, so the binding keeps the spelling typed while the tensor carries the mark.
    let BAR = string(Alakazam.COMBINING_MACRON), DOT = string(Alakazam.COMBINING_DOT)
        @Tensor Wbar_ab
        @test name(Wbar) == BAR * "W"
        @test [i[1].name for i in indices(Wbar)] == ["a", "b"]

        # a decorated index: the token is looked up as the variable the caller declared
        @Tensor V_αdot
        @test first(indices(V))[1].name == DOT * "α"

        # a name with no indices at all
        @Tensor Ubar
        @test name(Ubar) == BAR * "U"
        @test isempty(indices(Ubar))

        # and the undecorated form is unchanged
        @Tensor Q_ab
        @test name(Q) == "Q"
        @test [i[1].name for i in indices(Q)] == ["a", "b"]
    end

end
