@testset "Serialisation Import" begin
    #temp store global stuff
    aa=deepcopy(get_global_values())
    for (k, v) in get_global_values()
        @test !isempty(v)
    end
    # @test !isempty(get_global_values())
    clear_global!()
    for (k, v) in get_global_values()
        @test isempty(v)
    end
    @test isempty(SET_REGISTRY)
    for (k, v) in aa
        @test !isempty(v)
    end
    # expr2=load("test_expression_headers")
    _TESTFILE2 = joinpath(tempdir(), "test_expression_headers")
    expr2 = @load _TESTFILE2
    @test !isempty(get_global_values())
    #matches previous load
    @test aa==get_global_values()

    # @test get_global_values()==aa
    a=get_index("a", "flavour_indices")
    b=get_index("b", "flavour_indices")
    d=SymmetricTensor("T", [a=>upper, b=>upper])
    d2=SymmetricTensor("T", [a=>lower, b=>lower])
    @syms ħ
    expr=(2*ħ)*d*d2
    @test expr==expr2 #check laoding symbolic exp
    #even reconstructs indices in user scope from symbols (dont need get_index)
    @test (e; true)
end