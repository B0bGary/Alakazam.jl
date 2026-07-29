@testset "Basic Tensor Manipulation" begin
    #setup
    if 1 == 1
        flavour_indices = IndexSet("flavour", 2, EpsilonTensor, upper)
        spacetime_indices = IndexSet("spacetime", 4, upper)
        @indices a:f ∈ spacetime_indices
        @indices i:z ∈ flavour_indices
        x = Coordinate("x", spacetime_indices)
        θ = Coordinate("θ", flavour_indices)
    end
    ###############################################
    ### begin test cases
    ###############################################
    #consistent  zero() behaviour
    @test Alakazam. zero()==TensorExpression(0)
    @test !(Alakazam. zero()==TensorExpression(1))
    @test Alakazam. zero()==TensorExpression()
    @test Alakazam. zero()==0
    @test Alakazam. zero()==0.0
    ##############
    #basic creation
    T1=Tensor("T")
    T2=Tensor("T", a=>lower)
    #should throw error since free indices don't match
    @test_throws FreeIndexMismatch() T1+T2
    ##############
    #fetching properties
    @test indices(T1)==[]
    @test indices(T2)==[a=>lower]

    #test metric popping
    texp1=pop_metric(T2)
    @test contract_metrics(texp1)==T2 #should contract properly


    #test kronecker elimination
    exp=KroneckerDelta([a=>upper, c=>lower])*Tensor("T", a=>lower)
    @test eliminate_kronecker(exp)==Tensor("T", c=>lower)

    #multiple kroneckers
    exp2=KroneckerDelta([a=>lower, b=>upper])KroneckerDelta([a=>upper, c=>lower])*Tensor("T", c=>upper)
    @test eliminate_kronecker(exp2)==Tensor("T", b=>upper)

    #default pos is upper, and su(2) indices so use epsilon as metric
    T1=Tensor("T", i=>lower, j=>lower)
    @test pop_metric(T1)==EpsilonTensor([i=>lower, k=>lower])*EpsilonTensor([j=>lower, l=>lower])*Tensor("T", [k=>upper, l=>upper])

    #contracting an epsilon needs to respect conventions for raising/lowering indices, so should give a minus sign
    tt=EpsilonTensor([i=>upper, k=>upper])*Tensor("T", i=>lower, j=>lower)
    @test eliminate_epsilon(tt)==-Tensor("T", k=>upper, j=>lower)

    @test (update_indices!(1);true)
    @test (rename_dummies(Tensor("T", [c=>lower,c=>upper]));true)

    @test (to_latex(tt);true)
    @test (to_latex((1+im/2)*tt);true)

    @test (fresh_indices(T1,tt);true)
end