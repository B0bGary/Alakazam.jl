@testset "Basic Tensor Arithmetic" begin
    #setup
    if 1 == 1
        flavour_indices = IndexSet("flavour_indices", 2, upper)
        @indices a:z ∈ flavour_indices

        x = Coordinate("x", flavour_indices)
        θ = Coordinate("θ", flavour_indices)
        θbar = Coordinate("Θ", flavour_indices)
    end
    ###############################################
    ### begin test cases
    ###############################################


    ###############################################
    ### BASIC ARITHMETIC
    ###############################################

    #TENSORS

    #Tensor addition basics
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    S1 = Tensor("S", [a => lower, b => lower, c => lower], x)
    @test (T1 + S1; true)
    @test (-T1; true) #check negation


    T2 = Tensor("T", [a => lower, b => lower, d => lower], x)
    T3 = Tensor("T", [a => lower, b => upper, c => lower], x)
    @test_throws Alakazam.FreeIndexMismatch T1 + T2
    @test_throws Alakazam.FreeIndexMismatch T1 + T3

    #add consts to funcs
    F = Tensor("f", x)
    @test (F + 1; true)
    @syms α
    @test (F + α; true)
    #cant add a const to tensor
    @test_throws Alakazam.FreeIndexMismatch T1 + α

    #division by numbers
    @test (T1 / 2; true)
    @test (T1 / α; true)

    #multiplication
    @test_throws Alakazam.SameLevelDummyError T1 * T3 # indices a and c are on same level, incorrect Einstein notation error.
    T4 = Tensor("T", [a => upper, b => upper, c => upper], x)
    @test (T4 * T1; true) #einstein summation is correct
    T5 = Tensor("T", [a => upper, b => upper, b => lower], x)
    @test_throws Alakazam.OverusedIndices T5 * T1 #too many repeated indices


    ### TENSOR TERMS
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    T2 = Tensor("T", [d => lower, e => lower, f => lower], x)
    T3 = Tensor("T", [i => lower, j => lower, k => lower], x)
    @test ((T1 * T2) * T3) == (T1 * (T2 * T3)) #multiplication should be associative

    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    T2 = Tensor("T", [d => lower, e => lower, f => lower], x)
    T3 = Tensor("T", [d => lower, a => lower, k => lower], x)
    @test_throws Alakazam.SameLevelDummyError ((T1 * T2) * T3)#adding last term gives incorrect einstein not
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    T2 = Tensor("T", [d => lower, a => upper, f => lower], x)
    T3 = Tensor("T", [e => lower, a => lower, k => lower], x)
    @test_throws Alakazam.OverusedIndices ((T1 * T2) * T3)#adding last term has index already used in summation


    #### TERM TIMES TERM
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    T2 = Tensor("T", [d => lower, e => lower, f => lower], x)
    T3 = Tensor("T", [i => lower, j => lower, k => lower], x)
    T4 = Tensor("T", [s => lower, t => lower, u => lower], x)
    #5 ways of bracketing:
    @test ((T1 * T2) * (T3 * T4)) == ((T1 * (T2 * T3)) * T4) #multiplication should be associative
    @test ((T1 * T2) * (T3 * T4)) == (((T1 * T2) * T3) * T4) #multiplication should be associative
    @test (T1 * (T2 * (T3 * T4))) == ((T1 * T2) * (T3 * T4)) #multiplication should be associative
    @test (T1 * ((T2 * T3) * T4)) == ((T1 * T2) * (T3 * T4)) #multiplication should be associative

    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    T2 = Tensor("T", [d => lower, a => upper, f => lower], x)
    T3 = Tensor("T", [a => lower, j => lower, k => lower], x)
    T4 = Tensor("T", [s => lower, a => upper, u => lower], x)
    @test_throws Alakazam.OverusedIndices ((T1 * T2) * (T3 * T4)) #the two sums use same dummt indices
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    T2 = Tensor("T", [d => lower, a => upper, f => lower], x)
    T3 = Tensor("T", [b => lower, j => lower, k => lower], x)
    T4 = Tensor("T", [s => lower, b => upper, u => lower], x)
    @test_throws Alakazam.OverusedIndices ((T1 * T2) * (T3 * T4)) #the b is used too many times
    T1 = Tensor("T", [a => lower, i => lower, c => lower], x)
    T2 = Tensor("T", [d => lower, a => upper, f => lower], x)
    T3 = Tensor("T", [b => lower, j => lower, k => lower], x)
    T4 = Tensor("T", [s => lower, b => upper, c => lower], x)
    @test_throws Alakazam.SameLevelDummyError ((T1 * T2) * (T3 * T4)) #the c is on wrong level for einstein not

    @test (-(T1 * T2); true) #test negation of TT

    #tensor exp
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    R1 = Tensor("R", [a => lower, b => lower, c => lower], x)
    T4 = Tensor("T", [s => lower, t => upper, u => lower], x)
    T5 = Tensor("T", [d => lower, e => upper, t => lower], x)
    #exp and tens
    @test ((T1 + R1) * T4) == (T1 * T4 + R1 * T4) #associativity
    @test T4 * ((T1 + R1)) == (T4 * T1 + T4 * R1) #associativity
    #exp and term
    T45 = T4 * T5
    @test (T1 + R1) * T45 == (T1 * T45 + R1 * T45)#associativity
    @test T45 * (T1 + R1) == (T45 * T1 + T45 * R1)#associativity
    T1 = Tensor("T", [a => lower, b => lower, c => lower], x)
    R1 = Tensor("R", [a => lower, b => lower, c => lower], x)
    T4 = Tensor("T", [c => lower, t => upper, u => lower], x)
    T5 = Tensor("T", [d => lower, e => upper, t => lower], x)
    T45 = T4 * T5
    @test_throws Alakazam.SameLevelDummyError ((T1 + R1) * T45) #c is repeated
    T1 = Tensor("T", [a => lower, b => lower, t => lower], x)
    R1 = Tensor("R", [a => lower, b => lower, t => lower], x)
    T4 = Tensor("T", [e => lower, t => upper, u => lower], x)
    T5 = Tensor("T", [d => lower, e => upper, t => lower], x)
    T45 = T4 * T5
    @test_throws Alakazam.OverusedIndices ((T1 + R1) * T45) #t already used in summation

    ##tex/tex
    T = Tensor("T", [a => lower, b => lower, c => lower], x)
    S = Tensor("S", [a => lower, b => lower, c => lower], x)
    R = Tensor("R", [d => lower, e => lower, f => lower], x)
    P = Tensor("P", [d => lower, e => lower, f => lower], x)
    @test (T + S) * (R + P) == T * R + T * P + S * P + S * R#associativity, even for terms out of "FOIL" order

    @test (-((T + S) * (R + P)); true) #test negation of whole expressions defined

    #test  ==(a::TensorTerm, b::TensorExpression) = b == a
    @test T4 * T5 == 1*T4 * T5
    @test !(1*T4*T5 == T4)
    @test !(1==T4)

end