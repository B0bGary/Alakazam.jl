@testset "Young Tableaux" begin

    #setup
    if 1 == 1
        flavour_indices = IndexSet("flavour_indices", 2, upper)
        @indices a:z ∈ flavour_indices
        x = Coordinate("x", [flavour_indices])
        θ = Coordinate("θ", [flavour_indices])
        θbar = Coordinate("Θ", [flavour_indices])
    end
    ###############################################
    ### begin test cases
    ###############################################
    ##
    #stop double-ups of tableaux
    index_assignment = [1, 2]
    row_lengths = [1, 1]
    asym_tabelau = YoungTableau(row_lengths, index_assignment)
    T1 = AntisymmetricTensor("T", [a => lower, b => lower, c => lower], asym_tabelau)
    @test length(tableaux(T1)) == 1
    index_assignment = [1, 2, 3]
    row_lengths = [1, 1, 1]
    asym_tabelau = YoungTableau(row_lengths, index_assignment)
    T1 = AntisymmetricTensor("T", [a => lower, b => lower, c => lower], asym_tabelau)
    @test length(tableaux(T1)) == 1

    index_assignment = [1, 2]
    row_lengths = [2]
    asym_tabelau = YoungTableau(row_lengths, index_assignment)
    T1 = SymmetricTensor("T", [a => lower, b => lower, c => lower], asym_tabelau)
    @test length(tableaux(T1)) == 1
    index_assignment = [1, 2, 3]
    row_lengths = [3]
    asym_tabelau = YoungTableau(row_lengths, index_assignment)
    T1 = SymmetricTensor("T", [a => lower, b => lower, c => lower], asym_tabelau)
    @test length(tableaux(T1)) == 1

    ################
    ##### irreps
    ###
    #simplest nontrivial tensor symmetry, irrep of S_3
    index_assignment = [1, 2, 3]
    row_lengths = [2, 1]
    t_tabelau = YoungTableau(row_lengths, index_assignment)
    T1 = Tensor("T", [a => lower, b => lower, c => lower], [t_tabelau])
    T2 = Tensor("T", [c => lower, b => lower, a => lower], [t_tabelau])
    T3 = Tensor("T", [b => lower, a => lower, c => lower], [t_tabelau])
    # @show Alakazam.young_projection(T1)+Alakazam.young_projection(T2)
    T4 = Tensor("T", [c => lower, a => lower, b => lower], [t_tabelau])
    T5 = Tensor("T", [b => lower, c => lower, a => lower], [t_tabelau])
    T6 = Tensor("T", [a => lower, c => lower, b => lower], [t_tabelau])
    @test Alakazam.young_projection(T1) - Alakazam.young_projection(Alakazam.young_projection(T1)) == Alakazam.zero() #idempotency of projector  
    @test Alakazam.young_projection(T1 + T4 + T5) == Alakazam.zero() #bianchi like identity




    ###
    #Riemann tensor
    index_assignment = [1, 3, 2, 4]
    row_lengths = [2, 2]
    riemann_tabelau = YoungTableau(row_lengths, index_assignment)
    

    R1 = Tensor("R", [a => lower, b => lower, c => lower, d => lower], [riemann_tabelau])
    # @show Alakazam.young_projection(R1)
    @test Alakazam.young_projection(Alakazam.young_projection(R1)) - Alakazam.young_projection(R1) == Alakazam.zero() #idempotency of projector

    #test the bianchi id
    R2 = Tensor("R", [a => lower, c => lower, d => lower, b => lower], [riemann_tabelau])
    R3 = Tensor("R", [a => lower, d => lower, b => lower, c => lower], [riemann_tabelau])
    @test Alakazam.young_projection(R1) + Alakazam.young_projection(R2) + Alakazam.young_projection(R3) == Alakazam.zero() #Bianchi identity is manifest

    #swapping syms
    R4 = Tensor("R", [a => lower, b => lower, d => lower, c => lower], [riemann_tabelau]) #antisym in last
    R5 = Tensor("R", [d => lower, c => lower, a => lower, b => lower], [riemann_tabelau]) #pair swap and antisym in first two 
    # @show tableaux(R4)
    @test R1 + R4 ≅ Alakazam.zero()
    @test R1 + R5 ≅ Alakazam.zero()

    #check that built in Riemann tensor has correct symmetries
    Riem = RiemannTensor([d => lower, c => lower, b => lower, a => lower])
    @test Riem - R1 ≅ Alakazam.zero()

    ###
    # totally antisymmetric tensor
    A1 = AntisymmetricTensor("A", [a => lower, b => lower, c => lower, d => lower], [x])
    A2 = AntisymmetricTensor("A", [d => lower, a => lower, b => lower, c => lower], [x])
    @test A1+A2≅Alakazam.zero()
    ###
    # totally symmetric tensor
    S1 = SymmetricTensor("S", [a => lower, b => lower, c => lower, d => lower], [x])
    S2 = SymmetricTensor("S", [d => lower, a => lower, b => lower, c => lower], [x])
    @test S1-S2≅Alakazam.zero()
    
    
    A = AntisymmetricTensor("A", [a => upper, b => upper, c => upper], [x])
    #1st bianchi identity of Riemann
    @test Alakazam.hard_simplify(Alakazam.young_projection(A) * Alakazam.young_projection(R1)) == Alakazam.zero()


    @test -1== to_dumless_same_tensor_mod_syms(AntisymmetricTensor("A", [a => upper, b => upper, c => upper], [x]),AntisymmetricTensor("A", [a => upper, c => upper, b => upper], [x]))
    @test 1==to_dumless_same_tensor_mod_syms(RiemannTensor("R", [a => upper, b => upper, c => upper,d=>upper], [x]),RiemannTensor("R", [c => upper, d => upper, a => upper,b=>upper], [x]))

    #see if we can compute large tableau w/o hanging
    @test (YoungTableau([4,4,4,4],[1:16;]);true)

    ###
    # mono-term vs multi-term symmetry.
    # Mono-term symmetries are those expressible as a signed permutation of a single monomial's
    # slots; they are what Butler-Portugal (Cadabra's canonicalise, xAct's ToCanonical) solves,
    # and what hard_simplify's graph isomorphism solves independently.
    Ri(i,j,k,l) = Tensor("R", [i=>lower, j=>lower, k=>lower, l=>lower], [riemann_tabelau])
    Rup(i,j,k,l) = Tensor("R", [i=>upper, j=>upper, k=>upper, l=>upper], [riemann_tabelau])

    @test Alakazam.hard_simplify(Ri(a,b,c,d) + Ri(b,a,c,d)) == Alakazam.zero()   # antisym in 1,2
    @test Alakazam.hard_simplify(Ri(a,b,c,d) - Ri(c,d,a,b)) == Alakazam.zero()   # pair exchange
    # the antisym-against-sym zero, found by the 4-cycle search rather than by pairwise comparison
    @test Alakazam.hard_simplify(AntisymmetricTensor("A",[a=>upper,b=>upper])
                                *SymmetricTensor("B",[a=>lower,b=>lower])) == Alakazam.zero()

    # Multi-term symmetries are linear relations across several monomials, each of which is
    # separately non-zero and separately canonical. Butler-Portugal cannot see these, and neither
    # can graph isomorphism -- they need the Young projector, which hard_simplify runs only on
    # request. Note neither identity is a pairwise t == ±u relation: Bianchi ties three monomials
    # together at once, and the contracted form relates two by a factor of 1/2.
    bianchi = Ri(a,b,c,d) + Ri(a,c,d,b) + Ri(a,d,b,c)
    @test Alakazam.hard_simplify(bianchi) != Alakazam.zero()                      # not mono-term
    @test Alakazam.hard_simplify(bianchi; multiterm=true) == Alakazam.zero()

    # a contracted consequence of Bianchi: R_{abcd}R^{acbd} = 1/2 R_{abcd}R^{abcd}
    contracted = Ri(a,b,c,d)*Rup(a,c,b,d) - (1//2)*Ri(a,b,c,d)*Rup(a,b,c,d)
    @test Alakazam.hard_simplify(contracted) != Alakazam.zero()
    @test Alakazam.hard_simplify(contracted; multiterm=true) == Alakazam.zero()
    @test cong(Ri(a,b,c,d)*Rup(a,c,b,d), (1//2)*Ri(a,b,c,d)*Rup(a,b,c,d))         # ≅ opts in

    # the budget is evaluated from tableau shape alone, before any projection is performed
    @test young_projector_size(riemann_tabelau) == 16
    @test young_projector_size(Ri(a,b,c,d)*Rup(a,c,b,d)) == 256
    @test young_projector_size(Tensor("V", [a=>lower])) == 1                      # no tableau
    @test Alakazam.hard_simplify(contracted; multiterm=true, multiterm_budget=8) == contracted

    @test canonicalise(S1-S2)==0
    @test canonicalise(S2)==S1
    @test (sort(S2)==S1);
    @test (S2<S1;true)

    r1 = Tensor("r", [a => lower, b => lower, c => lower, d => lower], [riemann_tabelau])
    r2 = Tensor("h", [a => lower, b => lower, c => lower, d => lower], [riemann_tabelau])
    r3 = Tensor("m", [a => lower, b => lower, c => lower, d => lower], [riemann_tabelau])
    #make sure large sum doesnt crash
    young_projection(R1)+young_projection(r2)+young_projection(r1)+young_projection(r3)
end