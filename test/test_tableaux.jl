@testset "Young Tableaux" begin

    #setup
    if 1==1
        flavour_indices = IndexSet("flavour_indices", 2, upper)
        a = Index("a", flavour_indices, 1)
        b = Index("b", flavour_indices, 1)
        c = Index("c", flavour_indices, 1)
        d = Index("d", flavour_indices, 1)
        e = Index("e", flavour_indices, 1)
        f = Index("f", flavour_indices, 1)
        i = Index("i", flavour_indices, 1)
        j = Index("j", flavour_indices, 1)
        k = Index("k", flavour_indices, 1)
        l = Index("l", flavour_indices, 1)
        q = Index("q", flavour_indices, 1)
        n = Index("n", flavour_indices, 1)
        r = Index("r", flavour_indices, 1)
        s = Index("s", flavour_indices, 1)
        t = Index("t", flavour_indices, 1)
        u = Index("u", flavour_indices, 1)
        v = Index("v", flavour_indices, 1)
        w = Index("w", flavour_indices, 1)
        y = Index("y", flavour_indices, 1)
        z = Index("z", flavour_indices, 1)
        x = Coordinate("x", [flavour_indices])
        θ = Coordinate("θ", [flavour_indices])
        θbar = Coordinate("Θ", [flavour_indices])
    end
    ###############################################
    ### begin test cases
    ###############################################
    
    ###
    ###
    #simplest nontrivial tensor symmetry, irrep of S_3
    ###
    ###
    labels = [1, 2, 3]
    row_lengths = [2, 1]
    t_tabelau = YoungTableau(row_lengths, labels)
    T1 = Tensor("T", [a => lower, b => lower, c => lower], [t_tabelau])
    T2 = Tensor("T", [c => lower, b => lower, a => lower], [t_tabelau])
    T3 = Tensor("T", [b => lower, a => lower, c => lower], [t_tabelau])
    # @show Alakazam.young_projection(T1)+Alakazam.young_projection(T2)
    T4 = Tensor("T", [c => lower, a => lower, b => lower], [t_tabelau])
    T5 = Tensor("T", [b => lower, c => lower, a => lower], [t_tabelau])
    T6 = Tensor("T", [a => lower, c => lower, b => lower], [t_tabelau])
    @test Alakazam.young_projection(T1) - Alakazam.young_projection(Alakazam.young_projection(T1)) == Alakazam.zero #idempotency of projector  
    @test Alakazam.young_projection(T1 + T4 + T5) == Alakazam.zero #bianchi like identity

    ###
    ###
    #Riemann tensor
    ###
    ###
    labels = [1, 3, 2, 4]
    row_lengths = [2, 2]
    riemann_tabelau = YoungTableau(row_lengths, labels)
    @show riemann_tabelau


    R1 = Tensor("R", [a => lower, b => lower, c => lower, d => lower], [riemann_tabelau])
    @show Alakazam.young_projection(R1)
    @test Alakazam.young_projection(Alakazam.young_projection(R1)) - Alakazam.young_projection(R1) == Alakazam.zero #idempotency of projector

    #test the bianchi id
    R2 = Tensor("R", [a => lower, c => lower, d => lower, b => lower], [riemann_tabelau])
    R3 = Tensor("R", [a => lower, d => lower, b => lower, c => lower], [riemann_tabelau])
    @test Alakazam.young_projection(R1) + Alakazam.young_projection(R2) + Alakazam.young_projection(R3) == Alakazam.zero #Bianchi identity is manifest

    #swapping syms
    R4 = Tensor("R", [a => lower, b => lower, d => lower, c => lower], [riemann_tabelau]) #antisym in last
    R5 = Tensor("R", [d => lower, c => lower, a => lower, b => lower], [riemann_tabelau]) #pair swap and antisym in first two 
    @test Alakazam.young_projection(R1) + Alakazam.young_projection(R4) == Alakazam.zero
    @test Alakazam.young_projection(R1) + Alakazam.young_projection(R5) == Alakazam.zero

    #check that built in Riemann tensor has correct symmetries
    Riem = RiemannTensor([d => lower, c=> lower, b => lower, a => lower])
    @test Alakazam.young_projection(Riem) - Alakazam.young_projection(R1) == Alakazam.zero 



    ###
    ###
    # totally antisymmetric tensor
    ### 
    ###
    A1 = AntisymmetricTensor("A", [a => lower, b => lower, c => lower, d => lower], [x])
    @show Alakazam.young_projection(A1)
    A2 = AntisymmetricTensor("A", [d => lower, a => lower, b => lower, c => lower], [x])
    @test Alakazam.young_projection(A1) + Alakazam.young_projection(A2) == Alakazam.zero #check that antisymmetric tensor is antisymmetric

    ###
    ###
    # totally symmetric tensor
    ### 
    ###
    S1 =SymmetricTensor("S", [a => lower, b => lower, c => lower, d => lower], [x])
    S2 =SymmetricTensor("S", [d => lower, a => lower, b => lower, c => lower], [x])
    @show S1.tableaux[1]
    @show S2.tableaux[1]
    @show Alakazam.compute_young_projector(S1.tableaux[1])
    @show Alakazam.young_projection(S1)
    @show Alakazam.young_projection(S2)
    @test Alakazam.young_projection(S1) - Alakazam.young_projection(S2) == Alakazam.zero #check that symmetric tensor is symmetric

    # A= AntisymmetricTensor("A", [a => upper, b => upper, c => upper], [x])
    # @show Alakazam.simplify(Alakazam.young_projection(A)*Alakazam.young_projection(R1))
    #want to show a nontrivial identity here following from ricci identity...
    # can't tell if simplify is wrong, or whether this is an edge case where something doesn't eval to zero
end