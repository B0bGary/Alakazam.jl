@testset "Tensor Macros" begin

    ####setup

    flavour_indices = IndexSet("flavour_indices", 2,upper)

    a = Index("a", flavour_indices, 1)
    b = Index("b", flavour_indices,1)
    c = Index("c", flavour_indices,1)
    d = Index("d", flavour_indices,1)
    e = Index("e", flavour_indices,1)
    f = Index("f", flavour_indices,1)
    i = Index("i", flavour_indices,1)
    j = Index("j", flavour_indices,1)
    k = Index("k", flavour_indices,1)
    l = Index("l", flavour_indices,1)
    q = Index("q", flavour_indices,1)
    n = Index("n", flavour_indices,1)
    r = Index("r", flavour_indices,1)
    s = Index("s", flavour_indices,1)
    t = Index("t", flavour_indices,1)
    u = Index("u", flavour_indices,1)
    v = Index("v", flavour_indices,1)
    w = Index("w", flavour_indices,1)
    y = Index("y", flavour_indices,1)
    z = Index("z", flavour_indices,1)

    x = Coordinate("x", [flavour_indices])
    θ = Coordinate("θ", [flavour_indices])
    θbar = Coordinate("Θ", [flavour_indices])


    # Underscore expr::Symbol cases
    @test Tensor("A", [a => lower, b => lower]) == @Tensor A_ab
    @test Tensor("A", [a => upper, b => upper]) == @Tensor A^ab
    @test Tensor("A", [a => upper, b => upper,c=>lower,d=>upper]) == @Tensor A^ab_c^d
    @test Tensor("A", [a => lower, b => upper,c=>lower,d=>upper]) == @Tensor A_a^b_c^d


     @test Tensor("A", [a => lower, b => lower]) ≅ @Tensor "A_ab"

     A=@Tensor "A_[ab]c(d|e|f)" #check that tableau is parsed faithfully
     A2=Tensor("A", [a => lower, b => lower,c=>lower,d => lower, e => lower, f => lower],[YoungTableau([1,1], [1,2]),YoungTableau([2], [4,6])])
     println(A.tableaux)
     println(A2.tableaux)
    #  @test A.tableaux==A2.tableaux


     #to do: equality checking for an array of tableaux, need to define == for them, and somehow check if the arrays are the same even if the tableaux are in different orders etc.
        #need isomorphism to check if the tableaux are the same, current it only checks labels or if one of them can be rearranged using symmetries (even if the other doesnt have)
        #probably encode everything into a single adj matrix??
    B=@Tensor "B_ab{x,θ,θbar}"
    B2=Tensor("B", [a => lower, b => lower],[θbar])
    println(B2)
    println(B.function_of)
    println(x)
    println(θbar)
    println(B2.function_of)
    @test Tensor("B", [a => lower, b => lower],x) ≅ B
    #isoequal should check function_of also...
    
end

@testset "Young Tableaux" begin
    flavour_indices = IndexSet("flavour_indices", 2,upper)

    a = Index("a", flavour_indices, 1)
    b = Index("b", flavour_indices,1)
    c = Index("c", flavour_indices,1)
    d = Index("d", flavour_indices,1)
    e = Index("e", flavour_indices,1)
    f = Index("f", flavour_indices,1)
    i = Index("i", flavour_indices,1)
    j = Index("j", flavour_indices,1)
    k = Index("k", flavour_indices,1)
    l = Index("l", flavour_indices,1)
    q = Index("q", flavour_indices,1)
    n = Index("n", flavour_indices,1)
    r = Index("r", flavour_indices,1)
    s = Index("s", flavour_indices,1)
    t = Index("t", flavour_indices,1)
    u = Index("u", flavour_indices,1)
    v = Index("v", flavour_indices,1)
    w = Index("w", flavour_indices,1)
    y = Index("y", flavour_indices,1)
    z = Index("z", flavour_indices,1)

    x = Coordinate("x", [flavour_indices])
    θ = Coordinate("θ", [flavour_indices])
    θbar = Coordinate("Θ", [flavour_indices])

    ###############################################
    ### begin test cases
    ###############################################

    ##
    #simplest nontrivial tensor
    labels=[1,2,3]
    row_lengths=[2,1]
    t_tabelau=YoungTableau(row_lengths,labels)
    # println(riemann_tabelau)
    T1=Tensor("T", [a => lower, b => lower, c => lower], [t_tabelau])
    T2=Tensor("T", [c => lower, b => lower, a => lower], [t_tabelau]) #should be minus
    T3=Tensor("T", [c => lower, a => lower, b => lower], [t_tabelau]) # plus perm
    @test Alakazam.young_projection(T1)- Alakazam.young_projection(Alakazam.young_projection(T1))==Alakazam.zero #idempotency of projector  
    @test Alakazam.young_projection(T1)+Alakazam.young_projection(T2) ==Alakazam.zero #correct antisym behaviour


    ##
    #Riemann tensor
    labels=[1,3,2,4]
    row_lengths=[2,2]
    riemann_tabelau=YoungTableau(row_lengths,labels)

    #test the bianchi id
    R1=Tensor("R", [a => lower, b => lower, c => lower, d => lower], [riemann_tabelau])
    R2=Tensor("R", [a => lower, c => lower, d => lower, b => lower], [riemann_tabelau]) 
    R3=Tensor("R", [a => lower, d => lower, b => lower, c => lower], [riemann_tabelau])

    @test Alakazam.young_projection(Alakazam.young_projection(R1))-Alakazam.young_projection(R1)==Alakazam.zero #idempotency of projector
    @test Alakazam.young_projection(R1)+Alakazam.young_projection(R2)+Alakazam.young_projection(R3)==Alakazam.zero #Bianchi identity is manifest
    
    
    R4=Tensor("R", [a => lower, b => lower, d => lower, c => lower], [riemann_tabelau]) #antisym in last
    R5=Tensor("R", [d => lower, c => lower, a => lower, b => lower], [riemann_tabelau]) #pair swap and antisym in first two 
    @test Alakazam.young_projection(R1)+Alakazam.young_projection(R4)==Alakazam.zero
    @test Alakazam.young_projection(R1)+Alakazam.young_projection(R5)==Alakazam.zero


end