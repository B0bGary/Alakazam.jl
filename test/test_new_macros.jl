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

