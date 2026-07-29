@testset "Tensor Equalities" begin


    ####setup

    spacetime_indices = IndexSet("spacetimelc", 2, lower)
    flavour_indices = IndexSet("flavour_indices", 2, upper,fermionic)
    lightcone_indices = IndexSet("lightcone", 2, lower)
    spinor_indices=IndexSet("spinor_indices", 2, lower)
    spinor_conj_indices=IndexSet("spinor_conj_indices", 2, lower)
    @indices a:z ∈ flavour_indices
    x = Coordinate("x", [spacetime_indices, spacetime_indices])
    θ = Coordinate("θ", [spacetime_indices])
    θbar = Coordinate("Θ", [spacetime_indices])

    #############

    ####
    #simple tensors
    @test Tensor("A", [i=>upper, j=>upper])≅Tensor("A", [i=>upper, j=>upper]) #should be exactly equal
    @test ! (Tensor("A", [i=>upper, j=>upper])≅Tensor("A", [j=>upper, i=>upper])) #not equal since no syms
    @test SymmetricTensor("A", [i=>upper, j=>upper])≅SymmetricTensor("A", [j=>upper, i=>upper]) # equal since syms
    @test !(SymmetricTensor("A", [i=>upper, j=>upper])≅SymmetricTensor("B", [j=>upper, i=>upper])) # not equal since different names
    @test !(AntisymmetricTensor("A", [i=>upper, j=>upper]) ≅ AntisymmetricTensor("A", [j=>upper, i=>upper])) # not equal since antisym
    @test !(AntisymmetricTensor("A", [i=>upper, i=>lower]) ≅ AntisymmetricTensor("A", [i=>lower, i=>upper])) # not equal since antisym

    ####
    #simple products and contractions
    #should be not equal by antisym
    @test !((AntisymmetricTensor("A", [i=>upper, j=>upper], x)*SymmetricTensor("B", [i=>lower, k=>lower], x)) ≅ (AntisymmetricTensor("A", [j=>upper, i=>upper], x)*SymmetricTensor("B", [i=>lower, k=>lower], x)))
    #should be equal since contractions imply equal to  zero(), even tho they look unequal by antisymmetry.
    @test (AntisymmetricTensor("A", [i=>upper, j=>upper], x)*SymmetricTensor("B", [i=>lower, j=>lower], x))≅(AntisymmetricTensor("A", [j=>upper, i=>upper], x)*SymmetricTensor("B", [i=>lower, j=>lower], x))
    #equal by double antisym
    @test (AntisymmetricTensor("A", [i=>upper, j=>upper], x)*AntisymmetricTensor("B", [i=>lower, j=>lower], x)) ≅ (AntisymmetricTensor("A", [j=>upper, i=>upper], x)*AntisymmetricTensor("B", [j=>lower, i=>lower], x))


    ####
    #testing the grading and swaps 
    #should be equal by grading (even*odd=even)
    @test (Tensor("A", [i=>upper], x)*Tensor("B", [j=>lower, k=>lower], x))≅(Tensor("B", [j=>lower, k=>lower], x) * Tensor("A", [i=>upper], x))
    #should be equal by grading (even*even=even)
    @test (Tensor("A", [i=>upper, t=>upper], x)*Tensor("B", [j=>lower, k=>lower], x))≅(Tensor("B", [j=>lower, k=>lower], x) * Tensor("A", [i=>upper, t=>upper], x))
    #should be not equal by grading (odd*odd=odd) 
    @test ! ((Tensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [j=>lower, k=>lower, l=>lower], x))≅(Tensor("B", [j=>lower, k=>lower, l=>lower], x) * Tensor("A", [i=>upper, t=>upper, s=>upper], x)))

    #should be equal by grading plus antisymm (odd*odd=odd) 
    @test (AntisymmetricTensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [j=>lower, k=>lower, l=>lower], x))≅(Tensor("B", [j=>lower, k=>lower, l=>lower], x) * AntisymmetricTensor("A", [i=>upper, s=>upper, t=>upper], x))
    #should be equal by two graded swaps
    @test (Tensor("A", [i=>upper], x)*Tensor("B", [j=>lower], x)*Tensor("C", [k=>lower], x))≅(Tensor("C", [k=>lower], x)*Tensor("A", [i=>upper], x)*Tensor("B", [j=>lower], x))
    #should be not equal
    @test ! ((Tensor("A", [i=>upper], x)*Tensor("B", [j=>lower], x)*Tensor("C", [k=>lower], x))≅(Tensor("C", [k=>lower], x)*Tensor("B", [j=>lower], x)*Tensor("A", [i=>upper], x)))

    #should be equal by grading plus antisymm thanks to the contraction relabelling
    @test (AntisymmetricTensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [t=>lower], x)*Tensor("B", [s=>lower], x))≅(AntisymmetricTensor("A", [i=>upper, s=>upper, t=>upper], x)*Tensor("B", [s=>lower], x)*Tensor("B", [t=>lower], x))

    #should be equal since equal to  zero() via contraction
    #sort of have to use sort here... i think this is a bad case, not sure how to get around  zero() both sides are really equal to  zero()
    @test (sort(SymmetricTensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [t=>lower], x)*Tensor("B", [s=>lower], x)))≅(SymmetricTensor("A", [i=>upper, t=>upper, s=>upper], x)*Tensor("B", [s=>lower], x)*Tensor("B", [t=>lower], x))

    #should be equal since the tableaux encodes totally symmetric
    @test (Tensor("A", [i=>upper, j => upper, k=>upper, l=>upper], x, YoungTableau([4], [1, 2, 3, 4]))≅Tensor("A", [j=>upper, i => upper, l=>upper, k=>upper], x, YoungTableau([4], [1, 2, 3, 4])))

    @test (1*Tensor("A", [i=>upper, j => upper, k=>upper, l=>upper], x, YoungTableau([1, 1, 1, 1], [1, 2, 3, 4]))≅-Tensor("A", [j=>upper, i => upper, k=>upper, l=>upper], x, YoungTableau([1, 1, 1, 1], [1, 2, 3, 4])))

end