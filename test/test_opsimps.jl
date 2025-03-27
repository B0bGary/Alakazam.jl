@testset "Operand Simplifications" begin
        

    ####setup
    
spacetime_indices = IndexSet("spacetimelc", 2,lower)
flavour_indices = IndexSet("flavour_indices", 2,upper)
lightcone_indices = IndexSet("lightcone", 2,lower)
spinor_indices=IndexSet("spinor_indices", 2,lower)
spinor_conj_indices=IndexSet("spinor_conj_indices", 2,lower)

a = Index("a", flavour_indices, 1)
b = Index("b", flavour_indices,1)
c = Index("c", flavour_indices,1)
d = Index("d", flavour_indices,1)
e = Index("e", flavour_indices,1)
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
x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices])
θbar = Coordinate("Θ", [spacetime_indices])
# grassman_deriv=Derivative("𝛛",zero, [θ], [i => lower])
# grassman_deriv_bar=Derivative("∂",zero, [θbar], [i => lower])
# spatial_deriv=Derivative("∇",zero, [x])
minus_deriv=Derivative("δ",zero, [x])
# add_derivative_action(grassman_deriv_bar,Tensor("Θ",[j=>lower],θbar),EpsilonTensor([i =>lower, j=>lower]))
# add_derivative_action(grassman_deriv,Tensor("θ",[j=>lower],θ),EpsilonTensor([i =>lower, j=>lower]))
# add_derivative_action(grassman_deriv_bar,Tensor("Θ",[j=>upper],θbar),KroneckerDelta([i =>lower, j=>upper]))
# add_derivative_action(grassman_deriv,Tensor("θ",[j=>upper],θ),KroneckerDelta([i =>lower, j=>upper]))
# covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[i=>lower],θ)*spatial_deriv
# covar_unbar_expansion=grassman_deriv-im/2*Tensor("Θ",[ i=>lower],θbar)*spatial_deriv
# covar=Derivative("̄D",zero, [x,θbar], [ i => lower])
# covar_unbar=Derivative("D",zero, [x,θ], [ i => lower])


#############

    ####
    #simple tensors
    a1=Tensor("A",[i=>upper,j=>upper])
    a2=Tensor("A",[k=>upper,l=>upper])

    # exp=EpsilonTensor([k=>upper,l=>upper])*Derivative("δ",a1, [x])+EpsilonTensor([i=>upper,j=>upper])*Derivative("δ",a2, [x])
    e1=EpsilonTensor([i=>lower,j=>lower])*Derivative("δ",a1, [x])
    e1=simplify(e1)
    for smd in e1.summands
        println(smd)
        for term in smd[2].terms
            if typeof(term)<:OperatorSuperType
                println(term.indices)
                println(term.operator_indices)
                println(term.operand_indices)
                println(term.operand.summands[1][2].indices)
            end
        end
    end

    e2=EpsilonTensor([l=>lower,k=>lower])*Derivative("δ",a2, [x])
    e2=simplify(e2)
    for smd in e2.summands
        println(smd)
        for term in smd[2].terms
            if typeof(term)<:OperatorSuperType
                println(term.indices)
                println(term.operator_indices)
                println(term.operand_indices)
                println(term.operand.summands[1][2].indices)
            end
        end
    end


    res=e1-e2
    res=simplify(res)
    println(res)
    # println(e1.operator_indices)
    # println(e1.operand_indices)
    # @test ==Tensor("A",[i=>upper,j=>upper]) #should be exactly equal
    a1=Tensor("A",[i=>lower,j=>lower])
    a2=Tensor("A",[k=>lower,l=>lower])
    e1=EpsilonTensor([i=>lower,j=>lower])*Derivative("δ",a2, [x])+EpsilonTensor([k=>lower,l=>lower])*Derivative("δ",a1, [x])
    e2=simplify(EpsilonTensor([j=>upper,i=>upper])*(EpsilonTensor([l=>upper,k=>upper])*e1))
    println(e2)
end