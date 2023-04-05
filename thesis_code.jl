include("alakazam.jl")



#define index sets
flavour_indices = IndexSet("flavour", 2)
spinor_indices = IndexSet("spinor", 2)
spacetime_indices = IndexSet("spacetime", 4)
spinor_conj_indices = IndexSet("spinor_conj", 2)
flavour_conj_indices = IndexSet("flavour_conj", 2)


#define indices
i = Index("i", flavour_indices, 0)
j = Index("j", flavour_indices)
k = Index("k", flavour_indices, 0)
l = Index("l", flavour_indices)
m = Index("m", flavour_indices)
n = Index("n", flavour_indices)
p = Index("p", flavour_indices)
q = Index("q", flavour_indices)
r = Index("r", flavour_indices)
s = Index("s", flavour_indices)

ic = Index("̄i", flavour_conj_indices, 0)
jc = Index("̄j", flavour_conj_indices)
kc = Index("̄k", flavour_conj_indices, 0)
lc = Index("̄l", flavour_conj_indices)
mc = Index("̄m", flavour_conj_indices)
nc = Index("̄n", flavour_conj_indices)
pc = Index("̄p", flavour_conj_indices)
qc = Index("̄q", flavour_conj_indices)

α = Index("α", spinor_indices, 1)
β = Index("β", spinor_indices, 1)
γ = Index("γ", spinor_indices, 1)
μ = Index("μ", spinor_indices, 1)
ν = Index("ν", spinor_indices, 1)

αdot = Index("̇α", spinor_conj_indices, 1)
βdot = Index("̇β", spinor_conj_indices, 1)
γdot = Index("̇γ", spinor_conj_indices, 1)
μdot = Index("̇μ", spinor_conj_indices, 1)
νdot = Index("̇ν", spinor_conj_indices, 1)

a = Index("a", spacetime_indices, 0)
b = Index("b", spacetime_indices)
c = Index("c", spacetime_indices)
d = Index("d", spacetime_indices)
e = Index("e", spacetime_indices)
f = Index("f", spacetime_indices)
g = Index("g", spacetime_indices)


#define coordinates, tensors will be functions of these coordinates (for use in derivatives)
x = Coordinate("x", spacetime_indices)
θ = Coordinate("θ", [flavour_indices, spinor_indices])
θbar = Coordinate("̄θ", [flavour_conj_indices, spinor_conj_indices])

# define our grassman basis
θiα=Tensor("θ",[i=>lower, α=>upper],θ)
barθαi=Tensor("θ",[αdot=>lower, ic=>upper],θbar)
θij=Tensor("θ",[i=>lower, α=>upper],θ)*Tensor("θ",[j=>lower, α=>lower],θ)


# define a basic proof of concept derivative
grassman_deriv_normalrep=Derivative("∂",zero, [θ], [j => upper, β => lower])
#note that the derivative isn't defined on literally θiα, but on tensors with the same index patterns
add_derivative_action(grassman_deriv_normalrep,θiα,KroneckerDelta([i=>lower,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>lower, α=>lower],θ),KroneckerDelta([i=>lower,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))


#define our superfield expansions tensors

SUPERFIELD=Tensor("j", x)+Tensor("ψ",[i=>upper, α=>lower],x)*θiα+Tensor("F",[i=>upper, j=>upper],x)*θij
println("our superfield expansions is: ", SUPERFIELD)

#apply our derivative
D_SUPERFIELD=Derivative("∂",SUPERFIELD, [θ], [k => upper, γ => lower])

apply_derivative(D_SUPERFIELD)