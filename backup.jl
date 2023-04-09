include("alakazam.jl")



#define index sets
flavour_indices = IndexSet("flavour", 2,lower)
spinor_indices = IndexSet("spinor", 2)
spacetime_indices = IndexSet("spacetime", 4,lower)
spinor_conj_indices = IndexSet("spinor_conj", 2,lower)
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

# ic = Index("̄i", flavour_conj_indices, 0)
# jc = Index("̄j", flavour_conj_indices)
# kc = Index("̄k", flavour_conj_indices, 0)
# lc = Index("̄l", flavour_conj_indices)
# mc = Index("̄m", flavour_conj_indices)
# nc = Index("̄n", flavour_conj_indices)
# pc = Index("̄p", flavour_conj_indices)
# qc = Index("̄q", flavour_conj_indices)

α = Index("α", spinor_indices, 1)
β = Index("β", spinor_indices, 1)
γ = Index("γ", spinor_indices, 1)
μ = Index("μ", spinor_indices, 1)
ν = Index("ν", spinor_indices, 1)
ρ = Index("ρ", spinor_indices, 1)
τ = Index("τ", spinor_indices, 1)

ω = Index("ω", spinor_conj_indices, 1)
χ = Index("χ", spinor_conj_indices, 1)
υ = Index("υ", spinor_conj_indices, 1)
ξ = Index("ξ", spinor_indices, 1)
π = Index("π", spinor_conj_indices, 1)
# νdot = Index("̇ν", spinor_conj_indices, 1)
# τdot = Index("̇τ", spinor_indices, 1)

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

# define our grassman basis ̄θ 
θiα=Tensor("θ",[i=>lower, α=>upper],θ)
barθωi=Tensor("̄θ",[ω=>lower, i=>upper],θbar)
barθij=Tensor("̄θ",[i=>lower, ξ=>upper],θbar)*Tensor("θ",[j=>lower, ξ=>lower],θbar)
barθupij=Tensor("̄θ",[i=>upper, υ=>upper],θbar)*Tensor("θ",[j=>upper, υ=>lower],θbar)
θij=Tensor("θ",[i=>lower, ν=>upper],θ)*Tensor("θ",[j=>lower, ν=>lower],θ)
θupij=Tensor("θ",[i=>upper, μ=>upper],θ)*Tensor("θ",[j=>upper, μ=>lower],θ)
θjα=Tensor("θ",[j=>lower, α=>upper],θ)

# define a basic proof of concept derivative
grassman_deriv_normalrep=Derivative("∂",zero, [θ], [j => upper, β => lower])
grassman_deriv_conjrep=Derivative("∂",zero, [θbar], [j => lower, χ => upper])
spatial_deriv=Derivative("∂",zero, [x], [b => lower])
#note that the derivative isn't defined on literally θiα, but on tensors with the same index patterns
add_derivative_action(grassman_deriv_normalrep,θiα,KroneckerDelta([i=>lower,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>lower, α=>lower],θ),KroneckerDelta([i=>lower,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>upper, α=>upper],θ),EpsilonTensor([i=>upper,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>upper, α=>lower],θ),EpsilonTensor([i=>upper,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))

add_derivative_action(grassman_deriv_conjrep,barθαi,KroneckerDelta([i=>upper,j=>lower])*KroneckerDelta([χ=>upper,ω=>lower]))


#define our superfield expansions tensors
# exp=(θjα*θupij)*Tensor("G",[i=>lower, α=>lower],x)
# println(exp)
# SUPERFIELD=(Tensor("j", x)+θiα*Tensor("ψ",[i=>upper, α=>lower],x)+θij*SymmetricTensor("F",[i=>upper, j=>upper],x)
# +(θjα*θupij)*Tensor("G",[i=>lower, α=>lower],x)+θij*θupij*Tensor("H",x)
# +θjα*barθαi*Tensor("K",[j=>upper, α=>lower,i => lower, αdot =>upper],x)
# +Tensor("θ",[k=>lower, α=>upper],θ)*barθupij*Tensor("L",[k=>upper, α =>lower, i=>lower,j => lower],x)
# )
SUPERFIELD=(Tensor("j", x))
covar_expansion=grassman_deriv_normalrep-im*Tensor("σ",[b =>upper, π=>upper, β =>lower])*Tensor("̄θ",[π=>lower, j=>upper],θbar)*spatial_deriv
println("our superfield expansions is: ", SUPERFIELD)

#apply our derivative
D_SUPERFIELD=Derivative("D",SUPERFIELD, [x,θ], [s => upper, τ => lower])
println("We want to compute: ", D_SUPERFIELD)

println("Appyling the derivative:")
out=expand_derivative(D_SUPERFIELD,D_SUPERFIELD,covar_expansion)
out=apply_derivative(out)
println(out)
println("Getting rid of epsilon tensors")
out=contract_metrics(out,EpsilonTensor)
println(out)
println("Eliminating deltas")
out=eliminate_kronecker(out)
println(out)
println("Renaming dummies")
rename_dummies(out)
println(out)
DD_SUPERFIELD=Derivative("∂",out, [θ], [l => upper, ρ => lower])
println("We want to compute: ", DD_SUPERFIELD)
println("Appyling the derivative:")
out=apply_derivative(DD_SUPERFIELD)
println(out)
println("Getting rid of epsilon tensors")
out=contract_metrics(out,EpsilonTensor)
println(out)
println("Eliminating deltas")
out=eliminate_kronecker(out)
println(out)
println("Renaming dummies")
rename_dummies(out)
println(out)
println("Simplifying")
out=simplify(out)
println(out)