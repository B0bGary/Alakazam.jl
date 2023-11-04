using Alakazam

######################################
######################################
## setup and definitions
######################################

#define index sets
flavour_indices = IndexSet("flavour", 2,lower)
spinor_indices = IndexSet("spinor", 2)
spacetime_indices = IndexSet("spacetime", 4,lower)
spinor_conj_indices = IndexSet("spinor_conj", 2,lower)

#define indices

# these indices are of grading 0 (bosonic)
i = Index("i", flavour_indices,0)
j = Index("j", flavour_indices)
k = Index("k", flavour_indices)
l = Index("l", flavour_indices)
m = Index("m", flavour_indices)
n = Index("n", flavour_indices)
p = Index("p", flavour_indices)
q = Index("q", flavour_indices)
r = Index("r", flavour_indices)
s = Index("s", flavour_indices)
t = Index("t", flavour_indices)

# these indices are of grading 1 (ie ibjects with an odd number of these are fermionic)
α = Index("α", spinor_indices, 1)
β = Index("β", spinor_indices, 1)
γ = Index("γ", spinor_indices, 1)
μ = Index("μ", spinor_indices, 1)
ν = Index("ν", spinor_indices, 1)
ρ = Index("ρ", spinor_indices, 1)
τ = Index("τ", spinor_indices, 1)
σ = Index("σ", spinor_indices, 1)

αdot = Index("̇α", spinor_conj_indices, 1)
βdot = Index("̇β", spinor_conj_indices, 1)
γdot = Index("̇γ", spinor_conj_indices, 1)
ρdot = Index("̇ρ", spinor_conj_indices, 1)
μdot = Index("̇μ", spinor_conj_indices, 1)
νdot = Index("̇ν", spinor_conj_indices, 1)
τdot = Index("̇τ", spinor_conj_indices, 1)
σdot = Index("̇σ", spinor_conj_indices, 1)


a = Index("a", spacetime_indices, 0)
b = Index("b", spacetime_indices)
c = Index("c", spacetime_indices)
d = Index("d", spacetime_indices)
e = Index("e", spacetime_indices)
f = Index("f", spacetime_indices)
g = Index("g", spacetime_indices)

#define coordinates, tensors will be functions of these, ie tensor fields. these are useful for taking derivatives
x = Coordinate("x", spacetime_indices)
θ = Coordinate("θ", [flavour_indices, spinor_indices])
θbar = Coordinate("Θ", [flavour_indices, spinor_conj_indices])

# define our grassman basis ̄θ for use in superfield expansions 
θiα=Tensor("θ",[i=>lower, α=>upper],θ)
barθαi=Tensor("Θ",[i=>upper,αdot=>lower],θbar)
barθij=Tensor("Θ",[i=>lower, νdot=>upper],θbar)*Tensor("Θ",[j=>lower,
νdot=>lower],θbar)
barθupij=Tensor("Θ",[i=>upper, μdot=>upper],θbar)*Tensor("Θ",
[j=>upper, μdot=>lower],θbar)
θij=Tensor("θ",[i=>lower, ν=>upper],θ)*Tensor("θ",[j=>lower,
ν=>lower],θ)
θupij=Tensor("θ",[i=>upper, μ=>upper],θ)*Tensor("θ",[j=>upper,
μ=>lower],θ)



# θjα=Tensor("θ",[j=>lower, α=>upper],θ)
# θ4=Tensor("θ",[m=>lower, ν=>upper],θ)*Tensor("θ",[n=>lower, 
# ν=>lower],θ)*Tensor("θ",[m=>upper, μ=>upper],θ)*Tensor("θ",[n=>upper,
# μ=>lower],θ)
# θbar4=Tensor("Θ",[p=>lower, νdot=>upper],θbar)*Tensor("Θ",[q=>lower,
# νdot=>lower],θbar)*Tensor("Θ",[p=>upper, μdot=>upper],θbar)*Tensor("Θ",
# [q=>upper, μdot=>lower],θbar)
# θᵅᵝ=Tensor("θ",[y=>upper, α=>upper],θ)*Tensor("θ",[y=>lower,
# β=>upper],θ)
# θᵦᵧ =Tensor("θ",[y=>upper, β=>lower],θ)*Tensor("θ",[y=>lower,
# γ=>lower],θ)
# barθᵅᵝ=Tensor("Θ",[z=>upper, αdot=>upper],θbar)*Tensor("Θ",
# [z=>lower, βdot=>upper],θbar)
# barθᵦᵧ =Tensor("Θ",[z=>upper, βdot=>lower],θbar)*Tensor("Θ",
# [z=>lower, γdot=>lower],θbar)
# θiβ=Tensor("θ",[i=>lower, β=>upper],θ)
# barθjα=Tensor("Θ",[j=>upper,αdot=>lower],θbar)
# barθjk=Tensor("Θ",[j=>lower, νdot=>upper],θbar)*Tensor("Θ",[k=>lower,
# νdot=>lower],θbar)
# θjk=Tensor("θ",[j=>lower, ν=>upper],θ)*Tensor("θ",[k=>lower, ν=>lower],θ)
# θʲᵏ=Tensor("θ",[j=>upper, ν=>upper],θ)*Tensor("θ",[k=>upper, ν=>lower],θ)





# defining the derivative operators of our theory
# each type of derivative should be given a different symbol
# note that zero here is just the placeholder 0 tensor expression and will be replaced when
# the derivative is evaluated later
# the 3rd parameter is which coordinates the derivative is w/r to (it will annihilate tensors that are not
# functions of these coordinates)
# the 4th parameter are the indices associated with the derivative
grassman_deriv_normalrep=Derivative("∂",Alakazam.zero, [θ], [j => upper, β => lower])
grassman_deriv_conjrep=Derivative("𝛛",Alakazam.zero, [θbar], [j => lower, βdot => upper])
spatial_deriv=Derivative("∇",Alakazam.zero, [x], [b => lower])

#define the action of our derivatives on some tensors
#note that the derivative isn't defined on literally "θiα", but on tensors, with the same name and index structure (mod index names)
# its important that the action is defined on a tensor with different indices to the deriv indices, else the 
# software will interpret the action to only be for a contraction
add_derivative_action(grassman_deriv_normalrep,θiα,KroneckerDelta([i=>lower,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>lower, α=>lower],θ),KroneckerDelta([i=>lower,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>upper, α=>upper],θ),EpsilonTensor([i=>upper,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>upper, α=>lower],θ),EpsilonTensor([i=>upper,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))
add_derivative_action(grassman_deriv_conjrep,barθαi,KroneckerDelta([i=>upper,j=>lower])*KroneckerDelta([βdot=>upper,αdot=>lower]))
add_derivative_action(grassman_deriv_conjrep,Tensor("Θ",[i=>lower,αdot=>lower],θbar),EpsilonTensor([i=>lower,j=>lower])*KroneckerDelta([βdot=>upper,αdot=>lower]))
add_derivative_action(grassman_deriv_conjrep,Tensor("Θ",[i=>upper,αdot=>upper],θbar),KroneckerDelta([i=>upper,j=>lower])*EpsilonTensor([αdot=>upper,βdot=>upper]))
add_derivative_action(grassman_deriv_conjrep,Tensor("Θ",[i=>lower,αdot=>upper],θbar),EpsilonTensor([i=>lower,j=>lower])*EpsilonTensor([αdot=>upper,βdot=>upper]))


#create our superfield expansion (just a few terms to illustrate)

SUPERFIELD=(Tensor("j", x)+θiα*Tensor("ψ",[i=>upper, α=>lower],x)+θij*Tensor("F",[i=>upper, j=>upper],x,[YoungTableau([2],[1,2])]))

# note that here, F is a symmetric tensor. To encode this information, we use a vector of Young tableaux (since our tensor might be reducible rep).
# a tableau is a tiling pattern of squares. Squares in the same row denote indices that are symmetric,
# squares in the same column are antisymmetric. 
# When creating a tableau, the first vector is the length of each row. The 2nd vector is labeling the squares with indices,
# left to right, top to bottom. 

# For example, a totally symmetric tesnor, eg S^{μνστ} would have a Young tableau
# YoungTableau([4],[1,2,3,4]):
# [1][2][3][4]

# a totally antisym tensor A^{μνστ} can be encoded in a tableau YoungTableau([1,1,1,1],[1,2,3,4]):
# [1]
# [2]
# [3]
# [4]

# tensor M^{μνστ}  that is symmetric in the first two, symmetric in the last two, but antisymmetric under swaps 1<->3, and 2<->4
# can be encoded in a tableau YoungTableau([2,2],[1,2,3,4]):
# [1][2]
# [3][4]


custom_sort_order(["θ","Θ","ϵ"])


######################################
######################################
## computing the constraint equations
######################################

# to make expressions simple, we can have a cov. derivative eg D, but then expand it out in terms of the full expression
covar=Derivative("D",Alakazam.zero, [x,θ], [j => upper, β => lower])
covar_expansion=grassman_deriv_normalrep-im*Tensor("σ",[b =>upper, γdot=>upper, β =>lower])*Tensor("Θ",[j=>upper,γdot=>lower],θbar)*spatial_deriv

covar_bar=Derivative("̄D",Alakazam.zero, [x,θbar], [j => lower, βdot => upper])
covar_expansion_bar=grassman_deriv_conjrep-im*Tensor("σ",[b =>upper,βdot=>upper, β =>lower])*Tensor("θ",[j=>lower,β=>upper],θ)*spatial_deriv

println("our superfield expansions is: ", SUPERFIELD)

#apply our derivative
# in this example, we'll compute the 2nd derivative, but symmetrised in indices
 println("Computing Symmetrised deriv: ")

D_SUPERFIELD=Derivative("D",SUPERFIELD, [x,θbar], [s => upper, τ => lower])
out=expand_derivative(D_SUPERFIELD,covar,covar_expansion)
out=apply_derivative(out)
out=simplify(out)
DD_SUPERFIELD=Derivative("D",out,[x,θbar], [r => upper, ρ => lower])
out=expand_derivative(DD_SUPERFIELD,covar,covar_expansion)
out=apply_derivative(out)
out=EpsilonTensor([τ=>upper,ρ=>upper])*out
out=eliminate_epsilon(out)

D_SUPERFIELD2=Derivative("D",SUPERFIELD, [x,θbar], [r => upper, τ => lower])
out2=expand_derivative(D_SUPERFIELD2,covar,covar_expansion)
out2=apply_derivative(out2)
out2=simplify(out2)
DD_SUPERFIELD2=Derivative("D",out2,[x,θbar], [s => upper, ρ => lower])
out2=expand_derivative(DD_SUPERFIELD2,covar,covar_expansion)
out2=apply_derivative(out2)
out2=EpsilonTensor([τ=>upper,ρ=>upper])*out2
out2=eliminate_epsilon(out2)


println("Result:")
res2=simplify(out+out2)/2
println(res2)


#example extraction of a constraint

#theta=2, thetabar=0 term in expansion:
tt=2
tb=0
constraint=rename_dummies(get_weighted_terms(get_weighted_terms(res2,
"θ",tt),"Θ",tb))
#get the constraint on the 2nd independant structure associated with θ^2
constraint=project_component(constraint,θ,[i=>upper, α=>lower],tt,tb,θbar,
[i=>lower, αdot=>upper],2,1)
constraint=rename_dummies(contract_metrics(simplify(
eliminate_epsilon(eliminate_kronecker(constraint)))))


######################################
######################################
## compute the supercurrent squared
######################################


SC²=SUPERFIELD × SUPERFIELD
out=get_weighted_terms(get_weighted_terms(SC²,"θ",4),"Θ",4)
#integrate over superspace to obtain the supercurrent-squared term
ttb=project_component(out,θ,[i=>upper, α=>lower],4,4,θbar,[i=>lower,
αdot=>upper],1,1)
println(ttb)