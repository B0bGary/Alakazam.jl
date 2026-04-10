using Alakazam
zero=Alakazam.zero

#define index sets
println("beginning...")

spacetime_indices = IndexSet("spacetimelc", 2,lower)
flavour_indices = IndexSet("flavour_indices", 2,upper)
lightcone_indices = IndexSet("lightcone", 2,lower)
spinor_indices=IndexSet("spinor_indices", 2,lower)
spinor_conj_indices=IndexSet("spinor_conj_indices", 2,lower)

p = Index("+", lightcone_indices, 1)
m = Index("-", lightcone_indices, 1)

a = Index("a", flavour_indices, 1)
b = Index("b", flavour_indices,1)
c = Index("c", flavour_indices,1)
d = Index("d", flavour_indices,1)
e = Index("e", flavour_indices,1)
f = Index("f", flavour_indices,1)
g = Index("g", flavour_indices,1)
h = Index("h", flavour_indices,1)
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
α = Index("α", flavour_indices,1)
β = Index("β", flavour_indices,1)
γ = Index("γ", flavour_indices,1)   

x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices])
θbar = Coordinate("Θ", [spacetime_indices])

grassman_deriv=Derivative("d",zero, [θ], [i => lower])
grassman_deriv_bar=Derivative("∂",zero, [θbar], [i => lower])
spatial_deriv=Derivative("∇",zero, [x])
minus_deriv=Derivative("δ",zero, [x])


add_derivative_action(grassman_deriv_bar,Tensor("Θ",[j=>lower],θbar),EpsilonTensor([j =>lower, i=>lower]))
add_derivative_action(grassman_deriv,Tensor("θ",[j=>lower],θ),EpsilonTensor([j =>lower, i=>lower]))
add_derivative_action(grassman_deriv_bar,Tensor("Θ",[j=>upper],θbar),KroneckerDelta([i =>lower, j=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[j=>upper],θ),KroneckerDelta([i =>lower, j=>upper]))

println("def derivs...")

covar_expansion=grassman_deriv_bar-im/2*Tensor("θ",[i=>lower],θ)*spatial_deriv
covar_unbar_expansion=-1*grassman_deriv-im/2*Tensor("Θ",[ i=>lower],θbar)*spatial_deriv
covar=Derivative("̄D",zero, [x,θbar], [ i => lower])
covar_unbar=Derivative("D",zero, [x,θ], [ i => lower])


custom_sort_order(["θ","Θ","ϵ"])

f=Tensor("f", [x,θ,θbar])


dbarf=Derivative("̄D",f, [x,θbar], [ k=>lower])
out=expand_derivative(dbarf,covar,covar_expansion)
# println("out:")
# println(to_latex(out))
ddbarf=Derivative("D",out, [x,θ], [ j=>lower])
out=expand_derivative(ddbarf,covar_unbar,covar_unbar_expansion)
println(out)
println(to_latex(out))
out1=apply_derivative(out)
# println("out1:")
# println(to_latex(out1))


df=Derivative("D",f, [x,θ], [ j=>lower])
out=expand_derivative(df,covar_unbar,covar_unbar_expansion)
dbardf=Derivative("̄D",out, [x,θbar], [ k=>lower])
out=expand_derivative(dbardf,covar,covar_expansion)
out2=apply_derivative(out)

println("out1:")
println(to_latex(out1))
println("out2:")    
println(to_latex(out2))

println("anticommutator:")
com=hard_simplify(out1+out2)

println(to_latex(com))

println("commutator:")
acom=hard_simplify(out1-out2)
println(to_latex(acom))

println("D^2")
df=Derivative("D",f, [x,θ], [ j=>lower])
out=expand_derivative(df,covar_unbar,covar_unbar_expansion)
dbardf=Derivative("D",out, [x,θ], [ k=>lower])
out=expand_derivative(dbardf,covar_unbar,covar_unbar_expansion)
out2=apply_derivative(out)
println("out2:")
println(to_latex(out2))
oft=simplify(EpsilonTensor([j=>upper,k=>upper])*out2)
println(oft)
println(to_latex(oft))