using Alakazam
using Test
zero=Alakazam.zero

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

grassman_deriv=Derivative("𝛛",zero, [θ], [i => lower])
grassman_deriv_bar=Derivative("∂",zero, [θbar], [i => lower])
spatial_deriv=Derivative("∇",zero, [x])
minus_deriv=Derivative("δ",zero, [x])


add_derivative_action(grassman_deriv_bar,Tensor("Θ",[j=>lower],θbar),EpsilonTensor([i =>lower, j=>lower]))
add_derivative_action(grassman_deriv,Tensor("θ",[j=>lower],θ),EpsilonTensor([i =>lower, j=>lower]))
add_derivative_action(grassman_deriv_bar,Tensor("Θ",[j=>upper],θbar),KroneckerDelta([i =>lower, j=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[j=>upper],θ),KroneckerDelta([i =>lower, j=>upper]))


covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[i=>lower],θ)*spatial_deriv
covar_unbar_expansion=grassman_deriv-im/2*Tensor("Θ",[ i=>lower],θbar)*spatial_deriv
covar=Derivative("̄D",zero, [x,θbar], [ j => lower])
covar_unbar=Derivative("D",zero, [x,θ], [ j => lower])

custom_sort_order(["θ","Θ","ϵ"])

SUPERFIELD_S=(Tensor("Θ", x)
+Tensor("θ",[ i=> lower],θ)*Tensor("A",[ i=>upper],x)
+Tensor("Θ",[ i=> lower],θbar)*Tensor("̄A",[ i=>upper],x)
+Tensor("θ",[ i=> lower],θ)*Tensor("θ",[ i=> upper],θ)*Tensor("B",x)
+Tensor("Θ",[ i=> lower],θbar)*Tensor("Θ",[i=> upper],θbar)*Tensor("̄B",x)
+Tensor("θ",[ i=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("C",[  i=>upper,j=>upper],x)
+Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ i=> lower],θbar)*Tensor("S",[   i=>upper],x)
+Tensor("θ",[ i=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[ j=> upper],θbar)*Tensor("̄S",[   i=>upper],x)
+Tensor("θ",[ i=> lower],θ)*Tensor("θ",[ i=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("Θ",[ k=> upper],θbar)*Tensor("T",x)
)


write_to_file("test",SUPERFIELD_S)

deser=read_from_file("test")


@test simplify(SUPERFIELD_S)== deser