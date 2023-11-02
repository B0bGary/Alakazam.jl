

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
t = Index("t", flavour_indices)
u = Index("u", flavour_indices)
v = Index("v", flavour_indices)
y = Index("y", flavour_indices)


 S=SymmetricTensor("S", [n=> lower, m => upper, p => upper, q=> upper, r=> upper])
# A=AntisymmetricTensor("A", [n=> lower, m => upper])
# A.tableaux[1].tableau
S.tableaux[1].tableau
S.tableaux[1].adj_matrix

S=Tensor("T", [n=> lower, m => upper, p => upper, q=> upper, r=> upper,s=>upper,r=>upper,s=>upper], YoungTableau([3,2,2],[1,2,4,5,6,7,8]))
# A=AntisymmetricTensor("A", [n=> lower, m => upper])
# A.tableaux[1].tableau
S.tableaux[1].tableau
S.tableaux[1].adj_matrix


T1=Tensor("T", [n=> lower, p => upper, q => upper,n=>upper,t=>lower,t=>upper,r=>upper,j=>upper], YoungTableau([1,1,1,1,1,1,1,1],[1,2,3,4,5,6,7,8]))
T2=Tensor("T", [t=>upper,n=> lower, q => upper,r=>upper,j=>upper, p => upper,t=>lower,n=>upper], YoungTableau([1,1,1,1,1,1,1,1],[1,2,3,4,5,6,7,8]))

# same_tensor_mod_tableaux(T1,T2)


T1=Tensor("T", [i=>upper, j=>upper], YoungTableau([1,1],[1,2]))
T2=Tensor("T", [i=>upper, j=>upper], YoungTableau([1,1],[1,2]))

println(T1.tableaux[1].adj_matrix)


exp=2*T1+5*T2
# rename_dummies(exp)
# same_tensor_mod_syms(T1,T2)
exp


# T3=Tensor("T", [s=>upper,n=> upper, s => lower,r=>upper,p=>upper, m => upper], YoungTableau([1,1,1,1,1,1],[1,2,3,4,5,6]))

# get_dummies_with_pos_index(T3)