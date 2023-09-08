
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


T1=Tensor("T", [n=> upper, m => upper, p => upper, q=> upper], YoungTableau([1,1,1,1],[1,2,3,4]))
T2=Tensor("T", [m=> upper, n => upper, p => upper, q=> upper], YoungTableau([1,1,1,1],[1,2,3,4]))

same_tensor_mod_tableaux(T1,T2)

2*T1+5*T2