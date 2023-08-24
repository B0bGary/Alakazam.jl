
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