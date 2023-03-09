using Symbolics
using OrderedCollections
using Revise
using StatsBase
import SymbolicUtils
import Base.:*
import DataStructures
include("indices.jl")
include("tensor_term.jl")

flavour_indices=IndexSet("flavour")
spacetime_indices=IndexSet("spacetime")

i=Index("i",flavour_indices)
j=Index("j",flavour_indices)
k=Index("k",flavour_indices,0)
l=Index("l",flavour_indices)
m=Index("m",flavour_indices)
n=Index("n",flavour_indices)
p=Index("p",flavour_indices)
q=Index("q",flavour_indices)


a=Index("a",spacetime_indices)
b=Index("b",spacetime_indices)
c=Index("c",spacetime_indices)
d=Index("d",spacetime_indices)
e=Index("e",spacetime_indices)
f=Index("f",spacetime_indices)
g=Index("g",spacetime_indices)

println(i)
println(j)
println(flavour_indices)

D=Tensor("D")
println(D)
# A=Tensor("A",OrderedDict(i=>upper, j=>upper))
# println(A)
# B=Tensor("B",OrderedDict(j=>lower, k=>lower))
# println(B)
# K=Tensor("K", OrderedDict(i=>lower, j=>lower, k=>lower, l=>upper))
# grading(K)


J=Tensor("J",[n=>upper, j=>upper, j=>lower, g=>upper, p=> upper, p=> lower, e=>upper, f=>upper, q=> lower, q=>upper, g=> lower, f=>lower])
K=Tensor("K",[n=>lower, m=>upper, m=>lower])

F=TensorTerm([J,K])

println(F)
println(F.indices)
println(J.indices)

# println(get_dummies(J))
# println(J)
# rename_dummies(J)
# println(J)

println(F)
rename_dummies(F)
println(F)

#eps=Tensor("eps",[k l],[])
# C=A*B

# println(C)
# println(getFreeIndices(C))