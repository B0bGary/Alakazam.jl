# import Pkg; 
# Pkg.add("Symbolics")
# Pkg.add("OrderedCollections")
# Pkg.add("Revise")
# Pkg.add("StatsBase")
# Pkg.add("SymbolicUtils")
# Pkg.add("DataStructures")
using Symbolics
using OrderedCollections
using Revise
using StatsBase
import SymbolicUtils
import Base.:*
import Base.:+
import Base.:/
import Base.:(==)
import Base.:-
import DataStructures
include("indices.jl")
include("custom_exceptions.jl")
include("tensor_term.jl")
include("tensor_ops.jl")
include("metric.jl")
include("operator.jl")
include("operator_ops.jl")
include("derivative.jl")
include("derivative_ops.jl")
include("kronecker.jl")
include("epsilon_tensor.jl")

flavour_indices = IndexSet("flavour", 2)
spinor_indices = IndexSet("spinor", 2)
spacetime_indices = IndexSet("spacetime", 4)
spinor_conj_indices = IndexSet("spinor_conj", 2)
flavour_conj_indices = IndexSet("flavour_conj", 2)


i = Index("i", flavour_indices, 0)
j = Index("j", flavour_indices)
k = Index("k", flavour_indices, 0)
l = Index("l", flavour_indices)
m = Index("m", flavour_indices)
n = Index("n", flavour_indices)
p = Index("p", flavour_indices)
q = Index("q", flavour_indices)

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


x = Coordinate("x", spacetime_indices)
θ = Coordinate("θ", [flavour_indices, spinor_indices])
θbar = Coordinate("̄θ", [flavour_conj_indices, spinor_conj_indices])

# println(i)
# println(j)
# println(flavour_indices)

# D=Tensor("D")
# println(D)
# A=Tensor("A",OrderedDict(i=>upper, j=>upper))
# println(A)
# B=Tensor("B",OrderedDict(j=>lower, k=>lower))
# println(B)
# K=Tensor("K", OrderedDict(i=>lower, j=>lower, k=>lower, l=>upper))
# grading(K)


J = Tensor("J", [n => upper, j => upper, j => lower, g => upper, p => upper, p => lower, e => upper, f => upper, q => lower, q => upper, g => lower, f => lower])
K = Tensor("K", [n => lower, m => upper, m => lower])

F = TensorTerm([J, K])

# println(F)
# println(F.indices)
# println(J.indices)

# println(get_dummies(J))
# println(J)
# rename_dummies(J)
# println(J)

# println(F)
# rename_dummies(F)
# println(F)

J = Tensor("J", [i => upper, k => upper, k => lower, β => lower], x)
K = Tensor("K", [i => upper, j => upper, j => lower, β => lower], θ)
bb = 2 * K - 1 * K - 1 / 2 * K
# println(bb)

# get_free_indices(J)


test_eps = 0
test_delta = 0
test_met = 1
test_deriv = 0


exp= EpsilonTensor([α => upper, k => upper, k => lower, β => lower])

if test_eps == 1

    println("---------------------------------")
    println("Epsilon")
    eps = EpsilonTensor([i => lower, j => lower])
    println("Have an epsilon tensor")
    println(eps)

    println("lets eliminate it from expression")
    H = Tensor("H", [j => upper, k => upper, k => lower, β => lower], x)
    exp = eps * H
    println(exp)
    exp = eliminate_epsilon(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [i => upper, k => upper, k => lower, β => lower], x)
    exp = eps * H
    println(exp)
    exp = eliminate_epsilon(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [i => lower, k => upper, k => lower, β => lower], x)
    exp = eps * H
    println(exp)
    exp = eliminate_epsilon(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [j => lower, k => upper, k => lower, β => lower], x)
    exp = eps * H
    println(exp)
    exp = eliminate_epsilon(exp)
    println(exp)

    eps = EpsilonTensor([i => upper, j => upper])
    println("Have an epsilon tensor")
    println(eps)
    println("lets eliminate it from expression")
    H = Tensor("H", [j => upper, k => upper, k => lower, β => lower], x)
    exp = eps * H
    println(exp)
    exp = eliminate_epsilon(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [i => upper, k => upper, k => lower, β => lower], x)
    exp = eps * H
    println(exp)
    exp = eliminate_epsilon(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [i => lower, k => upper, k => lower, β => lower], x)
    exp = eps * H
    println(exp)
    exp = eliminate_epsilon(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [j => lower, k => upper, k => lower, β => lower], x)
    exp = eps * H
    println(exp)
    exp = eliminate_epsilon(exp)
    println(exp)


end






if test_met == 1

    println("---------------------------------")
    println("Metrics")


    println("Here's a tensor")
    println(J)
    println("Here's me popping out the metric")
    B = pop_metric(J)
    println(B)
    println("Here's me contracting the metric")
    B=contract_metrics(B)
    println(B)
    B=eliminate_epsilon(B)
    println(B)
    println("cash money")


end

# O=Operator("O",J,[a=>upper, b=>lower])
# println("first op is ")
# println(O)

# O=Operator("O",J*K,[a=>upper, b=>lower])
# println("2nd op is ")
# println(O)
# println("doing prod ")
# ex=product_rule(O)
# typeof(ex)
# println(ex)
# println(bb.summands[1][2])
# println(bb.summands[2][2])

#eps=Tensor("eps",[k l],[])
# C=A*B


# println(C)
# println(getFreeIndices(C))



if test_deriv == 1
    println("---------------------------------")
    println("Derivatives")
    D = Derivative("D", J * K, θ, [i => upper, α => lower])
    println("Want to find derivative:")
    println(D)
    println("D is with respect to")
    println(D.coord)
    println("J is a function of")
    println(J.function_of)
    println("K is a function of")
    println(K.function_of)
    println("Applying derivative")
    D = apply_derivative(D)
    println(D)

end





if test_delta == 1

    println("---------------------------------")
    println("Kroncker")
    δ = KroneckerDelta([i => lower, j => upper])
    println("Have a delta")
    println(δ)
    println("lets eliminate it from expression")
    H = Tensor("H", [j => upper, k => upper, k => lower, β => lower], x)
    exp = δ * H
    println(exp)
    exp = eliminate_kronecker(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [i => upper, k => upper, k => lower, β => lower], x)
    exp = δ * H
    println(exp)
    exp = eliminate_kronecker(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [i => lower, k => upper, k => lower, β => lower], x)
    exp = δ * H
    println(exp)
    exp = eliminate_kronecker(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [j => lower, k => upper, k => lower, β => lower], x)
    exp = δ * H
    println(exp)
    exp = eliminate_kronecker(exp)
    println(exp)

    δ = KroneckerDelta([i => upper, j => lower])
    println("Have a delta")
    println(δ)
    println("lets eliminate it from expression")
    H = Tensor("H", [j => upper, k => upper, k => lower, β => lower], x)
    exp = δ * H
    println(exp)
    exp = eliminate_kronecker(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [i => upper, k => upper, k => lower, β => lower], x)
    exp = δ * H
    println(exp)
    exp = eliminate_kronecker(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [i => lower, k => upper, k => lower, β => lower], x)
    exp = δ * H
    println(exp)
    exp = eliminate_kronecker(exp)
    println(exp)
    println("lets eliminate it from expression")
    H = Tensor("H", [j => lower, k => upper, k => lower, β => lower], x)
    exp = δ * H
    println(exp)
    exp = eliminate_kronecker(exp)
    println(exp)


end

