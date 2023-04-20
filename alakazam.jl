# import Pkg; 
# Pkg.add("Symbolics")
# Pkg.add("OrderedCollections")
# Pkg.add("Revise")
# Pkg.add("StatsBase")
# Pkg.add("SymbolicUtils")
# Pkg.add("DataStructures")
using Symbolics
using OrderedCollections
# using Revise
using StatsBase
import SymbolicUtils
import Base.:*
import Base.:+
import Base.:copy
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
include("tensor_replacements.jl")


sort_order=Dict()

function custom_sort_order(v::Vector{String})
    global sort_order
    sort_order = Dict(zip(v, 1:length(v)))
end

# flavour_indices = IndexSet("flavour", 2)
# spinor_indices = IndexSet("spinor", 2)
# spacetime_indices = IndexSet("spacetime", 4)
# spinor_conj_indices = IndexSet("spinor_conj", 2)
# flavour_conj_indices = IndexSet("flavour_conj", 2)


# i = Index("i", flavour_indices, 0)
# j = Index("j", flavour_indices)
# k = Index("k", flavour_indices, 0)
# l = Index("l", flavour_indices)
# m = Index("m", flavour_indices)
# n = Index("n", flavour_indices)
# p = Index("p", flavour_indices)
# q = Index("q", flavour_indices)
# r = Index("r", flavour_indices)
# s = Index("s", flavour_indices)

# ic = Index("̄i", flavour_conj_indices, 0)
# jc = Index("̄j", flavour_conj_indices)
# kc = Index("̄k", flavour_conj_indices, 0)
# lc = Index("̄l", flavour_conj_indices)
# mc = Index("̄m", flavour_conj_indices)
# nc = Index("̄n", flavour_conj_indices)
# pc = Index("̄p", flavour_conj_indices)
# qc = Index("̄q", flavour_conj_indices)

# α = Index("α", spinor_indices, 1)
# β = Index("β", spinor_indices, 1)
# γ = Index("γ", spinor_indices, 1)
# μ = Index("μ", spinor_indices, 1)
# ν = Index("ν", spinor_indices, 1)

# αdot = Index("̇α", spinor_conj_indices, 1)
# βdot = Index("̇β", spinor_conj_indices, 1)
# γdot = Index("̇γ", spinor_conj_indices, 1)
# μdot = Index("̇μ", spinor_conj_indices, 1)
# νdot = Index("̇ν", spinor_conj_indices, 1)

# a = Index("a", spacetime_indices, 0)
# b = Index("b", spacetime_indices)
# c = Index("c", spacetime_indices)
# d = Index("d", spacetime_indices)
# e = Index("e", spacetime_indices)
# f = Index("f", spacetime_indices)
# g = Index("g", spacetime_indices)


# x = Coordinate("x", spacetime_indices)
# θ = Coordinate("θ", [flavour_indices, spinor_indices])
# θbar = Coordinate("̄θ", [flavour_conj_indices, spinor_conj_indices])

# # println(i)
# # println(j)
# # println(flavour_indices)

# # D=Tensor("D")
# # println(D)
# # A=Tensor("A",OrderedDict(i=>upper, j=>upper))
# # println(A)
# # B=Tensor("B",OrderedDict(j=>lower, k=>lower))
# # println(B)
# # K=Tensor("K", OrderedDict(i=>lower, j=>lower, k=>lower, l=>upper))
# # grading(K)


# J = Tensor("J", [n => upper, j => upper, j => lower, g => upper, p => upper])
# K = Tensor("K", [n => lower, m => upper, m => lower])

# F = TensorTerm([J, K])

# # println(F)
# # println(F.indices)
# # println(J.indices)

# # println(get_dummies(J))
# # println(J)
# # rename_dummies(J)
# # println(J)

# # println(F)
# # rename_dummies(F)
# # println(F)

# J = Tensor("J", [i => upper, k => upper, j => lower, β => lower], x)
# K = Tensor("K", [n => upper, m => upper, m => lower, α => lower], θ)
# bb = 2 * K - 1 * K - 1 / 2 * K
# # println(bb)

# # get_free_indices(J)


test_eps = 0
test_delta = 0
test_met = 0
test_deriv = 0
test_contract =0 
test_replace =0
test_newtimes =0
test_strings=0
test_rep2=0

test_deriv_exp=0


test_deriv_registry=0

println("executing")

if test_deriv_registry==1
    # F = 2*Tensor("F", [i => upper, k => upper, j => lower, β => lower], θ)*Tensor("G", [p => upper, q => upper, γ => lower], θ)

    delta=KroneckerDelta([n => upper, k=>lower])
    delta2=KroneckerDelta([γ => upper, μ=>lower])
    theta=Tensor("θ", [n => upper, μ=>lower],θ)
    new_deriv_D = Derivative("D",zero, [θ,x], [k => lower, γ => upper])
    add_derivative_action(new_deriv_D,theta,delta*delta2)
    println(derivative_actions)

    act_theta=theta=Tensor("θ", [i => upper, α=>lower],θ)
    D=Derivative("D",act_theta, [θ,x], [j => lower, β => upper])
    res=apply_derivative(D)
    println(res)
    δ2=KroneckerDelta([αdot => upper, βdot=>lower])
    println(δ2)
end





if test_deriv_exp==1

    F = 2*Tensor("F", [i => upper, k => upper, j => lower, β => lower], θ)*Tensor("G", [p => upper, q => upper, γ => lower], θ)
    D = Derivative("D",F, θ, [i => upper, α => lower])
    K =Tensor("K", [ α => lower, β => upper], θ)
    expression=2*D
    println(expression)
    new_deriv=Derivative("∇",zero, θ, [i => upper, α => lower])+K*Derivative("∂",zero, θ, [i => upper, β => lower])
    out=expand_derivative(expression,D,new_deriv)
    println(out)
    println(product_rule(out))

end





if test_rep2==1
    # defining patterns to swap indices
    H=DummyPatternTensor("?H")
    # println(term2histogram(F*H))
    # println(term2histogram(F*H× F))

    J = Tensor("F", [i => upper, k => upper, j => lower, β => lower], x)
    Z = Tensor("Z", [i => upper, k => upper, j => lower, β => lower], x)
    K = Tensor("K", [n => upper, m => upper, m => lower, α => lower], θ)
    L = Tensor("L", [q => upper, p => upper], θ)
    # println(L.indices)
    H=DummyPatternTensor("?H",K.indices)
    J2 = DummyPatternTensor("?D", [i => upper, k => upper, j => lower, DummyPatternIndex("γ") => lower], x)
    
    J3 = DummyPatternTensor("?D", [i => upper, k => upper, j => lower, μ => lower], x)

    println("Consider the tensor expression:")
    # println(L)
    exp1=K*J*Z*L
    println(exp1)
    println("Now, suppose we want to replace tensors of the form")
    
    pat1=J*J2
    println(pat1)
    println("where '?' denotes that the next character is a wildcard")
    println("Then, suppose that we want to replace such matching tensors with an expression")
    rep= TensorTerm(J2)
    println(rep)
    # a=term2histogram(J*J*K) #exp
    # b=term2histogram(J*K*J) #pat
    # println(term2histogram(pat1))
    # println(term2histogram(exp1))
    # println(patt_match(term2histogram(pat1),term2histogram(exp1),pat1,exp1))
    println("Then, our pattern matching algorithm yields:")

    outt=substitute_terms(exp1,pat1,rep,true)
    println(outt)
    # tens=outt.summands[1][2].terms[1]
    # println(tens)
    # println(typeof(tens))

    # dum=parse_pattern2tensor("?A^{i**}_{*?l?m}?B_{mn}")
    # rep=parse_pattern2tensor("?A^{i**}_{*?m?l}")

    # println(typeof(dum))
    # println(dum)

    # println(dum)
    # println(typeof(dum))
    # println(dum.indices)

    # K=Tensor("K", [i => upper, j => lower, k => lower])
    # J=Tensor("J", [i => upper, j => lower, k => lower])
    # println(K)

    # out=sub_tensor(K,dum,rep) # should map to Kikj
    # println(out)
    # println(J)
    # out=sub_tensor(J,DummyTensor(),K)
    # println(out)


end






order=[]


if test_strings==1

    println(term2strpos(K × J))

end


if test_newtimes==1

    println(J)
    println(J × J)
     ex2=(J + K)
    # println(J)
    # println("sum:")
    # println(ex2)
    #  println(J × ex2)
     println(ex2 × ex2)
    # println(J  × (J × (J × J)))

    # println((J  × J) × (J × J))
    # println(1*(J  × J × J) × (J × J))
    # println(1*(J  × J) × (1*(J × J)))
    # println((J  × J) × (1*(J × J)))
end


if test_replace==1

    ex2=EpsilonTensor([ γ => upper, β => upper]) * KroneckerDelta([ γ => lower, α => upper])*J
    println(ex2)
   # ex2=replace(ex2, J,K)
    println(ex2)
end



if test_contract==1

    ex2=EpsilonTensor([ γ => upper, β => upper]) * KroneckerDelta([ γ => lower, α => upper])*J
    # println(ex2)
    # println(get_free_indices(ex2))
     exp= EpsilonTensor([ μ => upper, ν => upper]) * (EpsilonTensor([α => upper, β => upper])*J) +(ex2)*EpsilonTensor([ μ => upper, ν => upper])
    # ######## something is broken...
    println(exp)
    #  println(get_free_indices(exp))
    # println(ex2)
    # println(get_free_indices(ex2))
     exp=eliminate_epsilon(exp) ############### its editing in place!!!!! noooooooo
     println(exp)
     exp=eliminate_kronecker(exp)
    
    # println(exp)
    # exp=eliminate_epsilon(exp)
    println(exp)
    
end

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

