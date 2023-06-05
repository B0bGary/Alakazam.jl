include("superspace_old.jl")

tt=0
tb=0
#println(res)
t11=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
println(t11)
#     # println("t:",w1,", tbar:",w2)
#     merge_sum(t11)

#     t11=merge_sum(sort(rename_dummies(eliminate_epsilon(rename_dummies(sort(pop_metric(t11)))))))

#     t11=merge_sum(eliminate_kronecker(t11))
#     # for summand in t11.summands
#     #     println(summand)
#     #     println(get_free_indices(summand[2]))
#     # end
#     t11=rename_dummies(sort(merge_sum(sort(pop_metric(t11)))))
#t11=θab*Tensor("N",[r=>upper, f=>upper,  a=>upper, b=>upper],x)+t11+θij*Tensor("P",[r=>upper, f=>upper,  i=>upper,j=>upper],x)
# t11=t11+θ3*Tensor("Q",[r=>upper, f=>upper,  b=>upper, j=>upper],x)
# t11=t11+θ4*Tensor("R",[r=>upper, f=>upper],x)
#   t11=t11+θ3*Tensor("S",[  b=>upper, j=>upper],x)
 #t11+=θia*Tensor("ψ",[   i=>upper, a=>upper],x)
# t11+=θ4*Tensor("Tp",x)
# # # println(sort(rename_dummies(t11)))

 t11=project_component(t11,θ,[i=>upper, a=>upper],tt,tb)[1]



#thet2::Coordinate=nothing,deriv2inds::Vector{Pair{Index,IndexPosition}}=nothing, theta1lev::Integer=1, theta2lev::Integer=1)
    # println(exp)

println(t11)
# t11=t11[1]
# # # println("\$\$",to_latex(eliminate_epsilon(t11)),"=0\$\$")
#  t11=rename_dummies(pop_metric(t11))
#   t11=rename_dummies(sort(t11))
# t11=rename_dummies(sort(t11))
# #  println(t11)
#    t11=strip_thetas(t11)
#    t11=eliminate_kronecker(t11)
#  t11=eliminate_epsilon(t11)
#  println(t11)
#  t11=rename_dummies(t11)
 #t11=merge_sum(eliminate_epsilon(eliminate_kronecker(KroneckerDelta([c=>upper, g=>lower])*(KroneckerDelta([d=>upper, f=>lower])*((EpsilonTensor([r=>lower, s=>lower])*(EpsilonTensor([k=>upper, n=>upper])*t11)))))))
#t11=simplify(merge_sum(eliminate_epsilon(EpsilonTensor([πdot =>lower, ωdot=> lower])*(EpsilonTensor([C =>lower, y=> lower])*(EpsilonTensor([A =>lower, u=> lower])*(EpsilonTensor([τdot =>lower, σdot=> lower])*(EpsilonTensor([s =>lower, t=> lower])*(EpsilonTensor([σ =>upper, τ=> upper])*(EpsilonTensor([r =>lower, v=> lower])*eliminate_kronecker(merge_sum(t11)))))))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(EpsilonTensor([ πdot=>lower, ωdot=> lower])*(EpsilonTensor([ B=>lower, v=> lower])*(EpsilonTensor([ τdot=>lower, σdot=> lower])*(EpsilonTensor([ s=>lower, q=> lower])*(EpsilonTensor([ z=>lower, t=> lower])*eliminate_kronecker(merge_sum(t11))))))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker(EpsilonTensor([πdot=>lower, ωdot=>lower])*(EpsilonTensor([y=>lower, q=>lower])*(EpsilonTensor([τdot=>lower, σdot=>lower])*(EpsilonTensor([u=>lower, A=>lower])*merge_sum(t11)))))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(EpsilonTensor([r=>lower, t=>lower])*(EpsilonTensor([s=>lower, p=>lower])*(EpsilonTensor([ρ=>upper, ν=>upper])*(EpsilonTensor([q=>lower, u=>lower])*merge_sum(t11))))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker(EpsilonTensor([t=>lower, z=>lower])*(EpsilonTensor([πdot=>lower, ωdot=>lower])*(EpsilonTensor([s=>lower, q=>lower])*(EpsilonTensor([τdot=>lower, σdot=>lower])*(EpsilonTensor([v=>lower, B=>lower])*merge_sum(t11))))))))))
# println(t11)
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(EpsilonTensor([C=>lower, y=>lower])*(EpsilonTensor([ωdot=>lower, πdot=>lower])*(EpsilonTensor([τdot=>lower, σdot=>lower])*(EpsilonTensor([u=>lower, A=>lower])*(EpsilonTensor([s=>lower, t=>lower])*(EpsilonTensor([ρ=>upper, ν=>upper])*(EpsilonTensor([v=>lower, r=>lower])*merge_sum(t11)))))))))))
# t11=merge_sum(simplify(merge_sum(eliminate_kronecker((KroneckerDelta([s=>upper, l=>lower])*merge_sum(t11))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon((EpsilonTensor([m=>lower, p=>lower])*merge_sum(t11))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker((KroneckerDelta([s=>upper, n=>lower])*(KroneckerDelta([r=>upper, q=>lower])*(EpsilonTensor([m=>lower, p=>lower])*(EpsilonTensor([μdot=>lower, νdot=>lower])*merge_sum(t11))))))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker(EpsilonTensor([r=>lower, t=>lower])*(EpsilonTensor([u=>lower, q=>lower])*((EpsilonTensor([s=>lower, p=>lower])*(EpsilonTensor([ν=>upper, ρ=>upper])*merge_sum(t11))))))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker((EpsilonTensor([u=>lower, q=>lower])*((EpsilonTensor([s=>lower, p=>lower])*merge_sum(t11)))))))))
#t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker((KroneckerDelta([s=>upper, q=>lower])*(KroneckerDelta([r=>upper, u=>lower])*(EpsilonTensor([v=>lower, p=>lower])*(EpsilonTensor([σ=>upper, ω=>upper])*(EpsilonTensor([μdot=>lower, νdot=>lower])*merge_sum(t11)))))))))))

# t11=merge_sum(simplify(merge_sum(eliminate_epsilon((eliminate_kronecker(KroneckerDelta([n=>upper,s =>lower])*(KroneckerDelta([e=>upper, g=>lower])*merge_sum(t11))))))))
# t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker(EpsilonTensor([a=>upper, e=>upper])*((EpsilonTensor([j=>upper, n=>upper])*(EpsilonTensor([b=>upper, g=>upper])*(EpsilonTensor([i=>upper, t=>upper])*merge_sum(t11))))))))))


# t11=contract_metrics(t11)
# t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker(EpsilonTensor([v=>upper, t=>upper])*((EpsilonTensor([h=>upper, e=>upper])*(EpsilonTensor([g=>upper, a=>upper])*(EpsilonTensor([i=>upper, n=>upper])*merge_sum(t11))))))))))

#  t11=merge_sum(simplify(merge_sum(eliminate_epsilon(eliminate_kronecker((EpsilonTensor([F=>upper, D=>upper])*(EpsilonTensor([X2=>upper, H=>upper])*(EpsilonTensor([B=>upper, v=>upper])*(EpsilonTensor([z=>upper, Y=>upper])*merge_sum(t11))))))))))
# t11=1/16*t11

# #n=4,0 2,0
# # t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(EpsilonTensor([b=>upper,c=>upper])*(EpsilonTensor([j=>upper,l=>upper])*eliminate_kronecker(t11)))))
# t11=(merge_sum(eliminate_kronecker(rename_dummies((t11)))))
#t11=merge_sum(simplify(eliminate_kronecker(KroneckerDelta([ n=>upper, r=>lower])*(eliminate_kronecker(rename_dummies((t11)))))))

#t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(EpsilonTensor([t=>upper,r=>upper])*(EpsilonTensor([B=>upper,v=>upper])*((EpsilonTensor([y=>upper,g=>upper])*(EpsilonTensor([X=>upper,h=>upper])*eliminate_kronecker(t11))))))))
#t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(EpsilonTensor([r=>lower,s=>lower])*eliminate_kronecker(t11))))

# t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(EpsilonTensor([y=>lower,z=>lower])*(EpsilonTensor([t=>lower,u=>lower])*eliminate_kronecker(t11)))))

t11=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(t11)))))



#  t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(EpsilonTensor([s=>upper, F=>upper])*(EpsilonTensor([E=>upper, r=>upper])*eliminate_kronecker(t11)))))


# t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(KroneckerDelta([s=>upper, F=>lower])*(EpsilonTensor([r=>upper, D=>upper])*eliminate_kronecker(t11)))))





#    t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(EpsilonTensor([Z=>lower, d=>lower])*(EpsilonTensor([C=>lower, l=>lower])*eliminate_kronecker(t11)))))

# t11=eliminate_epsilon(rename_dummies(eliminate_epsilon((EpsilonTensor([Z=>lower, d=>lower])*eliminate_kronecker(t11)))))



#   t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(eliminate_kronecker(KroneckerDelta([s=>upper, E=>upper])*res))))
# temp=eliminate_epsilon(rename_dummies(eliminate_epsilon(EpsilonTensor([r=>lower, E=>lower])*(EpsilonTensor([s=>lower, F=>lower])*eliminate_kronecker(bt242)))))

println("\$\$ \\theta=",tt,"\\quad \\overline{\\theta}=",tb,":\$\$")
println("\\begin{dmath}",to_latex(t11),"=0\\end{dmath}")
println(t11)

