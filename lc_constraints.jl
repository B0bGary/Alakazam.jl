
t=3
tb=0
println(res)
t11=rename_dummies(get_weighted_terms(get_weighted_terms(res,"θ",t),"Θ",tb))
    # println("t:",w1,", tbar:",w2)
    merge_sum(t11)

    t11=merge_sum(sort(rename_dummies(eliminate_epsilon(rename_dummies(sort(pop_metric(t11)))))))

    t11=merge_sum(eliminate_kronecker(t11))
    # for summand in t11.summands
    #     println(summand)
    #     println(get_free_indices(summand[2]))
    # end
    t11=rename_dummies(sort(merge_sum(sort(pop_metric(t11)))))
println(t11)
println("\$\$",to_latex(eliminate_epsilon(t11)),"=0\$\$")
t11=strip_thetas(t11)
println(t11)
t11=eliminate_kronecker(merge_sum(t11))
println(t11)
#n=4,0 2,0
# t11=eliminate_epsilon(rename_dummies(eliminate_epsilon(EpsilonTensor([b=>upper,c=>upper])*(EpsilonTensor([j=>upper,l=>upper])*eliminate_kronecker(t11)))))
t11=(merge_sum(eliminate_kronecker(rename_dummies((t11)))))
# t11=eliminate_epsilon(EpsilonTensor([d=>upper,c=>upper])*eliminate_kronecker(rename_dummies((t11))))

println("\$\$ \\theta=",t,"\\quad \\overline{\\theta}=",tb,":\$\$")
println("\$\$",to_latex(t11),"=0\$\$")
println(t11)

