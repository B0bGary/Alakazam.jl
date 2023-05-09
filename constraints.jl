

w1=1
w2=2

# println(to_latex(res))
# for w1 in range(1,4)
#     for w2 in range(1,4)
    t11=get_weighted_terms(get_weighted_terms(res,"θ",w1),"Θ",w2)
    # println("t:",w1,", tbar:",w2)
    merge_sum(t11)

    t11=merge_sum(sort(rename_dummies(eliminate_epsilon(rename_dummies(sort(pop_metric(t11)))))))

    t11=merge_sum(eliminate_kronecker(t11))
    # for summand in t11.summands
    #     println(summand)
    #     println(get_free_indices(summand[2]))
    # end
    out=merge_sum(sort(pop_metric(t11)))
println(out)
    out=strip_thetas(out)

    # for summand in out.summands
    #         println(get_free_indices(summand[2]))
    # end
    # println("constraint:")
    out=rename_dummies(eliminate_epsilon(eliminate_epsilon(EpsilonTensor([k=>lower,r=>lower])*(EpsilonTensor([γdot=>lower,ρdot=>lower])*eliminate_kronecker(out)))))
    println(out)
    # end
# end