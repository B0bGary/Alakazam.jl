using Latexify
using LaTeXStrings
include("tensor_ops.jl")

w1=3
w2=4

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
    # println("this wokr?")
    # println(copy(out))
    out=merge_sum(sort(pop_metric(t11)))
    # println("this ??")
    # println(copy(out))
    out=strip_thetas(out)
    # println("this hbu??")
    # println(copy(out))
    # for summand in out.summands
    #         println(get_free_indices(summand[2]))
    # end
    # println("constraint:")
    out=eliminate_epsilon(eliminate_kronecker(out))
 out=eliminate_epsilon(eliminate_epsilon(EpsilonTensor([s=>lower,p=>lower])*(EpsilonTensor([t=>lower,v=>lower])*(EpsilonTensor([r=>lower,m=>lower])*(EpsilonTensor([l=>lower,q=>lower])*out)))))
  #  out=rename_dummies(eliminate_epsilon(eliminate_epsilon(EpsilonTensor([k=>lower,r=>lower])*(EpsilonTensor([γdot=>lower,ρdot=>lower])*out))))
    println("\$\$\\theta :",w1,",\\quad \\bar{\\theta}:",w2,"\\quad :\\quad")
    println(to_latex(out))
    println("=0\$\$")
    println(out)
    # end
# end