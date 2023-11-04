
# θ = Coordinate("θ", [flavour_indices, flavour2_indices])
# θbar = Coordinate("Θ", [flavour_conj_indices, spinor_conj_indices])

θMAX = 4
θbarMAX = 4

function project_component(exp::TensorExpression, thet::Coordinate, derivinds::Vector{Pair{Index,IndexPosition}}, theta_weight::Integer, theta_bar_weight::Integer, thet2::Coordinate=nothing,deriv2inds::Vector{Pair{Index,IndexPosition}}=nothing, theta1lev::Integer=1, theta2lev::Integer=1)
    # println(exp)
    # println("here")

    exp = get_weighted_terms(exp, "θ", theta_weight)
    # println(exp)
    exp = get_weighted_terms(exp, "Θ", theta_bar_weight)

    
    final = exp
    # println(exp)
    if theta_weight == 1
        res = partially_project_component(exp, thet, derivinds)
        out = res[1]
        # out = get_weighted_terms(out, "θ", 0)
        # out = get_weighted_terms(out, "Θ", 0)
        final = out
    elseif theta_weight == 2
        # println("aqui")
        res = partially_project_component(exp, thet, derivinds)
        out = res[1]
        usind = res[2]
        res2 = partially_project_component(out, thet, derivinds)
        out = res2[1]
        usind2 = res2[2]
        # println("eps ma")  
        # println(out)

        exp1 = -1 / 2 * simplify(EpsilonTensor([usind[1][1] => op_pos(usind[1][2]), usind2[1][1] => op_pos(usind2[1][2])]) * out)
        exp2 = -1 / 2 * simplify(EpsilonTensor([usind[2][1] => op_pos(usind[2][2]), usind2[2][1] => op_pos(usind2[2][2])]) * out)
        # println("ici")
        # o1= get_weighted_terms(exp1, "θ", 0)
        # o2= get_weighted_terms(exp2, "θ", 0)
        final = (exp1=>exp2)[theta1lev]
    elseif theta_weight == 3
        # println("aqui")
        res = partially_project_component(exp, thet, derivinds)
        out = res[1]
        usind = res[2]
        res2 = partially_project_component(out, thet, derivinds)
        out = res2[1]
        usind2 = res2[2]
        res3 = partially_project_component(out, thet, derivinds)
        out = res3[1]
        usind3 = res3[2]
        # println("eps ma")  
        # println(out)

        exp1 = simplify(EpsilonTensor([usind2[2][1] => op_pos(usind2[2][2]), usind3[2][1] => op_pos(usind3[2][2])]) * (EpsilonTensor([usind2[1][1] => op_pos(usind2[1][2]), usind[1][1] => op_pos(usind[1][2])]) * out))
        #  exp1=(out)
        # println("ici")
        # # o1= get_weighted_terms(exp1, "θ", 0)
        # # o2= get_weighted_terms(exp2, "θ", 0)
        final = -1 * exp1 / 9
    elseif theta_weight == 4

        res = partially_project_component(exp, thet, derivinds)
        println("1st deriv done")
        out = res[1]
        usind = res[2]
        res2 = partially_project_component(simplify(out), thet, derivinds)
        println("2nd deriv done")
        out = res2[1]
        usind2 = res2[2]
        res3 = partially_project_component(simplify(out), thet, derivinds)
        println("3rd deriv done")
        out = res3[1]
        usind3 = res3[2]
        res4 = partially_project_component(simplify(out), thet, derivinds)
        println("4th deriv done")
        out = res4[1]
        # out=simplify(res4[1])
        usind4 = res4[2]
        # println("eps ma")  
        # println(out)
        println("starting eps contraction")
        exp1 = simplify(EpsilonTensor([usind4[2][1] => op_pos(usind4[2][2]), usind3[2][1] => op_pos(usind3[2][2])]) *
                        (EpsilonTensor([usind2[2][1] => op_pos(usind2[2][2]), usind[2][1] => op_pos(usind[2][2])])
                         *
                         (EpsilonTensor([usind4[1][1] => op_pos(usind4[1][2]), usind2[1][1] => op_pos(usind2[1][2])])
                          *
                          (EpsilonTensor([usind3[1][1] => op_pos(usind3[1][2]), usind[1][1] => op_pos(usind[1][2])]) * out))))
        final = exp1 / 144
    end

    if theta_bar_weight != 0
        if theta_bar_weight == 1
            res = partially_project_component(final, thet2, deriv2inds)
            out = res[1]
            # out = get_weighted_terms(out, "θ", 0)
            # out = get_weighted_terms(out, "Θ", 0)
            final = out
        elseif theta_bar_weight == 2
            # println("aqui")
            res = partially_project_component(final, thet2, deriv2inds)
            out = res[1]
            usind = res[2]
            res2 = partially_project_component(out, thet2, deriv2inds)
            out = res2[1]
            usind2 = res2[2]
            # println("eps ma")  
            # println(out)

            exp1 = -1 / 2 * simplify(EpsilonTensor([usind[1][1] => op_pos(usind[1][2]), usind2[1][1] => op_pos(usind2[1][2])]) * out)
            exp2 = -1 / 2 * simplify(EpsilonTensor([usind[2][1] => op_pos(usind[2][2]), usind2[2][1] => op_pos(usind2[2][2])]) * out)
            # println("ici")
            # o1= get_weighted_terms(exp1, "θ", 0)
            # o2= get_weighted_terms(exp2, "θ", 0)
            final = (exp1=>exp2)[theta2lev]
        elseif theta_bar_weight == 3
            # println("aqui")
            res = partially_project_component(final, thet2, deriv2inds)
            out = res[1]
            usind = res[2]
            res2 = partially_project_component(out, thet2, deriv2inds)
            out = res2[1]
            usind2 = res2[2]
            res3 = partially_project_component(out, thet2, deriv2inds)
            out = res3[1]
            usind3 = res3[2]
            # println("eps ma")  
            # println(out)

            exp1 = simplify(EpsilonTensor([usind2[2][1] => op_pos(usind2[2][2]), usind3[2][1] => op_pos(usind3[2][2])]) * (EpsilonTensor([usind2[1][1] => op_pos(usind2[1][2]), usind[1][1] => op_pos(usind[1][2])]) * out))
            #  exp1=(out)
            # println("ici")
            # # o1= get_weighted_terms(exp1, "θ", 0)
            # # o2= get_weighted_terms(exp2, "θ", 0)
            final = -1 * exp1 / 9
        elseif theta_bar_weight == 4

            res = partially_project_component(final, thet2, deriv2inds)
            println("1st deriv done")
            out = res[1]
            usind = res[2]
            res2 = partially_project_component(simplify(out), thet2, deriv2inds)
            println("2nd deriv done")
            out = res2[1]
            usind2 = res2[2]
            res3 = partially_project_component(simplify(out), thet2, deriv2inds)
            println("3rd deriv done")
            out = res3[1]
            usind3 = res3[2]
            res4 = partially_project_component(simplify(out), thet2, deriv2inds)
            println("4th deriv done")
            out = res4[1]
            # out=simplify(res4[1])
            usind4 = res4[2]
            # println("eps ma")  
            # println(out)
            println("starting eps contraction")
            exp1 = simplify(EpsilonTensor([usind4[2][1] => op_pos(usind4[2][2]), usind3[2][1] => op_pos(usind3[2][2])]) *
                            (EpsilonTensor([usind2[2][1] => op_pos(usind2[2][2]), usind[2][1] => op_pos(usind[2][2])])
                             *
                             (EpsilonTensor([usind4[1][1] => op_pos(usind4[1][2]), usind2[1][1] => op_pos(usind2[1][2])])
                              *
                              (EpsilonTensor([usind3[1][1] => op_pos(usind3[1][2]), usind[1][1] => op_pos(usind[1][2])]) * out))))
            final = exp1 / 144
        end
    end

    return final

end



function project_component(exp::TensorExpression, thet::Coordinate, derivinds::Vector{Pair{Index,IndexPosition}}, theta_weight::Integer, theta_bar_weight::Integer, theta1lev::Integer=1)
    # println(exp)
    # println("here")

    exp = get_weighted_terms(exp, "θ", theta_weight)
    # println(exp)
    exp = get_weighted_terms(exp, "Θ", theta_bar_weight)

    
    final = exp
    # println(exp)
    if theta_weight == 1
        res = partially_project_component(exp, thet, derivinds)
        out = res[1]
        # out = get_weighted_terms(out, "θ", 0)
        # out = get_weighted_terms(out, "Θ", 0)
        final = out
    elseif theta_weight == 2
        # println("aqui")
        res = partially_project_component(exp, thet, derivinds)
        out = res[1]
        usind = res[2]
        res2 = partially_project_component(out, thet, derivinds)
        out = res2[1]
        usind2 = res2[2]
        # println("eps ma")  
        # println(out)

        exp1 = -1 / 2 * simplify(EpsilonTensor([usind[1][1] => op_pos(usind[1][2]), usind2[1][1] => op_pos(usind2[1][2])]) * out)
        exp2 = -1 / 2 * simplify(EpsilonTensor([usind[2][1] => op_pos(usind[2][2]), usind2[2][1] => op_pos(usind2[2][2])]) * out)
        # println("ici")
        # o1= get_weighted_terms(exp1, "θ", 0)
        # o2= get_weighted_terms(exp2, "θ", 0)
        final = (exp1=>exp2)[theta1lev]
    elseif theta_weight == 3
        # println("aqui")
        res = partially_project_component(exp, thet, derivinds)
        out = res[1]
        usind = res[2]
        res2 = partially_project_component(out, thet, derivinds)
        out = res2[1]
        usind2 = res2[2]
        res3 = partially_project_component(out, thet, derivinds)
        out = res3[1]
        usind3 = res3[2]
        # println("eps ma")  
        # println(out)

        exp1 = simplify(EpsilonTensor([usind2[2][1] => op_pos(usind2[2][2]), usind3[2][1] => op_pos(usind3[2][2])]) * (EpsilonTensor([usind2[1][1] => op_pos(usind2[1][2]), usind[1][1] => op_pos(usind[1][2])]) * out))
        #  exp1=(out)
        # println("ici")
        # # o1= get_weighted_terms(exp1, "θ", 0)
        # # o2= get_weighted_terms(exp2, "θ", 0)
        final = -1 * exp1 / 9
    elseif theta_weight == 4

        res = partially_project_component(exp, thet, derivinds)
        println("1st deriv done")
        out = res[1]
        usind = res[2]
        res2 = partially_project_component(simplify(out), thet, derivinds)
        println("2nd deriv done")
        out = res2[1]
        usind2 = res2[2]
        res3 = partially_project_component(simplify(out), thet, derivinds)
        println("3rd deriv done")
        out = res3[1]
        usind3 = res3[2]
        res4 = partially_project_component(simplify(out), thet, derivinds)
        println("4th deriv done")
        out = res4[1]
        # out=simplify(res4[1])
        usind4 = res4[2]
        # println("eps ma")  
        # println(out)
        println("starting eps contraction")
        exp1 = simplify(EpsilonTensor([usind4[2][1] => op_pos(usind4[2][2]), usind3[2][1] => op_pos(usind3[2][2])]) *
                        (EpsilonTensor([usind2[2][1] => op_pos(usind2[2][2]), usind[2][1] => op_pos(usind[2][2])])
                         *
                         (EpsilonTensor([usind4[1][1] => op_pos(usind4[1][2]), usind2[1][1] => op_pos(usind2[1][2])])
                          *
                          (EpsilonTensor([usind3[1][1] => op_pos(usind3[1][2]), usind[1][1] => op_pos(usind[1][2])]) * out))))
        final = exp1 / 144
    end
    return final

end



function partially_project_component(exp::TensorExpression, thet::Coordinate, derivinds::Vector{Pair{Index,IndexPosition}})
    println("doing this")
    d_indices = Vector{Pair{Index,IndexPosition}}()
    for pr in derivinds
        ty = pr[1].index_set
        ind = pop!(get_unused_indices(exp, ty))
        push!(d_indices, ind => pr[2])
    end

    if thet.name == "θ"
        n_deriv = Derivative("𝛛", exp, [thet], d_indices)
    elseif thet.name == "Θ"
        n_deriv = Derivative("∂", exp, [thet], d_indices)
    end
    println("dont drink and derive")
    out = apply_derivative(n_deriv)
    # out = simplify(out)
    println("got here")
    return (out, d_indices)
end