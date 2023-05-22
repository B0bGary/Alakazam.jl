
# θ = Coordinate("θ", [flavour_indices, flavour2_indices])
# θbar = Coordinate("Θ", [flavour_conj_indices, spinor_conj_indices])


function project_component(exp::TensorExpression, thet::Coordinate, derivinds::Vector{Pair{Index,IndexPosition}}, theta_weight::Integer,theta_bar_weight::Integer)
    # println(exp)
    # println("here")
    exp= get_weighted_terms(exp, "θ", theta_weight)
    # println(exp)
    exp= get_weighted_terms(exp, "Θ", theta_bar_weight)
    # println(exp)
    if theta_weight == 1
        res=partially_project_component(exp,thet,derivinds)
        out=res[1]
        # out = get_weighted_terms(out, "θ", 0)
        # out = get_weighted_terms(out, "Θ", 0)
        return out
    elseif theta_weight==2
        println("aqui")
        res=partially_project_component(exp,thet,derivinds)
        out=res[1]
        usind=res[2]
        res2=partially_project_component(out,thet,derivinds)
        out=res2[1]
        usind2=res2[2]
        println("eps ma")  
        println(out)

        exp1=simplify(EpsilonTensor([usind[1][1]=>op_pos(usind[1][2]), usind2[1][1]=>op_pos(usind2[1][2]) ])*out)
        exp2=simplify(EpsilonTensor([usind[2][1]=>op_pos(usind[2][2]), usind2[2][1]=>op_pos(usind2[2][2]) ])*out)
        println("ici")
        # o1= get_weighted_terms(exp1, "θ", 0)
        # o2= get_weighted_terms(exp2, "θ", 0)
        return exp1=>exp2
elseif theta_weight==3
    # println("aqui")
    res=partially_project_component(exp,thet,derivinds)
    out=res[1]
    usind=res[2]
    res2=partially_project_component(out,thet,derivinds)
    out=res2[1]
    usind2=res2[2]
    res3=partially_project_component(out,thet,derivinds)
    out=res3[1]
    usind3=res3[2]
    # println("eps ma")  
    # println(out)

    exp1=simplify(EpsilonTensor([usind2[2][1]=>op_pos(usind2[2][2]), usind3[2][1]=>op_pos(usind3[2][2]) ])*(EpsilonTensor([usind2[1][1]=>op_pos(usind2[1][2]), usind[1][1]=>op_pos(usind[1][2]) ])*out))

    # println("ici")
    # # o1= get_weighted_terms(exp1, "θ", 0)
    # # o2= get_weighted_terms(exp2, "θ", 0)
    return exp1
elseif theta_weight==4
    res=partially_project_component(exp,thet,derivinds)
    out=res[1]
    usind=res[2]
    res2=partially_project_component(simplify(out),thet,derivinds)
    out=res2[1]
    usind2=res2[2]
    res3=partially_project_component(simplify(out),thet,derivinds)
    out=res3[1]
    usind3=res3[2]
    res4=partially_project_component(simplify(out),thet,derivinds)
    out=simplify(res4[1])
    usind4=res4[2]
    # println("eps ma")  
    # println(out)

    exp1=simplify(EpsilonTensor([usind4[2][1]=>op_pos(usind4[2][2]), usind3[2][1]=>op_pos(usind3[2][2]) ])*
    (EpsilonTensor([usind2[2][1]=>op_pos(usind2[2][2]), usind[2][1]=>op_pos(usind[2][2]) ])
    *(EpsilonTensor([usind4[1][1]=>op_pos(usind4[1][2]), usind2[1][1]=>op_pos(usind2[1][2]) ])
    *(EpsilonTensor([usind3[1][1]=>op_pos(usind3[1][2]), usind[1][1]=>op_pos(usind[1][2]) ])*out))))
    return exp1
end
    
end

function partially_project_component(exp::TensorExpression, thet::Coordinate, derivinds::Vector{Pair{Index,IndexPosition}})
    println("doing this")  
        d_indices = Vector{Pair{Index,IndexPosition}}()
        for pr in derivinds
            ty = pr[1].index_set
            ind = pop!(get_unused_indices(exp, ty))
            push!(d_indices, ind => pr[2])
        end

        if thet.name == θ.name
            n_deriv = Derivative("𝛛", exp, [θ], d_indices)
        elseif thet.name == θbar.name
            n_deriv = Derivative("∂", exp, [θbar], d_indices)
        end
        println("dont drink and derive")  
        out = apply_derivative(n_deriv)
        out = simplify(out)
        println("got here")  
        return (out, d_indices)
end