
function specify_indices(TEx::TensorExpression, reps::Dict{Index,Index})
    new_exp=zero
    w=copy(TEx)
    for summand in w.summands
        term=1
        for tens in summand[2].terms
            if typeof(tens)<:OperatorSuperType
                tens.operand=specify_indices(tens.operand,reps)
                tens.operand_indices=tens.operand.summands[1][2].indices
                # println(reps)
                # test=tens.operator_indices[1][1]
                # println(test )
                # println(haskey(reps,test))
                tens.operator_indices=[ haskey(reps,ind[1]) ? reps[ind[1]]=>ind[2] : ind[1]=>ind[2]  for ind in tens.operator_indices]
                # println("op should inclde", tens.operator_indices)
                tens.indices=vcat(tens.operator_indices,tens.operand_indices)
                # println("op has", tens.indices)
            else
                if !isempty(tens.indices)
                    tens.indices=[ haskey(reps,ind[1]) ? reps[ind[1]]=>ind[2] : ind[1]=>ind[2]  for ind in tens.indices]
                end
            end
            if !(typeof(term)<:SymorNum) && ((typeof(term)<:TensorTerm && tens in term.terms) || (typeof(term)<:TensorExpression && tens in term.summands[1][2].terms) )
                term=zero
                break
            else
                term=term*tens
            end
        end
        new_exp+=summand[1]*term
    end

    return new_exp
end

function to_lc_expression(TEx::TensorExpression, reps::Dict{Index,Index})
    #replace non-dummies
    exp=specify_indices(TEx, reps)
    # println(exp)
    #cookery for dummies....
    new_exp=zero

    for summand in exp.summands
        dums=get_dummies(summand[2])
        tempterm=TensorExpression(summand[2])
        # println("reset",tempterm)
        while !isempty(dums)
            dum=pop!(dums)
            if dum.index_set.name=="spacetimelc"
            tempterm=specify_indices(tempterm, Dict(dum => p))#+specify_indices(tempterm, Dict(dum => m))
            end
            # println(tempterm)
        end

        new_exp+=summand[1]*tempterm
    end
    # res_ex=zero
    # for summand in new_exp.summands
    #     flag=0
    #         for tens in summand[2].terms
    #             if typeof(tens)<:OperatorSuperType
    #                 if tens.operator_indices[1][1].name!=tens.operator_indices[2][1].name
    #                     flag=1
    #                  #   println("ignoring term", summand[2])
    #                     break
    #                 end
    #             end
    #         end
    #         if flag==0 &&term_weight(summand[2], "θ")<=2 && term_weight(summand[2], "Θ")<=2
    #             res_ex+=summand[1]*summand[2]
    #         end
    # end

    # end
    out_exp=zero
    flag=0
    for summand in new_exp.summands
        for tens in summand[2].terms
            if typeof(tens)<:OperatorSuperType
                if tens.operator_indices[1][1].name!=tens.operator_indices[2][1].name
                    
                    flag=1
                    break
                end
            end
        end
        if flag==0
            # println("got here, ",summand[2] )
        out_exp+=summand[1]*summand[2]
        end
        flag=0
    end



    return out_exp
end


function to_chiral_part(TEx::TensorExpression, reps::Index)
    new_exp=zero

    for summand in TEx.summands
        term=1

        for tens in summand[2].terms
            if tens.name=="θ" || tens.name=="Θ"
                if tens.indices[1][1].name!=reps.name
                    term=zero
                    break
                end
            end
            term*=tens
        end

        new_exp+=summand[1]*term
    end

    return new_exp
end



#define index sets

spacetime_indices = IndexSet("spacetimelc", 2,lower)

lightcone_indices = IndexSet("lightcone", 2,lower)
± = Index("±", lightcone_indices, 1)
∓ = Index("∓", lightcone_indices, 1)
p = Index("+", lightcone_indices, 1)
m = Index("-", lightcone_indices, 1)