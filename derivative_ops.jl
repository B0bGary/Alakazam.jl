include("operator.jl")

##############################################################################################################################
##############################################################################################################################

#apply the derivative and cancel terms that are annihilated

function apply_derivative(D::DerivativeSuperType)
    D=product_rule(D)
    new_summands=D.summands
    for summandind in eachindex(D.summands) # for each summand (coefficient, tensorterm)
        summand=new_summands[summandind]
        for tensor in summand[2].terms # for each individual tensor in tensorterm
            if typeof(tensor)<:DerivativeSuperType # if the tensor is an operator
                operands=tensor.operand #get the operand (should only be one thanks to prod rule)
                tens=operands.summands[1][2].terms[1] #get the only item in the operand
                if isempty(intersect(Set(tens.function_of),Set(tensor.coord))) #see if tensor is a function of the derivative coords
                    # println(summand[1])
                    new_summands[summandind]= 0 => summand[2] #annihilate the term
                end
            end
        end
    end
    
        D.summands=new_summands
    return merge_sum(D)
end

function expand_derivative(TEX::TensorExpression)


end

# graded leibniz

function product_rule(O::DerivativeSuperType)
    new_summand=Nothing
    # println("operator is")
    # println(O)
    for summandind= eachindex(O.operand.summands) #for each sum term in the expression
        summand=O.operand.summands[summandind] #get the summand (tensor term)
        grad=0 #initialise grading
        for tensorind in eachindex(summand[2].terms)# for each tensor in the summand (tensorterm)
            grd=(-1)^(grad*grading(O))
            operator_on_tensor=typeof(O)(O.name,summand[2].terms[tensorind],O.coord, O.operator_indices) #operator on nth term
            # println("this is the op on a tens")
            # println(operator_on_tensor)
            # println("help")
            new_prod=1
            for indsss in eachindex(summand[2].terms)
                if indsss!= tensorind
               new_prod=new_prod*summand[2].terms[indsss]
                else
                    new_prod=new_prod*operator_on_tensor
                end
            end
            # new_prod=[t  where t=summand[2].terms]
            # println(new_prod)
            # new_prod.terms[tensorind]=operator_on_tensor
            if new_summand==Nothing
                new_summand=grd*summand[1]*new_prod #adding the new term to the leibniz expanded version of the summand
            else
                new_summand=new_summand+grd*summand[1]*new_prod #adding the new term to the leibniz expanded version of the summand
            end
            
            grad=(grad+grading(summand[2].terms[tensorind]))%2
        end
        
    end
    new_summand
end
