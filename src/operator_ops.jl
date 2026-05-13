# include("operator.jl")

##############################################################################################################################
##############################################################################################################################

# graded leibniz

# function product_rule(O::LeibnizOperatorSuperType)
#     new_summand=Nothing
#     println("operator is")
#     println(O)
#     for summandind= eachindex(O.operand.summands) #for each sum term in the expression
#         summand=O.operand.summands[summandind] #get the summand (tensor term)
#         grad=0 #initialise grading
#         for tensorind in eachindex(summand[2].terms)# for each tensor in the summand (tensorterm)
#             grd=(-1)^(grad*grading(O))
#             println(grd)
#             operator_on_tensor=typeof(O)(O.name,summand[2].terms[tensorind],O.operator_indices) #operator on nth term
#             # println("this is the op on a tens")
#             # println(operator_on_tensor)
#             # println("help")
#             new_prod=1
#             for indsss in eachindex(summand[2].terms)
#                 if indsss!= tensorind
#                new_prod=new_prod*summand[2].terms[indsss]
#                 else
#                     new_prod=new_prod*operator_on_tensor
#                 end
#             end
#             # new_prod=[t  where t=summand[2].terms]
#             # println(new_prod)
#             # new_prod.terms[tensorind]=operator_on_tensor
#             if new_summand==Nothing
#                 new_summand=grd*summand[1]*new_prod #adding the new term to the leibniz expanded version of the summand
#             else
#                 new_summand=new_summand+grd*summand[1]*new_prod #adding the new term to the leibniz expanded version of the summand
#             end
            
#             grad=(grad+grading(summand[2].terms[tensorind]))%2
#         end
        
#     end
#     return new_summand
# end

function op_grading(T::OperatorSuperType)
    grade = 0
    for (ind, pos) in T.operator_indices
        grade = grade + ind.grading
    end
    return grade % 2
end


function update_indices(O::OperatorSuperType)
    update_indices(O.operand)
    O.operand_indices=O.operand.summands[1][2].indices
    O.indices=vcat(O.operator_indices,O.operand_indices)
    return 0
end



##############################################################################################################################
##############################################################################################################################

# printing and display
function to_latex(T::OperatorSuperType)
    out=""*uni_to_TeX(T.name)
    # println(io, T.operator_indices)
    # println(io, T.operand_indices)
    # println(io, T.indices)
    if !isempty(T.operator_indices)
        lowers=""
        uppers=""
        lowersempt::Bool=1
        uppersempt::Bool=1
        # println("oper inds")
        # println(T.operator_indices)
        for (ind, pos) in T.operator_indices
            if pos.position==0
                uppers=uppers*uni_to_TeX(ind.name)
                uppersempt=0
                lowers=lowers * ""
            else
                uppers=uppers * ""
                lowersempt=0
                lowers=lowers * uni_to_TeX(ind.name) 
            end
        end
        if !uppersempt
            out*="^{"*uppers*"}"
        end
        if !lowersempt
            out*= "_{"*lowers*"}"
        end
    end
    out*="\\Big("
    out*=to_latex(T.operand)
    out*="\\Big)"
    return Base.replace(out * " ", "\$" => "")
end

function Base.show(io::IO, T::OperatorSuperType)
    print(io, T.name)
    # println(io, T.operator_indices)
    # println(io, T.operand_indices)
    # println(io, T.indices)
    if !isempty(T.operator_indices)
        lowers=""
        uppers=""
        lowersempt::Bool=1
        uppersempt::Bool=1
        # println("oper inds")
        # println(T.operator_indices)
        for (ind, pos) in T.operator_indices
            if pos.position==0
                uppers=uppers*ind.name
                uppersempt=0
                lowers=lowers * "*"
            else
                uppers=uppers * "*"
                lowersempt=0
                lowers=lowers * ind.name 
            end
        end
        if !uppersempt
            print(io, "^{"*uppers*"}")
        end
        if !lowersempt
            print(io, "_{"*lowers*"}")
        end
    end
    print(io,"[[")
    print(io, T.operand)
    print(io,"]]")
end

function operator_weight(O::OperatorSuperType,w::String)
    if singleterm(O.operand)
        return term_weight(O.operand.summands[1][2],w)
    else
        throw(MixedWeightOperandTerms())
    end
end

function ==(a::OperatorSuperType, b::OperatorSuperType)
    # println(a,b,"a.operand==b.operand",a.operand==b.operand)
    if a.name == b.name && a.indices == b.indices && a.operand==b.operand
        # println("got here ever???")
        # println(a)
        # println(b)
        return true
    else

        return false
    end
end

function copy(T::OperatorSuperType)
    return typeof(T)(T.name, copy(T.operand), copy(T.operator_indices), copy(T.function_of), deepcopy(T.weights), deepcopy(T.tableaux))
end
