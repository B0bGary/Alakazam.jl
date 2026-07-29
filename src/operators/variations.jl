
"""Computes the Euler-Lagrange equations for a given lagrangian and tensor field."""
function EL(lag::TensorExpression, A::Tensor, DA::DerivativeSuperType)
    #compute the variational derivative of the lagrangian with respect to the tensor field A
    v=VariationalDerivative(lag, A)
    set_operand!(DA, A)
    #apply the derivative operator to the variational derivative
    dv=apply_derivative(v)
    v∂=VariationalDerivative(lag, DA)
    dv∂=apply_derivative(v∂)
    dn=apply_derivative(simplify(PD(name(DA), simplify(dv∂), coord(DA), op_indices(DA))))
    #simplify the result
    return simplify(dv-dn)
end

copy(V::VariationalSuperType) = typeof(V)(symcopy(V.data), symcopy(V.var))

"""
Variational Derivative. Shouldn't be created by user under normal operation.
"""
mutable struct VariationalDerivative <: VariationalSuperType
    data::OperatorData
    var::TensorSuperType
    function VariationalDerivative(operand::Tensorial, var::TensorSuperType, args...)
        new(OperatorData("δ", TensorExpression(operand), args...), var)
    end
    function VariationalDerivative(data::OperatorData, var)
        new(data, var)
    end
end

"""Custom grading for variational derivatives to account for the w/r variable"""
function grading(D::VariationalSuperType)
    grade = (grading(operand(D))+grading(D.var))%2
    return grade
end
"""Custom grading for variational derivatives to account for the w/r variable"""
function op_grading(D::VariationalSuperType)
    grade =grading(D.var)
    return grade % 2
end

coord(V::VariationalSuperType)=get_operand_function_of(operand(V))

function apply_derivative(D::VariationalSuperType)
    D = product_rule(D)
    out =  zero()
    for summandind in eachindex(summands(D)) # for each summand (coefficient, tensorterm)
        summand = summands(D)[summandind]
        term = summand[1]
        derived=false
        for tensor in summand[2].terms # for each individual tensor in tensorterm
            if typeof(tensor)<:VariationalSuperType
                operands = operand(tensor) #get the operand (should only be one thanks to prod rule)
                tens = summands(operands)[1][2].terms[1] #get the only item in the operand
                #compare the operand to see if it is the tensor that the variation is with respect to
                if tens isa DerivativeSuperType && tensor.var isa DerivativeSuperType #take variation w/r to derivative piece
                    if name(tens) == name(tensor.var)
                        a=isTensor(operand(tens))
                        b=isTensor(operand(tensor.var))
                        if a isa TensorSuperType && b isa TensorSuperType
                            if name(a) == name(b)
                                for inds in zip(indices(tens), indices(tensor.var))
                                    if inds[1][2] != inds[2][2] #check if the positions are the same
                                        #positions not same, use metric...
                                        term*=Metric(inds[1][1] => inds[1][2], inds[2][1] => inds[1][2])
                                    else
                                        term *= KroneckerDelta(inds[1][1]=>inds[1][2], inds[2][1] => op_pos(inds[2][2]))
                                    end
                                    derived=true
                                end
                            end


                        end

                    end
                elseif name(tens) == name(tensor.var) #same name
                    #replace with delta functions for the indices
                    for inds in zip(indices(tens), indices(tensor.var))
                        if inds[1][2] != inds[2][2] #check if the positions are the same
                            term*=Metric(inds[1][1] => inds[1][2], inds[2][1] => inds[1][2])
                        else
                            term *= KroneckerDelta(inds[1][1]=>inds[1][2], inds[2][1] => op_pos(inds[2][2]))
                        end
                        derived=true
                    end
                else
                    term *=  zero() #variation of the tensor with respect to something else is 0
                end
            else
                term *= tensor
            end
        end
        if isZeroTens(term)
            term =  zero()
        end
        if isZeroTens(out)
            out =  zero()
        end
        if derived == true
            out += term
        end

    end
    out = TensorExpression(out)
    return merge_sum(out)
end
