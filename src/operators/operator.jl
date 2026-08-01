
mutable struct OperatorData <: TensorDataSuperType #?
    name::String
    operand::TensorExpression
    operator_data::TensorData
    #these are the combined properties of op and oprnd
    indices::Vector{Pair{IndexSuperType,IndexPosition}}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}
    tableaux::Set{YoungTableau}
    labels::Vector{Pair{String,IndexPosition}}

    function OperatorData(name::String, operand::TensorExpression, args...)
        given_data=false
        for arg in args
            if arg isa TensorData
                op_data=arg
                top=Tensor(op_data)
                valid_index_placement(top, operand)
                given_data=true
                break
            end
        end
        if !given_data
            top = Tensor(name, args...)#change eventually to not have to do full construction
            op_data = data(top)
            valid_index_placement(top, operand)
        end
        # println(top)
        # @show operand
        # @show indices(op_data)
        # @show indices(operand)
        return new(name, operand, op_data,
            vcat(indices(op_data), indices(operand)),
            vcat(function_of(op_data), function_of(operand)),
            weights(operand),#merge(weights(operand),weights(op_data)), ### need to add weights... need to think about best way, since not homo across summands in texp. see operator_weight func
            tableaux(op_data),###### actually.... not right... need to shift....  vcat(tableaux(op_data),tableaux(operand)). complicated, would need checks to see if all terms have tabs
            labels(op_data)#vcat(labels(op_data),labels(operand))
        )
    end
    #special case for copy
    function OperatorData(::Symbol, args...)
        return new(args...)
    end
    #special case for copy
    function OperatorData(op_data::TensorData, operand::TensorExpression)
        return new(name(op_data), operand, op_data,
            vcat(indices(op_data), indices(operand)),
            vcat(function_of(op_data), function_of(operand)),
            weights(operand),#merge(weights(operand),weights(op_data)), ### need to add weights... need to think about best way, since not homo across summands in texp. see operator_weight func
            tableaux(op_data),###### actually.... not right... need to shift....  vcat(tableaux(op_data),tableaux(operand)). complicated, would need checks to see if all terms have tabs
            labels(op_data)#vcat(labels(op_data),labels(operand))
        )
    end
end

OperatorData(name::String, operand::TensorTerm, args...) = OperatorData(name, TensorExpression(operand), args...)
OperatorData(name::String, operand::TensorSuperType, args...) = OperatorData(name, TensorExpression(operand), args...)


mutable struct Operator <: OperatorSuperType
    data::OperatorData
    function Operator(name::String, operand, args...)
        dt = OperatorData(name, operand, args...)
        new(dt)
    end
    function Operator(operd::OperatorData)
        new(operd)
    end
end

Base.hash(d::OperatorSuperType, h::UInt) = begin
    return hash(
        (
            name(d),
            op_indices(d),
        ),
        h
    )
end