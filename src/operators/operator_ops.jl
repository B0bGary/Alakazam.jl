##############################################################################################################################
##############################################################################################################################
## get properties
##############################################################################################################################
##############################################################################################################################
operand(O::OperatorSuperType) = data(O).operand
op_data(O::OperatorSuperType) = data(O).operator_data
op_indices(O::OperatorSuperType) = data(O).operator_indices
op_function_of(O::OperatorSuperType) = function_of(op_data(O))
op_weights(O::OperatorSuperType) = weights(op_data(O))
op_tableaux(O::OperatorSuperType) = tableaux(op_data(O))
op_labels(O::OperatorSuperType) = labels(op_data(O))
# @inline uid(t::OperatorSuperType)=op_data(t).uid
operand_indices(O::OperatorSuperType) = indices(operand(O))

# function set_operand!(O::OperatorSuperType, oprd::TensorExpression)
#     data(O).operand = oprd
#     update_indices!(O)#make sure everything is updated when we do these sort of things
#     return O
# end


coord(O::OperatorSuperType) = get_operand_function_of(operand(O))


function get_operand_function_of(expr::TensorExpression)
    coords = Coordinate[]
    for (_, term) in expr.summands
        for t in term.terms
            append!(coords, function_of(t))
        end
    end
    return coords
end




function set_operand!(O::OperatorSuperType, oprd::Tensorial)
    data(O).operand = TensorExpression(oprd)
    update_indices!(O)#make sure everything is updated when we do these sort of things
    update_functionof!(O)
    return O
end

function set_operand_deep!(O::TensorExpression, oprd::Tensorial)
    # @show O, oprd
    for summand in summands(O)
        for tterm in terms(summand)
            n = 0 #keep track of new length of tterm w/ weird operator exchanges...
            for term in terms(tterm)
                if term isa OperatorSuperType
                    if operand(term) == zero()
                        set_operand!(term, oprd)
                    else
                        set_operand_deep!(term, oprd)
                    end
                end
                n += length(indices(term))
            end
            update_indices!(tterm, n)
        end

    end
    update_functionof!(O)
    return O
end

function update_functionof!(O::OperatorSuperType)
    a = op_function_of(O)
    b = function_of(operand(O))
    isempty(b) && (data(O).function_of=a; return O)
    isempty(a) && (data(O).function_of=b; return O)

    contained = true
    for x in b
        if !(x in a)
            ;
            contained = false;
            break;
        end
    end
    contained && (data(O).function_of=a; return O)

    out = copy(a)
    for x in b
        _insert_coord!(out, x)
    end
    data(O).function_of = out
    return O
end
function update_functionof!(Tex::TensorExpression)
    for summand in summands(Tex)
        for term in terms(summand)
            if term isa OperatorSuperType
                update_functionof!(term)
            end
        end
    end
end

"""Only updates the external indices"""
function set_indices!(O::OperatorSuperType, inds::Vector{<:Pair{<:IndexSuperType,IndexPosition}}) #where {M<:IndexSuperType}
    data(O).indices = inds
    return O
end

function set_operand_indices!(O::OperatorSuperType, inds::Vector{<:Pair{<:IndexSuperType,IndexPosition}}) #where {M<:IndexSuperType}
    set_indices!(O, vcat(op_indices(O), inds))
    return O
end
function set_op_indices!(O::OperatorSuperType, inds::Vector{<:Pair{<:IndexSuperType,IndexPosition}})# where {M<:IndexSuperType}
    data(O).operator_indices = inds
    return O
end

function set_indices_deep!(O::OperatorSuperType, rep_guide::Dict{<:Pair{<:IndexSuperType,IndexPosition},<:Pair{<:IndexSuperType,IndexPosition}})# where {M<:IndexSuperType}
    op_repl_inds=[k in keys(rep_guide) ? rep_guide[k] : k for k in op_indices(O)]
    set_op_indices!(O, op_repl_inds)
    # @show op_repl_inds
    set_indices_deep!(operand(O), rep_guide)
    update_indices!(O)
    return O
end


function update_indices!(O::OperatorSuperType)#good
    update_indices!(operand(O))
    return set_indices!(O, vcat(op_indices(O), indices(operand(O))))
end

"""
Grading of the operator, not including the operand. Used for product rule.
"""
function op_grading(T::OperatorSuperType)
    grade = 0
    for (ind, pos) in op_indices(T)
        grade = grade + grading(ind)
    end
    return grade % 2
end

"""
Grading of an operator: its own grading composed with that of what it acts on.

Without this an operator falls through to the generic `grading(::TensorSuperType)`, which sums over
`indices`. This would accidentally ignore operands with non-standard grading.
"""
function grading(T::OperatorSuperType)
    opd = operand(T)
    g = (opd === nothing || isZeroTens(opd)) ? 0 : grading(opd)
    return (op_grading(T) + g) % 2
end


###############################################################
###############################################################
#multiple dispatch for operators.
function (O::OperatorSuperType)(T::Tensorial)
    texp = TensorExpression(T)
    nO = copy(O)
    set_operand!(nO, texp)
    return nO
end

"""If a single term, pull out the global factor."""
function simplify(T::OperatorSuperType)#incase of operators..
    (coefin, O)=factor_coef(T)
    simplify(coefin*TensorExpression(O))
end


function factor_coef(O::OperatorSuperType)
    O=copy(O)
    coefin=1
    if singleterm(operand(O))
        coefin=summands(operand(O))[1][1]
        set_operand!(O, TensorExpression(summands(operand(O))[1][2]))
    end
    return (coefin, O)
end

function *(N::SymorNum, O::OperatorSuperType)
    (coefin, O)=factor_coef(O)
    TensorExpression(N*coefin, O)
end
*(O::OperatorSuperType, N::SymorNum) = N*O




##############################################################################################################################
##############################################################################################################################

# printing and display
function to_latex(T::OperatorSuperType)
    out = "" * uni_to_TeX(name(T))
    if !isempty(op_indices(T))
        lowers = ""
        uppers = ""
        lowersempt::Bool = 1
        uppersempt::Bool = 1
        # println("oper inds")
        # println(T.operator_indices)
        for (ind, pos) in op_indices(T)
            if pos.position == 0
                uppers = uppers * uni_to_TeX(ind.name)
                uppersempt = 0
                lowers = lowers * ""
            else
                uppers = uppers * ""
                lowersempt = 0
                lowers = lowers * uni_to_TeX(ind.name)
            end
        end
        if !uppersempt
            out *= "^{" * uppers * "}"
        end
        if !lowersempt
            out *= "_{" * lowers * "}"
        end
    end
    out *= "\\left("
    out *= to_latex(operand(T))
    out *= "\\right)"
    return Base.replace(out * " ", "\$" => "")
end

function show_op_base(io::IO, T::OperatorSuperType)
    print(io, name(T))
    if !isempty(op_indices(T))
        lowers = ""
        uppers = ""
        lowersempt::Bool = 1
        uppersempt::Bool = 1
        for (ind, pos) in op_indices(T)
            if pos.position == 0
                if typeof(ind) == DummyPatternIndex
                    uppers *= "?"
                end
                uppers = uppers * ind.name
                uppersempt = 0
                lowers = lowers * "∘"
            else
                uppers = uppers * "∘"
                lowersempt = 0
                if typeof(ind) == DummyPatternIndex
                    lowers *= "?"
                end
                lowers = lowers * ind.name
            end
        end
        if !uppersempt
            print(io, "^{" * uppers * "}")
        end
        if !lowersempt
            print(io, "_{" * lowers * "}")
        end
    end
    lab = labels(T)
    if !isempty(lab)
        lowers = ""
        uppers = ""
        lowersempts::Bool = 1
        uppersempts::Bool = 1
        for (ind, pos) in lab
            if pos.position == 0
                uppers = uppers * ind
                uppersempts = 0
                lowers = lowers * "∘"
            else
                uppers = uppers * "∘"
                lowersempts = 0
                lowers = lowers * ind
            end
        end
        if !uppersempts
            out *= "^{" * uppers * "}"
        end
        if !lowersempts
            out *= "_{" * lowers * "}"
        end
    end

end

function Base.show(io::IO, T::OperatorSuperType)
    show_op_base(io, T)
    print(io, "⟦")
    print(io, operand(T))
    print(io, "⟧")
end

function operator_weight(O::OperatorSuperType, w::String)
    if singleterm(operand(O))
        return term_weight(operand(O).summands[1][2], w)
    else
        throw(MixedWeightOperandTerms())
    end
end

function ==(a::OperatorSuperType, b::OperatorSuperType)
    # @show a
    # @show b
    # println("operator eqlt")
    if name(a) == name(b) && indices(a) == indices(b) && operand(a) == operand(b)#may want better notion at some pt?
        return true
    else

        return false
    end
end

# function copy(T::OperatorSuperType)
#     return typeof(T)(name(T), copy(operand(T)), copy(indices(op_data(T))), copy(function_of(op_data(T))), deepcopy(weights(op_data(T))), deepcopy(tableaux(op_data(T))), deepcopy(labels(op_data(T))))
# end
# function copy(T::OperatorSuperType)
#     return typeof(T)((symcopy(getfield(T, f)) for f in fieldnames(typeof(T)))...)
# end
copy(O::OperatorSuperType) = typeof(O)(symcopy(O.data))
function symcopy(d::OperatorData)
    return OperatorData(
        :copy,                          # dispatch to the raw new(args...) constructor
        d.name,
        symcopy(d.operand),             # TensorExpression operand
        d.operator_data,                # immutable and shared now
        d.operator_indices,             # copy on write
        d.indices,                # copy on write now
        d.function_of,            # might change if operand is updated -> but with current implementation, i think its ok
        d.weights,                      # shared 
        d.tableaux,                     # shared
        d.labels,                       # shared
    )
end
symcopy(o::OperatorSuperType) = copy(o)



copy_owned(d::OperatorData) = OperatorData(
    :copy,
    d.name,
    copy_owned(d.operand),          # operand tree owned all the way down
    d.operator_data,                # immutable and shared
    copy(d.operator_indices),       # OWNED
    copy(d.indices),                # OWNED
    copy(d.function_of),
    d.weights,
    d.tableaux,
    d.labels,
)
copy_owned(o::OperatorSuperType) =
    typeof(o)((copy_owned(getfield(o, f)) for f in fieldnames(typeof(o)))...)
