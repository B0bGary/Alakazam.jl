copy(O::LieBracketSuperType) = typeof(O)(symcopy(O.A),O.grading,symcopy(O.data))
"""
Lie Bracket Operator. Obeys derivation rule thanks to subtyping LeibnizOperatorSuperType. 
Antisymmetry is encoded in LieBracketSuperType.
For [A,B], B is the operand, allowing for graded Leibniz naturally.
    """
mutable struct Lie <: LieBracketSuperType
    A::TensorExpression
    grading::Int
    data::OperatorData
    function Lie(A::Tensorial, B::Tensorial, grading::Int=BOSONIC)
        @assert grading in (0, 1) "grading must be 0 (bosonic) or 1 (fermionic)"
        op_data = A isa OperatorSuperType ?
                  OperatorData("Lie"*sprint(print, A), TensorExpression(B)) :
                  OperatorData("Lie"*sprint(print, A), TensorExpression(B), data(A), indices(A))
        #note that the name here is just for distinction in op data
        return new(TensorExpression(A), grading, op_data)
    end
    function Lie(v::Vector, grading::Int)
        @assert length(v) == 2 "Lie vector form requires exactly two elements; nest vectors for deeper brackets"
        left = v[1] isa Tensorial ? v[1] : error("Lie: left element must be Tensorial or nested bracket, got $(typeof(v[1]))")
        right = v[2] isa Tensorial ? v[2] : error("Lie: right element must be Tensorial or nested bracket, got $(typeof(v[2]))")
        return Lie(left, right, grading)
    end
    function Lie(T::TensorExpression, grade::Int, data::OperatorData)
        return new(T, grade, data)
    end
end

""" ── @Lie macro ────────────────────────────────────────────────────────────────
Recursively rewrites nested [] / {} into Lie(A, B, grading) calls,
with grading determined per-node from bracket type:   @Lie [A, {B, {C, D}}]
→ Lie(A, Lie(B, Lie(C, D, FERMIONIC), FERMIONIC), BOSONIC)
"""
macro Lie(expr)
    return esc(_expand_lie(expr))
end

# ── recursive macro helper ────────────────────────────────────────────────────
#
#   Each call handles one bracket node.
#   Returns a Julia Expr that constructs the appropriate Lie or passes through
#   a bare symbol.
#
function _expand_lie(expr)
    # ── bracket node: [] or {} ────────────────────────────────────────────────
    if expr isa Expr && expr.head in (:vect, :braces, :bracescat)
        grading = expr.head === :vect ? :0 : :1

        args = expr.args
        @assert length(args) == 2 """
            @Lie: each bracket must have exactly two elements; nest brackets for deeper structure.
            Got $(length(args)) elements in: $expr
            """

        left = _expand_lie(args[1])   # recurse: may itself be a bracket node
        right = _expand_lie(args[2])   # recurse: may itself be a bracket node

        return :(Lie($left, $right, $grading))

        # ── leaf: bare symbol or pre-built Tensorial expression ──────────────────
    elseif expr isa Symbol || expr isa Expr
        return expr   # passed through as-is; runtime will check <: Tensorial

    else
        error("@Lie: unexpected expression type $(typeof(expr)): $expr")
    end
end

function Base.show(io::IO, T::LieBracketSuperType)
    T.grading==0 ? print(io, "[") : print(io, "{")
    print(io, T.A)
    print(io, ",")
    print(io, operand(T))
    T.grading==0 ? print(io, "]") : print(io, "}")
end


function to_latex(T::LieBracketSuperType)
    out=""
    T.grading==0 ? out*="\\left[" : out*="\\left\\{"
    out*=to_latex(T.A)
    out*=","
    out*=to_latex(operand(T))
    T.grading==0 ? out*="\\right]" : out*="\\right\\}"
    return out
end

"""Gets the first argument of a Lie bracket."""
function Base.first(bracket::LieBracketSuperType)
    return bracket.A
end
function second(bracket::LieBracketSuperType)
    return operand(bracket)
end

function set_first!(bracket::LieBracketSuperType, Anew::Tensorial)
    bracket.A=TensorExpression(Anew)
    set_op_indices!(bracket, indices(Anew))
    update_indices!(bracket)
    return bracket
end
function set_second!(bracket::LieBracketSuperType, Anew::Tensorial)
    set_operand!(bracket, Anew)
    return bracket
end
#### add copy no bang version.


function swap(exp::TensorExpression)
    out= zero()
    for summand in summands(exp)
        out+=summand[1]*swap(summand[2])
    end
    return out
end
function swap(tt::TensorTerm)
    out=1
    for term in terms(tt)
        out*=swap(term)
    end
    return out
end
swap(t::Tensorial) = t


"""Applies (anti)symmetric of the Lie bracket."""
function swap(bracket::LieBracketSuperType)
    A1=first(bracket)
    B1=second(bracket)
    b=copy(bracket)
    set_second!(b, A1)
    set_first!(b, B1)
    return (-1)^(b.grading+1)*b #[,] vs {,}
end

# function expand_commutator(bracket::LieBracketSuperType)
#     return bracket.A*operand(bracket)-(-1)^(bracket.grading)*operand(bracket)*bracket.A
# end
function expand_commutators(TEX::TensorExpression)
    _expand_commutators(product_rule(TEX))
end

function _expand_commutators(TEX::TensorExpression)
    new_summand =  zero()
    for summand in summands(TEX)
        nterm = summand[1]
        for term in summand[2].terms
            if typeof(term) <: LieBracketSuperType
                nterm *= _expand_commutators(term)
            else
                nterm *= term
            end
        end
        new_summand += nterm
    end
    return new_summand
end

function _expand_commutators(bracket::LieBracketSuperType)
    A = bracket.A
    B = operand(bracket)
    # ── check for a known commutation relation ────────────────────────────────
    if isTensor(A) isa MatrixSuperType && isTensor(B) isa MatrixSuperType
        grading, result = get_commutator(isTensor(A), isTensor(B))
        if result !=  zero()
            return result
        end
    end
    # ── no known relation: expand explicitly ──────────────────────────────────
    A_expanded = _expand_commutators(A isa TensorExpression ? A : TensorExpression(A))
    B_expanded = _expand_commutators(B isa TensorExpression ? B : TensorExpression(B))
    return A_expanded * B_expanded - (-1)^(bracket.grading) * B_expanded * A_expanded
end
# passthrough for bare tensors
expand_commutators(T::Union{TensorSuperType,TensorTerm}) = expand_commutators(TensorExpression(T))



function signed_equals(a::LieBracketSuperType, b::LieBracketSuperType)#check to see if ops commute
    if a.grading!=b.grading
        return (false, 0)
    elseif first(a)==first(b) && second(a)==second(b)
        return (true, 1)
    elseif first(a)==second(b) && second(a)==first(b)
        return (true, (-1)^(a.grading+1))
    end
    return (false, 0)
end