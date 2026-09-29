struct RegisteredExpansion
    canonical_deriv::CovDerivSuperType
    expansion::TensorExpression
end
const _NO_MATCH = Dict{Index,Index}()
const DERIVATIVE_EXPANSIONS = Dict{String,RegisteredExpansion}()

Base.:(==)(a::RegisteredExpansion, b::RegisteredExpansion) =
    all(getfield(a, f) == getfield(b, f) for f in fieldnames(RegisteredExpansion))

clear_covariant!() = empty(DERIVATIVE_EXPANSIONS)

function register_expansion!(D::CovDerivSuperType, expansion::Tensorial)
    DERIVATIVE_EXPANSIONS[D.data.name] = RegisteredExpansion(D, TensorExpression(expansion))
end

macro covariant(expr)
    @assert expr.head == :(=) "Expected `@register_expansion D[μ] = expansion`"
    lhs, rhs = expr.args
    return quote
        let _D = $(esc(lhs)), _rhs = $(esc(rhs))
            register_expansion!(_D, _rhs)
        end
    end
end

function get_expansion_for(D::CovDerivSuperType)::TensorExpression
    haskey(DERIVATIVE_EXPANSIONS, D.data.name) || throw(KeyError("No expansion registered for $(D.data.name)"))
    reg = DERIVATIVE_EXPANSIONS[D.data.name]

    canon_op_indices = op_indices(reg.canonical_deriv)
    D_op_indices = op_indices(D)
    if length(canon_op_indices) != length(D_op_indices)
        println("Registered expansion for $(D.data.name) has $(length(canon_op_indices)) operator indices, but found $(length(D_op_indices))")
        throw(FreeIndexMismatch())
    end

    # @show reg.canonical_deriv
    mapping = construct_index_rep_map(reg.canonical_deriv, D)
    # @show mapping
    # convert Dict{String,Index} → Dict{Pair{Index,IndexPosition}, Pair{Index,IndexPosition}}
    # using the canonical expansion's index positions as keys
    ref = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
    for (canon_name, actual_index) in mapping
        # find the canonical index+position pair from the registered expansion
        for summand in reg.expansion.summands
            for tensor in summand[2].terms
                for pair in indices(tensor)
                    if pair[1].name == canon_name
                        ref[pair] = actual_index => pair[2]
                    end
                end
            end
        end
    end

    expansion = copy(reg.expansion)
    for summand in expansion.summands
        for tensor in summand[2].terms
            set_indices_deep!(tensor, ref)
        end
        update_indices!(summand[2])
    end
    return expansion
end
function _expand_derivative(TEX::Union{TensorExpression,<:CovDerivSuperType}, D::CovDerivSuperType, outertex=nothing)
    return expand_derivative(TEX, D, get_expansion_for(D), outertex)
end


function expand_derivative(TEX::Union{TensorExpression,<:CovDerivSuperType}, outertex::Union{TensorExpression,Nothing}=nothing)
    if TEX isa CovDerivSuperType
        # @show TEX
        # println("is cov deriv")
        D = set_operand!(copy(TEX), zero())
        expanded = _expand_derivative(TEX, D, outertex)
        return expand_derivative(expanded, outertex)   # sweep again to capture any other derivs
    end

    newtex = zero()
    for summand in summands(TEX)
        new_summand_terms = summand[1]
        for term in terms(summand[2])
            if term isa CovDerivSuperType
                # @show term ,"cov"
                D = set_operand!(copy(term), zero())
                new_summand_terms *= expand_derivative(_expand_derivative(term, D, outertex), outertex)
            elseif term isa OperatorSuperType
                # @show term
                if outertex===nothing
                    new_summand_terms *= set_operand!(copy(term), expand_derivative(operand(term), TEX))
                else
                    new_summand_terms *= set_operand!(copy(term), expand_derivative(operand(term), outertex))
                end
            else
                new_summand_terms *= term
            end
        end
        newtex += new_summand_terms
    end

    if any(term isa CovDerivSuperType for summand in summands(newtex) for term in terms(summand[2]))
        return expand_derivative(newtex, outertex)
    end
    return newtex
end

"""Fully expand, apply derivatives and simplify"""
function derive(TEX::Union{TensorExpression,<:CovDerivSuperType})
    #remove prod rule here
    return simplify(apply_derivative(expand_derivative(TEX)))
end
"""Fully expand, apply derivatives and simplify"""
function derive(TEX::T) where T<:GeneratorDerivSuperType
    return simplify(apply_derivative(product_rule(TEX)))
end

"""
Apply derivatives and simplify.

A partial derivative carries no covariant expansion, so `expand_derivative` has nothing to rewrite
and is skipped.
"""
function derive(TEX::T) where T<:PartialSuperType
    return simplify(apply_derivative(TEX))
end

#####
#props
######
"""
Returns which coordinates a derivative is with respect to.
"""
coord(D::DerivativeSuperType) = D.coord


###############################################################
###############################################################

#user-level macro for adding deriv actions
macro action(expr)
    @assert expr.head == :(=)

    lhs, rhs = expr.args
    @assert lhs.head == :call

    D = lhs.args[1]
    arg = lhs.args[2]

    return esc(:(add_derivative_action($D, $arg, $rhs)))
end


#####add a macro to register covariant deriv... or something

macro covariant(deriv::CovDerivSuperType, expansion::Tensorial)
    @assert expr.head == :(=)

    lhs, rhs = expr.args


end


###############################################################
###############################################################

#derivation actions and patterns


derivative_actions = Dict{DerivativeSuperType,Dict{TensorSuperType,TensorExpression}}() #where {D<:}
# derivative_expansions = Dict{CovDerivSuperType,TensorExpression}()

"""
Same as == currently.
"""
function hard_equals(a::DerivativeSuperType, b::DerivativeSuperType)#check to see if derivatives commute
    (a isa PartialDerivative || b isa PartialDerivative) && return a == b
    if name(a) == name(b) && op_indices(a) == op_indices(b) && operand(a) == operand(b) #basic check
        return true
    end
    opra = isTensor(operand(a))
    oprab = isTensor(operand(b))
    if oprab !== nothing && opra !== nothing && opra isa DerivativeSuperType && oprab isa DerivativeSuperType #nested ops
        #check to see if the operand of one derivative is the other, and if the operators are the same (if they commute then this is true)
        #note that i think this is not deep!! if many partials etc
        if op_indices(opra) == op_indices(b) && op_indices(oprab) == op_indices(a) && name(opra) == name(b) && name(oprab) == name(a) && operand(opra) == operand(oprab)
            return derivatives_commute(a::DerivativeSuperType, b::DerivativeSuperType)
        else
            return false
        end
    end
    return false
end


"""
    metric_compatible(D::Type, g::Type) -> Bool

Whether derivatives of type `D` annihilate metrics of type `g` whatever those depend on, as a
Levi-Civita connection does its own metric. False unless declared with [`@compatible`](@ref).

Such a derivative carries the metric out past itself when popping, and `apply_derivative` sets its
action on the metric to zero.
"""
metric_compatible(::Type, ::Type) = false
metric_compatible(O::DerivativeSuperType, g::TensorSuperType) = metric_compatible(typeof(O), typeof(g))

"""
    @compatible D g

Declare derivatives of type `D` compatible with metrics of type `g`, and with any of their
subtypes: it defines [`metric_compatible`](@ref) to be true for them.

    @compatible Derivative Metric          # ∇g = 0
    @compatible CovDerivSuperType Metric   # every covariant derivative
"""
macro compatible(args...)
    # `@compatible D, g` arrives as a single tuple expression rather than as two arguments
    if length(args) == 1 && args[1] isa Expr && args[1].head === :tuple
        args = Tuple(args[1].args)
    end
    length(args) == 2 || error("@compatible takes a derivative type and a metric type; got $(length(args)) arguments.")
    D, g = esc(args[1]), esc(args[2])
    return quote
        $(_check_compatible)($D, $g)
        (::typeof($metric_compatible))(::Type{<:$D}, ::Type{<:$g}) = true
        nothing
    end
end

function _check_compatible(D, g)
    (D isa Type && D <: DerivativeSuperType) ||
        throw(ArgumentError("@compatible: $D is not a derivative type."))
    (g isa Type && g <: MetricSuperType) ||
        throw(ArgumentError("@compatible: $g is not a metric type."))
    return nothing
end

"""
    commutes_with_metric(O, g) -> Bool

Whether the metric `g` may be carried out past the derivative `O`: whether `O(g·T) = g·O(T)`, which
holds exactly when `O` annihilates `g`: `O` is `metric_compatible` with it, or `g` does not depend
on what `O` differentiates. `pop_metric` consults this before moving a metric it has popped off an
operand out through the derivative.
"""
commutes_with_metric(O::DerivativeSuperType, g) =
    metric_compatible(O, g) || isdisjoint(function_of(g), coord(O))
# commutes_with_metric(::DerivativeSuperType, g) = false

function ==(a::DerivativeSuperType, b::DerivativeSuperType)#check to see if derivatives commute
    #partials: the whole chain up to allowed swaps, so that equality is transitive (and hash agrees)
    if a isa PartialDerivative || b isa PartialDerivative
        return a isa PartialDerivative && b isa PartialDerivative && _chain_signed_equals(a, b, false) == (true, 1)
    end
    # @show "using dervi"
    if name(a) == name(b) && op_indices(a) == op_indices(b) && operand(a) == operand(b) #basic check
        # @show name(a) == name(b)
        # @show  op_indices(a) == op_indices(b)
        # @show operand(a) == operand(b) 
        # println("found equal...")
        return true
    end
    # @show "not equal yet"
    opra = isTensor(operand(a))
    oprab = isTensor(operand(b))
    if oprab !== nothing && opra !== nothing && opra isa DerivativeSuperType && oprab isa DerivativeSuperType #nested ops
        #check to see if the operand of one derivative is the other, and if the operators are the same (if they commute then this is true)
        #note that i think this is not deep!! if many partials etc
        if op_indices(opra) == op_indices(b) && op_indices(oprab) == op_indices(a) && name(opra) == name(b) && name(oprab) == name(a) && operand(opra) == operand(oprab)
            return false#derivatives_commute(a::DerivativeSuperType, b::DerivativeSuperType)
        else
            return false
        end
    end
    return false
end
"""
A partial is hashed by its whole chain, whose order equality ignores, so the hashes are summed.

"""
function Base.hash(d::DerivativeSuperType, h::UInt)
    chn = d isa PartialDerivative ? get_operator_chain(d, PartialDerivative) : (d,)
    acc = zero(UInt)
    for op in chn
        hh = hash(name(op))
        @inbounds for p in op_indices(op)
            hh = hash(p[1].uid, hash(p[2].position, hh))
        end
        acc += hh
    end
    return hash(acc, h)
end

"""
Checks if two derivative expressions are the same, up to an overall sign.
"""
function signed_equals(a::DerivativeSuperType, b::DerivativeSuperType)#check to see if derivatives commute
    if a isa PartialDerivative || b isa PartialDerivative
        (a isa PartialDerivative && b isa PartialDerivative) || return (false, 0)
        return _chain_signed_equals(a, b)
    end
    if name(a) == name(b) && op_indices(a) == op_indices(b) && operand(a) == operand(b) #basic check
        return (true, 1)
    end
    opra = isTensor(operand(a))
    oprab = isTensor(operand(b))
    if oprab !== nothing && opra !== nothing
        if opra isa DerivativeSuperType && oprab isa DerivativeSuperType #nested ops
            #check to see if the operand of one derivative is the other, and if the operators are the same (if they commute then this is true)
            #note that i think this is not deep!! if many partials etc
            if op_indices(opra) == op_indices(b) && op_indices(oprab) == op_indices(a) && name(opra) == name(b) && name(oprab) == name(a) && operand(opra) == operand(oprab)
                return derivatives_commute_signed(a::DerivativeSuperType, b::DerivativeSuperType)
            else
                return (false, 0)
            end
        else
            !(opra isa DerivativeSuperType) && !(oprab isa DerivativeSuperType)
            val = same_tensor_mod_dummies(opra, oprab)
            return (Bool(abs(val)), val)
        end
    end
    return (false, 0)
end

function signed_equals(a::OperatorSuperType, b::OperatorSuperType)#check to see if ops commute
    if name(a) == name(b) && op_indices(a) == op_indices(b)
        if operand(a) == operand(b) #basic check
            return (true, 1)
        elseif isTensor(operand(a))!==nothing && isTensor(operand(b))!==nothing
            val = same_tensor_mod_dummies(isTensor(operand(a)), isTensor(operand(b)))
            return (Bool(abs(val)), val)
        end
    end
    # println("checking operatos")
    return (false, 0)
end

"""
Commutation of derivative operators. While partials automatically commute, they do not nec. commute with covariant derivatives.
"""
mutable struct CovariantCommutator
    grading::Int #±1 for commutator or anticommutator
    A::DerivativeSuperType
    B::DerivativeSuperType
    C::TensorExpression
end

derivative_algebra = Dict{Tuple{CommutatorKey,CommutatorKey},CovariantCommutator}() #integer is 1 for comm, -1 for anticommut, expression is the result of the commutator/anticommutator (eg for fermionic derivatives)

# ==(C1::Commutator,C2::Commutator)= f1==f2 for fields in C1,C2
==(a::CovariantCommutator, b::CovariantCommutator) = begin
    a.grading == b.grading &&
        a.A == b.A &&
        a.B == b.B &&
        a.C == b.C
end

"""
Non-commuting derivatives.
"""
function add_commutator!(grading::Integer, A::DerivativeSuperType, B::DerivativeSuperType, C::TensorExpression)
    if free_indices(A * B) != free_indices(C) && !isZeroTens(C)
        throw(InvalidCommutatorDefinition(A, B, C))
    end
    derivative_algebra[minmax(_ckey(A), _ckey(B))] = CovariantCommutator(grading, A, B, C)
end
add_commutator!(grading::Integer, A::DerivativeSuperType, B::DerivativeSuperType, C::Int64) = add_commutator!(grading, A, B, TensorExpression(C))
function add_commutator!(va::CovariantCommutator)
    k=minmax(_ckey(va.A), _ckey(va.B))
    derivative_algebra[k]=va
end

function get_commutator(A::DerivativeSuperType, B::DerivativeSuperType)
    key=minmax(_ckey(A), _ckey(B))
    if haskey(derivative_algebra, key)
        comrel=derivative_algebra[key]
        dictout=Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
        sign=1
        if _ckey(comrel.A)==_ckey(A)
            merge!(dictout, get_pattern(comrel.A, A))
            merge!(dictout, get_pattern(comrel.B, B))
        else
            sign*=(-comrel.grading)#fermions: {B,A}={A,B}, bosons [A,B]=-[B,A]
            merge!(dictout, get_pattern(comrel.A, B))
            merge!(dictout, get_pattern(comrel.B, A))
        end

        out=sign*set_from_pattern(comrel.C, dictout)
        return (comrel.grading, out)

    else #no commutation rels defined, no swap allowed! (different to matrix behaviour)
        return (0, zero())
    end
end

function get_pattern(Asym::DerivativeSuperType, Aact::DerivativeSuperType)
    dictout=Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
    for (i, j) in zip(indices(Asym), indices(Aact))
        dictout[i]=j
    end
    return dictout
end

function commute(D::DerivativeSuperType)
    oprd=isTensor(operand(D))
    if #=oprd!==nothing &&=#oprd isa DerivativeSuperType
        # keyy=minmax(name(D),name(oprd))
        (sign, swout)=get_commutator(D, oprd)
        if sign!=0
            oprdoprd=operand(oprd)
            Dp=copy(D)
            set_operand!(Dp, oprdoprd)
            oprdout=copy(oprd)
            set_operand!(oprdout, Dp)

            return sign*oprdout+set_operand_deep!(swout, oprdoprd)
        end
    end
    return D #do nothing
end

function commute_right(D::DerivativeSuperType, N)
    D=product_rule(D)
    if isTensor(D)!==nothing
        D=isTensor(D)
        if D isa DerivativeSuperType
            oprd=isTensor(operand(D))
            if #=oprd!=false &&=#oprd isa DerivativeSuperType
                # keyy=minmax(name(D),name(oprd))
                (sign, swout)=get_commutator(D, oprd)
                if sign!=0
                    oprdoprd=operand(oprd)
                    Dp=copy(D)
                    set_operand!(Dp, oprdoprd)
                    oprdout=copy(oprd)
                    # @show oprdout
                    # Dpo=Dp
                    Dp=commute_right(Dp, N)
                    # @show Dpo
                    # @show Dp
                    # Dp===nothing && (Dp =  zero())
                    set_operand!(oprdout, Dp)
                    return sign*oprdout+set_operand_deep!(swout, oprdoprd)
                end
            end
        end
    else
        return commute_right(D, N)
    end
    return D
end


function commute_right(tt::Tensorial, N)
    transform(t -> commute_right(t, N), tt)
end

commute_right(tt::TensorSuperType, N) = tt



"""
Checks to see if derivatives commute, up to a sign.
"""
function derivatives_commute_signed(D1::DerivativeSuperType, D2::DerivativeSuperType)
    if typeof(D1) <: PartialDerivative && typeof(D2) <: PartialDerivative
        _partials_swap(D1, D2) || return (false, 0)
        #use grading of the derivative itself, not the grading of full object since we commute ops not expressions
        if (op_grading(D1) * op_grading(D2) == 0) #commute if bb or fb
            return (true, 1) #regular partial derivatives commute
        else
            return (true, -1) #fermionic derivatives anticommute... add anti behaviour later
        end
    elseif haskey(derivative_algebra, D1 => D2) #non partials
        return derivative_algebra[D1=>D2][1] == (1, zero())
    else
        return (false, 0)
    end
end

function derivatives_commute(D1::DerivativeSuperType, D2::DerivativeSuperType)
    return derivatives_commute_signed(D1, D2) == (true, 1)
end

"""
Whether two partials may be exchanged. An index at the default position hides a metric,
∂^a = g^{ab}∂_b, so the partials swap only if each hidden metric passes through the other
derivative, as `commutes_with_metric` decides for the metric the set builds.
"""
function _partials_swap(D1::DerivativeSuperType, D2::DerivativeSuperType)
    for (D, E) in ((D1, D2), (D2, D1)), (i, pos) in op_indices(D)
        S = i.index_set
        (pos == S.default_position && S.metric_type !== nothing) || continue   #nothing hidden
        (metric_compatible(typeof(E), S.metric_type) || isdisjoint(metric_function(S), coord(E))) || return false
    end
    return true
end

"""
    _chain_signed_equals(a, b, mod_dummies=true) -> (Bool, sign)

Whether two partials are equal up to a sign. The chains are the same up to allowed swaps
(`_partials_swap`) when they hold the same derivatives and every pair that may not swap keeps its
relative order
"""
function _chain_signed_equals(a::PartialDerivative, b::PartialDerivative, mod_dummies::Bool=true)
    u, v = get_operator_chain(a, PartialDerivative), get_operator_chain(b, PartialDerivative)
    n = length(u)
    n == length(v) || return (false, 0)
    #the k-th occurrence of a derivative in u is the k-th in v: repeats then never cross
    p = zeros(Int, n)
    used = falses(n)
    for i in 1:n
        j = findfirst(k -> !used[k] && name(u[i]) == name(v[k]) && op_indices(u[i]) == op_indices(v[k]), 1:n)
        j === nothing && return (false, 0)
        used[j] = true
        p[i] = j
    end
    sign = 1
    for i in 1:n, j in (i + 1):n
        p[i] < p[j] && continue
        _partials_swap(u[i], u[j]) || return (false, 0)
        sign *= (-1)^(op_grading(u[i]) * op_grading(u[j]))
    end
    ra, rb = operand(u[end]), operand(v[end])
    ra == rb && return (true, sign)
    mod_dummies || return (false, 0)
    ta, tb = isTensor(ra), isTensor(rb)
    if ta !== nothing && tb !== nothing && !(ta isa DerivativeSuperType) && !(tb isa DerivativeSuperType)
        val = same_tensor_mod_dummies(ta, tb)
        return isZero(val) ? (false, 0) : (true, sign * val)
    end
    return (false, 0)
end


function same_operator(D1::DerivativeSuperType, D2::DerivativeSuperType)
    name(D1) == name(D2) || return (false => _NO_MATCH)
    ind1 = op_indices(D1);
    ind2 = op_indices(D2)
    length(ind1) == length(ind2) || return (false => _NO_MATCH)

    # cheap structural pass: positions and index sets must agree
    @inbounds for i in eachindex(ind1)
        (ind1[i][2].position == ind2[i][2].position &&
         ind1[i][1].index_set === ind2[i][1].index_set) || return (false => _NO_MATCH)
    end

    # only now build the map
    doubles = Dict{Index,Index}()
    @inbounds for i in eachindex(ind1)
        a, b = ind1[i][1], ind2[i][1]
        if haskey(doubles, a)
            doubles[a] == b || return (false => _NO_MATCH)
        else
            doubles[a] = b
        end
    end
    return (true => doubles)
end

function add_derivative_action(D::T, exp, res::E) where {T<:DerivativeSuperType,E<:SymorNum}
    add_derivative_action(D, exp, TensorExpression(res))
end
"""
Define action of a derivative operator on some basis of tensors.
"""
function add_derivative_action(D::T, exp::TensorTerm, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}
    length(terms(exp)) == 1 || throw(InvalidDerivativeBasis(D, exp))
    return add_derivative_action(D, terms(exp)[1], res)
end

function add_derivative_action(D::T, exp::TensorExpression, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}
    (length(summands(exp)) == 1 && length(terms(term(summands(exp)[1]))) == 1) ||
        throw(InvalidDerivativeBasis(D, exp))
    return add_derivative_action(D, terms(term(summands(exp)[1]))[1], res)
end
function add_derivative_action(D::T, exp::TensorSuperType, res::E) where {T<:DerivativeSuperType,E<:Union{TensorSuperType,TensorExpression,TensorTerm}}
    key = nothing
    in_dict = false
    pattern_ref = Dict()
    res = TensorExpression(res)
    for ky in keys(derivative_actions)
        out = same_operator(D, ky)
        B = out[1]
        o = out[2]
        if B == true
            key = ky
            in_dict = true
            pattern_ref = o
            break
        end
    end

    #check to make sure action not already defined, else add

    #now we have index pattern if we found a match
    if in_dict
        action_dict = derivative_actions[key] #all registered derivatives to act on
        for k in keys(action_dict)#looks through all tensors
            if same_tensor_mod_pattern(exp, k, pattern_ref)[1]
                return false #already in the dict
            end
        end
        #no matches found.... proceed

        #need to make sure correct indices
        pat_adj_tensor = index_to_pattern(exp, pattern_ref)
        pat_adj_out = index_to_pattern(res, pattern_ref)
        derivative_actions[key][pat_adj_tensor] = pat_adj_out #adds the expression
        return true
    else
        derivative_actions[D] = Dict()
        derivative_actions[D][exp] = res #adds the expression
        return true
    end
end


#given an expression/tensor rename/relabel indices according to a particular pattern

function index_to_pattern(T1::TensorExpression, pattern_ref::Dict{Index,Index})
    new_exp = zero()
    for summand in summands(T1)
        new_summand = summand[1]
        for term in terms(summand[2])
            new_summand *= index_to_pattern(term, pattern_ref)
        end
        new_exp += new_summand
    end
    return new_exp
end


function index_to_pattern(T1::TensorSuperType, pattern_ref::Dict{Index,Index})
    newT = copy_owned(T1)
    ind1 = indices(newT)
    for ind in eachindex(ind1)
        if haskey(pattern_ref, ind1[ind][1])
            new_ind = pattern_ref[ind1[ind][1]]
            ind1[ind] = new_ind => ind1[ind][2]
        end
    end
    return newT
end


# checks if tensors are the same in terms of structure of indices, ignoring the specfic name
# eg differentiating between contracted and non contracted
function same_tensor_mod_pattern(T1::T, T2::T_2, index_dict::Dict{Index,Index}) where {T<:TensorSuperType,T_2<:TensorSuperType}
    name(T1) == name(T2) || return (false, _NO_MATCH)

    ind1 = indices(T1)
    ind2 = indices(T2)
    n = length(ind1)
    n == length(ind2) || return (false, _NO_MATCH)

    # cheap pass: positions and index sets must agree, no allocation
    @inbounds for i in 1:n
        ind1[i][2] == ind2[i][2] || return (false, _NO_MATCH)
        ind1[i][1].index_set === ind2[i][1].index_set || return (false, _NO_MATCH)
    end

    # second pass: build the free-index correspondence
    temp_dict = Dict{Index,Index}()
    sizehint!(temp_dict, n)
    @inbounds for i in 1:n
        a, b = ind1[i][1], ind2[i][1]
        if haskey(index_dict, a)
            index_dict[a] == b || return (false, _NO_MATCH)
        elseif haskey(temp_dict, a)
            temp_dict[a] == b || return (false, _NO_MATCH)
        else
            temp_dict[a] = b
        end
    end
    return (true, temp_dict)
end


##############################################################################################################################
##############################################################################################################################

#applying derivatives

#leibniz out the operands of all derivatives in an expression

function apply_derivative(TEX::TensorExpression)
    s = summands(TEX)
    isempty(s) && return TEX

    found = false
    for (_, tt) in s, t in terms(tt)
        if t isa DerivativeSuperType
            found = true;
            break
        end
    end
    found || return TEX

    out = CoeffPair[]
    sizehint!(out, length(s))
    for (coef, tt) in s
        nt = coef
        for t in terms(tt)
            if t isa DerivativeSuperType
                d = copy(t)
                set_operand!(d, apply_derivative(operand(d)))
                nt *= apply_derivative(d)
            else
                nt *= t
            end
        end
        isZeroTens(nt) && continue
        append!(out, summands(TensorExpression(nt)))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end
# A single product, or a lone tensor, is not a `TensorExpression`, so neither
# matched a method and an identity written as one product raised a MethodError
# while the same identity written as a difference went through. Both delegate,
# so the semantics are the ones the expression method already has: a term with
# no derivative in it comes back unchanged.
apply_derivative(TT::TensorTerm) = apply_derivative(TensorExpression(TT))
apply_derivative(T::TensorSuperType) = apply_derivative(TensorExpression(T))

#apply the derivative and cancel terms that are annihilated
function apply_derivative(D::DerivativeSuperType)
    Dx = product_rule(D)
    s = summands(Dx)
    isempty(s) && return zero()

    out = CoeffPair[]
    sizehint!(out, length(s))
    factors = TensorSuperType[]
    for (coef, tt) in s
        ts = terms(tt)
        empty!(factors)
        sizehint!(factors, length(ts))
        c = coef
        simple = true
        annihilated = false
        for tensor in ts
            if tensor isa DerivativeSuperType
                tens = terms(term(summands(operand(tensor))[1]))[1]
                fo = function_of(tens)
                hit = false
                for cc in coord(tensor)
                    if cc in fo
                        ;
                        hit = true;
                        break;
                    end
                end
                #∂ of what it does not depend on, or a compatible metric: the whole summand vanishes
                if !hit || metric_compatible(tensor, tens)
                    annihilated = true
                    break
                end
                res = perform_derivative_action(tensor)
                if res isa TensorSuperType
                    push!(factors, res)
                else
                    sm = summands(TensorExpression(res))
                    if length(sm) == 1
                        c = c * sm[1][1]
                        append!(factors, terms(sm[1][2]))
                    else
                        simple = false
                        break
                    end
                end
            else
                push!(factors, tensor)
            end
        end
        annihilated && continue

        local oterm
        if simple
            term_out = TensorTerm(copy(factors))
            _check_term_indices(term_out)
            oterm = c * term_out
        else
            oterm = coef                       # rare: rebuild the general (quadratic) way
            for tensor in ts
                oterm *= (tensor isa DerivativeSuperType ?
                          perform_derivative_action(tensor) : tensor)
            end
        end
        isZeroTens(oterm) && continue
        append!(out, summands(oterm))
    end
    isempty(out) && return zero()
    return TensorExpression(_RAW, out)
end

"""
Performs the actual derivative operation on tensors for which it is defined.
Assumes already has used product rule.
If there are rules (as defined by add_derivative_action), then the derivative will be performed, else expressions will be left alone eg ∇B
"""
function perform_derivative_action(D::DerivativeSuperType)
    oprd = terms(summands(operand(D))[1][2])[1]
    dname = name(D)
    oname = name(oprd)

    for (ky, action_dict) in derivative_actions
        name(ky) == dname || continue                 # cheap reject before same_operator
        matched, o = same_operator(D, ky)
        matched || continue
        for (k, replacement) in action_dict
            name(k) == oname || continue              # cheap reject before the full comparison
            ok, om = same_tensor_mod_pattern(oprd, k, o)
            ok || continue
            # build the inverse map once, in place
            all_invs = Dict{Index,Index}()
            sizehint!(all_invs, length(o) + length(om))
            for (a, b) in o
                ;
                all_invs[b] = a;
            end
            for (a, b) in om
                ;
                all_invs[b] = a;
            end     # om entries win, matching merge(inv_o, inv_om)
            return index_to_pattern(replacement, all_invs)
        end
    end
    return D
end

"""
expand_derivative(TEX,D,DExpansion)

Replaces a symbolic derivative with its full form, eg substituting a covariant derivative D with all of its connection terms.
'TEX' is the original expression. 'D' is the derivative to replace, and 'DExpansion' is a tensor expression that will be a sum of tensors/operators.
For each operator in 'DExpansion', the operand will be replaced by the matching operand in 'TEX'.
"""
function expand_derivative(TEX::Union{TensorExpression,<:DerivativeSuperType}, D::DerivativeSuperType, DExpansion::TensorExpression, outertex=nothing)
    #### keeping all the comments here in case something breaks in future
    if (free_indices(D) != free_indices(DExpansion))
        @show free_indices(D)
        @show free_indices(DExpansion)
        println(D)
        println(DExpansion)
        println(TEX)
        println("It appears that your expansion is not compatible.")
        throw(FreeIndexMismatch())
    end

    #think about: hellish dummie indice relabeling//////

    new_TEX = zero() #construct new expression term by term
    TEX = 1 * TEX #ensures tensor exp
    avail = Dict{IndexSet,Vector{IndexSuperType}}() #pairs of index sets, with the unused indices
    dummymap = Dict{String,Index}()#dict of what to replace pattern indices with
    for summand in TEX.summands
        new_summand = summand[1]
        for tensor in summand[2].terms #look at each tensor in the tensor term
            if typeof(tensor) <: DerivativeSuperType #apply to operands first recursively
                if outertex === nothing
                    set_operand!(tensor, expand_derivative(operand(tensor), D, DExpansion, TEX))
                else
                    set_operand!(tensor, expand_derivative(operand(tensor), D, DExpansion, outertex))
                end
                #feed in term... need to make sure we dont reuse outer indices
                #recursively expand the operand if it has the same derivative in it
                update_ext_indices!(tensor)
            end
            if typeof(tensor) == typeof(D) && name(tensor) == name(D) #operator found
                oprd = operand(tensor)
                mapping = construct_index_rep_map(D, tensor) #index map so we use correct indices in expansion
                replace = copy(DExpansion) #new copy of expansion to replace in our exp
                #figure out what indices we're allowed to use in the replacement
                sets = get_index_sets(tensor) #Set{IndexSet}()
                if outertex === nothing
                    ref = TEX
                else
                    ref = outertex
                end
                for set in sets

                    a1 = get_unused_indices(ref, set)              # already a Vector now
                    b = get_unused_indices(operand(tensor), set)
                    a = filter(in(b), a1)                          # vector intersection
                    push!(avail, set => a) #ones we can replace with
                end
                #construct the replacement
                new_rep_summand = zero()
                for repsummand in replace.summands
                    new_rep_term = repsummand[1]
                    ###############fix index logic for multiple tensors in the expressions
                    ## currently, K^{*β}_{α*}∂^{i*}_{*β}[[0*]]   => goes to K^{*β}_{α*}∂^{i*}_{*β}[[F^{ik**}_{**jβ}]], where F^{ik**}_{**jβ}
                    found_op = false
                    for reptensor in repsummand[2].terms ##add somewhere to handle a term with no op (like a YM covar).
                        # @show reptensor

                        ##println(reptensor)
                        if typeof(reptensor) <: OperatorSuperType
                            found_op = true
                            # #if we found an operator,
                            set_operand!(reptensor, oprd)#set to correct operand
                            opinds = copy(op_indices(reptensor))          # one copy
                            #dums=dummies(reptensor)#gets the placeholder dummy indices from the pattern
                            for pair in eachindex(opinds) # replace dummies with correct etc
                                if haskey(mapping, opinds[pair][1].name) #if the index is fixed already from before
                                    opinds[pair] = mapping[opinds[pair][1].name] => opinds[pair][2]#should replace with the right thing
                                elseif haskey(dummymap, opinds[pair][1].name)
                                    opinds[pair] = dummymap[opinds[pair][1].name] => opinds[pair][2]

                                else #if just a placeholder index in the pattern, take from the list
                                    if !haskey(avail, opinds[pair][1].index_set)
                                        #if we're using a dummy index type thats not in the orignal expression
                                        #changed to ref!!
                                        a11 = get_unused_indices(ref, opinds[pair][1].index_set)
                                        b1 = get_unused_indices(operand(tensor), opinds[pair][1].index_set)
                                        aa = filter(in(b1), a11)
                                        new_index = safe_pop!(aa, opinds[pair][1].index_set)
                                    else
                                        new_index = safe_pop!(avail[opinds[pair][1].index_set], opinds[pair][1].index_set) # get an unused index
                                    end
                                    dummymap[opinds[pair][1].name] = new_index #makes sure summation dummies do the right thing
                                    opinds[pair] = new_index => opinds[pair][2]
                                end
                                # for repindex in reptensor.
                            end
                            set_op_indices!(reptensor, opinds)
                        else
                            inds = copy(indices(reptensor))
                            for pair in eachindex(inds) # replace dummies with correct etc
                                if haskey(mapping, inds[pair][1].name) #if the index is fixed already from before
                                    inds[pair] = mapping[inds[pair][1].name] => inds[pair][2]#should replace with the right thing
                                elseif haskey(dummymap, inds[pair][1].name)
                                    inds[pair] = dummymap[inds[pair][1].name] => inds[pair][2]

                                else #if just a placeholder index in the pattern, take from the list
                                    if !haskey(avail, inds[pair][1].index_set)
                                        #if we're using a dummy index type thats not in the orignal expression
                                        #changed to ref!!!
                                        a11 = OrderedSet{Index}(get_unused_indices(ref, inds[pair][1].index_set))
                                        b1 = OrderedSet{Index}(get_unused_indices(operand(tensor), inds[pair][1].index_set)) ##will check for dummies inside...
                                        aa = intersect(a11, b1)
                                        new_index = safe_pop!(aa, inds[pair][1].index_set)
                                        # new_index = pop!(get_unused_indices(ref, reptensor.indices[pair][1].index_set))
                                    else
                                        new_index = safe_pop!(avail[inds[pair][1].index_set], inds[pair][1].index_set) # get an unused index
                                    end
                                    #new_index = pop!(avail[reptensor.indices[pair][1].index_set]) # get an unused index
                                    dummymap[inds[pair][1].name] = new_index #makes sure summation dummies do the right thing
                                    inds[pair] = new_index => inds[pair][2]


                                end

                            end
                            set_indices!(reptensor, inds)
                        end
                        update_indices!(reptensor)
                        new_rep_term *= reptensor
                    end
                    # update_ext_indices!(new_rep_term)
                    # @show new_rep_term
                    # update_ext_indices!(new_rep_summand)
                    if found_op == false
                        new_rep_term *= oprd #if there was no operator, we need to add the operand at the end, since the pattern is just for the operator part
                    end
                    new_rep_summand = new_rep_summand ⨥ new_rep_term
                end
                new_summand *= new_rep_summand
            else#normal tensor
                new_summand *= tensor
            end
        end
        new_TEX += new_summand #*summand[1] * #might be squaring...
    end
    return new_TEX
end

# graded leibniz
"""
Perform the graded Leibniz product rule for all derivatives in the expression.
"""
function product_rule(TEX::TensorExpression)
    s = summands(TEX)
    isempty(s) && return TEX

    # nothing to expand → return unchanged
    found = false
    for (_, tt) in s, t in terms(tt)
        if t isa LeibnizOperatorSuperType
            found = true;
            break
        end
    end
    found || return TEX

    out = CoeffPair[]
    sizehint!(out, length(s))
    for (coef, tt) in s
        nterm = coef
        for t in terms(tt)
            nterm *= (t isa LeibnizOperatorSuperType ? product_rule(t) : t)
        end
        isZeroTens(nterm) && continue
        append!(out, summands(TensorExpression(nterm)))
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end

function product_rule(O::LeibnizOperatorSuperType)
    set_operand!(O, product_rule(operand(O)))
    isderiv = O isa DerivativeSuperType
    cO = coord(O)
    s = summands(operand(O))

    # out = CoeffPair[]
    out = CoeffPair[]
    sizehint!(out, sum(length(terms(t)) for (_, t) in s; init=0))   # operand may be zero
    for (coef, tterm) in s
        ts = terms(tterm)
        grad = 0
        for tensorind in eachindex(ts)
            t = ts[tensorind]

            annihilates = false
            if isderiv
                fo = function_of(t)
                annihilates = true
                for c in cO
                    if c in fo
                        annihilates = false
                        break
                    end
                end
                annihilates |= metric_compatible(O, t)
            end
            if !annihilates
                grd = (-1)^grad
                new_terms = copy(ts)
                new_terms[tensorind] = _with_operand(O, TensorExpression(t))
                contribution = (grd * coef) * TensorTerm(new_terms)
                isempty(summands(contribution)) || append!(out, summands(contribution))
            end
            grad = (grad + grading(t) * op_grading(O)) % 2
        end
    end
    isempty(out) && return zero()
    return TensorExpression(out)
end

# fast: construct with the new operand directly with no copy
_with_operand(O::DerivativeSuperType, newop::TensorExpression) =
    typeof(O)(OperatorData(op_data(O), op_indices(O), newop), coord(O))

_with_operand(O::OperatorSuperType, newop::TensorExpression) =
    typeof(O)(OperatorData(op_data(O), op_indices(O), newop))

_with_operand(O::LieBracketSuperType, newop::TensorExpression) =
    typeof(O)(O.A, O.grading, OperatorData(op_data(O), op_indices(O), newop))

_with_operand(O::VariationalSuperType, newop::TensorExpression) =
    typeof(O)(OperatorData(op_data(O), op_indices(O), newop), O.var)
"""
Put every chain of nested partial derivatives into its canonical order.

"""
function symmetrise_partials(TEX::TensorExpression)
    s = summands(TEX)
    isempty(s) && return TEX

    # nothing to reorder, return unchanged, no rebuild
    found = false
    @inbounds for k in eachindex(s)
        for t in terms(s[k][2])
            if t isa PartialDerivative
                found = true;
                break
            end
        end
        found && break
    end
    found || return TEX

    out = zero()
    for (c, tt) in s
        out += c * symmetrise_partials(tt)
    end
    return out
end

function symmetrise_partials(TT::TensorTerm)
    ts = terms(TT)
    any(t -> t isa PartialDerivative, ts) || return TensorExpression(TT)#so we always have some output type to help the compiler

    new_term = 1
    for term in ts
        new_term *= (term isa PartialDerivative ? symmetrise_partials(term) : term)
    end
    return TensorExpression(new_term)
end

"""
`D` with its chain of partials in canonical order, times the graded sign of the reordering, and the
chains inside what it acts on likewise.
"""
function symmetrise_partials(D::PartialDerivative)
    chn = get_operator_chain(D, PartialDerivative)
    root = symmetrise_partials(operand(chn[end]))
    frees = Set(p[1] for p in free_indices_ordered(D))
    keys = [(name(op), [(i.index_set.name, pos.position, i in frees ? i.name : "") for (i, pos) in op_indices(op)])
            for op in chn]
    remaining = collect(eachindex(chn))
    order = Int[]
    while !isempty(remaining)
        best = 0
        for (r, k) in enumerate(remaining)
            #it may move to the front only past partials it swaps with
            all(j -> _partials_swap(chn[j], chn[k]), @view remaining[1:(r - 1)]) || continue
            (best == 0 || isless(keys[k], keys[best])) && (best = k)
        end
        push!(order, best)
        deleteat!(remaining, findfirst(==(best), remaining))
    end
    return graded_sign(order, chn) * TensorExpression(rebuild_nested(chn[order], root))
end

"""Innermost operand of a nest of derivatives, partial or covariant."""
function root_oprnd(D::DerivativeSuperType)
    opr = isTensor(operand(D))
    if opr isa DerivativeSuperType
        return root_oprnd(opr)
    else
        return operand(D)
    end
end


function update_operator_data!(T::TensorExpression)
    for summand in summands(T)
        for term in terms(summand[2])
            if typeof(term) <: OperatorSuperType
                update_operator_data!(term)
            end
        end
    end
end

#yeah deep right?? didnt finish... TODO?
function update_operator_data!(O::OperatorSuperType)
    for summand in operand(O).summands
        for term in terms(summand[2])
            if typeof(term) <: OperatorSuperType
                update_operator_data!(term)
            end
        end
    end
end



"""
    rebuild_nested(ops, root)

Rebuild a nested chain of operators over `root`, `ops` given outermost-first. Inverse of
`strip_derivatives`, and the way back from a permuted chain in `symmetrise_partials`.

`root` may be a tensor or an expression.
"""
function rebuild_nested(ops, root)
    out = root
    for k in length(ops):-1:1
        out = _with_operand(ops[k], out isa TensorExpression ? out : TensorExpression(out))
    end
    return out
end

function graded_sign(perm, ops)
    sign = 1
    for i in 1:length(perm)
        for j in (i+1):length(perm)
            # inversion in permutation
            if perm[i] > perm[j]
                gi = grading(op_indices(ops[perm[i]]))
                gj = grading(op_indices(ops[perm[j]]))
                sign *= (-1)^(gi * gj)
            end
        end
    end
    return sign
end

function swap_to_perm(N, swap)
    perm = collect(1:N)
    i, j = swap
    perm[i], perm[j] = perm[j], perm[i]
    return perm
end

function graded_swap(ops, perm)
    perm2 = swap_to_perm(length(ops), perm)
    perm_full = ops[perm2]
    sign = graded_sign(perm, ops)
    return (perm_full, sign)
end

"""
    get_operator_chain(T, through=OperatorSuperType)

Collect the nested operators of `T`, descending only through operators that are instances of
`through`. Used to find chains of nested partials derivs over which to symmetrise.
"""
function get_operator_chain(T::OperatorSuperType, through::Type=OperatorSuperType)
    chain = Any[T]                        # the head is always included, so the chain is never empty
    current = T
    while current isa through
        nxt = isTensor(operand(current))  #get the operand tensor
        (nxt isa OperatorSuperType && nxt isa through) || break
        push!(chain, nxt)
        current = nxt
    end
    return chain
end
function update_ext_indices!(O::OperatorSuperType) #probably needs to include dummies...
    if operand(O) isa TensorSuperType || operand(O) isa TensorTerm
        set_indices!(O, vcat(operator_indices(O), indices(operand(O))))
    elseif operand(O) isa TensorExpression
        if !isempty(operand(O).summands)
            set_indices!(O, vcat(op_indices(O), indices(operand(O).summands[1][2])))
        else
            set_indices!(O, op_indices(O))
        end
    end
end

function *(TEx::TensorExpression, T::OperatorSuperType)
    # println("doing this sort of multiplcation5")
    return update_indices!(TEx * TensorExpression(T))
end

function *(T::OperatorSuperType, TEx::TensorExpression)
    # println("doing this sort of multiplcation6")
    return update_indices!(TensorExpression(T) * TEx)
end

# function copy(T::DerivativeSuperType)
#     # return typeof(T)((isa(getfield(T, f), TensorDataSuperType) ? copy(getfield(T, f),:copy) : deepcopy(getfield(T, f)) for f in fieldnames(typeof(T)))...)
#     return typeof(T)((deepcopy(getfield(T, f)) for f in fieldnames(typeof(T)))...)
# end
# function copy(T::DerivativeSuperType)
#     return typeof(T)((symcopy(getfield(T, f)) for f in fieldnames(typeof(T)))...)
# end
#faster w/ concrete fields
copy(T::DerivativeSuperType) = typeof(T)(symcopy(T.data), symcopy(T.coord))
