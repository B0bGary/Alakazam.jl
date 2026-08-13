#########
## functions on tensors that are diff. forms
#########

function isForm(T::TensorSuperType, index_type::Union{IndexSet,Nothing}=nothing)
    if index_type===nothing
        return isAntisymmetric(T,"spacetime") 
    else
        return isAntisymmetric(T,index_type.name) 
    end
end

"""Take the Hodge star of a form.
This operation will automatically be taken with respect to spacetime indices, if they are defined on the tensor.
"""
star(T::TensorSuperType, index_type::Union{IndexSet,Nothing}=nothing) = ⋆(T, index_type)

"""Take the Hodge star of a form.
This operation will automatically be taken with respect to spacetime indices, if they are defined on the tensor.
"""
function ⋆(T::TensorSuperType, index_type::Union{IndexSet,Nothing}=nothing)
    if index_type === nothing && length(indices(T)) == 0
        throw(InvalidDifferentialForm("The tensor " * tens2str(T) * " has no indices, please specify the spacetime indices for the Hodge star."))
    end
    if !isForm(T,index_type) ### actually... need to do a different check because we may have flavour indices
        throw(InvalidDifferentialForm("The tensor " * tens2str(T) * " should be totally antisymmetric."))
    end
    has_spacetime = false
    mismatch = false
    tinds = nothing
    if index_type === nothing
        for (ind, pos) in indices(T)
            is = ind.index_set
            if is.name == "spacetime" #automatically detect if spacetime indices exist
                has_spacetime = true
                tinds = is
                break
            end
            if tinds === nothing
                tinds = is
            elseif tinds != ind.index_set
                mismatch = true
            end
        end
    else
        tinds = index_type
    end
    if has_spacetime == false && mismatch == true
        throw(InvalidDifferentialForm("The tensor " * tens2str(T) * "has multiple index types, none of which are 'spacetime'. Indicate which indices should be dualised."))
    end
    dimension = tinds.range
    #add logic to raise/lower indices... pop to canonical...
    filtered = [x for x in indices(T) if x[1].index_set == tinds] #only contract with the spacetime indices (incase we have flavour etc)
    p = min(length(filtered), dimension) #form degree (ignoring eg flavour inds)
    new_deg = dimension - p
    if new_deg < 0 #can't have negative degree forms
        return  zero()
    end
    # println(T)
    #start the epsilon with the free inds
    epsindlist = Vector{Pair{IndexSuperType,Alakazam.IndexPosition}}([tinds.indices[i] => op_pos(filtered[i][2]) for i in 1:p] ∪ [tinds.indices[i] => op_pos(tinds.default_position) for i in (p+1):dimension])#will work after popping the tensor
    # println(typeof(epsindlist))
    # println(T)
    eps = EpsilonTensor(epsindlist) #create an n dim epsilon
    Tdual = eps × T #append the epsilon first -> this order should ensure we keep (atleast) the same free inds...
    forminds = [x for x in indices(terms(Tdual)[2]) if x[1].index_set == tinds]
    epsinds = indices(terms(Tdual)[1])
    contraction = 1 // factorial(p) * prod(KroneckerDelta([ti[1] => op_pos(ti[2]), ei[1] => op_pos(ei[2])]) for (ti, ei) in zip(forminds, epsinds[dimension-p+1:end]); init=1)#form the contraction with epsilon. we go backwards to preserve
    Tdual *= contraction
    Tdual = eliminate_kronecker(Tdual)
    return Tdual
end

"""At least `n` indices of `I` outside `used`, generating fresh ones when the declared pool runs short."""
function _spare_indices(I::IndexSet, used, n::Integer)
    out = IndexSuperType[x for x in I.indices if x ∉ used]
    while length(out) < n
        nid = generate_indices(I)
        nid === nothing && throw(OutOfIndices(I))
        push!(out, nid)
    end
    out
end

"""Copy of `T` with its `tinds` indices replaced by `newinds` in order, every other index kept."""
function _set_form_indices(T::TensorSuperType, tinds::IndexSet, newinds)
    out = Vector{Pair{IndexSuperType,IndexPosition}}()
    k = 0
    for (ind, pos) in indices(T)
        if ind.index_set == tinds
            k += 1
            push!(out, newinds[k])
        else
            #a spinor or flavour index rides along untouched
            push!(out, ind => pos)
        end
    end
    Tc = copy(T)
    set_indices!(Tc, out)
    Tc
end

"""Take the wedge exterior product between two forms"""
wedge(T1::TensorSuperType, T2::TensorSuperType,index_type::Union{IndexSet,Nothing}=nothing) = ∧(T1::TensorSuperType, T2::TensorSuperType,index_type)

"""Take the wedge exterior product between two forms"""
function ∧(T1::TensorSuperType, T2::TensorSuperType, index_type::Union{IndexSet,Nothing}=nothing)

    if !isForm(T1,index_type) || !isForm(T2,index_type)
        throw(InvalidDifferentialForm("Both tensors must be antisymmetric forms."))
    end

    #similar to star
    if index_type===nothing
        tinds = nothing
        for (ind, _) in indices(T1)
            tinds = ind.index_set
            break
        end
    else
        tinds=index_type
    end

    # extra spacetime
    inds1 = [x for x in indices(T1) if x[1].index_set == tinds]
    inds2 = [x for x in indices(T2) if x[1].index_set == tinds]
    # combinds1=vcat(inds1,inds2)
    T = T1 × T2
    combinds = [x for x in indices(T) if x[1].index_set == tinds]
    # println(combinds)
    # println([x for x in tinds.indices])
    # println([y for y in combinds])
    # println([y[1] for y in combinds])

    unused = _spare_indices(tinds, [y[1] for y in combinds], tinds.range)
    p = length(inds1)
    q = length(inds2)
    r = p + q
    dim = tinds.range
    poss = combinds[1][2]#find the pos of free
    rem = dim - r

    if r > dim
        return  zero()
    end

    #fill out the outermost eps
    # println(combinds)
    epsindlist = Vector{Pair{IndexSuperType,Alakazam.IndexPosition}}(
        vcat(
            [combinds[i][1] => poss for i in 1:r],
            [unused[i] => poss for i in (r+1):dim]
        )
    )
    # println(epsindlist)
    eps1 = EpsilonTensor(epsindlist) #the leadings inds will not change

    unused = _spare_indices(tinds, [y[1] for y in vcat(combinds, indices(eps1))], r)
    epsindlist2 = Vector{Pair{IndexSuperType,Alakazam.IndexPosition}}(
        vcat(
            [unused[i] => op_pos(poss) for i in 1:r],
            [indices(eps1)[i][1] => op_pos(poss) for i in (r+1):dim] #contract with end of other eps
        )
    )
    eps2 = EpsilonTensor(epsindlist2)
    #only the form indices are renamed; a spinor or flavour index the factor
    #carries stays where it is, else a gravitino bilinear comes back bosonic
    T1c = _set_form_indices(T1, tinds, [epsindlist2[i][1] => poss for i in 1:p])
    T2c = _set_form_indices(T2, tinds, [epsindlist2[i][1] => poss for i in p+1:r])
 
    prefact = 1 // (factorial(p) * factorial(q) * factorial(dim - r))
    return prefact * eps1 * eps2 * T1c * T2c
end

"""
The index set a form's slots belong to, chosen from the slots it actually has.
"""
_form_index_set(T::TensorSuperType, x::Coordinate, index_type::Union{IndexSet,Nothing}) =
    _form_index_set(indices(T), x, index_type)

function _form_index_set(pairs, x::Coordinate, index_type::Union{IndexSet,Nothing})
    index_type === nothing || return index_type
    sets = IndexSet[]
    for (ind, _) in pairs
        ind.index_set in sets || push!(sets, ind.index_set)
    end
    #a zero-form carries no slot to read the set from, so the coordinate says
    isempty(sets) && return first(x.index_set)
    cands = IndexSet[S for S in sets if S in x.index_set]
    #if the coordinate declares none of them, fall back to every set rather than refusing outright
    isempty(cands) && (cands = sets)
    for S in cands
        S.name == "spacetime" && return S
    end
    first(cands)
end

"""Exterior derivative of a form."""
𝖽(T::TensorSuperType, x::Coordinate, args...; kwargs...) = exterior_derivative(T, x, args...; kwargs...)

"""
Exterior derivative of a p-form, as the antisymmetrised partial derivative

    (dω)_{a₀a₁…a_p} = Σᵢ (-1)ⁱ ∂_{aᵢ} ω_{a₀…âᵢ…a_p}

The new slot is the first. Its index is generated unless one is given, and the
derivative is named `dname` so several may coexist.
"""
function exterior_derivative(T::TensorSuperType, x::Coordinate,
                             newidx::Union{IndexSuperType,Nothing}=nothing,
                             index_type::Union{IndexSet,Nothing}=nothing;
                             dname::String="∂")
    tinds = _form_index_set(T, x, index_type)
    isForm(T, tinds) || throw(InvalidDifferentialForm(
        "The tensor " * tens2str(T) * " should be totally antisymmetric."))

    formpairs = [pr for pr in indices(T) if pr[1].index_set == tinds]
    p = length(formpairs)
    pos = p == 0 ? tinds.default_position : formpairs[1][2]
    p + 1 > tinds.range && return zero()   #no form of degree above the dimension

    if newidx === nothing
        newidx = _spare_indices(tinds, [pr[1] for pr in indices(T)], 1)[1]
    end

    slots = vcat(Pair{IndexSuperType,IndexPosition}[newidx => pos], formpairs)
    out = zero()
    for i in 1:(p+1)
        #drop the slot the derivative carries, the rest stay on the form in order
        remaining = [slots[k] for k in 1:(p+1) if k != i]
        Tc = _set_form_indices(T, tinds, remaining)
        out += (-1)^(i - 1) * PartialDerivative(dname, x, slots[i])(Tc)
    end
    out
end

"""Copy of `T` with its `tinds` slots relabelled from `from` to `to`, in order."""
function _relabel_slots(T::Tensorial, from, to)
    ref = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
    for (f, t) in zip(from, to)
        ref[f] = t
    end
    Tc = copy(T)
    set_indices_deep!(Tc, ref)
    Tc
end

"""Copy of an expression with its slots relabelled, every term at once."""
function _relabel_expression(T::TensorExpression, ref::Dict)
    out = zero()
    for (coef, term) in summands(T)
        new_term = 1
        for t in terms(term)
            tc = copy(t)
            set_indices_deep!(tc, ref)
            new_term *= tc
        end
        out += coef * new_term
    end
    out
end

"""
Exterior derivative of a form carried by an expression, so that d∘d composes.

The slots are the expression's free indices, fixed once and used for every
term. Reading them off each term separately would order them differently in
∂_a⟦A_b⟧ and ∂_b⟦A_a⟧, and the antisymmetrisation would not close.
"""
function exterior_derivative(T::TensorExpression, x::Coordinate,
                             newidx::Union{IndexSuperType,Nothing}=nothing,
                             index_type::Union{IndexSet,Nothing}=nothing;
                             dname::String="∂")
    #free_indices is a multiset, so take its keys
    free = Pair{IndexSuperType,IndexPosition}[k for (k, _) in collect(free_indices(T))]
    #same choice as the tensor method: not free[1]'s set, which depends both on declaration order
    #and on how free_indices happens to enumerate
    tinds = _form_index_set(free, x, index_type)
    formpairs = sort([pr for pr in free if pr[1].index_set == tinds], by=pr -> name(pr[1]))
    p = length(formpairs)
    pos = p == 0 ? tinds.default_position : formpairs[1][2]
    p + 1 > tinds.range && return zero()

    if newidx === nothing
        newidx = _spare_indices(tinds, [pr[1] for pr in free], 1)[1]
    end

    slots = vcat(Pair{IndexSuperType,IndexPosition}[newidx => pos], formpairs)
    out = zero()
    for i in 1:(p+1)
        remaining = [slots[k] for k in 1:(p+1) if k != i]
        ref = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
        for (f, t) in zip(formpairs, remaining)
            ref[f] = t
        end
        out += (-1)^(i - 1) * PartialDerivative(dname, x, slots[i])(_relabel_expression(T, ref))
    end
    out
end

exterior_derivative(T::TensorTerm, x::Coordinate, args...; kwargs...) =
    exterior_derivative(TensorExpression(T), x, args...; kwargs...)

"""
Interior product (contraction) of a vector with a differential form.
"""
ι(T1::TensorSuperType, T2::TensorSuperType) = ⨼(T1, T2)
contraction(T1::TensorSuperType, T2::TensorSuperType) = ⨼(T1, T2)

"""
Interior product (contraction) of a vector with a differential form.
"""
function ⨼(T1::TensorSuperType, T2::TensorSuperType)
    if length(indices(T1)) != 1 #will need to update this, really
        throw(ArgumentError("Vector must have one index of the contraction type."))
    end
    if !isForm(T2)
        throw(InvalidDifferentialForm("The tensor " * tens2str(T2) * " should be totally antisymmetric."))
    end

    T2cp = copy_owned(T2)
    set_index!(T2cp,1,(indices(T1)[1][1] => op_pos(indices(T1)[1][2])))

    return T1 * T2cp
end


