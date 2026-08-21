#########
## functions on tensors that are diff. forms
#########

"""The positions of `T`'s slots belonging to `iset`."""
_form_slots(T::TensorSuperType, iset::IndexSet) =
    [i for (i, ind) in enumerate(indices(T)) if ind[1].index_set == iset]

"""
Whether `T` is a differential form in the given index set, that is, whether the
slots belonging to that set are totally antisymmetric among themselves.
"""
function isForm(T::TensorSuperType, index_type::Union{IndexSet,Nothing}=nothing)
    isempty(indices(T)) && return true            #a 0-form
    iset = _form_index_set(indices(T), nothing, index_type)
    iset === nothing && return length(indices(T)) <= 1
    slots = _form_slots(T, iset)
    length(slots) <= 1 && return true   #a 0- or 1-form is trivially antisymmetric
    for tab in tableaux(T)
        #a column tableau totally antisymmetrises the slots it covers
        length(tab.row_lengths) == length(tab.index_assignment) || continue
        issubset(slots, tab.index_assignment) && return true
    end
    return false
end

"""Take the Hodge star of a form.
This operation will automatically be taken with respect to spacetime indices, if they are defined on the tensor.
"""
star(T::Tensorial, index_type::Union{IndexSet,Nothing}=nothing) = ⋆(T, index_type)

"""
Hodge star of a summand that is a product rather than a single tensor.

"""
function _star_term(T, tinds::IndexSet)
    #Structural order here, not `_form_free`'s canonical one. The star keeps the names it
    #was given by reusing the argument's slots first, and a second star recovers them by
    #reading that leading position back 
    formpairs = Pair{IndexSuperType,IndexPosition}[pr for pr in free_indices_ordered(T)
                                                   if pr[1].index_set == tinds]
    dimension = tinds.range
    p = length(formpairs)
    dimension - p < 0 && return zero()
    if length(unique(pos for (_, pos) in formpairs)) > 1
        throw(InvalidDifferentialForm("The " * tinds.name * " slots of this summand sit at " *
            "mixed levels. Raise or lower them to a common level before taking the Hodge star."))
    end
    epspos = isempty(formpairs) ? op_pos(tinds.default_position) : op_pos(formpairs[1][2])
    nfree = dimension - p
    nreuse = min(nfree, p)
    avoid = OrderedSet{IndexSuperType}()
    for i in _all_indices(T)
        push!(avoid, i)
    end
    freeslots = IndexSuperType[formpairs[k][1] for k in 1:nreuse]
    if nfree > nreuse
        for idx in _spare_indices(tinds, avoid, nfree - nreuse)
            push!(freeslots, idx)
            push!(avoid, idx)
        end
    end
    #The argument's slots move aside to fresh names, which then become the epsilon's contracted ones
    temps = _spare_indices(tinds, avoid, p)
    ref = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
    for (k, pr) in enumerate(formpairs)
        ref[pr[1] => upper] = temps[k] => upper
        ref[pr[1] => lower] = temps[k] => lower
    end
    Tr = copy(T)
    isempty(ref) || set_indices_deep!(Tr, ref)

    epsindlist = Vector{Pair{IndexSuperType,IndexPosition}}(
        vcat([f => epspos for f in freeslots], [t => epspos for t in temps]))
    return (1 // factorial(p)) * EpsilonTensor(epsindlist) * Tr
end

"""
Hodge star of an expression, summand by summand.
"""
function ⋆(T::TensorExpression, index_type::Union{IndexSet,Nothing}=nothing)
    out = zero()
    for (coef, term) in summands(T)
        f = _lone_factor(term)
        if f isa TensorSuperType
            out += coef * ⋆(f, index_type)
        else
            tinds = index_type === nothing ?
                _form_index_set(free_indices_ordered(f), nothing, nothing) : index_type
            out += coef * _star_term(f, tinds)
        end
    end
    out
end

function ⋆(T::TensorTerm, index_type::Union{IndexSet,Nothing}=nothing)
    f = _lone_factor(T)
    f isa TensorSuperType && return ⋆(f, index_type)
    tinds = index_type === nothing ?
        _form_index_set(free_indices_ordered(f), nothing, nothing) : index_type
    _star_term(f, tinds)
end

"""Take the Hodge star of a form.
This operation will automatically be taken with respect to spacetime indices, if they are defined on the tensor.
"""
function ⋆(T::TensorSuperType, index_type::Union{IndexSet,Nothing}=nothing)
    if index_type === nothing && length(indices(T)) == 0
        throw(InvalidDifferentialForm("The tensor " * tens2str(T) * " has no indices, please specify the spacetime indices for the Hodge star."))
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
    #the form check needs the set the slots live in, so it runs once that is known
    if !isForm(T, tinds)
        throw(InvalidDifferentialForm("The tensor " * tens2str(T) * " should be totally antisymmetric in its " * tinds.name * " indices."))
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
    if length(unique(pos for (_, pos) in filtered)) > 1
        throw(InvalidDifferentialForm("The " * tinds.name * " slots of " * tens2str(T) *
            " sit at mixed levels. Raise or lower them to a common level before taking the Hodge star."))
    end
    return _star_term(T, tinds)
end

"""
`n` indices of `I` outside `used`, generating fresh ones when the declared pool runs short.

Pass an empty `used` for simply the first `n` indices of the set.
"""
function _spare_indices(I::IndexSet, used, n::Integer)
    pool = IndexSuperType[x for x in I.indices if x ∉ used]
    [safe_pop!(pool, I) for _ in 1:n]
end

"""
The slots a form-like object offers to the form machinery: its free indices in
that set.

For a plain tensor these are just its slots, so nothing changes. For a product
such as `⋆F = ϵ^{ab}{}_{ef}F^{ef}` the contracted pairs are not slots of the form,
and only the free `a`,`b` are, which is what lets `⋆` and `∧` consume their own
output.
"""
function _form_free(T, tinds::IndexSet)
    out = Pair{IndexSuperType,IndexPosition}[]
    for pr in free_indices_ordered(T)
        pr[1].index_set == tinds && push!(out, pr)
    end
    _slot_sort(out)
end

"""Every index in play, dummies included, as a set."""
_all_indices(T) = collect_indices!(OrderedSet{IndexSuperType}(), T)

"""
Rename `idxs` throughout `T` to names outside `avoid`, both levels together so a
dummy pair stays a pair.
"""
function _rename_away(T, tinds::IndexSet, idxs, avoid)
    isempty(idxs) && return T
    ref = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
    for idx in idxs
        fresh = _spare_indices(tinds, avoid, 1)[1]
        push!(avoid, fresh)
        ref[idx => upper] = fresh => upper
        ref[idx => lower] = fresh => lower
    end
    Tc = copy(T)
    set_indices_deep!(Tc, ref)
    Tc
end

"""
Copy of `T` with its `tinds` slots brought to `pos`, a metric contracted onto each
slot that sits at the other level. Metrics are appended so the new slots are read
back in the order the old ones were.
"""
function _to_position(T, tinds::IndexSet, pos::IndexPosition, avoid)
    out = T
    for pr in _form_free(T, tinds)
        pr[2] == pos && continue
        tmp = _spare_indices(tinds, avoid, 1)[1]
        push!(avoid, tmp)
        ref = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}(
            (pr[1] => pr[2]) => (tmp => pr[2]))
        Tc = copy(out)
        set_indices_deep!(Tc, ref)
        out = Tc * Metric(pr[1] => pos, tmp => pos)
    end
    out
end

"""Copy of a product with its `tinds` free slots relabelled from `from` to `to`, in order."""
function _set_form_indices(T::TensorTerm, tinds::IndexSet, newinds)
    ref = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
    for (old, new) in zip(_form_free(T, tinds), newinds)
        ref[old] = new
    end
    Tc = copy(T)
    set_indices_deep!(Tc, ref)
    Tc
end

"""
Copy of `T` with its `tinds` indices replaced by `newinds` in order, every other index kept.

`newinds[k]` goes to the k-th slot in canonical order, matching `_form_free`, so a tensor
whose indices were written out of that order takes the new names by rank rather than by
the position they sit in.
"""
function _set_form_indices(T::TensorSuperType, tinds::IndexSet, newinds)
    formpairs = Pair{IndexSuperType,IndexPosition}[pr for pr in indices(T)
                                                   if pr[1].index_set == tinds]
    rank = Dict(pr => k for (k, pr) in enumerate(_slot_sort(formpairs)))
    out = Vector{Pair{IndexSuperType,IndexPosition}}()
    for (ind, pos) in indices(T)
        if ind.index_set == tinds
            push!(out, newinds[rank[ind => pos]])
        else
            #a spinor or flavour index rides along untouched
            push!(out, ind => pos)
        end
    end
    Tc = copy(T)
    if Tc isa OperatorSuperType
        ref = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()
        for (old, new) in zip(indices(T), out)
            ref[old] = new
        end
        set_indices_deep!(Tc, ref)
    else
        set_indices!(Tc, out)
    end
    Tc
end

"""Take the wedge exterior product between two forms"""
wedge(T1::TensorSuperType, T2::TensorSuperType,index_type::Union{IndexSet,Nothing}=nothing; expand::Bool=false) = ∧(T1, T2, index_type; expand=expand)

"""Take the wedge exterior product between two forms"""
function ∧(T1::TensorSuperType, T2::TensorSuperType, index_type::Union{IndexSet,Nothing}=nothing; expand::Bool=false)
    if !isForm(T1,index_type) || !isForm(T2,index_type)
        throw(InvalidDifferentialForm("Both tensors must be antisymmetric forms."))
    end
    _wedge_core(T1, T2, index_type; expand=expand)
end

"""
The canonical slot order for a form: its slots sorted by index name.
    d(α∧β) = dα∧β + (-1)^{deg α} α∧dβ
"""
_slot_sort(prs) = sort(collect(prs), by = pr -> name(pr[1]))

"""
The wedge itself, without the form check.
"""
function _wedge_core(T1::Tensorial, T2::Tensorial, index_type::Union{IndexSet,Nothing}=nothing; expand::Bool=false)
    if index_type===nothing
        tinds = _form_index_set(vcat(free_indices_ordered(T1), free_indices_ordered(T2)),
                                nothing, nothing)
    else
        tinds=index_type
    end
    #Neither side has a slot at all, so both are scalars and the wedge is their product.
    tinds === nothing && return T1 * T2

    #the slots are the free indices, so a factor may itself be a contracted product
    inds1 = _form_free(T1, tinds)
    inds2 = _form_free(T2, tinds)
    avoid = OrderedSet{IndexSuperType}()
    for i in union(_all_indices(T1), _all_indices(T2))
        push!(avoid, i)
    end

    t1_free = OrderedSet{IndexSuperType}(pr[1] for pr in inds1)
    t2_free = OrderedSet{IndexSuperType}(pr[1] for pr in inds2)
    squatters = IndexSuperType[i for i in _all_indices(T1)
                               if i in t2_free && !(i in t1_free)]
    if !isempty(squatters)
        T1 = _rename_away(T1, tinds, squatters, avoid)
        inds1 = _form_free(T1, tinds)
    end
    first_all = _all_indices(T1)
    clash = IndexSuperType[i for i in _all_indices(T2) if i in first_all]
    if !isempty(clash)
        T2 = _rename_away(T2, tinds, clash, avoid)
        inds2 = _form_free(T2, tinds)
    end
    #A 0-form is a scalar, and wedging with one is plain multiplication: f∧ω = fω, with no antisymmetrisation to do
    (isempty(inds1) || isempty(inds2)) && return T1 * T2

    poss = isempty(inds1) ? (isempty(inds2) ? tinds.default_position : inds2[1][2]) : inds1[1][2]

    T1 = _to_position(T1, tinds, poss, avoid); inds1 = _form_free(T1, tinds)
    T2 = _to_position(T2, tinds, poss, avoid); inds2 = _form_free(T2, tinds)

    combinds = _slot_sort(vcat(inds1, inds2))

    #the fresh names must dodge the dummies inside either factor, not just the slots
    unused = _spare_indices(tinds, union(_all_indices(T1), _all_indices(T2)), tinds.range)
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

    unused = _spare_indices(tinds,
        union(_all_indices(T1), _all_indices(T2),
              OrderedSet{IndexSuperType}(y[1] for y in indices(eps1))), r)
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
 
    #The wedge is algebraic on the cotangent space and so cannot depend on the
    #metric, but contracting two epsilons carries a factor of det(η), which is
    #(-1)^t for t timelike directions. Dividing it out keeps the result equal to
    #((p+q)!/p!q!) α_[..β..] in every signature rather than flipping sign between
    #Lorentzian and Euclidean.
    detsign = (-1)^(tinds.signature[1])
    prefact = detsign // (factorial(p) * factorial(q) * factorial(dim - r))
    
    expand || return prefact * eps1 * eps2 * T1c * T2c
    return eliminate_kronecker(prefact * eliminate_epsilon(eps1 * eps2) * T1c * T2c)
end

"""
The (coefficient, factor) pairs a wedge argument splits into.

A single form is its own only part; a sum contributes one part per summand. A
summand holding more than one factor is a tensor product, not a form, so it has
no wedge.
"""
_wedge_parts(T::TensorSuperType) = Tuple{Any,Tensorial}[(1, T)]
_wedge_parts(T::TensorTerm) = Tuple{Any,Tensorial}[(1, _lone_factor(T))]
_wedge_parts(T::TensorExpression) =
    Tuple{Any,Tensorial}[(coef, _lone_factor(term)) for (coef, term) in summands(T)]

"""
A one-factor product is just that factor, so the tensor path with its form check
still applies to it. A genuine product is handed on whole, its slots being its
free indices.
"""
function _lone_factor(T::TensorTerm)
    fs = terms(T)
    length(fs) == 1 ? fs[1] : T
end
_lone_factor(T) = T

"""
Wedge product of expressions, by bilinearity.
"""
function ∧(T1::Tensorial, T2::Tensorial, index_type::Union{IndexSet,Nothing}=nothing; expand::Bool=false)
    out = zero()
    for (c1, f1) in _wedge_parts(T1)
        for (c2, f2) in _wedge_parts(T2)
            out += c1 * c2 * _wedge_core(f1, f2, index_type; expand=expand)
        end
    end
    out
end

"""
The index set a form's slots belong to, chosen from the slots it actually has.
"""
_form_index_set(T::TensorSuperType, x::Coordinate, index_type::Union{IndexSet,Nothing}) =
    _form_index_set(indices(T), x, index_type)

function _form_index_set(pairs, x::Union{Coordinate,Nothing}, index_type::Union{IndexSet,Nothing})
    index_type === nothing || return index_type
    sets = IndexSet[]
    for (ind, _) in pairs
        ind.index_set in sets || push!(sets, ind.index_set)
    end
    #a zero-form carries no slot to read the set from, so the coordinate says
    isempty(sets) && return x === nothing ? nothing : first(x.index_set)
    #isForm has no coordinate to narrow the candidates with, so it keeps them all
    cands = x === nothing ? sets : IndexSet[S for S in sets if S in x.index_set]
    #if the coordinate declares none of them, fall back to every set rather than refusing outright
    isempty(cands) && (cands = sets)
    for S in cands
        S.name == "spacetime" && return S
    end
    first(cands)
end

"""Exterior derivative of a form."""
𝖽(T::Tensorial, x::Coordinate, args...; kwargs...) = exterior_derivative(T, x, args...; kwargs...)

"""
The coordinate to differentiate along when the caller did not name one.
"""
function _infer_coordinate(T::Tensorial, tinds::IndexSet)
    own = Coordinate[]
    for c in get_operand_function_of(TensorExpression(T))
        tinds in c.index_set && !(c in own) && push!(own, c)
    end
    length(own) == 1 && return own[1]
    length(own) > 1 && throw(InvalidDifferentialForm(
        "This form is a function of several coordinates over '" * tinds.name * "' (" *
        join(sort([c.name for c in own]), ", ") * "), so the one to differentiate along " *
        "is ambiguous. Pass it explicitly."))
    reg = Coordinate[c for c in values(COORD_REGISTRY) if tinds in c.index_set]
    isempty(reg) && throw(InvalidDifferentialForm(
        "No coordinate is declared over '" * tinds.name * "', so there is nothing to " *
        "differentiate along. Declare one, or pass it explicitly."))
    sort!(reg, by = c -> (length(c.index_set), c.name))
    return reg[1]
end

"""
    exterior_derivative(T)
    𝖽(T)

Exterior derivative with the coordinate left implicit
"""
function exterior_derivative(T::Tensorial,
                             newidx::Union{IndexSuperType,Nothing}=nothing,
                             index_type::Union{IndexSet,Nothing}=nothing; kwargs...)
    tinds = _form_index_set(free_indices_ordered(T), nothing, index_type)
    tinds === nothing && throw(InvalidDifferentialForm(
        "This form carries no slots to read an index set from, so its coordinate cannot " *
        "be inferred. Pass one explicitly."))
    exterior_derivative(T, _infer_coordinate(T, tinds), newidx, index_type; kwargs...)
end

𝖽(T::Tensorial, newidx::Union{IndexSuperType,Nothing}=nothing,
  index_type::Union{IndexSet,Nothing}=nothing; kwargs...) =
    exterior_derivative(T, newidx, index_type; kwargs...)

"""
Exterior derivative of a p-form, as the antisymmetrised partial derivative

    (dω)_{a₀a₁…a_p} = Σᵢ (-1)ⁱ ∂_{aᵢ} ω_{a₀…âᵢ…a_p}

Its index is generated unless one is given, and
the derivative is named `dname` so several may coexist.
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
        #every index already in play, dummies included, not just the free slots
        newidx = _spare_indices(tinds, collect_indices!(OrderedSet{IndexSuperType}(), T), 1)[1]
    end

    #canonical order
    slots = _slot_sort(vcat(Pair{IndexSuperType,IndexPosition}[newidx => pos], formpairs))
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

        newidx = _spare_indices(tinds, collect_indices!(OrderedSet{IndexSuperType}(), T), 1)[1]
    end

    slots = _slot_sort(vcat(Pair{IndexSuperType,IndexPosition}[newidx => pos], formpairs))
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


