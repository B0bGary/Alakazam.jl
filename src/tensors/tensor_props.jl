######################
## Functions for extracting basic data from tensors.
## Functions are grouped into different files for readability.
######################
"""Vector of tensors inside TensorTerm"""
@inline terms(T::TensorTerm) = T.terms
"""
Get just the TensorTerm parts
"""
@inline terms(T::Vector{Pair{SymorNum,TensorTerm}}) = [t[2] for t in T]
"""
Get just the TensorTerm part
"""
@inline terms(T::Pair{SymorNum,TensorTerm}) = [T[2]]
"""Get just the tensor terms that compose an expression."""
@inline terms(T::TensorExpression) = terms(summands(T))
@inline coefficients(T::Pair{SymorNum,TensorTerm}) = [T[1]]
@inline coefficients(T::Vector{Pair{SymorNum,TensorTerm}}) = [t[1] for t in T]
@inline coef(t::Pair{SymorNum,TensorTerm}) = t[1]
@inline term(t::Pair{SymorNum,TensorTerm}) = t[2]

function set_terms!(T::TensorTerm, newterms::Vector{TensorSuperType})
    T.terms = newterms
    return T
end
@inline summands(T::TensorExpression) = T.summands
@inline function isconst(T::TensorExpression)
    length(terms(T)) == 1 || return false
    tt = T.summands[1][2]              # the single TensorTerm
    n = length(terms(tt))
    return n == 0 || (n == 1 && isempty(name(terms(tt)[1])))
end
Base.length(T::TensorExpression) = length(summands(T))
Base.length(T::TensorTerm) = length(terms(T))

"""First TensorTerm in the expression."""
function first_summand(T::TensorExpression)
    if isempty(summands(T))
        return nothing
    end
    return summands(T)[1][2]
end

"""First numerical/symbolic coefficient in the expression."""
function first_coef(T::TensorExpression)
    return summands(T)[1][1]
end


first_summands(T::TensorExpression) = summands(T)[1][2]



"""Get the data associated to a tensor object.

Plain tensors return their shared [`TensorData`](@ref), operators return their `OperatorData`.
"""
data(T::TensorSuperType) = T.data
data(T::OperatorSuperType) = T.data

"""Re-bundle a tensor's shared props and its indices into a `TensorData`.
"""
tensor_data(T::TensorSuperType) = data(T)
tensor_data(T::OperatorSuperType) = op_data(T)


"""Return the name String of a tensor."""
name(T::TensorDataSuperType)= T.name
"""Return the name String of a tensor."""
name(T::TensorSuperType)=name(data(T))

"""Rename a tensor. Data is immutable, so this returns a new tensor with the copied data but a new name.
"""
function set_name!(T::TensorSuperType, newname::String)
    T.data = _with_name(T.data, newname)
    return newname
end
function set_name!(T::OperatorSuperType, newname::String)
    return T.data.name = newname
end

indices(::SymorNum) = Vector{<:Pair{<:IndexSuperType,IndexPosition}}()
indices(T::OperatorData) = T.indices
indices(T::TensorSuperType) = T.indices
indices(T::OperatorSuperType) = indices(data(T))
indices(T::TensorTerm) = T.indices
"""
Gets the vector of index/position pairs attached to an expression.
This is based off the first term in the expression, so note that other terms may have different dummy index names.
"""
function indices(T::TensorExpression)
    if isZeroTens(T)
        return Vector{Pair{IndexTypes,IndexPosition}}()
    end
    return indices(first_summand(T))
end

function set_indices!(T::TensorSuperType, inds::Vector{<:Pair{<:IndexSuperType,IndexPosition}})# where {M<:IndexSuperType}
    T.indices = inds
    return T
end
"""Only updates the external indices, not those of the underlying tensors. Use set_indices_deep! otherwise."""
function set_indices!(T::TensorTerm, inds::Vector{<:Pair{<:IndexSuperType,IndexPosition}}) #where {M<:IndexSuperType}
    T.indices = inds
    return T
end


set_indices_deep!(T::TensorSuperType, inds::Vector{<:Pair{<:IndexSuperType,IndexPosition}}) = set_indices!(T, inds)#where {M<:IndexSuperType}
function set_indices_deep!(T::TensorSuperType, ref::Dict{<:Pair{<:IndexSuperType,IndexPosition},<:Pair{<:IndexSuperType,IndexPosition}}) #where {I<:IndexSuperType,M<:IndexSuperType}
    # @show indices(T)
    # @show [ref[k] for k in indices(T)]
    # @show T
    set_indices!(T, Pair{IndexTypes,IndexPosition}[get(ref, k, k) for k in indices(T)])
    # @show T
    update_indices!(T)
end


"""Overwrite the indices in a tensorterm according to ref dict"""
function set_indices_deep!(T::TensorTerm, ref::Dict{<:Pair{<:IndexSuperType,IndexPosition},<:Pair{<:IndexSuperType,IndexPosition}})# where {M<:IndexSuperType}
    for tterm in terms(T)
        set_indices_deep!(tterm, ref)
    end
    set_indices!(T, Pair{IndexTypes,IndexPosition}[get(ref, k, k) for k in indices(T)])
    update_indices!(T)
end


"""Overwrite the indices in a tensor expression - using a reference set"""
function set_indices_deep!(T::TensorExpression, rep_guide::Dict{<:Pair{<:IndexSuperType,IndexPosition},<:Pair{<:IndexSuperType,IndexPosition}}) #where {M<:IndexSuperType}
    #convert to dict
    for tterm in terms(T)
        set_indices_deep!(tterm, rep_guide)
    end
    update_indices!(T)
end



"""Sets the index on a vector in place. Assumes user has supplied an unshared index vector."""
function set_index!(T::TensorSuperType, vin::Integer, new_ind::Pair{I,IndexPosition}) where {I<:IndexSuperType}
    indices(T)[vin] = new_ind
end


function function_of(T::TensorDataSuperType)
    return T.function_of
end

function function_of(T::TensorSuperType)
    return function_of(data(T))
end
function function_of(T::TensorTerm)
    fnc = coordfree()
    for t in terms(T)
        f = function_of(t)
        isempty(f) && continue
        fnc === _EMPTY_FUNCTION_OF && (fnc = Coordinate[])
        for c in f
            _insert_coord!(fnc, c)
        end
    end
    return fnc
end

function function_of(T::TensorExpression)
    fnc = coordfree()
    for s in summands(T)
        f = function_of(term(s))
        isempty(f) && continue
        fnc === _EMPTY_FUNCTION_OF && (fnc = Coordinate[])
        for c in f
            _insert_coord!(fnc, c)
        end
    end
    return fnc
end


function weights(T::TensorDataSuperType)
    return T.weights
end

function weights(T::TensorSuperType)
    return weights(data(T))
end

function weights(T::TensorTerm)
    w = Dict{String,Number}()
    for t in terms(T)
        for (k, v) in weights(t)
            w[k] = get(w, k, 0) + v
        end
    end
    return w
end


function weights(T::TensorExpression)
    w = Dict{String,Number}()
    for t in summands(T)
        for (k, v) in weights(term(t))
            w[k] = get(w, k, 0) + v
        end
    end
    return w
end


function tableaux(T::TensorDataSuperType)
    return T.tableaux
end
function tableaux(T::TensorSuperType)
    return tableaux(data(T))
end


function labels(T::TensorDataSuperType)
    return T.labels
end
function labels(T::TensorSuperType)
    return labels(data(T))
end


"""
Check if a TensorExpression is a single TensorTerm. Returns true or false.
"""
function singleterm(TEx::TensorExpression)
    if length(summands(TEx)) == 1 && !coeff_equal(first_coef(TEx), 0)
        return true
    end
    return false
end
"""
Check if a TensorTerm is a single Tensor. Returns true or false.
"""
function singleterm(TT::TensorTerm)
    if length(terms(TT)) == 1
        return true
    end
    return false
end



################################
# Get indices logic
################################
"""
Gets the set of free index/↑↓position pairs attached to an expression. 
"""
function free_indices(T::Tensorial)
    counter(free_indices_ordered(T))
end

# function free_indices_ordered(T::Tensorial)
#     dp = dummies_positional(T)
#     filter(i -> i ∉ dp, indices(T))
# end
function free_indices(T::Number)
    return Vector{Pair{IndexTypes,IndexPosition}}()
end

"""
Gets the set of free indices, without upper/lower positions, attached to an expression.
"""
# function free_indices_nopos(T::Tensorial)
#     out = IndexTypes[]
#     for d in free_indices(T)
#         push!(out, first(first(d)))
#     end
#     return out
# end
function free_indices_nopos(T::Tensorial, dp::Vector)#shortcut version
    out = IndexTypes[]
    for (idx, pos) in indices(T)
        (idx => pos) in dp && continue   # skip dummies
        f = idx
        f in out || push!(out, f)
    end
    return out
end

"""
Gets the set of IndexSets associated to a tensorial object.
"""
function get_index_sets(T::Tensorial)
    out = IndexSet[]
    for ind in indices(T)
        s = first(ind).index_set
        s in out || push!(out, s)
    end
    return out
end

"""
Gets the set of repeated summation indices in a Tensor.
"""
dummies(T::Tensorial) = _dummies(T, Val(false))
"""
Gets the set of repeated summation index pairs in a Tensor.
"""
dummies_positional(T::Tensorial) = _dummies(T, Val(true))

"""
Gets the set of all dummies across terms in a TensorExpression.
"""
function _dummies(T::TensorExpression, ::Val{positional}) where {positional}
    out = positional ?
          Vector{Pair{IndexTypes,IndexPosition}}() :
          Vector{IndexTypes}()
    v = Val(positional)
    for summand in summands(T)
        for d in _dummies(summand[2], v)
            d in out || push!(out, d)      # union, vector-style
        end
    end
    return out
end

# function _dummies(T::Union{TensorSuperType,TensorTerm}, v::Val{positional}) where {positional}
#     dummies, seen = _dummies_raw(T, v)
#     isempty(dummies) && return dummies
#     if positional
#         sort!(dummies, by=p -> findfirst(q -> first(q) == first(p), seen))
#     else
#         sort!(dummies, by=idx -> findfirst(q -> first(q) == idx, seen))
#     end
#     return dummies
# end
function _dummies_raw(T::Union{TensorSuperType,TensorTerm}, ::Val{positional}) where {positional}
    inds = indices(T)
    n = length(inds)
    seen = Vector{Pair{IndexTypes,IndexPosition}}()
    if n <= 1
        return positional ? (Pair{IndexTypes,IndexPosition}[], seen) : (IndexTypes[], seen)
    end

    # small collections → plain vectors, linear scan. No hash tables.

    sizehint!(seen, n)
    dummies = positional ?
              Vector{Pair{IndexTypes,IndexPosition}}() :
              Vector{IndexTypes}()

    for (idx, pos) in inds
        isLightcone(idx) && continue

        if (idx => pos) in seen        # linear, fast for small n
            println(T)
            println(idx => pos)
            @show indices(T)
            throw(SameLevelDummyError(T))
        elseif (idx => op_pos(pos)) in seen
            if positional
                push!(dummies, idx => op_pos(pos))
                push!(dummies, idx => pos)
            else
                idx in dummies || push!(dummies, idx) 
            end
        else
            push!(seen, idx => pos)
        end
    end

    return dummies, seen    # return both, unsorted
end
"""
Gets pairs of vector position/dummy indices from tensorial objects. (try to make faster)
"""
function dummies_with_pos_index(T::Tensorial)
    dp = dummies_positional(T)           # single scan
    dset = first.(dp)                    # the dummy indices
    return [i => idx for (i, (idx, _)) in enumerate(indices(T)) if idx in dset]
end

"""
Gets vector of pairs, first is the vector index, second is the pair-index/position for dummy indices from tensorial objects. (try to make faster)
"""
function dummies_with_pos_index_pos_ul(T::Tensorial)
    dp = dummies_positional(T)
    dset = first.(dp)
    return [i => pair for (i, pair) in enumerate(indices(T)) if first(pair) in dset]
end

# (upper-mask, lower-mask, non-lightcone-mask)
@inline function _index_masks(T)
    up = UInt64(0);
    dn = UInt64(0);
    nl = UInt64(0)
    inds = indices(T)
    @inbounds for k in eachindex(inds)
        idx, pos = inds[k]
        b = _bit(idx)
        m = -UInt64(pos.position)        # 0x0000… if upper, 0xFFFF… if lower
        dn |= b & m
        up |= b & ~m
        nl |= b & ~(-UInt64(idx.index_set.islightcone))
    end
    return (up, dn, nl)
end
function valid_index_placement(T1, T2)
    up1, dn1, nl1 = _index_masks(T1)
    up2, dn2, nl2 = _index_masks(T2)

    # dummy of one appearing in the other (lightcone excluded)
    d1 = up1 & dn1 & nl1
    d2 = up2 & dn2 & nl2
    if ((d1 & (up2 | dn2)) | (d2 & (up1 | dn1))) != 0
        _ccd_exact(T1, T2) && throw(OverusedIndices(T1, T2))
    end

    # same index at same position in both, and not lightcone
    if (((up1 & up2) | (dn1 & dn2)) & nl1) != 0
        _ccl_exact(T1, T2) || throw(SameLevelDummyError(T1, T2))
    end

    return true
end
"""
One-pass index validity check for a finished factor list.

Equivalent to the pairwise screening `valid_index_placement` does while a term is built up, but
O(k) instead of O(k²): an index may appear at most twice, and if it appears twice it must be once
raised and once lowered. Used where a term is assembled in one go rather than by repeated `*`.
"""
function _check_term_indices(T::TensorTerm)
    inds = indices(T)
    length(inds) < 2 && return T
    seen = Dict{UInt64,Tuple{Int,Bool}}()      # uid => (occurrences, position of first sighting)
    @inbounds for (idx, pos) in inds
        isLightcone(idx) && continue
        cnt, firstpos = get(seen, idx.uid, (0, false))
        cnt += 1
        cnt > 2 && throw(OverusedIndices(T))
        cnt == 2 && firstpos == pos.position && throw(SameLevelDummyError(T))
        seen[idx.uid] = (cnt, cnt == 1 ? pos.position : firstpos)
    end
    return T
end

function _ccl_exact(T1::Tensorial, T2::Tensorial)
    i2 = indices(T2)
    @inbounds for p1 in indices(T1)
        isLightcone(p1[1]) && continue
        for p2 in i2
            p2 === p1 && return false
        end
    end
    return true
end

function _uniq_first(dp)
    out = IndexTypes[]
    for p in dp
        f = first(p)
        f in out || push!(out, f)
    end
    out
end

function free_indices_ordered(T::Tensorial)
    inds = indices(T)
    out = Pair{IndexTypes,IndexPosition}[]
    @inbounds for i in eachindex(inds)
        _is_dummy_at(inds, i) || push!(out, inds[i])
    end
    return out
end

function free_indices_nopos(T::Tensorial)
    inds = indices(T)
    out = IndexTypes[]
    @inbounds for i in eachindex(inds)
        _is_dummy_at(inds, i) && continue
        idx = inds[i][1]
        idx in out || push!(out, idx)
    end
    return out
end
function _dummies(T::Union{TensorSuperType,TensorTerm}, ::Val{positional}) where {positional}
    inds = indices(T)
    n = length(inds)
    out = positional ? Pair{IndexTypes,IndexPosition}[] : IndexTypes[]
    n <= 1 && return out
    @inbounds for i in 1:n
        idx, pos = inds[i]
        isLightcone(idx) && continue
        isfirst = true                          # only emit at first occurrence
        for k in 1:(i-1)
            inds[k][1] === idx && (isfirst=false; break)
        end
        isfirst || continue
        for j in (i+1):n
            inds[j][1] === idx || continue
            inds[j][2] === pos && throw(SameLevelDummyError(T))
            if positional
                push!(out, inds[i]);
                push!(out, inds[j])
            else
                push!(out, idx)
            end
            break
        end
    end
    return out
end

@inline _bit(idx) = UInt64(1) << (idx.uid & 63)

function _ccd_exact(T1::Tensorial, T2::Tensorial)
    i1 = indices(T1);
    i2 = indices(T2)
    @inbounds for a in eachindex(i1)
        _is_dummy_at(i1, a) || continue
        idx = i1[a][1]
        for b in eachindex(i2)
            i2[b][1] === idx && return true
        end
    end
    @inbounds for b in eachindex(i2)
        _is_dummy_at(i2, b) || continue
        idx = i2[b][1]
        for a in eachindex(i1)
            i1[a][1] === idx && return true
        end
    end
    return false
end
"""Is indices(T)[i] a dummy, same index appearing at the opposite position elsewhere in T?"""
@inline function _is_dummy_at(inds, i)
    @inbounds idx, pos = inds[i]
    isLightcone(idx) && return false
    opp = op_pos(pos)
    @inbounds for k in eachindex(inds)
        k == i && continue
        p = inds[k]
        p[1] === idx && p[2] === opp && return true
    end
    return false
end
##############################################################################################################################
##############################################################################################################################

#gradings
"""Returns  zero if the object is bosonic, 1 if fermionic. Assignment is based on the grading of the indices."""
function grading(inds::Vector{Pair{I,IndexPosition}}) where {I<:IndexSuperType}
    grade = 0
    for (ind, pos) in inds
        grade = grade + grading(ind)
    end
    return grade % 2
end

function grading(T::TensorSuperType)
    grade = 0
    for (ind, pos) in indices(T)
        grade = grade + grading(ind)
    end
    return grade % 2
end

"""Hard override for the grading of an object (for example, maybe a Ghost/Grassmann with no indices)"""
grading(T::Fermion) = 1
grading(T::Boson) = 0

function grading(T::TensorTerm)
    grade = 0
    for tens in terms(T)
        grade = grade + grading(tens)
    end
    return grade % 2
end

function grading(T::TensorExpression)
    grade = 0
    for tens in terms(first_summands(T)) # should be same grading for everything in sum since they have matching free indices
        grade = grade + grading(tens)
    end
    return grade % 2
end


########################
## Symmetry Properties
########################

"""
Checks a tensor's Young Tableau to see if it is a totally antisymmetric tensor.
Tensors with 0 or 1 index will be considered antisymmetric, so that they can be interpretted as differential forms.
"""
function isAntisymmetric(T::TensorSuperType, iset_name::Union{Nothing,String}=nothing)
    if iset_name === nothing
        inl = length(indices(T))
        if isempty(tableaux(T))
            return inl <= 1 #no tableau -> still a 0 or 1-form
        end
        #should add some version to check total asym on certain inds (eg H^a_μνσ)
        tab = first(tableaux(T)) #totally antisymm so only should have 1 anyway (in most cases...)
        # println(tab.row_lengths)
        return length(tab.row_lengths) == inl#should be a column of n rows
    else
        filtered_pos = [
            i for (i, ind) in enumerate(indices(T))
                  if ind[1].index_set.name == iset_name
        ]
        isempty(tableaux(T)) && return length(filtered_pos) <= 1 #vectors and scalars are deg 1 and 0 forms, so "antisym"
        # Check whether any tableau corresponds to the spacetime indices
        for tab in tableaux(T)
            # println(tab)
            if tab.index_assignment == filtered_pos #same assignments
                return length(tab.row_lengths) == length(filtered_pos) #antisym assignment
            end
        end
        return false
    end
end

"""
Checks a tensor's Young Tableau to see if it is a totally symmetric tensor.
Tensors with 0 or 1 index will be considered symmetric.
"""
function isSymmetric(T::TensorSuperType, iset_name::Union{Nothing,String}=nothing)
    if iset_name === nothing
        inl = length(indices(T))

        isempty(tableaux(T)) && return inl <= 1

        tab = first(tableaux(T)) #totally symm so only should have 1 anyway
        return length(tab.row_lengths) == 1 &&
               tab.row_lengths[1] == inl  # single row of length n
    else
        filtered_pos = [
            i for (i, ind) in enumerate(indices(T))
                  if ind[1].index_set.name == iset_name
        ]

        isempty(tableaux(T)) && return length(filtered_pos) <= 1

        for tab in tableaux(T)
            if tab.index_assignment == filtered_pos
                return length(tab.row_lengths) == 1 &&
                       tab.row_lengths[1] == length(filtered_pos)
            end
        end

        return false
    end
end




function hasIndices(T::TensorSuperType)
    return length(indices(T)) != 0
end

function hasIndices(T::TensorTerm)
    return length(indices(T)) != 0
end


function isZeroTens(T::TensorExpression)
    if length(summands(T)) == 0 || (length(summands(T)) == 1 && isZero(summands(T)[1][1]))
        return true
    else
        return false
    end
end
function isZeroTens(T::SymorNum)
    return isequal(T, 0)
end

"""Returns false, or the TensorTerm. Casts Tensors to a TensorTerm."""
function isTensorTerm(T::TensorExpression)
    if length(summands(T)) == 1 && isequal(summands(T)[1][1], 1)
        return summands(T)[1][2]
    else
        return nothing
    end
end
isTensorTerm(T::TensorSuperType) = TensorTerm(T)
isTensorTerm(T::TensorTerm) = T


"""
Check if a TensorExpression is a single summand of one tensor with a 1 coefficient. Returns the tensor or false.
"""
function isTensor(T::TensorExpression)
    term = isTensorTerm(T)
    if term !== nothing && length(term.terms) == 1
        return term.terms[1]
    else
        return nothing
    end

end
"""
Check if a TensorTerm is a single tensor. Returns the tensor or false.
"""
function isTensor(T::TensorTerm)
    if length(terms(T)) == 1
        return terms(T)[1]
    else
        return nothing
    end

end
isTensor(T::TensorSuperType) = T



# function isZero(x)#1.11 version
#     return ((typeof(x) <: SymorNum) && isequal(x, 0))
# end

#1.12
isZero(x::Number) = iszero(x)
function isZero(x)
    x isa SymorNum || return false
    SymbolicUtils.isconst(x) && return iszero(SymbolicUtils.unwrap_const(x))
    return false
end

function _same_coords(a::Vector{Coordinate}, b::Vector{Coordinate})
    length(a) == length(b) || return false
    for x in a
        x in b || return false
    end
    return true
end