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

    unused = [x for x in tinds.indices if x ∉ [y[1] for y in combinds]]
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

    unused = [x for x in unused if x ∉ (y[1] for y in vcat(combinds, indices(eps1)))]
    epsindlist2 = Vector{Pair{IndexSuperType,Alakazam.IndexPosition}}(
        vcat(
            [unused[i] => op_pos(poss) for i in 1:r],
            [indices(eps1)[i][1] => op_pos(poss) for i in (r+1):dim] #contract with end of other eps
        )
    )
    eps2 = EpsilonTensor(epsindlist2)
    T1c = copy(T1)
    set_indices!(T1c, [epsindlist2[i][1] => poss for i in 1:p])
    T2c = copy(T2)
    set_indices!(T2c,[epsindlist2[i][1] => poss for i in p+1:r])
 
    prefact = 1 // (factorial(p) * factorial(q) * factorial(dim - r))
    return prefact * eps1 * eps2 * T1c * T2c
end

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


