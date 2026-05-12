######################
## Functions for extracting basic data from tensors
######################

function get_data(T<:TensorSuperType)
    return T.data
end

function name(T<:TensorSuperType)
    return T.data.name
end

function indices(T<:TensorSuperType)
    return T.data.indices
end

function function_of(T<:TensorSuperType)
    return T.data.function_of
end

function weights(T<:TensorSuperType)
    return T.data.weights
end

function tableaux(T<:TensorSuperType)
    return T.data.tableaux
end





"""
Checks a tensor's Young Tableau to see if it is a totally antisymmetric tensor.
Tensors with 0 or 1 index will be considered antisymmetric, so that they can be interpretted as differential forms.
"""
function isAntisymmetric(T::TensorSuperType)
    inl=length(indices(T))
    if isempty(tableaux(T))
        return  inl<= 1 #no tableau -> still a 0 or 1-form
    end
    #should add some version to check total asym on certain inds (eg H^a_μνσ)
    tab = tableaux(T)[1] #totally antisymm so only should have 1 anyway
    # println(tab.row_lengths)
    return length(tab.row_lengths) == inl#should be a column of n rows
end

"""
Checks a tensor's Young Tableau to see if it is a totally symmetric tensor.
Tensors with 0 or 1 index will be considered symmetric.
"""
function isSymmetric(T::TensorSuperType)
    if isempty(T.tableaux)
        return length(T.indices) <= 1 #no tableau 
    end
    tab = T.tableaux[1] #totally symm so only should have 1 anyway
    return tab.row_lengths[1] == length(T.indices)#should only have a single row of length n
end
