

function same_tensor_mod_tableaux(T1::TensorSuperType,T2::TensorSuperType)
    
    for ind in 1:size(T1.indices)[1]

        indx1=T1.indices[ind][1]
        indx2=T1.indices[ind][1]
        if indx1.name==indx2.name
            continue
        end
    end

    return false

end