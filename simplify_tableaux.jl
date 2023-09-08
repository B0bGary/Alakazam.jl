

function same_tensor_mod_tableaux(T1::TensorSuperType,T2::TensorSuperType)
   
    #names don't match
    if T1.name!=T2.name
        return (false,0)
    end

    inds1=get_free_indices_nopos(T1)
    inds2=get_free_indices_nopos(T2)

    #not the same free indices
    if inds1!=inds2
        return (false,0)
    end

    # free indices are a permutation
    inds1=get_index_index(T1)
    inds2=get_index_index(T2)
    println(inds1)
    println("these are the indind paris")
    dums1=get_dummies(T1)
    dums2=get_dummies(T2)
    # 
    needperm1=setdiff(inds1,inds2)
    needperm2=setdiff(inds2,inds1)
    #everything that is not in the same position in the two tensors (we'll ingnore dummies)


    grading=1
    found=false



    for indpair1 in needperm1
        if indpair1[1] in dums1
           continue #is a dummy, ignore 
        end
        for indpair2 in needperm2
            # println(indpair1)
            # println(indpair1[1])
            # println(indpair1[2])
            if indpair1[1].name!=indpair2[1].name
                continue #not the same index
            else
                # same index name, see if exists permuate
                println(T1.tableaux)
                for tab in T1.tableaux
                    if tab.adj_matrix[indpair1[2],indpair2[2]]==0
                        println("not in tableau")
                        continue
                    else
                        grading*=tab.adj_matrix[indpair1[2],indpair2[2]]
                        println("in tabelau")
                        found=true
                        break#found a perm
                    end
                end
                if found==false
                    println("NO MATCH FOUND")
                return (false,0) #no match found in any tableaux
                end
                found=false
            end
        end

    end
    println("SAME TENSORS")
    return (true,grading)


end


#broken
function combine_same_tensors(TEx1::TensorExpression,TEx2::TensorExpression)
    println("running this 22222")
    if length(TEx1.summands)==1 && length(TEx2.summands)==1
        (res,val)=same_tensor_mod_tableaux(TEx1.summands[1][2], TEx2.summands[1][2])
        if res
            new_co=val*TEx1.summands[1][1]+ TEx2.summands[1][1]
            return (true,new_co)
        end
    end
    return (false,0)
end

function combine_same_tensors(TEx1::TensorTerm,TEx2::TensorTerm)
    println("running this 11111")
    if length(TEx1.terms)==1 && length(TEx2.terms)==1
        println(TEx1.terms[1])
        println(TEx2.terms[1])
        (res,val)=same_tensor_mod_tableaux(TEx1.terms[1], TEx2.terms[1])
        println("got here")
        if res
            new_co=val
            return (true,new_co)
        end
    end
    return (false,0)
end