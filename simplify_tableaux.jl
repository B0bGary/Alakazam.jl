







#since i remember nothing, lets start fresh

function same_tensor_mod_syms(T1::TensorSuperType,T2::TensorSuperType)

      #names don't match
      if T1.name!=T2.name
        return (false,0)
    end


    dums1=get_dummies(T1)
    dums2=get_dummies(T2)

    #if we dont need to worry about dummies, do simple algo
    if isempty(dums1) && isempty(dums2)
        return dumless_same_tensor_mod_syms(T1, T2)
    end

end


# the "simple algo" in question
function dumless_same_tensor_mod_syms(T1::TensorSuperType,T2::TensorSuperType)
 # free indices are a permutation
 inds1=get_index_index(T1)
 inds2=get_index_index(T2)
 println(inds1)
 println("these are the indind paris")
#  dums1=get_dummies(T1)
#  dums2=get_dummies(T2)
 # 
 needperm1=setdiff(inds1,inds2) #indices of T1 that don't match T2
 needperm2=setdiff(inds2,inds1) #indices of T2 that don't match T1
 #everything that is not in the same position in the two tensors (we'll ingnore dummies)

 #we need to see if there's a map from T2 to T1 via the tableaux


end









# check if two tensors are the same, mod some symmetry dictated by young tableaux
# note that we must account for dummy indices. individually, there must be some
# check to see if all pairs of dummy indices can be permuted so that the configuration
# is "reachable". Then, take this new tensor with same dummies, and see if the rest of 
# the indices work. note that dummy indices may have different names which makes things hard

function same_tensor_mod_syms(T1::TensorSuperType,T2::TensorSuperType)

    #names don't match
    if T1.name!=T2.name
        return (false,0)
    end

    dums1=get_dummies(T1)
    dums2=get_dummies(T2)

    # take the dummies of one tensor, and see if the dummies of the other can be moved to make them match

    # maybe a while loop is better? keep going till all pairings checked
    grading=1

    for dum1 in dums1
        #find the position of this pair of dummies in the first tensor
        occurrences1 = [i for (i, pair) in enumerate(T1.indices) if pair[1] == dum1] 

        for dum2 in dums2
            #find the position of this pair of dummies in the 2nd tensor
            occurrences2 = [i for (i, pair) in enumerate(T2.indices) if pair[1] == dum2] 
            #...
            println(occurrences1)
            println(occurrences2)

            #if there exists a swap so that these dummies are on same indices do the swap and then 
            #continue to the next pair dum1 (keep track of grading!)
            yts=T1.tableaux # same tensor thus same tab

            for occ1 in occurrences1
                for occ2 in occurrences2
                    for yt1 in yts
                        if ty1[occ1,occ2]!=0
                            # match found
                            break
                        end
                    end
                end
                #some sort of no match found termination idk
            end


        end
    end

    # if all dummies matched, good sign. now, check remaining indices, leaving the dummies untouched



end












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

    permed= Set{IndexSuperType}()

    while needperm1 !=needperm2
        indpair1=pop!(needperm1)

        for indpair1 in needperm1
            if indpair1[1] in dums1 || indpair1[1] in permed
                continue #is a dummy, ignore 
            end
            for indpair2 in needperm2
                # println(indpair1)
                # println(indpair1[1])
                # println(indpair1[2])
                if indpair1[1].name!=indpair2[1].name
                    continue #not the same index
                else
                    if indpair1[2]==indpair2[2]
                        println("found a thingo")
                        break #already swapped
                    end
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
                            println(indpair1[1])
                            println(indpair2[1])
                            for ind3 in needperm2
                                if ind3[2]==indpair1[2]
                                    println("doing a swap")
                                    ind3=ind3[1]=>indpair2[2]
                                    println(ind3)
                                    println(needperm1)
                                    println(needperm2)
                                    break
                                end
                            end
                            indpair2=indpair2[1]=>indpair1[2]
                            push!(permed, indpair1[1])
                            println(permed)
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

function decompose2tabelaux(T1::TensorSuperType)

    indnum=length(T1.indices)
    permutations = permutations(indnum,T1.yo)
    permutation_matrices=[]
    for yt in T1.tableaux
        push!(permutation_matrices,yt.adj_matrix)
    end

    for perm_matrix in permutation_matrices
        for perm in permutations
            sign = sign_of_permutation(perm, perm_matrix)
            
            # Apply the permutation and sign to the tensor
            permuted_tensor = permutedims(tensor, perm)       ########write the function to create all perms
                                                            ### need to compute coefficnets in front
            
            # Update the result tensor with the permuted tensor and sign
            result += sign * permuted_tensor
        end
    end


end


function sign_of_permutation(perm, perm_matrix)
    sign = 1
    for i in 1:length(perm)
        for j in (i+1):length(perm)
            if perm[i] > perm[j]
                sign *= perm_matrix[indices[i], indices[j]]
            end
        end
    end
    return sign
end
