
#since i remember nothing, lets start fresh

function same_tensor_mod_syms(T1::TensorSuperType, T2::TensorSuperType)

    #names don't match
    if T1.name != T2.name
        return (false, 0)
    end
    dums1 = get_dummies(T1)
    dums2 = get_dummies(T2)
    #if we dont need to worry about dummies, do simple algo
    if isempty(dums1) && isempty(dums2)
        out = dumless_same_tensor_mod_syms(T1, T2)
        if out != 0
            return (true, out)
        end
    else
        out = same_tensor_mod_dummies(T1, T2)
        if out != 0
            return (true, out)
        end
    end

    return (false, 0)
end


# this works for simple tensors with no dummies

function do_swaps(T1::TensorSuperType, patt::OrderedSet, movers::Dict) #where {I<:IndexSuperType}
    #elemetns in patt: indiceobj, intpos
    #elemtns in mover: intpos, indiceobj
    if length(patt) ==1 || length(patt)==0
       return 1 #no swap needed
    else
        #loop over movers, see if a move exists
        for (intpos, indiceobj) in movers
            if patt[1][1] == indiceobj #same index, different pos
                for tab in T1.tableaux
                    # println(patt[1][2])
                    if tab.adj_matrix[patt[1][2], intpos] != 0
                        grad = tab.adj_matrix[patt[1][2], intpos]
                        movers[intpos] = movers[patt[1][2]]
                        delete!(movers, patt[1][2])
                        popfirst!(patt)
                        return grad * do_swaps(T1, patt, movers)
                        break
                    end

                end
            else
                continue
            end
        end
        #no matches
        return 0
    end
end


# the "simple algo" in question
function dumless_same_tensor_mod_syms(T1::TensorSuperType, T2::TensorSuperType)
    inds1 = get_index_index(T1)
    inds2 = get_index_index(T2)
    needperm1 = setdiff(inds1, inds2) #indices of T1 that don't match T2
    needperm2 = setdiff(inds2, inds1) #indices of T2 that don't match T1
    inds_to_move = index_index_to_dict(needperm2)

    
    # println("do dum by default")
    # println(T1)
    # println(needperm1)
    # println(inds_to_move)

    out = do_swaps(T1, needperm1, inds_to_move)
    # println(out)

    return out
end




# the "simple algo" in question (again)
function to_dumless_same_tensor_mod_syms(T1::TensorSuperType, T2::TensorSuperType)
    #name,indpos
    # println(T1)
    # println(T2)

    inds1 = get_index_index(T1)
    dums1=get_dummies(T1)
    dums2=get_dummies(T2)
    inds2 = get_index_index(T2)
    # #remove dummies
    # println("inds before")
    # println(inds1)
    # println(dums1)
    inds1=OrderedSet([pair for pair in (inds1) if !(pair[1] in dums1)])
    # println("inds after")
    # println(inds1)
    inds2=OrderedSet([pair for pair in (inds2) if !(pair[1] in dums2)])

    # println(inds2)

    needperm1 = setdiff(inds1, inds2) #indices of T1 that don't match T2
    needperm2 = setdiff(inds2, inds1) #indices of T2 that don't match T1
    inds_to_move = index_index_to_dict(needperm2)



    # println("rem dums verios")
    # println(T1)
    # println(needperm1)
    # println(inds_to_move)

    out = do_swaps(T1, needperm1, inds_to_move)

    # println("todumless")
    # println(out)
    return out
end



#now lets try and get it working if there are dummy indices
#(so beautiful)

function same_tensor_mod_dummies(T1::TensorSuperType, T2::TensorSuperType)

    #best approach (maybe): can we swap all the dummies into correct positions first using modified version of above?
    dums11 = get_dummies_with_pos_index_and_pos_ul(T1)
    #elements are (i, pair), where pair is name=>upper/lower
    contractions1 = [(x, y) for x in dums11 for y in dums11 if x[2][1] == y[2][1] && x != y]
    #prev, pairs in contr were just ints

    dums22 = get_dummies_with_pos_index_and_pos_ul(T2)
    contractions2 = [(x, y) for x in dums22 for y in dums22 if x[2][1] == y[2][1] && x != y]

    # create a multigraph representing the tensor indices. weights 1 for sym, weights -1 for antisym
    # and a 2nd edge if there is a contraction between two indices. (with weight 2, just to make sure things play nice)
    # then, if a graph isomorphism exists, the dummies can be permuated so that theyre in the same positions.
    # in fact, the isomorphism will tell us exactly which ones to swap, which we'll do, and then ignore the dummies thereafter

    # a lil bit of python to hold me tight
    g1 = nx.MultiGraph()
    g1.add_nodes_from(range(1, length(T1.indices)))
    #add data from adj matrices
    for tab in T1.tableaux
        adm = tab.adj_matrix
        dims = size(adm)[1]
        for i in 1:dims
            for j in (i+1):dims
                g1.add_edge(i, j, weight=adm[i, j])
            end
        end
    end
    g2 = g1.copy()
    #add an extra edge for contractions

    # for simplicity, we should add an extra self loop for upper indices to make it easier later!

    for pair in contractions1
        g1.add_edge(pair[1][1], pair[2][1], weight=2)
        # println(pair)
        # println(pair[1][2][2])
        if pair[1][2][2]==upper
            g1.add_edge(pair[1][1], pair[1][1], weight=-2)
        else
            g1.add_edge(pair[2][1], pair[2][1], weight=-2)
        end
    end
    for pair in contractions2
        g2.add_edge(pair[1][1], pair[2][1], weight=2)
        if pair[1][2][2]==upper
            g2.add_edge(pair[1][1], pair[1][1], weight=-2)
        else
            g2.add_edge(pair[2][1], pair[2][1], weight=-2)
        end
    end
    # println("our graphs are")
    # #p# Draw the graph





    # pos = nx.spring_layout(g1)  # Layout algorithm, you can choose other layouts too
    # nx.draw(g1, pos, with_labels=true, node_color="skyblue", node_size=800, font_size=10, font_color="black")
    # # Extract edge weights as labels
    # edge_labels = Dict{Tuple{Int,Int},Float64}()
    # for (u, v, d) in g1.edges(data=true)
    #     edge_labels[(u, v)] = d["weight"]
    # end
    # # Add edge weight labels to the graph
    # nx.draw_networkx_edge_labels(g1, pos, edge_labels=edge_labels, font_size=10)
    # plt.title("Tensor")
    # plt.show()

    # # plot_graph(g1)

    # pos = nx.spring_layout(g2)  # Layout algorithm, you can choose other layouts too
    # nx.draw(g2, pos, with_labels=true, node_color="skyblue", node_size=800, font_size=10, font_color="black")
    # # Extract edge weights as labels
    # edge_labels = Dict{Tuple{Int,Int},Float64}()
    # for (u, v, d) in g2.edges(data=true)
    #     edge_labels[(u, v)] = d["weight"]
    # end
    # # Add edge weight labels to the graph
    # nx.draw_networkx_edge_labels(g2, pos, edge_labels=edge_labels, font_size=10)
    # plt.title("Tensor 2")
    # plt.show()

    ISO = nx.vf2pp_isomorphism(g2, g1)
    # println(ISO)

    #so tldr from all this:
    # we now have an iso of graphs ie we can swap indices.
    # the iso gives a map from g2 to g1 that guarentees that dummy indices are in the same place (in the right pairs) under the mapping
    # (but other indices are just some arb. allowed permutation)
    # we should thus make a permutation on T2 which brings it to the form of T1. manage the grading of this
    #also tho pay attention to which indices are upper and lower in the result since may be important?
    #regardless, once this is done, we can just ignore dummies
    # println(ISO)

    if ISO == nothing
        return 0
    end

    T3=copy(T2)
    #we'll map dummies of T2 to dummies of T1, resulting in a new tensor with some stuff in the wrong place
    #we'll next compute the grading of this process, then palm off T3 do do patt matching with T1 ignoring the dummies

    for (oldpos, newpos) in pairs(ISO)
        T3.indices[newpos]=T2.indices[oldpos]
    end

    # println(T2)
    # println(T3)

    out=bubble_sort_index_with_grading!(T2, ISO)
    grad=out[1]
    # println("ubblesort grad")
    # println(grad)

    grad*=to_dumless_same_tensor_mod_syms(T1,T3)
    #now, check T1, T3
    # println("dumlessgrad")
    # println(grad)

    return grad


    # # need to check this iso to find grading

    # ######

    # if ISO != nothing
    #     println("ladies and gentlemen, we have isomorphism")
    #     println(ISO)
    #     return 1
    # else
    #     return 0
    # end

end








# function plot_graph(g1)
#     @pywith g1 begin
#         pos = nx.spring_layout(g1)  # Layout algorithm, you can choose other layouts too
#         nx.draw(g1, pos, with_labels=true, node_color="skyblue", node_size=800, font_size=10, font_color="black")
#         # Extract edge weights as labels
#         edge_labels = Dict{Tuple{Int,Int},Float64}()
#         for (u, v, d) in g1.edges(data=true)
#             edge_labels[(u, v)] = d["weight"]
#         end
#         # Add edge weight labels to the graph
#         nx.draw_networkx_edge_labels(g1, pos, edge_labels=edge_labels, font_size=10)
#         plt.title("Tensor")
#         plt.show()
#     end
# end





















# # check if two tensors are the same, mod some symmetry dictated by young tableaux
# # note that we must account for dummy indices. individually, there must be some
# # check to see if all pairs of dummy indices can be permuted so that the configuration
# # is "reachable". Then, take this new tensor with same dummies, and see if the rest of 
# # the indices work. note that dummy indices may have different names which makes things hard

# function same_tensor_mod_syms(T1::TensorSuperType,T2::TensorSuperType)

#     #names don't match
#     if T1.name!=T2.name
#         return (false,0)
#     end

#     dums1=get_dummies(T1)
#     dums2=get_dummies(T2)

#     # take the dummies of one tensor, and see if the dummies of the other can be moved to make them match

#     # maybe a while loop is better? keep going till all pairings checked
#     grading=1

#     for dum1 in dums1
#         #find the position of this pair of dummies in the first tensor
#         occurrences1 = [i for (i, pair) in enumerate(T1.indices) if pair[1] == dum1] 

#         for dum2 in dums2
#             #find the position of this pair of dummies in the 2nd tensor
#             occurrences2 = [i for (i, pair) in enumerate(T2.indices) if pair[1] == dum2] 
#             #...
#             println(occurrences1)
#             println(occurrences2)

#             #if there exists a swap so that these dummies are on same indices do the swap and then 
#             #continue to the next pair dum1 (keep track of grading!)
#             yts=T1.tableaux # same tensor thus same tab

#             for occ1 in occurrences1
#                 for occ2 in occurrences2
#                     for yt1 in yts
#                         if ty1[occ1,occ2]!=0
#                             # match found
#                             break
#                         end
#                     end
#                 end
#                 #some sort of no match found termination idk
#             end


#         end
#     end

#     # if all dummies matched, good sign. now, check remaining indices, leaving the dummies untouched



# end






# function same_tensor_mod_tableaux(T1::TensorSuperType, T2::TensorSuperType)

#     #names don't match
#     if T1.name != T2.name
#         return (false, 0)
#     end

#     inds1 = get_free_indices_nopos(T1)
#     inds2 = get_free_indices_nopos(T2)

#     #not the same free indices
#     if inds1 != inds2
#         return (false, 0)
#     end

#     # free indices are a permutation
#     inds1 = get_index_index(T1)
#     inds2 = get_index_index(T2)
#     # println(inds1)
#     # println("these are the indind paris")
#     dums1 = get_dummies(T1)
#     dums2 = get_dummies(T2)
#     # 
#     needperm1 = setdiff(inds1, inds2)
#     needperm2 = setdiff(inds2, inds1)
#     #everything that is not in the same position in the two tensors (we'll ingnore dummies)


#     grading = 1
#     found = false

#     permed = Set{IndexSuperType}()

#     while needperm1 != needperm2
#         indpair1 = pop!(needperm1)

#         for indpair1 in needperm1
#             if indpair1[1] in dums1 || indpair1[1] in permed
#                 continue #is a dummy, ignore 
#             end
#             for indpair2 in needperm2
#                 # println(indpair1)
#                 # println(indpair1[1])
#                 # println(indpair1[2])
#                 if indpair1[1].name != indpair2[1].name
#                     continue #not the same index
#                 else
#                     if indpair1[2] == indpair2[2]
#                         println("found a thingo")
#                         break #already swapped
#                     end
#                     # same index name, see if exists permuate
#                     println(T1.tableaux)
#                     for tab in T1.tableaux
#                         if tab.adj_matrix[indpair1[2], indpair2[2]] == 0
#                             println("not in tableau")
#                             continue
#                         else
#                             grading *= tab.adj_matrix[indpair1[2], indpair2[2]]
#                             println("in tabelau")
#                             found = true
#                             println(indpair1[1])
#                             println(indpair2[1])
#                             for ind3 in needperm2
#                                 if ind3[2] == indpair1[2]
#                                     println("doing a swap")
#                                     ind3 = ind3[1] => indpair2[2]
#                                     println(ind3)
#                                     println(needperm1)
#                                     println(needperm2)
#                                     break
#                                 end
#                             end
#                             indpair2 = indpair2[1] => indpair1[2]
#                             push!(permed, indpair1[1])
#                             println(permed)
#                             break#found a perm
#                         end
#                     end
#                     if found == false
#                         println("NO MATCH FOUND")
#                         return (false, 0) #no match found in any tableaux
#                     end
#                     found = false
#                 end
#             end

#         end
#     end
#     println("SAME TENSORS")
#     return (true, grading)
# end


# #broken
# function combine_same_tensors(TEx1::TensorExpression, TEx2::TensorExpression)
#     println("running this 22222")
#     if length(TEx1.summands) == 1 && length(TEx2.summands) == 1
#         (res, val) = same_tensor_mod_tableaux(TEx1.summands[1][2], TEx2.summands[1][2])
#         if res
#             new_co = val * TEx1.summands[1][1] + TEx2.summands[1][1]
#             return (true, new_co)
#         end
#     end
#     return (false, 0)
# end

# function combine_same_tensors(TEx1::TensorTerm, TEx2::TensorTerm)
#     println("running this 11111")
#     if length(TEx1.terms) == 1 && length(TEx2.terms) == 1
#         println(TEx1.terms[1])
#         println(TEx2.terms[1])
#         (res, val) = same_tensor_mod_tableaux(TEx1.terms[1], TEx2.terms[1])
#         println("got here")
#         if res
#             new_co = val
#             return (true, new_co)
#         end
#     end
#     return (false, 0)
# end

# function decompose2tabelaux(T1::TensorSuperType)

#     indnum = length(T1.indices)
#     permutations = permutations(indnum, T1.yo)
#     permutation_matrices = []
#     for yt in T1.tableaux
#         push!(permutation_matrices, yt.adj_matrix)
#     end

#     for perm_matrix in permutation_matrices
#         for perm in permutations
#             sign = sign_of_permutation(perm, perm_matrix)

#             # Apply the permutation and sign to the tensor
#             permuted_tensor = permutedims(tensor, perm)       ########write the function to create all perms
#             ### need to compute coefficnets in front

#             # Update the result tensor with the permuted tensor and sign
#             result += sign * permuted_tensor
#         end
#     end


# end


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
