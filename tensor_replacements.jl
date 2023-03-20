


## O(n) runtime approx amortised
function term2histogram(term::TensorTerm)
    pat = Dict{String,Pair{Number,Vector{Number}}}() # tensor, (counts, [positions in expression]) #dont check indices yet
    #some sort of dictionary with repeated elements, constant time access etc after hashing

    for (i, tens) in enumerate(term.terms)

        old_val = get(pat, tens.name, nothing)
        if old_val === nothing
            push!(pat, tens.name => (1 => [i])) #if not in hist, add it 
        else

            push!(pat, tens.name => (old_val[1] + 1 => push!(old_val[2], i))) #increment the count

        end
    end
    return invert_histogram(pat)
end


#change the histogram to make the counts the keys
function invert_histogram(hist::Dict{String,Pair{Number,Vector{Number}}})

    #maybe n^2 time? less i think. n for going through all, Ordered dict is linked list, at worst O(n) insertion
    new_hist = OrderedDict{Number,Vector{Pair{String,Vector{Number}}}}() #counts => [tensorname =>[positions in term]]
    for (key, value) in zip(keys(hist), values(hist))
        # println(key)
        # println(value)
        old_val = get(new_hist, value[1], nothing)
        if old_val === nothing
            push!(new_hist, value[1] => [key => value[2]]) #if not in hist, add it 
        else

            push!(new_hist, value[1] => push!(old_val, key => value[2])) #increment the count

        end

        # push!(new_hist, value)

    end
    return new_hist

end



######################################################################################################
#######################################################################################################################
function sub_tensor(exp::TensorTerm, pattern::TensorTerm, rep::TensorTerm, repeat::Bool=true)

    #lets do no indices version first

    exp_hist = term2histogram(exp)
    pat_hist = term2histogram(pattern)


    dummy_matches = Vector{Pair{Number,Number}}() # dummy in pos x matches expression pos y

    exps = zip(keys(exp_hist), values(exp_hist))


    flag_no_match = true
    # for each term in pattern, try and find match in the expression. if can't find, then return exp
    for (pkey, pvalue) in zip(keys(pat_hist), values(pat_hist))

        for (ekey, evalue) in exps

            if ekey >= pkey
                #possible match, expression has more of this symbol than pattern
                for ptens in pvalue #loop through TensorSuperType pairs
                    for etens in evalue
                        if ptens[1][1] == "?"#if pattern is a dummy

                        end
                    end


                end

            end



        end
        #if loop through all with no matches, then break
        if flag_no_match
            break
        end



    end


    return exp

end
########################################################################################################################################


function match_dummy_logic(dummy_matching::Dict{String,String},dummy_name::String, exp_name::String)
    dummy=copy(dummy_matching)
    if get(dummy, dummy_name, nothing) == nothing #if dummy isnt already assigned a val in match
        dummy[dummy_name] = exp_name

    elseif dummy[dummy_name] != exp_name #dummy doesnt match

        return (false, dummy) #dummy tensor name pattern
    end

    return (true, dummy)

end







#dynamic prog attempt
function patt_match(pat::OrderedDict{N,Vector{Pair{String,Vector{Number}}}} where {N<:Number}, exp::OrderedDict{N,Vector{Pair{String,Vector{Number}}}} where {N<:Number},
    tensorpatt::TensorTerm, tensorexp::TensorTerm, dummy_matching::Dict{String,String}=Dict{String,String}())
    #list of potential dummies is included, by default its empty

    #dict are of form  [#occurances => [tensorname =>[postions],  tensorname =>[postions]], ]

    if length(pat) == 1 #if only one thing in the pattern   ---- No, this is only one type of frequency
        println(first(pat)[2])
        (pkeys, pvalues) = first(pat)
        #####all of this works!!!!!! if single tensor left in pattern!!   
        if length(first(pat)[2]) == 1 #if only one type of tensor in the pattern



            # if first(pat)[1]==1 #if only needs to match once only one thing left to match



            #handles the case where only one this etc
            if length(exp) == 1
                looper = [first(exp)]
                println("this one used")
            else
                looper = zip(keys(exp), values(exp))
            end
            println(looper)
            for (ekey, evalue) in looper
                if pkeys[1] <= ekey #if occ in pat is less than occ in exp

                    #need to loop over evals!
                    for eval in evalue
                        println("occurs less")
                        println(pvalues[1][1])
                        println(eval[1])

                        pdummy = false
                        println(pvalues[1][1][1])
                        if pvalues[1][1][1] == '?'
                            pdummy = true
                            println("dummy found")
                        end

                        if pdummy || pvalues[1][1] == eval[1]#names match
                            dummy_temp=dummy_matching
                            if pdummy
                
                                out = match_dummy_logic(dummy_matching, pvalues[1][1], eval[1])
                                if out[1] == false
                                    continue
                                else
                                    dummy_temp = out[2]
                                end
                            end
                            println(dummy_matching)

                            println("names match")
                            match = []

                            #check indices
                            #check_index_match

                            #construct the match pairings
                            # for i in 1:pkeys[1] #match with x many terms
                            #     push!(match, pvalues[1][2][i] => eval[2][i])
                            # end

                            #check if indices match

                            j = 1
                            println(ekey)
                            for i in 1:ekey #match with x many terms
                                # println(ekey)
                                println(tensorpatt.terms[pvalues[1][2][j]].indices)
                                println(tensorexp.terms[eval[2][i]].indices)
                                if check_index_match(pvalues[1][2][j], eval[2][i], tensorpatt, tensorexp)
                                    println("indices match")
                                    push!(match, pvalues[1][2][j] => eval[2][i])
                                    j += 1
                                    if j > pkeys[1] #if matched with all
                                        println("breaking")
                                        break
                                    end

                                else
                                    println("indices dont match")
                                end
                            end
                            if j <= pkeys[1]
                                println(j)
                                println(pkeys[1])
                                if !pdummy 
                                return (false, [-1 => -1])
                                else
                                    continue
                                end
                            end


                            # #feed new back into method with quantity decreased by counts
                            # if patt_match(exp2)

                            #needs to loop over all n matches
                            dummy_matching=dummy_temp
                            return (true, match) #found the match #found the match
                        end
                    end

                end
                continue #else, just keep going
            end
            return (false, [-1 => -1]) #no match


        else

            ##now for multiple things left in the pattern of this frequency

            if length(exp) == 1
                looper = [first(exp)]
                println("this one used")
            else
                looper = zip(keys(exp), values(exp))
            end
            println(looper)
            for (ekey, evalue) in looper
                #recall, ekey=#num occurances
                #eval= vector,  [tensorname =>[postions],  tensorname =>[postions], ...]

                ########## loop through evalues?

                for eval in evalue

                    # run the loop until we find a match. 
                    # take the first item of the pattern, try and match it with this expression term. 
                    # try and match remainder of pattern with remainder of expression term

                    #freqyency with first item in the pattern
                    check = patt_match(OrderedDict(first(pat)[1] => [first(pat)[2][1]]), OrderedDict(ekey => [eval]), tensorpatt, tensorexp) #single expression term

                    if check[1] #if it's a match
                        #construct remaining terms
                        exp2 = copy(exp)
                        pat2 = copy(pat)

                        #remove the first item from patt
                        rem = first(pat)[2]
                        popfirst!(rem) #removes first item, wont have errors since we know more than one thing
                        push!(pat2, first(pat)[1] => rem)

                        # remexp=first(exp)[ekey]
                        # popfirst!(remexp) #removes first item, wont have errors since we know more than one thing
                        # push!(pat2, first(pat)[1]=> rem)


                        #??????????????????????????????    pls check this lol idk anymore man
                        println("do ti")
                        println(eval)

                        println(evalue)
                        # println(filter!(p1=>p2 -> p1 != eval[1], evalue))

                        push!(exp2, ekey => filter!(x -> x[1] != eval[1], evalue))#remove tensor of this name
                        if ekey - pkeys[1] >= 1
                            # vnew=
                            if get(exp2,ekey-pkeys[1],nothing)==nothing
                                push!(exp2, (ekey - pkeys[1]) => [ eval[1] => eval[2][pkeys:end]])
                            else
                            println(exp2[ekey-pkeys[1]])
                            println(eval[1])
                            push!(exp2, (ekey - pkeys[1]) => push!(exp2[ekey-pkeys[1]], eval[1] => eval[2][pkeys:end])) #decrements occ of tensor 
                            end
                        end
                        #    push!(pat2, ekey => filter!(x=>y -> x != evalue[1], eval))#remove tensor of this name
                        ######################


                        #remove freqency of patt from exp


                        #return match of the remaining. effectively, we're peeling off one term at a time. we can then add mimo by saving things that dont match
                        out = patt_match(pat2, exp2, tensorpatt, tensorexp)
                        if out[1] == true
                            return (true, vcat(check[2], out[2]))
                        else
                            # 
                            continue
                        end

                    else
                        continue #else, look at the next term in the expression
                    end

                end


                return (false, [-1 => -1]) #no match


                #(will be more complex once we use wildcards)


                # decremement 



                #for each term in the expression


                #this is where we want to use recursion
                # return patt_match(first_term in_pattern, exp term) * patt_match(rest_of_terms_in_pattern, rest_of_exp)
                #can defs use mimo here


            end





        end


    else  #if multiple things to match in the pattern
        #use short circuit eval

        outs = []

        for (pkeys, pvalues) in pat #if multiple freqencies

            out = patt_match(OrderedDict(pkeys => pvalues), exp, tensorpatt, tensorexp)


            if out[1] == false
                return (false, [-1 => -1])
            else
                if isempty(outs)
                    outs = (true, out[2])
                else
                    outs = (true, vcat(outs[2], out[2]))
                end
            end

        end

        return outs

    end


end



function check_index_match(loc_pat::Number, loc_exp::Number, pat::TensorTerm, exp::TensorTerm)

    return pat.terms[loc_pat].indices == exp.terms[loc_exp].indices
end


function sub_tensor(exp::TensorTerm, pattern::TensorTerm, rep::TensorTerm, repeat::Bool=true)

    temp_dummy_tensor = Dict{String,TensorSuperType}() # keeps track of matches in a dummy pattern
    temp_dummy = Dict{String,Index}()

    exp_hist = term2histogram(exp)
    pat_hist = term2histogram(pattern)


    # go through the pattern and find mathces in the expression. if the names dont match, break and return exp no replacements

    #should construct a possible replacement. if it doesnt work when matching indices, should go to the next possible replacement etc 
    #else, if the names match, check the indices 

    #so should have a mapping ie ?A?B C and DCF => ?A - D, ?B - F, C-C, then check indicies. if indices work, use temp_dummy etc to construct
    #the replacement with the correct dummies

    #how can I loop through all possible matches? get max key, compare, slowly iterate through. once 



    matches = Dict{Number,Number}



    # pairs where pattern matched a tensor
    # i think that the keys should be exp, values are the expression. that way we can look through the keys and replace with the pattern replacements.
    # it will also make sure we dont have double-ups, can't replace the same tensor twice
    # use ints as key values so we can extract them directly from the original array


    # how to do this best?
    #suppose pattern.check first item of pattern to see if it matches anything. it may match with multiple things, in which case we need to store them all.
    #  -- i think fastest way is a depth first search. find first instance of a match, go to the 2nd item in pattern and look for match.
    # -- we can also then use mimo maybe to record if things have already matched??


    # should be a fast way eg mimo that makes it less than n^2 time
    # need a way to recognise  ?A ?A, same thing repeated but also ?A ?B


    # we need to keep track of all matches. else its n! time i think.
    # i think we need some sort of tree???
    # find all things that match first expression. add to an array. then first each match, check the 2nd pattern in a depth first way

    for tenpind in eachindex(pattern.terms)

        for tenseind in eachindex(exp.terms)

            # how do we deal with the rep? i think that we only know the correct rep after checking all the patterns
            # eg if we have $A C -> C B^i or something, only know how to replace if the entire pattern matches


            # if sub_tensor(exp.terms[tenseind], pattern.terms[tenpind], ,temp_dummy_tensor,temp_dummy)
            #     push!(matches, tenseind=>tenpind) # this term in the expression is related to the pattern
            #     break
            # end

            #so, what we should do
            # look at each item in the pattern. for each item, find a match. if no matches, break obvs. 
            # if we find a match for the entire pattern, construct a replacement, and then replace.
            # if we're repeating, we can probs skip over terms that we've already check?? nah because we might do a replacement, and then
            # a new valid replacement occurs
            # dont forget to keep track of gradings



        end

        # things could be made faster by using a LinkedList instead of a vector
        # then,we could remove anobject from the vector,and instead of shifting the indices, just change the
        # "next" item in the list by accessing previous/next of the removed item
        # idk about indexing tho
        # apparaently Julia built in doubley linked lsit from DataStructures has constant time access like a vector
        # this would be a good solution to eventually get around to implementing

        #to construct a new term at the end, use filter to get everything that's not being replaced, rather the removing one by one

    end


    #for each tensor term, see if it matches any of the replacement terms
    #add to a lsit of matches. if not repeat, stop once a single match found. else, find all matches
    #to replace go through each set of matches and replace

    # need to account for the fact that placeholders may be over multiple tensors
    # for example, need a way to recognise  ?A ?A, same thing repeated. 



    return exp
end


function build_replacement()
    rep_dum = typeof(rep) <: DummyPatternTensor
    rep = typeof(exp)(rep) #makes sure rep isnt a dummy if success


    if rep_dum
        # println("Namema")
        # println(exp.name)
        rep.name = exp.name
        # println(rep.name)
    end
    #needs to make sure dummy indices match pattern
    ri = rep.indices
    # println("about to replace indices with non dummies")
    for index in eachindex(ri)
        # println("indices exist")
        if typeof(ri[index][1]) <: DummyPatternIndex
            # println("dummy")
            # println(temp_dummy)
            # println(ri[index][1])
            if get(temp_dummy, ri[index][1].name, nothing) == nothing #problem with checking equality since technically different objects #use names instead
                throw(InvalidReplacementPattern())
            else
                ri[index] = temp_dummy[ri[index][1].name] => ri[index][2] # maps index from original expression to replacement
                # println("final assignment")##not happen
                # println(ri[index])
            end
        end
    end

end



# function parse_pattern2tenterm(pattern::String)


#returns a tensor term
function parse_pattern2tensor(pattern::String)
    pattern = strip(pattern)#remove spaces
    it = Iterators.Stateful(pattern)
    # # println(it)
    #add a way user can do custom patterns, eg \eps for eps tensor, \delta for delta func etc. some sort of vector containg possible presets (or a tree!)
    #else itll always output as a default tensor
    tensor_term = Vector{TensorSuperType}() #array of tensors
    tensor = []


    # need to : make many tensors, put them all in a tensor term


    #Construct the tensor
    for char in it #should only get here at the start of a tensor if valid syntax!!!!!!!!!!

        # # println("starting loop")
        # # println(char)
        if string(char) == "?"
            # # println("did this")
            s = string(iterate(it)[1])
            # # println(s)
            tensor = DummyPatternTensor(s) #?A  ie the name isn't important
        # # println(tensor)
        else
            # # println(string(char))
            tensor = Tensor(string(char)) #A     name is important #maybe something about coord dependance
            # # println(tensor)

        end
        next = peek(it)
        if next !== nothing
            char = iterate(it)[1]
            # # println(char)
            # # println("should be next^")
            # # println(string(char))
            #now, do indices
            #   while string(char) == "^" ||string(char)=="_"
            while string(char) == "^" || string(char) == "_"
                if string(char) == "^"
                    # # println("made it here, next is")
                    # # println(string(peek(it)[1]))
                    if string(peek(it)[1]) == "{"
                        iterate(it)
                        # # println("made it here,ffff next is")
                        # # println(string(peek(it)[1]))
                        ## do indices
                        ind = 1
                        while string(peek(it)[1]) != "}"
                            # # println("while")
                            if length(tensor.indices) < ind #if not already an index

                                if string(peek(it)[1]) == "?"
                                    iterate(it)
                                    push!(tensor.indices, DummyPatternIndex(string(iterate(it)[1])) => upper)
                                else

                                    push!(tensor.indices, Index(string(iterate(it)[1])) => upper)

                                end

                            elseif string(peek(it)[1]) != "*" #if already, and current value is not a *
                                if string(peek(it)[1]) == "?"
                                    iterate(it)
                                    tensor.indices[ind] = DummyPatternIndex(string(iterate(it)[1])) => upper
                                else

                                    tensor.indices[ind] = Index(string(iterate(it)[1])) => upper
                                end
                            else
                                iterate(it)
                            end
                            ind += 1
                        end
                        iterate(it)

                    else
                        push!(tensor.indices, Index(string(iterate(it)[1])) => upper) #single upper index
                    end
                end
                # char = iterate(it)[1]
                # # println("next show lowre")
                # # println(char)
                if string(char) == "_"
                    if string(peek(it)[1]) == "{"
                        iterate(it)
                        ## do indices
                        ind = 1
                        while string(peek(it)[1]) != "}"
                            # # println(tensor.indices)
                            # # println(length(tensor.indices))
                            if length(tensor.indices) < ind #if not already an index, stops out ofbounds error. add new index
                                if string(peek(it)[1]) == "?"
                                    iterate(it)
                                    push!(tensor.indices, DummyPatternIndex(string(iterate(it)[1])) => lower)
                                else

                                    push!(tensor.indices, Index(string(iterate(it)[1])) => lower)

                                end

                            elseif string(peek(it)[1]) != "*" #if already, and current value is not a *
                                if string(peek(it)[1]) == "?" #next index is a dummy, won't need exact match in replacement
                                    iterate(it)
                                    tensor.indices[ind] = DummyPatternIndex(string(iterate(it)[1])) => lower
                                    # # println( tensor.indices[ind])
                                else

                                    tensor.indices[ind] = Index(string(iterate(it)[1])) => lower
                                end
                            else
                                # # println(peek(it))
                                iterate(it) # removes * 
                            end
                            ind += 1
                        end
                        iterate(it) #remove }
                    else
                        push!(tensor.indices, Index(string(iterate(it)[1])) => lower) #single upper index
                    end

                end
                next = peek(it)
                # # println("next is")
                # # println(next)
                if next !== nothing
                    char = iterate(it)[1]
                else
                    break
                end
                # # println(char)
            end

        end
        push!(tensor_term, copy(tensor))#adds a copy
    end

    return TensorTerm(tensor_term) #should construct a term out of our set of patterns

end




function sub_tensor(exp::TensorSuperType, pattern::TensorSuperType, rep::TensorSuperType, temp_dummy_tensor::Dict{String,TensorSuperType}, temp_dummy::Dict{String,Index})


    #if replacement has dummy indices need logic to compute that
    rep_dum = typeof(rep) <: DummyPatternTensor
    rep = typeof(exp)(rep) #makes sure rep isnt a dummy if success


    if !hasIndices(pattern) && !hasIndices(exp) # if neither have indices, both scalars so do scalar replacement
        if typeof(pattern) <: DummyPatternTensor || exp.name == pattern.name
            # println("No indices, names match or dummy")
            if rep_dum
                if get(temp_dummy_tensor, pattern.name, nothing) == nothing #if dummy isnt already assigned a val in match
                    temp_dummy_tensor[pattern.name] = exp

                elseif temp_dummy_tensor[pattern.name] != exp #dummy doesnt match

                    return false #dummy tensor name pattern
                end

                rep.name = exp.name
            end
            return true
        end
        # println("No indices, names dont match")
        return false
    elseif hasIndices(pattern) && hasIndices(exp)
        if typeof(pattern) <: DummyPatternTensor || exp.name == pattern.name #both have indices, names match or one is a dummy

            #checks to make sure the dummy tensor names match
            if rep_dum
                if get(temp_dummy_tensor, pattern.name, nothing) == nothing #if dummy isnt already assigned a val in match
                    temp_dummy_tensor[pattern.name] = exp

                elseif temp_dummy_tensor[pattern.name] != exp #dummy doesnt match

                    return false #dummy tensor name pattern
                end
            end

            pat_ind = pattern.indices
            exp_ind = exp.indices
            if pat_ind == exp_ind #indices match, ez money  # no dummy pattern indices, dont have to do logic on replacement
                # println("Has indices, names match or dummy, and indices all match")
                if rep_dum
                    rep.name = exp.name
                end
                return true
            else
                # # println(pat_ind)
                # # println(exp_ind)
                pi_it = Iterators.Stateful(pat_ind)
                ei_it = Iterators.Stateful(exp_ind)
                # # println(pi_it)
                # # println(ei_it)

                # temp_dummy = Dict{String,Index}()

                for indp in pi_it
                    # # println(ei_it)
                    inde = iterate(ei_it)[1]
                    # println(inde[1])
                    # println(indp[1])
                    # println(typeof(indp[1]))
                    # println(typeof(inde[1]))
                    if indp[2] != inde[2]
                        # println("upper/lower dont match")
                        return false
                    elseif indp[1] == inde[1]
                        # println("Index matches")
                        #if names match, good, continue
                    elseif typeof(indp[1]) <: DummyPatternIndex
                        # println("dummy found")
                        if get(temp_dummy, indp[1].name, nothing) == nothing #if dummy isnt already assigned a val in match
                            temp_dummy[indp[1].name] = inde[1]
                            # println("checking index")
                            # println(inde[1])
                        elseif temp_dummy[indp[1].name] != inde[1] #dummy doesnt match
                            # println("Has indices, names match or dummy, dummy indices dont match the dummy pattern")
                            return false #dummy indice pattern didnt match
                        end

                    else
                        # println("Has indices, non-dummy index doesnt match")
                        return false #didnt match
                    end

                end
                #everything matched!
                # println("Has indices, names match or dummy, and indices all match or dummy pattern matches")


                #this whole section is for constructing the replacement right? so I could move it to a different function to construct the replacements once matches are found

                # # println(rep.name)
                return true
            end

        end
        # println("Has indices, names dont match")
        return false

    end

    return false

end







# else #needs to match more than once #frequency not equal to 1

#         (pkeys, pvalues)=first(pat)

#         #handles the case where only one this etc
#         if length(exp)==1
#             looper= [first(exp)]
#             println("this one used")
#         else
#             looper= zip(keys(exp), values(exp))
#         end
#         println(looper)
#         for (ekey , evalue) in looper
#             if pkeys[1]<=ekey #if occ in pat is less than occ in exp

#                 #need to loop over evals!
#                 for eval in evalue
#                     println("occurs less but many times")
#                     println(pvalues[1][1])
#                     println(eval[1])
#                     if pvalues[1][1]==eval[1]#names match
#                         println("names match")




#                         # return (true, [pvalues[1][2][1] => eval[2][1]]) #found the match
#                     end
#                 end

#             end
#             continue #else, just keep going
#         end
#         return (false, [-1=>-1]) #no match








# else #needs to match more than once #frequency not equal to 1
#     quant=pat[1][1] #number of times symbol occurs in the pattern
#     println("pattern is")
#     println(pat)
#     (pkeys, pvalues)=zip(keys(pat), values(pat))
#     for (ekey , evalue) in zip(keys(exp), values(exp))
#         if pkeys[1]<=ekey #if occ in pat is less than occ in exp
#             if pvalues[1][1]==evalue[1]#names match


#                 #only needed if not the last match

#                 #create copy of exp
#                 # exp2=copy(exp)
#                 # pat2=copy(pat)
#                 # # decremement 
#                 # push!(exp2, ekey => filter!(x=>y -> x != evalue[1], eval))#remove tensor of this name
#                 # push!(exp2, ekey-pkeys[1] => push!(exp2[ekey-pkeys[1]]), evalue[1]=> evalue[2][2:end]) #decrements occ of tensor

#                 # push!(pat2, ekey => filter!(x=>y -> x != evalue[1], eval))#remove tensor of this name

#                 match=[]
#                 #construct the match pairings
#                 for i in 1:pkeys[1]
#                     push!(match, pvalues[2][i] => evalue[2][i])
#                 end


#                 # #feed new back into method with quantity decreased by counts
#                 # if patt_match(exp2)

#                     #needs to loop over all n matches
#                 return (true, match) #found the match
#             end
#         end
#         continue #else, just keep going
#     end
#     return (false, [-1=>-1]) #no match

