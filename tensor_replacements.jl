


function sub_tensor(exp::TensorTerm,pattern::TensorTerm,rep::TensorTerm, repeat::Bool=true)

    temp_dummy_tensor = Dict{String,TensorSuperType}() # keeps track of matches in a dummy pattern
    temp_dummy:=Dict{String,Index}()

    matches:=Dict{Number,Number}; # pairs where pattern matched a tensor
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
            if sub_tensor(exp.terms[tenseind], pattern.terms[tenpind], ,temp_dummy_tensor,temp_dummy)
                push!(matches, tenseind=>tenpind) # this term in the expression is related to the pattern
                break
            end

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
    rep=typeof(exp)(rep) #makes sure rep isnt a dummy if success


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
    pattern=strip(pattern)#remove spaces
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
        push!(tensor_term,copy(tensor))#adds a copy
    end

    return TensorTerm(tensor_term) #should construct a term out of our set of patterns

end




function sub_tensor(exp::TensorSuperType, pattern::TensorSuperType, rep::TensorSuperType, temp_dummy_tensor::Dict{String,TensorSuperType}, temp_dummy::Dict{String,Index})


    #if replacement has dummy indices need logic to compute that
    rep_dum = typeof(rep) <: DummyPatternTensor
    rep=typeof(exp)(rep) #makes sure rep isnt a dummy if success


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