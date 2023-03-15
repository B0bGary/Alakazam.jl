

function parse_pattern2tensor(pattern::String)
    it = Iterators.Stateful(pattern)
    # # println(it)
    #add a way user can do custom patterns, eg \eps for eps tensor, \delta for delta func etc. some sort of vector containg possible presets (or a tree!)
    tensor = []

    #Construct the tensor
    for char in it
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
    end
    return tensor

end


function sub_tensor(exp::TensorSuperType, pattern::TensorSuperType, rep::TensorSuperType)


    #if replacement has dummy indices need logic to compute that
    rep_dum = typeof(rep) <: DummyPatternTensor
    rep=typeof(exp)(rep) #makes sure rep isnt a dummy if success


    if !hasIndices(pattern) && !hasIndices(exp) # if neither have indices, both scalars so do scalar replacement
        if typeof(pattern) <: DummyPatternTensor || exp.name == pattern.name
            # println("No indices, names match or dummy")
            if rep_dum
                rep.name = exp.name
            end
            return rep
        end
        # println("No indices, names dont match")
        return exp
    elseif hasIndices(pattern) && hasIndices(exp)
        if typeof(pattern) <: DummyPatternTensor || exp.name == pattern.name #both have indices, names match or one is a dummy
            pat_ind = pattern.indices
            exp_ind = exp.indices
            if pat_ind == exp_ind #indices match, ez money  # no dummy pattern indices, dont have to do logic on replacement
                # println("Has indices, names match or dummy, and indices all match")
                if rep_dum
                    rep.name = exp.name
                end
                return rep
            else
                # # println(pat_ind)
                # # println(exp_ind)
                pi_it = Iterators.Stateful(pat_ind)
                ei_it = Iterators.Stateful(exp_ind)
                # # println(pi_it)
                # # println(ei_it)

                temp_dummy = Dict{String,Index}()

                for indp in pi_it
                    # # println(ei_it)
                    inde = iterate(ei_it)[1]
                    # println(inde[1])
                    # println(indp[1])
                    # println(typeof(indp[1]))
                    # println(typeof(inde[1]))
                    if indp[2] != inde[2]
                        # println("upper/lower dont match")
                        return exp
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
                            return exp #dummy indice pattern didnt match
                        end

                    else
                        # println("Has indices, non-dummy index doesnt match")
                        return exp #didnt match
                    end

                end
                #everything matched!
                # println("Has indices, names match or dummy, and indices all match or dummy pattern matches")

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
                # # println(rep.name)
                return rep
            end

        end
        # println("Has indices, names dont match")
        return exp

    end

    return exp

end