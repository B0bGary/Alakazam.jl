

function parse_pattern2tensor(pattern::String)
    it = Iterators.Stateful(pattern)
    # println(it)
    #add a way user can do custom patterns, eg \eps for eps tensor, \delta for delta func etc. some sort of vector containg possible presets (or a tree!)
    tensor = []

    #Construct the tensor
    for char in it
        # println("starting loop")
        # println(char)
        if string(char) == "?"
            # println("did this")
            s = string(iterate(it)[1])
            # println(s)
            tensor = DummyPatternTensor(s) #?A  ie the name isn't important
            # println(tensor)
        else
            # println(string(char))
            tensor = Tensor(string(char)) #A     name is important #maybe something about coord dependance
            # println(tensor)

        end
        next = peek(it)
        if next !== nothing
            char = iterate(it)[1]
            # println(char)
            # println("should be next^")
            # println(string(char))
            #now, do indices
         #   while string(char) == "^" ||string(char)=="_"
            while string(char)=="^"||string(char)=="_"
            if string(char) == "^"
                # println("made it here, next is")
                # println(string(peek(it)[1]))
                if string(peek(it)[1]) == "{"
                    iterate(it)
                    # println("made it here,ffff next is")
                    # println(string(peek(it)[1]))
                    ## do indices
                    ind = 1
                    while string(peek(it)[1]) != "}"
                        # println("while")
                        if length(tensor.indices) < ind #if not already an index
                            push!(tensor.indices, Index(string(iterate(it)[1])) => upper)
                            
                        elseif string(peek(it)[1]) != "*" #if already, and current value is not a *
                            tensor.indices[ind]= Index(string(iterate(it)[1])) => upper
                        else
                            iterate(it)
                        end
                        ind+=1
                    end
                    iterate(it)
                    
                else
                    push!(tensor.indices, Index(string(iterate(it)[1])) => upper) #single upper index
                end
            end
            # char = iterate(it)[1]
            # println("next show lowre")
            # println(char)
            if string(char) == "_"
                if string(peek(it)[1]) == "{"
                    iterate(it)
                    ## do indices
                    ind = 1
                    while string(peek(it)[1]) != "}"
                        # println(tensor.indices)
                        # println(length(tensor.indices))
                        if length(tensor.indices) < ind #if not already an index, stops out ofbounds error. add new index
                            push!(tensor.indices, Index(string(iterate(it)[1])) => lower)
                            # println("kill ma")
                        elseif string(peek(it)[1]) != "*" #if there's already an index and the next thing isnt a placeholder
                            # println("poo")
                            # println(tensor.indices)
                            tensor.indices[ind]= Index(string(iterate(it)[1])) => lower
                        else
                            # println(peek(it))
                            iterate(it) # removes * 
                        end
                        ind+=1
                    end
                    iterate(it) #remove }
                else
                    push!(tensor.indices, Index(string(iterate(it)[1])) => lower) #single upper index
                end

            end
            next = peek(it)
            # println("next is")
            # println(next)
            if next !== nothing
                char = iterate(it)[1]
            else 
                break
            end
            # println(char)
        end
            
        end
    end
    return tensor

end


function sub_tensor(exp::TensorSuperType, pattern::TensorSuperType, rep::TensorSuperType)

    if !hasIndices(pattern) && !hasIndices(exp)
        if typeof(pattern) <: DummyPatternTensor || exp.name == pattern.name
            return rep
        end
        return exp
    elseif hasIndices(pattern) && hasIndices(exp)

    end

    return exp

end