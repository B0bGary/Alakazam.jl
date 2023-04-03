


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




function substitute_terms(exp::TensorTerm, pattern::TensorTerm, replace::TensorTerm, repeat::Bool=true)
    exp_hist = term2histogram(exp)
    pat_hist = term2histogram(pattern)
   # println(exp_hist)
    out = patt_match(pat_hist, exp_hist, pattern, exp) # output will be (true, [index replacements etc], dummy tensor name dict, dummy index dict)
    # replacemnts are pattern => exp index pairs
    grade = 0
    coef=1


        if out[1] #pattern matches
            # println("pattern matched")
            new_exp = TensorTerm(Vector{TensorSuperType}())
            exp2 = copy(exp.terms)
            # println("exp2")
            # println(exp2)
            #build new exp by getting all terms that need replacement
            
            println(values(out[2])) ######### may have to remove multiple things.............
            vlast=0
            for v in values(out[2]) #for each tensor that needs to be replaced 
                v = v[2]
                if vlast!=0&&vlast<v
                            #logic for grading, and removal
                    v-=1 #accounts for removal
                end
                
                # println(v)
                if v != 1
                    for i in 1:(v-1)
                        grade += grading(exp2[i])
                        grade %= 2
                    end
                    grade += grading(exp2[v]) #need to account for multiple terms being moved to the front
                end
                splice!(exp2, v) #removes the vth element, after calculating the grading to move it to the front
                vlast=v
            end

            # #construct the replacement
            # println(replace)
            rep2=copy(replace)
            for repind in eachindex(rep2.terms)
                rep = rep2.terms[repind]
                # println("replacement terms non empty")
                # println(rep)
                rep_dum = typeof(rep) <: DummyPatternTensor
                # rep = typeof(exp)(rep) #makes sure rep isnt a dummy if success


                if rep_dum
                    # println("Namema")
                    # println(exp.name)
                    println(out[3])
                    rep.name = out[3][rep.name]
                    # println(rep.name)
                end
                #needs to make sure dummy indices match pattern
                # println(rep)
                ri = rep.indices
                # println(ri)
                # println("indices")
                # println(ri)
                # println("about to replace indices with non dummies")
                for index in eachindex(ri)
                    # println("indices exist")
                    if typeof(ri[index][1]) <: DummyPatternIndex
                        # println("dummy index found")
                        # println("dummy")
                        # println(temp_dummy)
                        # println(ri[index][1])
                        if get(out[4], ri[index][1].name, nothing) === nothing #problem with checking equality since technically different objects #use names instead
                            throw(InvalidReplacementPattern())
                        else
                            println(out[4][ri[index][1].name] )
                            println( ri[index][2] )
                            ri[index] = out[4][ri[index][1].name] => ri[index][2] # maps index from original expression to replacement
                            # println("final assignment")##not happen
                            # println(ri[index])
                            # println(ri[index])
                        end
                    end

                end
                rep = typeof(exp.terms[1])(rep)
                push!(new_exp.terms, rep)
                # println(rep)
                # println(typeof(rep))
                # replace.terms[repind]=rep
            end


            out=(false, out[2],out[3],out[4])
            # println("new expressssss")
            # println(new_exp)
            # println(exp2)
           
            if repeat==true
               tex= substitute_terms(TensorTerm(exp2), pattern, replace, repeat)
            #    println("here")
               if typeof(tex)!=TensorExpression
                tex=TensorExpression(tex)
               end
            #    println( tex.summands[1][2].terms)
               new_exp.terms = vcat(new_exp.terms, tex.summands[1][2].terms)
                coef*=tex.summands[1][1]
            else
                new_exp.terms = vcat(new_exp.terms, exp2)
            end

            # println("final expression")
            # println(new_exp)
            return ((-1)^(grade))*coef *new_exp

            
            # for termind in eachindex(exp.terms)


            # end

            # #move all terms to be replaced to the front, accounting for grading

        else
            # out=(false, out[2],out[3],out[4])
            return exp #if no matches, return original expression
        end
end



function match_dummy_logic(dummy_matching::Dict{String,String}, dummy_name::String, exp_name::String)
    dummy = copy(dummy_matching)
    if get(dummy, dummy_name, nothing) == nothing #if dummy isnt already assigned a val in match
        dummy[dummy_name] = exp_name

    elseif dummy[dummy_name] != exp_name #dummy doesnt match

        return (false, dummy_matching) #dummy tensor name pattern
    end

    return (true, dummy)

end







#dynamic prog attempt
function patt_match(pat::OrderedDict{N,Vector{Pair{String,Vector{Number}}}} where {N<:Number}, exp::OrderedDict{N,Vector{Pair{String,Vector{Number}}}} where {N<:Number},
    tensorpatt::TensorTerm, tensorexp::TensorTerm, dummy_matching::Dict{String,String}=Dict{String,String}(), dummy_index_matching::Dict{String,String}=Dict{String,String}())
    #list of potential dummies is included, by default its empty
    println("exp is",exp)
    #dict are of form  [#occurances => [tensorname =>[postions],  tensorname =>[postions]], ]
    println(pat)
    println(length(pat))
    if length(pat) == 1 #if only one thing in the pattern   ---- No, this is only one type of frequency
        #println(first(pat)[2])
        println("exp is",exp)
        (pkeys, pvalues) = first(pat)
        #####all of this works!!!!!! if single tensor left in pattern!!   
        println("trying")
        if length(first(pat)[2]) == 1 #if only one type of tensor in the pattern
            


            # if first(pat)[1]==1 #if only needs to match once only one thing left to match



            #handles the case where only one this etc
            println(exp)
            if length(exp) == 1
                println("doing this for some reason")
                looper = [first(exp)]
                #println("this one used")
            else
                looper = zip(keys(exp), values(exp))
            end
            println(looper)
            for (ekey, evalue) in looper
                println(ekey, evalue)
                if pkeys[1] <= ekey #if occ in pat is less than occ in exp

                    #need to loop over evals!
                    for eval in evalue
                        #println("occurs less")
                        #println(pvalues[1][1])
                        #println(eval[1])

                        pdummy = false
                        #println(pvalues[1][1][1])
                        if pvalues[1][1][1] == '?'
                            pdummy = true
                            println("dummy found")
                        end

                        if pdummy || pvalues[1][1] == eval[1]#names match
                            dummy_temp = dummy_matching
                            dumindtemp = dummy_index_matching
                            if pdummy
                                println("dummy in patter")
                                out = match_dummy_logic(dummy_temp, pvalues[1][1], eval[1])
                                if out[1] == false
                                    continue
                                else
                                    dummy_temp = out[2]
                                    println(dummy_temp)
                                end
                            end
                            #println(dummy_matching)

                            #println("names match")
                            match = []

                            j = 1

                            #println(ekey)
                            for i in 1:ekey #match with x many terms
                                # #println(ekey)
                                #println(tensorpatt.terms[pvalues[1][2][j]].indices)
                                #println(tensorexp.terms[eval[2][i]].indices)
                                outt = check_index_match(pvalues[1][2][j], eval[2][i], tensorpatt, tensorexp, dummy_index_matching)
                                #  println(outt)
                                if outt[1]
                                    #println("indices match")
                                    push!(match, pvalues[1][2][j] => eval[2][i])
                                    j += 1
                                    if j > pkeys[1] #if matched with all
                                        #println("breaking")
                                        dumindtemp = outt[2]
                                        break
                                    end

                                else
                                    #println("indices dont match")
                                end

                            end
                            if j <= pkeys[1]
                                #println(j)
                                #println(pkeys[1])
                                if !pdummy
                                    println("ret1")
                                    return (false, [-1 => -1], dummy_matching, dummy_index_matching)
                                else
                                    continue
                                end
                            end


                            # #feed new back into method with quantity decreased by counts
                            # if patt_match(exp2)

                            #needs to loop over all n matches
                            # dummy_matching=dummy_temp
                            #          println(dummy_temp)
                            #         println(dumindtemp)
                            println("ret2 main found")
                            println(dummy_temp)
                            return (true, match, dummy_temp, dumindtemp) #found the match #found the match
                        end
                    end

                end
                continue #else, just keep going
            end
            println("ret3")
            return (false, [-1 => -1], dummy_matching, dummy_index_matching) #no match


        else
            println("multi pattern")
            ##now for multiple things left in the pattern of this frequency

            if length(exp) == 1
                println("small exp")
                looper = [first(exp)]
                #println("this one used")
            else
                println("big exp")
                looper = zip(keys(exp), values(exp))
            end
            #println(looper)
            for (ekey, evalue) in looper
                #recall, ekey=#num occurances
                #eval= vector,  [tensorname =>[postions],  tensorname =>[postions], ...]

                ########## loop through evalues?

                for eval in evalue

                    # run the loop until we find a match. 
                    # take the first item of the pattern, try and match it with this expression term. 
                    # try and match remainder of pattern with remainder of expression term

                    #freqyency with first item in the pattern
                    check = patt_match(OrderedDict(first(pat)[1] => [first(pat)[2][1]]), OrderedDict(ekey => [eval]), tensorpatt, tensorexp, dummy_matching, dummy_index_matching) #single expression term

                    if check[1] #if it's a match
                        println("found a match")
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
                        #println("do ti")
                        #println(eval)

                        #println(evalue)
                        # #println(filter!(p1=>p2 -> p1 != eval[1], evalue))

                        push!(exp2, ekey => filter!(x -> x[1] != eval[1], evalue))#remove tensor of this name
                        if ekey - pkeys[1] >= 1
                            # vnew=
                            if get(exp2, ekey - pkeys[1], nothing) == nothing
                                push!(exp2, (ekey - pkeys[1]) => [eval[1] => eval[2][pkeys:end]])
                            else
                                #println(exp2[ekey-pkeys[1]])
                                #println(eval[1])
                                push!(exp2, (ekey - pkeys[1]) => push!(exp2[ekey-pkeys[1]], eval[1] => eval[2][pkeys:end])) #decrements occ of tensor 
                            end
                        end
                        #    push!(pat2, ekey => filter!(x=>y -> x != evalue[1], eval))#remove tensor of this name
                        ######################


                        #remove freqency of patt from exp


                        #return match of the remaining. effectively, we're peeling off one term at a time. we can then add mimo by saving things that dont match
                        out = patt_match(pat2, exp2, tensorpatt, tensorexp, check[3], check[4])
                        if out[1] == true
                            println("ret4")
                            return (true, vcat(check[2], out[2]), out[3], out[4])
                        else
                            # 
                            continue
                        end

                    else
                        continue #else, look at the next term in the expression
                    end

                end

              
            end

            println("ret5")
            return (false, [-1 => -1], dummy_matching, dummy_index_matching) #no match




        end


    else  #if multiple things to match in the pattern
        #use short circuit eval

        println("multiple things in the pattern")

        outs = []

        for (pkeys, pvalues) in pat #if multiple freqencies

            out = patt_match(OrderedDict(pkeys => pvalues), exp, tensorpatt, tensorexp, dummy_matching, dummy_index_matching)


            if out[1] == false
                println("ret6")
                return (false, [-1 => -1], dummy_matching, dummy_index_matching)
            else
                if isempty(outs)
                    outs = (true, out[2], out[3], out[4])
                else
                    outs = (true, vcat(outs[2], out[2]), out[3], out[4])
                end
            end

        end

        return outs

    end
end



function check_index_match(loc_pat::Number, loc_exp::Number, pat::TensorTerm, exp::TensorTerm, dumind::Dict{String,String})
    dummy = copy(dumind)
    #println("checking indices")
    #println(typeof(pat.terms[loc_pat]))
    if typeof(pat.terms[loc_pat]) <: DummyPatternTensor
        #println("found dummy tensor")
        for ind = eachindex(pat.terms[loc_pat].indices)
            #println("index is of type")
            #println(typeof(pat.terms[loc_pat].indices[ind][1]))

            if typeof(pat.terms[loc_pat].indices[ind][1]) <: DummyPatternIndex
                if pat.terms[loc_pat].indices[ind][2] == exp.terms[loc_exp].indices[ind][2]#checks if same position
                    indexx = pat.terms[loc_pat].indices[ind][1]
                    ##println("checking dummy " * indexx.name)
                    if get(dummy, indexx.name, nothing) == nothing #if dummy isnt already assigned a val in match
                        dummy[indexx.name] = exp.terms[loc_exp].indices[ind][1].name
                        ##println("adding dummy " * indexx.name * " to dict with match:")
                        ##println(dummy[indexx.name])
                    elseif dummy[indexx.name] != exp.terms[loc_exp].indices[ind][1].name #dummy doesnt match
                        return (false, dumind) #dummy tensor name pattern
                    end
                else
                    return (false, dumind) #if different positions, then no match
                end
            else
                if pat.terms[loc_pat].indices[ind] != exp.terms[loc_exp].indices[ind]
                    return (false, dumind)
                end
            end

        end
        #       println(dummy)
        return (true, dummy)
    else
        return (pat.terms[loc_pat].indices == exp.terms[loc_exp].indices, dumind)
    end

    # if get(dummy, dummy_name, nothing) == nothing #if dummy isnt already assigned a val in match
    #     dummy[dummy_name] = exp_name

    # elseif dummy[dummy_name] != exp_name #dummy doesnt match

    #     return (false, dummy) #dummy tensor name pattern
    # end

    # return (true, dummy)

end




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