include("alakazam.jl")


#thesis specific code to get constraints
function strip_thetas(TEx::TensorExpression)
    used_inds=union(get_dummies(TEx),get_free_indices(TEx))
    unused_flav=Vector(setdiff(flavour_indices.indices,used_inds))
    unused_spin=Vector(setdiff(spinor_indices.indices,used_inds))
    unused_conj=Vector(setdiff(spinor_conj_indices.indices,used_inds))
    new_exp=zero
    for summand in TEx.summands
        new_term=summand[1]
        tcounter=1
        tbcounter=1
        for term in summand[2].terms
            if term.name=="θ"
                # println(term)
                for indexes in eachindex(term.indices)
                    if term.indices[indexes][2]==upper
                        if term.indices[indexes][1].index_set==flavour_indices
                            dl=KroneckerDelta([term.indices[indexes], unused_flav[2*tcounter-1]=>lower])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        else
                            dl=KroneckerDelta([term.indices[indexes], unused_spin[tcounter]=>lower])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        end
                    else
                        if term.indices[indexes][1].index_set==flavour_indices
                            dl=KroneckerDelta([term.indices[indexes], unused_flav[2*tcounter-1]=>upper])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        else
                            dl=KroneckerDelta([term.indices[indexes], unused_spin[tcounter]=>upper])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        end
                    end
                    
                end
                # println(new_term)
                tcounter+=1
            elseif term.name=="Θ"
                # println(term)
                for indexes in eachindex(term.indices)
                    if term.indices[indexes][2]==upper
                        if term.indices[indexes][1].index_set==flavour_indices
                            dl=KroneckerDelta([term.indices[indexes], unused_flav[2*tbcounter]=>lower])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        else
                            dl=KroneckerDelta([term.indices[indexes], unused_conj[tbcounter]=>lower])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        end
                    else
                        if term.indices[indexes][1].index_set==flavour_indices
                            dl=KroneckerDelta([term.indices[indexes], unused_flav[2*tbcounter]=>upper])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        else
                            dl=KroneckerDelta([term.indices[indexes], unused_conj[tbcounter]=>upper])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        end
                    end
                    
                end
                tbcounter+=1
                # println(new_term)
            else
                new_term*=term
            end
        end
        new_exp+=new_term
    end
    return new_exp
end

function specify_indices(TEx::TensorExpression, reps::Dict{Index,Index})
    new_exp=zero
    w=copy(TEx)
    for summand in w.summands
        term=1
        for tens in summand[2].terms
            if typeof(tens)<:OperatorSuperType
                tens.operand=specify_indices(tens.operand,reps)
                tens.operand_indices=tens.operand.summands[1][2].indices
                # println(reps)
                # test=tens.operator_indices[1][1]
                # println(test )
                # println(haskey(reps,test))
                tens.operator_indices=[ haskey(reps,ind[1]) ? reps[ind[1]]=>ind[2] : ind[1]=>ind[2]  for ind in tens.operator_indices]
                # println("op should inclde", tens.operator_indices)
                tens.indices=vcat(tens.operator_indices,tens.operand_indices)
                # println("op has", tens.indices)
            else
                if !isempty(tens.indices)
                tens.indices=[ haskey(reps,ind[1]) ? reps[ind[1]]=>ind[2] : ind[1]=>ind[2]  for ind in tens.indices]
                end
            end
            if !(typeof(term)<:SymorNum) && ((typeof(term)<:TensorTerm && tens in term.terms) || (typeof(term)<:TensorExpression && tens in term.summands[1][2].terms) )
                term=zero
                break
            else
                term=term*tens
            end
        end
        new_exp+=term
    end

    return new_exp
end

function to_lc_expression(TEx::TensorExpression, reps::Dict{Index,Index})
    #replace non-dummies
    exp=specify_indices(TEx, reps)
    # println(exp)
    #cookery for dummies....
    new_exp=zero

    for summand in exp.summands
        dums=intersect(get_dummies(summand[2]),Set(spacetime_indices.indices))
        tempterm=TensorExpression(summand[2])
        # println("reset",tempterm)
        while !isempty(dums)
            dum=pop!(dums)
            tempterm=specify_indices(tempterm, Dict(dum => p))+specify_indices(tempterm, Dict(dum => m))
            # println(tempterm)
        end

        new_exp+=summand[1]*tempterm
    end
    res_ex=zero
    for summand in new_exp.summands
        flag=0
            for tens in summand[2].terms
                if typeof(tens)<:OperatorSuperType
                    if tens.operator_indices[1][1].name!=tens.operator_indices[2][1].name
                        flag=1
                     #   println("ignoring term", summand[2])
                        break
                    end
                end
            end
            if flag==0 &&term_weight(summand[2], "θ")<=2 && term_weight(summand[2], "Θ")<=2
                res_ex+=summand[1]*summand[2]
            end
    end

    # end
    out_exp=zero
    flag=0
    for summand in new_exp.summands
        for tens in summand[2].terms
            if typeof(tens)<:OperatorSuperType
                if tens.operator_indices[1][1].name!=tens.operator_indices[2][1].name
                    
                    flag=1
                    break
                end
            end
        end
        if flag==0
            # println("got here, ",summand[2] )
        out_exp+=summand[1]*summand[2]
        end
        flag=0
    end



    return out_exp
end

function to_chiral_part(TEx::TensorExpression, reps::Index)
    new_exp=zero

    for summand in TEx.summands
        term=1

        for tens in summand[2].terms
            if tens.name=="θ" || tens.name=="Θ"
                if tens.indices[1][1].name!=reps.name
                    term=zero
                    break
                end
            end
            term*=tens
        end

        new_exp+=term
    end

    return new_exp
end

#define index sets

spacetime_indices = IndexSet("spacetimelc", 2,lower)
flavour_indices = IndexSet("flavour_indices", 4,upper)
lightcone_indices = IndexSet("lightcone", 2,lower)
± = Index("±", lightcone_indices, 1)
∓ = Index("∓", lightcone_indices, 1)
p = Index("+", lightcone_indices, 1)
m = Index("-", lightcone_indices, 1)

a = Index("a", spacetime_indices, 1)
b = Index("b", spacetime_indices,1)
c = Index("c", spacetime_indices,1)
d = Index("d", spacetime_indices,1)
e = Index("e", spacetime_indices,1)
f = Index("f", spacetime_indices,1)
g = Index("g", spacetime_indices,1)
h = Index("h", spacetime_indices,1)

i = Index("i", flavour_indices,0)
j = Index("j", flavour_indices,0)
k = Index("k", flavour_indices,0)
l = Index("l", flavour_indices,0)


x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices])
θbar = Coordinate("Θ", [spacetime_indices])

grassman_deriv=Derivative("𝛛",zero, [θ], [i => upper])
grassman_deriv_bar=Derivative("∂",zero, [θbar], [a => upper])
spatial_deriv=Derivative("∇",zero, [x], [a => lower, b => lower])

add_derivative_action(grassman_deriv_bar,Tensor("Θ",[c=>upper],θbar),KroneckerDelta([a =>lower, c=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[c=>upper],θ),KroneckerDelta([a =>lower, c=>upper]))


covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[b=>upper],θ)*spatial_deriv
covar_unbar_expansion=grassman_deriv-im/2*Tensor("Θ",[b=>upper],θbar)*spatial_deriv
covar=Derivative("̄D",zero, [x,θbar], [a => lower])


θgrass=Tensor("θ",[a=>upper],θ)
θbargrass=Tensor("Θ",[b=>upper],θbar)

custom_sort_order(["θ","Θ","ϵ"])
θci=Tensor("θ",[c=>upper, i=> lower],θ)
barθci=Tensor("Θ",[c=>upper, i=> lower],θbar)
#last two indices are the spacetime ones
SUPERFIELDpp=(Tensor("j", x)
+Tensor("θ",[i=> lower],θ)*Tensor("ψ",[c=>lower, i=>upper],x)
+Tensor("θ",[ i=> lower],θ)*Tensor("θ",[j=> lower],θ)*Tensor("F",[ i=>upper,j=>upper],x)
+Tensor("θ",[i=> lower],θ)*Tensor("θ",[j=> lower],θ)*Tensor("θ",[k=> lower],θ)*Tensor("S",[ i=>upper,j=>upper,k=>upper],x)
+Tensor("θ",[i=> lower],θ)*Tensor("θ",[j=> lower],θ)*Tensor("θ",[k=> lower],θ)*Tensor("θ",[ l=> lower],θ)*Tensor("T",[i=>upper,j=>upper,k=>upper,l=>upper],x)

# +Tensor("θ",[c=>upper],θ)*Tensor("S",[c => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("̄S",[c => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("Θ",[d=>upper],θbar)*Tensor("T",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("K",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("̄K",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("θ",[e=>upper],θ)*Tensor("̄H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("Θ",[f=>upper],θbar)*Tensor("F",[c => lower,d => lower,e => lower, f=>lower, a =>lower, b => lower], x)
)


D_SUPERFIELD=Derivative("̄D",SUPERFIELDpp, [x], [g => lower, h=>lower])

out=expand_derivative(D_SUPERFIELD,spatial_deriv,1*spatial_deriv)
o=out
out=apply_derivative(out)
println("Found deriv: ", out)
println("Elim delta: ")
out=eliminate_kronecker(out)
println(out)

# res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
res=to_lc_expression(out,Dict(g=>m, h=>m))
println("chiral part")
res=to_chiral_part(res,p)
println(res)

# t11=get_weighted_terms(get_weighted_terms(res,"θ",2),"Θ",1)

# SUPERFIELDpp=(Tensor("j",[p =>lower, p => lower], x)
# +Tensor("θ",[p=>upper],θ)*Tensor("S",[p => lower, p =>lower, p => lower], x)
# +Tensor("Θ",[p=>upper],θbar)*Tensor("̄s",[p => lower, p =>lower, p => lower], x)
# +Tensor("θ",[p=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("T",[p => lower,p => lower, p =>lower, p => lower], x))
# SUPERFIELDll=(Tensor("j",[p =>upper, p => upper], x)
# +Tensor("θ",[p=>lower],θ)*Tensor("S",[p => upper, p =>upper, p => upper], x)
# +Tensor("Θ",[p=>lower],θbar)*Tensor("̄s",[p => upper, p =>upper, p => upper], x)
# +Tensor("θ",[p=>lower],θ)*Tensor("Θ",[p=>lower],θbar)*Tensor("T",[p => upper,p => upper, p =>upper, p => upper], x))

# res2=SUPERFIELDpp*SUPERFIELDpp


# final=get_weighted_terms(get_weighted_terms(res2,"θ",1),"Θ",1)
# final=to_lc_expression(sort(final),Dict{Index, Index}())