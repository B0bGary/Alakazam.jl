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
            tempterm=specify_indices(tempterm, Dict(dum => plus))+specify_indices(tempterm, Dict(dum => minus))
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
flavour_indices = IndexSet("flavour_indices", 2,upper)
lightcone_indices = IndexSet("lightcone", 2,lower)
± = Index("±", lightcone_indices, 1)
∓ = Index("∓", lightcone_indices, 1)
plus = Index("+", lightcone_indices, 1)
minus = Index("-", lightcone_indices, 1)

a = Index("a", spacetime_indices, 1)
b = Index("b", spacetime_indices,1)
c = Index("c", spacetime_indices,1)
d = Index("d", spacetime_indices,1)
e = Index("e", spacetime_indices,1)
f = Index("f", spacetime_indices,1)
g = Index("g", spacetime_indices,1)
h = Index("h", spacetime_indices,1)
u = Index("u", spacetime_indices,1)
v = Index("v", spacetime_indices,1)

i = Index("i", flavour_indices,0)
j = Index("j", flavour_indices,0)
k = Index("k", flavour_indices,0)
l = Index("l", flavour_indices,0)
m = Index("m", flavour_indices,0)
n = Index("n", flavour_indices,0)
p = Index("p", flavour_indices,0)
r = Index("r", flavour_indices,0)
s = Index("s", flavour_indices,0)
t = Index("t", flavour_indices,0)

x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices])
θbar = Coordinate("Θ", [spacetime_indices])

grassman_deriv=Derivative("𝛛",zero, [θ], [a => lower,i=> upper])
grassman_deriv_bar=Derivative("∂",zero, [θbar], [a => lower, i=> lower])
spatial_deriv=Derivative("∇",zero, [x], [a => lower, b => lower])

add_derivative_action(grassman_deriv_bar,Tensor("Θ",[c=>upper, j=>upper],θbar),KroneckerDelta([a =>lower, c=>upper])*KroneckerDelta([j =>upper, i=>lower]))
add_derivative_action(grassman_deriv,Tensor("θ",[c=>upper, j=>lower],θ),KroneckerDelta([a =>lower, c=>upper])*KroneckerDelta([j =>lower, i=>upper]))


covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[b=>upper, i=>lower],θ)*spatial_deriv
covar_unbar_expansion=grassman_deriv-im/2*Tensor("Θ",[b=>upper, i=>upper],θbar)*spatial_deriv
covar=Derivative("̄D",zero, [x,θbar], [a => lower,i=>upper])


# θgrass=Tensor("θ",[a=>upper],θ)
# θbargrass=Tensor("Θ",[b=>upper],θbar)

custom_sort_order(["θ","Θ","ϵ"])
θci=Tensor("θ",[c=>upper, i=> lower],θ)
barθci=Tensor("Θ",[c=>upper, i=> upper],θbar)
θaj=Tensor("θ",[a=>upper, j=> lower],θ)
barθaj=Tensor("Θ",[a=>upper, j=> upper],θbar)
barθbk=Tensor("Θ",[b=>upper, k=> upper],θbar)
θbk=Tensor("θ",[b=>upper, k=> lower],θ)
barθdl=Tensor("Θ",[d=>upper, l=> upper],θbar)
θdl=Tensor("θ",[d=>upper, l=> lower],θ)
barθem=Tensor("Θ",[e=>upper, m=> upper],θbar)
θem=Tensor("θ",[e=>upper, m=> lower],θ)
barθfn=Tensor("Θ",[f=>upper, n=> upper],θbar)
θfn=Tensor("θ",[f=>upper, n=> lower],θ)
barθgr=Tensor("Θ",[g=>upper, r=> upper],θbar)
θgr=Tensor("θ",[g=>upper, r=> lower],θ)
barθhs=Tensor("Θ",[h=>upper, s=> upper],θbar)
#last two indices are the spacetime ones
SUPERFIELDpp=(Tensor("j", x)
# +θci*Tensor("ψ",[c=>lower, i=>upper],x)
# +barθci*Tensor("̄ψ",[c=>lower, i=>lower],x)
# +θci*θaj*Tensor("F",[c=>lower, a=>lower, i=>upper,j=>upper],x)
# +barθci*barθaj*Tensor("̄f",[c=>lower, a=>lower, i=>lower,j=>lower],x)
+θci*barθaj*Tensor("K",[c=>lower, a=>lower, i=>upper,j=>lower],x)
+θci*θaj*barθbk*Tensor("S",[c=>lower, a=>lower, b=>lower, i=>upper,j=>upper,k=>lower],x)
+barθci*barθaj*θbk*Tensor("̄s",[c=>lower, a=>lower, b=>lower, i=>lower,j=>lower,k=>upper],x)

+θci*θaj*barθbk*barθdl*Tensor("T",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>lower,l=>lower],x)

+θci*θaj*θbk*Tensor("H",[c=>lower, a=>lower, b=>lower, i=>upper,j=>upper,k=>upper],x)
+barθci*barθaj*barθbk*Tensor("̄h",[c=>lower, a=>lower, b=>lower, i=>lower,j=>lower,k=>lower],x)

+θci*θaj*θbk*θdl*Tensor("M",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>upper,l=>upper],x)
+barθci*barθaj*barθbk*barθdl*Tensor("̄m",[c=>lower, a=>lower, b=>lower, d=>lower, i=>lower,j=>lower,k=>lower,l=>lower],x)

+θci*θaj*θbk*barθdl*Tensor("N",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>upper, l=>lower],x)
+barθci*barθaj*barθbk*θdl*Tensor("̄n",[c=>lower, a=>lower, b=>lower, d=>lower, i=>lower,j=>lower,k=>lower, l=>upper],x)
+θci*θaj*θbk*barθdl*barθem*Tensor("P",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, i=>upper,j=>upper,k=>upper, l=>lower, m=>lower],x)
+barθci*barθaj*barθbk*θdl*θem*Tensor("̄p",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, i=>lower,j=>lower,k=>lower, l=>upper, m=>upper],x)

+θci*θaj*θbk*barθdl*barθem*barθfn*Tensor("Q",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, f=>lower, i=>upper,j=>upper,k=>upper, l=>lower, m=>lower, n=>lower],x)

+θci*θaj*θbk*θdl*barθem*Tensor("R",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, i=>upper,j=>upper,k=>upper, l=>upper, m=>lower],x)
+barθci*barθaj*barθbk*barθdl*θem*Tensor("̄r",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, i=>lower,j=>lower,k=>lower, l=>lower, m=>upper],x)

+θci*θaj*θbk*θdl*barθem*barθfn*Tensor("W",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, f=>lower, i=>upper,j=>upper,k=>upper, l=>upper, m=>lower, n=>lower],x)
+barθci*barθaj*barθbk*barθdl*θem*θfn*Tensor("̄w",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, f=>lower, i=>lower,j=>lower,k=>lower, l=>lower, m=>upper, n=>upper],x)

+θci*θaj*θbk*θdl*barθem*barθfn*barθgr*Tensor("U",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, f=>lower, g=>lower, i=>upper,j=>upper,k=>upper, l=>upper, m=>lower, n=>lower,r=>lower],x)
+barθci*barθaj*barθbk*barθdl*θem*θfn*θgr*Tensor("̄u",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, f=>lower, g=>lower, i=>lower,j=>lower,k=>lower, l=>lower, m=>upper, n=>upper,r=>upper],x)
+θci*θaj*θbk*θdl*barθem*barθfn*barθgr*barθhs*Tensor("V",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, f=>lower, g=>lower, h=>lower, i=>upper,j=>upper,k=>upper, l=>upper, m=>lower, n=>lower,r=>lower, s=>lower],x)

# +Tensor("θ",[c=>upper],θ)*Tensor("S",[c => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("̄S",[c => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("Θ",[d=>upper],θbar)*Tensor("T",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("K",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("̄K",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("θ",[e=>upper],θ)*Tensor("̄H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("Θ",[f=>upper],θbar)*Tensor("F",[c => lower,d => lower,e => lower, f=>lower, a =>lower, b => lower], x)
)


D_SUPERFIELD=Derivative("̄D",SUPERFIELDpp, [x], [u => lower, t=>upper])

out=expand_derivative(D_SUPERFIELD,covar,covar_expansion)
o=out
out=apply_derivative(out)
println("Found deriv: ", out)
println("Elim delta: ")
out=eliminate_kronecker(out)
println(out)

# res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
res=to_lc_expression(out,Dict(u=>minus))
# println("chiral part")
# res=to_chiral_part(res,p)
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