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

x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices])
θbar = Coordinate("Θ", [spacetime_indices])

grassman_deriv=Derivative("𝛛",zero, [θ], [a => lower])
grassman_deriv_bar=Derivative("∂",zero, [θbar], [a => lower])
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
θaj=Tensor("θ",[a=>upper, j=> lower],θ)
barθaj=Tensor("Θ",[a=>upper, j=> lower],θbar)
barθbk=Tensor("Θ",[b=>upper, k=> lower],θbar)
θbk=Tensor("θ",[b=>upper, k=> lower],θ)
barθdl=Tensor("Θ",[d=>upper, l=> lower],θbar)
θdl=Tensor("θ",[d=>upper, l=> lower],θ)
barθem=Tensor("Θ",[e=>upper, m=> lower],θbar)
θem=Tensor("θ",[e=>upper, m=> lower],θ)
barθfn=Tensor("Θ",[f=>upper, n=> lower],θbar)
θfn=Tensor("θ",[f=>upper, n=> lower],θ)
barθgr=Tensor("Θ",[g=>upper, r=> lower],θbar)
θgr=Tensor("θ",[g=>upper, r=> lower],θ)
barθhs=Tensor("Θ",[h=>upper, s=> lower],θbar)
#last two indices are the spacetime ones
SUPERFIELDpp=(Tensor("j", x)
+θci*Tensor("ψ",[c=>lower, i=>upper],x)
+barθci*Tensor("̄ψ",[c=>lower, i=>upper],x)
+θci*θaj*Tensor("F",[c=>lower, a=>lower, i=>upper,j=>upper],x)
+barθci*barθaj*Tensor("̄f",[c=>lower, a=>lower, i=>upper,j=>upper],x)
+θci*barθaj*Tensor("K",[c=>lower, a=>lower, i=>upper,j=>upper],x)
+θci*θaj*barθbk*Tensor("S",[c=>lower, a=>lower, b=>lower, i=>upper,j=>upper,k=>upper],x)
+barθci*barθaj*θbk*Tensor("̄s",[c=>lower, a=>lower, b=>lower, i=>upper,j=>upper,k=>upper],x)

+θci*θaj*barθbk*barθdl*Tensor("T",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>upper,l=>upper],x)

+θci*θaj*θbk*Tensor("H",[c=>lower, a=>lower, b=>lower, i=>upper,j=>upper,k=>upper],x)
+barθci*barθaj*barθbk*Tensor("̄h",[c=>lower, a=>lower, b=>lower, i=>upper,j=>upper,k=>upper],x)

+θci*θaj*θbk*θdl*Tensor("M",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>upper,l=>upper],x)
+barθci*barθaj*barθbk*barθdl*Tensor("̄m",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>upper,l=>upper],x)

+θci*θaj*θbk*barθdl*Tensor("N",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>upper, l=>lower],x)
+barθci*barθaj*barθbk*θdl*Tensor("̄n",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>upper, l=>lower],x)
+θci*θaj*θbk*barθdl*barθem*Tensor("P",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, i=>upper,j=>upper,k=>upper, l=>lower, m=>lower],x)
+barθci*barθaj*barθbk*θdl*θem*Tensor("̄p",[c=>lower, a=>lower, b=>lower, d=>lower, e=>lower, i=>upper,j=>upper,k=>upper, l=>lower, m=>lower],x)
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
res=to_lc_expression(out,Dict(g=>minus, h=>minus))
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