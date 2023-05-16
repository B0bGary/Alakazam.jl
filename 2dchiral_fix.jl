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


a = Index("a", spacetime_indices, 1)
b = Index("b", spacetime_indices,1)
c = Index("c", spacetime_indices,1)
d = Index("d", spacetime_indices,1)
e = Index("e", spacetime_indices,1)
f = Index("f", spacetime_indices,1)
g = Index("g", spacetime_indices,1)
h = Index("h", spacetime_indices,1)
i = Index("i", spacetime_indices,1)
j = Index("j", spacetime_indices,1)
k = Index("k", spacetime_indices,1)
l = Index("l", spacetime_indices,1)


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

#last two indices are the spacetime ones
SUPERFIELDpp=(Tensor("j",[a =>lower, b => lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("S",[c => lower, a =>lower, b => lower], x)
+Tensor("Θ",[c=>upper],θbar)*Tensor("̄s",[c => lower, a =>lower, b => lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("Θ",[d=>upper],θbar)*Tensor("T",[c => lower,d => lower, a =>lower, b => lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("K",[c => lower,d => lower, a =>lower, b => lower], x)
#+Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("̄K",[c => lower,d => lower, a =>lower, b => lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
#+Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("θ",[e=>upper],θ)*Tensor("̄H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("Θ",[f=>upper],θbar)*Tensor("F",[c => lower,d => lower,e => lower, f=>lower, a =>lower, b => lower], x)
)

SUPERFIELDmmmm=(Tensor("T",[a =>lower, b => lower,k => lower, l =>lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("B",[c => lower, a =>lower, b => lower,k => lower, l =>lower], x)
+Tensor("Θ",[c=>upper],θbar)*Tensor("̄b",[c => lower, a =>lower, b => lower,k => lower, l =>lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("Θ",[d=>upper],θbar)*Tensor("C",[c => lower,d => lower, a =>lower, b => lower,k => lower, l =>lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("K",[c => lower,d => lower, a =>lower, b => lower], x)
#+Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("̄K",[c => lower,d => lower, a =>lower, b => lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
#+Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("θ",[e=>upper],θ)*Tensor("̄H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("Θ",[f=>upper],θbar)*Tensor("F",[c => lower,d => lower,e => lower, f=>lower, a =>lower, b => lower], x)
)

SUPERFIELDm=(Tensor("w",[a =>lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("F",[c => lower, a =>lower], x)
+Tensor("Θ",[c=>upper],θbar)*Tensor("̄f",[c => lower, a =>lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("Θ",[d=>upper],θbar)*Tensor("G",[c => lower,d => lower, a =>lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("K",[c => lower,d => lower, a =>lower, b => lower], x)
#+Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("̄K",[c => lower,d => lower, a =>lower, b => lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
#+Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("θ",[e=>upper],θ)*Tensor("̄H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
#+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("Θ",[f=>upper],θbar)*Tensor("F",[c => lower,d => lower,e => lower, f=>lower, a =>lower, b => lower], x)
)


D_SUPERFIELD=Derivative("̄D",SUPERFIELDmmmm, [x], [i => lower])

out=expand_derivative(D_SUPERFIELD,covar,covar_expansion)
o=out
println("leibd: ",product_rule(o))
out=apply_derivative(out)
println("Found deriv: ", out)
println("Elim delta: ")
out=eliminate_kronecker(out)
println(out)

# res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
res=to_lc_expression(out,Dict(a=>m,b=>m,k=>m,l=>m,i=>p))

println(res)

# t11=get_weighted_terms(get_weighted_terms(res,"θ",2),"Θ",1)
res=to_chiral_part(res,p)
# res2=to_lc_expression(out2,Dict(a=>m,b=>m,k=>m,l=>m,i=> j=>p))
# res2=to_lc_expression(SUPERFIELDpp,Dict(a=>p,b=>p))*to_lc_expression(SUPERFIELDmmmm,Dict(a=>m,b=>m, k=>m, l=>m))


# final=get_weighted_terms(get_weighted_terms(res2,"θ",1),"Θ",1)
# final=to_chiral_part(final,p)
# # final=to_lc_expression(sort(final),Dict{Index, Index}())
