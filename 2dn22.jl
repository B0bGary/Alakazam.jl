include("alakazam.jl")


#thesis specific code to get constraints
function strip_thetas(TEx::TensorExpression)
    used_inds=union(get_dummies(TEx),get_free_indices_nopos(TEx))
    println(used_inds)
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


# function strip_thetas_lightcone(TEx::TensorExpression)
#     new_exp=zero
#     for summand in TEx.summands
#         term=1
#         for tens in summands[2].terms
#             if tens.name!=

#         end

#     end

# end


#define index sets
lightcone_indices = IndexSet("lightcone", 1,lower)
# spinor_indices = IndexSet("spinor", 2)
# spacetime_indices = IndexSet("spacetime", 4,lower)
# spinor_conj_indices = IndexSet("spinor_conj", 2,lower)
# flavour_conj_indices = IndexSet("flavour_conj", 2)


p = Index("±", lightcone_indices, 1)
m = Index("∓", lightcone_indices, 1)

x = Coordinate("x", [lightcone_indices, lightcone_indices])
θ = Coordinate("θ", [lightcone_indices])
θbar = Coordinate("Θ", [lightcone_indices])

A= Tensor("A",[p=>lower, p=>lower, p=>lower, p=>lower],x)
B= Tensor("B",[m=>upper, m=>upper],x)
exp= A*B
rename_dummies(exp)




grassman_deriv_bar=Derivative("∂",zero, [θbar], [m => lower])
spatial_deriv=Derivative("∇",zero, [x], [m => lower, m => lower])

add_derivative_action(grassman_deriv_bar,Tensor("Θ",[p=>upper],θbar),zero)
add_derivative_action(grassman_deriv_bar,Tensor("Θ",[m=>upper],θbar),TensorExpression(1))

covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[m=>upper],θ)*spatial_deriv
covar=Derivative("̄D",zero, [x,θbar], [m => lower])



grassman_deriv_barp=Derivative("∂",zero, [θbar], [p => lower])
spatial_derivp=Derivative("∇",zero, [x], [m => lower, p => lower])

add_derivative_action(grassman_deriv_barp,Tensor("Θ",[m=>upper],θbar),zero)
add_derivative_action(grassman_deriv_barp,Tensor("Θ",[p=>upper],θbar),TensorExpression(1))

covar_expansionp=-1*grassman_deriv_barp+im/2*Tensor("θ",[p=>upper],θ)*spatial_derivp
covarp=Derivative("d",zero, [x,θbar], [p => lower])

#basis tensors
# θp=Tensor("θ",[p=>upper],θ)
# θm=Tensor("θ",[m=>upper],θ)
# θbarp=Tensor("Θ",[p=>upper],θbar)
# θbarm=Tensor("Θ",[m=>upper],θbar)


# SUPERFIELDpp=(Tensor("j",[p =>lower], x)
# +θp*Tensor("ψ",[p =>lower],x)+θm*Tensor("ψ",[m =>lower],x)
# +θbarp*Tensor("χ",[p =>lower],x)+θbarm*Tensor("χ",[m =>lower],x)
# +θp*θm*Tensor("A",[p =>lower,m=>lower],x)+θbarp*θbarm*Tensor("B",[p =>lower,m=>lower],x)
# +θbarp*θp*θm*Tensor("S",[p =>lower,p=>lower,m=>lower],x)+θbarm*θp*θm*Tensor("S",[m =>lower,p=>lower,m=>lower],x)
# +θp*θbarp*θbarm*Tensor("̄S",[p =>lower,p =>lower,m=>lower],x)+θp*θbarp*θbarm*Tensor("̄S",[p =>lower,p =>lower,m=>lower],x)
# +θbarp*θp*θm*Tensor("T",[p =>lower,p=>lower,m=>lower],x)
# +θp*θbarp*θbarm*Tensor("̄T",[p =>lower,p =>lower,m=>lower],x)+θp*θbarp*θbarm*Tensor("̄T",[p =>lower,p =>lower,m=>lower],x)
# )

SUPERFIELDpp=(Tensor("j",[p =>lower, p => lower], x)
+Tensor("θ",[p=>upper],θ)*Tensor("S",[p => lower, p => lower, p => lower], x)
+Tensor("θ",[m=>upper],θ)*Tensor("S",[m => lower, p => lower, p => lower], x)

+Tensor("Θ",[p=>upper],θbar)*Tensor("̄S",[p => lower, p => lower, p => lower], x)
+Tensor("Θ",[m=>upper],θbar)*Tensor("̄S",[m => lower, p => lower, p => lower], x)

+Tensor("θ",[p=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("T",[p => lower,p => lower, p => lower, p => lower], x)
+Tensor("θ",[p=>upper],θ)*Tensor("Θ",[m=>upper],θbar)*Tensor("T",[p => lower,m => lower, p => lower, p => lower], x)
+Tensor("θ",[m=>upper],θ)*Tensor("Θ",[m=>upper],θbar)*Tensor("T",[m => lower,m => lower, p => lower, p => lower], x)
+Tensor("θ",[m=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("T",[m => lower,p => lower, p => lower, p => lower], x)

+Tensor("θ",[p=>upper],θ)*Tensor("θ",[m=>upper],θ)*Tensor("K",[p => lower,m => lower, p => lower, p => lower], x)

+Tensor("Θ",[p=>upper],θbar)*Tensor("Θ",[m=>upper],θbar)*Tensor("̄K",[p => lower,m => lower, p => lower, p => lower], x)

+Tensor("θ",[m=>upper],θ)*Tensor("θ",[p=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("H",[m => lower,p => lower,p => lower, p => lower, p => lower], x)
+Tensor("θ",[m=>upper],θ)*Tensor("θ",[p=>upper],θ)*Tensor("Θ",[m=>upper],θbar)*Tensor("H",[m => lower,p => lower,m => lower, p => lower, p => lower], x)

+Tensor("Θ",[m=>upper],θbar)*Tensor("Θ",[p=>upper],θbar)*Tensor("θ",[p=>upper],θ)*Tensor("̄H",[m => lower,p => lower,p => lower, p => lower, p => lower], x)
+Tensor("Θ",[m=>upper],θbar)*Tensor("Θ",[p=>upper],θbar)*Tensor("θ",[m=>upper],θ)*Tensor("̄H",[m=> lower,p => lower,m => lower, p => lower, p => lower], x)

+Tensor("θ",[p=>upper],θ)*Tensor("θ",[m=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("Θ",[m=>upper],θbar)*Tensor("F",[p => lower,m => lower,p => lower, m=>lower, p => lower, p => lower], x)

)

SUPERFIELDmm=(Tensor("j",[m =>lower, m => lower], x)
+Tensor("θ",[p=>upper],θ)*Tensor("S",[p => lower, m =>lower, m => lower], x)
+Tensor("θ",[m=>upper],θ)*Tensor("S",[m => lower, m =>lower, m => lower], x)

+Tensor("Θ",[p=>upper],θbar)*Tensor("̄S",[p => lower, m =>lower, m => lower], x)
+Tensor("Θ",[m=>upper],θbar)*Tensor("̄S",[m => lower, m =>lower, m => lower], x)

+Tensor("θ",[p=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("T",[p => lower,p => lower, m =>lower, m => lower], x)
+Tensor("θ",[p=>upper],θ)*Tensor("Θ",[m=>upper],θbar)*Tensor("T",[p => lower,m => lower, m =>lower, m => lower], x)
+Tensor("θ",[m=>upper],θ)*Tensor("Θ",[m=>upper],θbar)*Tensor("T",[m => lower,m => lower, m =>lower, m => lower], x)
+Tensor("θ",[m=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("T",[m => lower,p => lower, m =>lower, m => lower], x)

+Tensor("θ",[p=>upper],θ)*Tensor("θ",[m=>upper],θ)*Tensor("K",[p => lower,m => lower, m =>lower, m => lower], x)

+Tensor("Θ",[p=>upper],θbar)*Tensor("Θ",[m=>upper],θbar)*Tensor("̄K",[p => lower,m => lower, m =>lower, m => lower], x)

+Tensor("θ",[m=>upper],θ)*Tensor("θ",[p=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("H",[m => lower,p => lower,p => lower, m =>lower, m => lower], x)
+Tensor("θ",[m=>upper],θ)*Tensor("θ",[p=>upper],θ)*Tensor("Θ",[m=>upper],θbar)*Tensor("H",[m => lower,p => lower,m => lower, m =>lower, m => lower], x)

+Tensor("Θ",[m=>upper],θbar)*Tensor("Θ",[p=>upper],θbar)*Tensor("θ",[p=>upper],θ)*Tensor("̄H",[m => lower,p => lower,p => lower, m =>lower, m => lower], x)
+Tensor("Θ",[m=>upper],θbar)*Tensor("Θ",[p=>upper],θbar)*Tensor("θ",[m=>upper],θ)*Tensor("̄H",[m=> lower,p => lower,m => lower, m =>lower, m => lower], x)

+Tensor("θ",[p=>upper],θ)*Tensor("θ",[m=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("Θ",[m=>upper],θbar)*Tensor("F",[p => lower,m => lower,p => lower, m=>lower, m =>lower, m => lower], x)

)



println("our superfield expansions is: ", SUPERFIELDpp)

# #apply our derivative
println("Computing deriv: ")

D_SUPERFIELD=Derivative("̄D",SUPERFIELDpp, [x,θbar], [m => lower])

out=expand_derivative(D_SUPERFIELD,covar,covar_expansion)

out=apply_derivative(out)

# D_SUPERFIELD2=Derivative("d",SUPERFIELDmm, [x,θbar], [p => lower])
# out2=expand_derivative(D_SUPERFIELD2,covarp,covar_expansionp)
# out2=apply_derivative(out2)

res=out#+out2

println("Found deriv: ", res)
println("Elim delta: ")
res=eliminate_kronecker(res)
println(res)

t11=get_weighted_terms(get_weighted_terms(res,"θ",0),"Θ",0)
