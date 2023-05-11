include("alakazam.jl")

#define index sets
#thesis specific code to get constraints
function strip_thetas(TEx::TensorExpression)
    used_inds=union(get_dummies(TEx),get_free_indices(TEx))
    unused_flav=Vector(setdiff(flavour_indices.indices,used_inds))
    unused_spin=Vector(setdiff(flavour2_indices.indices,used_inds))
    # unused_conj=Vector(setdiff(spinor_conj_indices.indices,used_inds))
    new_exp=zero
    for summand in TEx.summands
        new_term=summand[1]
        tcounter=1
        tbcounter=1
        println("stripping summand")
        for term in summand[2].terms
            if term.name=="θ"
                println(term)
                for indexes in eachindex(term.indices)
                    if term.indices[indexes][2]==upper
                        if term.indices[indexes][1].index_set==flavour_indices
                            # println(indexes)
                            # println(2*tcounter-1)
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
                println(new_term)
            end
        end
        new_exp+=new_term
    end
    return new_exp
end

flavour_indices = IndexSet("flavour_indices", 4,upper)
flavour2_indices = IndexSet("flavour2_indices", 4,upper)


# μ=Index("μ", spacetime_indices, 1)
# ν=Index("ν", spacetime_indices, 1)
# ρ=Index("ρ", spacetime_indices, 1)
# σ=Index("σ", spacetime_indices, 1)
# τ=Index("τ", spacetime_indices, 1)
# λ=Index("λ", spacetime_indices, 1)


A = Index("A", flavour2_indices,1)
B = Index("B", flavour2_indices,1)
C = Index("C", flavour2_indices,1)
D = Index("D", flavour2_indices,1)
E = Index("E", flavour2_indices,1)
F = Index("F", flavour2_indices,1)
# g = Index("g", flavour2_indices,1)
# h = Index("h", flavour2_indices,1)
# x = Index("x", flavour2_indices,1)
# y = Index("y", flavour2_indices,1)
# z = Index("z", flavour2_indices,1)
# h = Index("h", flavour2_indices,1)

I = Index("I", flavour_indices,1)
J = Index("J", flavour_indices,1)
K = Index("K", flavour_indices,1)
L = Index("L", flavour_indices,1)
M = Index("M", flavour_indices,1)
N = Index("N", flavour_indices,1)
#  = Index("q", flavour_indices,0)
# r = Index("r", flavour_indices,0)
# s = Index("s", flavour_indices,0)
# t = Index("t", flavour_indices,0)
# u = Index("u", flavour_indices,0)
# v = Index("v", flavour_indices,0)

# x = Coordinate("x", [spacetime_indices])
θp = Coordinate("θ", [flavour_indices])
θm = Coordinate("θ", [flavour2_indices])
# θbar = Coordinate("Θ", [spacetime_indices])

grassman_derivp=Derivative("𝛛",zero, [θp], [ I =>upper])
grassman_derivm=Derivative("𝛛",zero, [θm], [ A =>upper])
# grassman_deriv_bar=Derivative("∂",zero, [θbar], [a => lower, i=>lower])
spatial_deriv=Derivative("∇",zero, [x])

# add_derivative_action(grassman_deriv_bar,Tensor("Θ",[c=>upper, j=>upper],θbar),KroneckerDelta([a =>lower, c=>upper])*KroneckerDelta([i =>lower, j=>upper]))
add_derivative_action(grassman_derivp,Tensor("θ",[ J=>lower],θp),KroneckerDelta([J =>lower, I=>upper]))
add_derivative_action(grassman_derivp,Tensor("θ",[ J=>upper],θp),Metric("δ",[J =>upper, I=>upper]))
# add_derivative_action(grassman_derivp,Tensor("θ",[ J=>lower],θ),KroneckerDelta([j =>lower, i=>upper])*EpsilonTensor([b =>upper, a=>upper]))
# add_derivative_action(grassman_derivp,Tensor("θ",[ J=>upper],θ),EpsilonTensor([j =>upper, i=>upper])*EpsilonTensor([b =>upper, a=>upper]))


# covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[b=>upper, i=>lower],θ)*spatial_deriv
covar_unbar_expansion=1*grassman_derivp#-im/2*Tensor("θ",[ i=>upper, a=>upper],θ)*spatial_deriv
covar=Derivative("D",zero, [θp], [ I => upper])
# covar_bar=Derivative("D",zero, [x,θbar], [a => lower, i=>lower])

# θgrass=Tensor("θ",[a=>upper],θ)
# θbargrass=Tensor("Θ",[b=>upper],θbar)

custom_sort_order(["θ","Θ","ϵ"])
# θci=Tensor("θ",[c=>upper, i=> lower],θ)
# barθci=Tensor("Θ",[c=>upper, i=> lower],θbar)

θI=Tensor("θ",[  I=> lower],θp)
θJ=Tensor("θ",[  J=> lower],θp)
θK=Tensor("θ",[  K=> lower],θp)
θL=Tensor("θ",[  L=> lower],θp)
# θab=Tensor("θ",[  i=> lower, a=>lower],θ)*Tensor("θ",[ i=> upper, b=>lower],θ)
# θij=Tensor("θ",[  i=> lower, a=>lower],θ)*Tensor("θ",[ j=> lower, a=>upper],θ)
# θ3=Tensor("θ",[  j=> lower, a=>upper],θ)*Tensor("θ",[  i=> lower, a=>lower],θ)*Tensor("θ",[ i=> upper, b=>lower],θ)
# θ4=Tensor("θ",[  j=> upper, b=>upper],θ)*Tensor("θ",[  j=> lower, a=>upper],θ)*Tensor("θ",[  i=> lower, a=>lower],θ)*Tensor("θ",[ i=> upper, b=>lower],θ)


#last two indices are the spacetime ones
Left=(Tensor("l")
+θI*Tensor("ψ",[   I=>upper])
+θI*θJ*EpsilonTensor([I=>upper, J=>upper])*Tensor("F")
+θI*θJ*θK*AntisymmetricTensor("ϵ",[I=>upper, J=>upper, K=>upper])*Tensor("S")
+θI*θJ*θK*θL*AntisymmetricTensor("ϵ",[I=>upper, J=>upper, K=>upper, L=>upper])*Tensor("T")
)

# Tens=(Tensor("T", x)
# +θia*Tensor("A",[   i=>upper, a=>upper],x)
# +θab*SymmetricTensor("B",[   a=>upper, b=>upper],x)
# +θij*SymmetricTensor("C",[   i=>upper, j=>upper],x)
# +θ3*Tensor("G",[b=>upper, j=>upper],x)
# +θ4*Tensor("H",x)
# )


# dScalar=Derivative("∇",Scalar, [x], [τ => lower, λ=>upper])

println(Left)


 D_SUPERFIELD=Derivative("D",Left, [θp], [M => upper])
# # println(covar)
# # println(covar_unbar_expansion)
# # println(covar_bar)
# # println(covar_expansion)
# # println(D_SUPERFIELD)
 out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
o=out
 println(o)
out=apply_derivative(out)
# println("Found deriv: ", out)
# println("Elim delta: ")
out=eliminate_kronecker(out)
println(out)
res=out
# # res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
# res=to_lc_expression(out,Dict(τ=>m, λ=>m))
# println(res)
# println("chiral part")
# res=to_chiral_part(res,p)
# println(res)




# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ r=>upper, f=>upper])
# # println(covar)
# # println(covar_unbar_expansion)
# # println(covar_bar)
# # println(covar_expansion)
# # println(D_SUPERFIELD)
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# println(o)
# out=apply_derivative(out)
# println("Found deriv: ", out)
# println("Elim delta: ")
# out=eliminate_kronecker(out)
# println(out)

# # res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
# # res=to_lc_expression(out,Dict(τ=>p,λ=>p))
# # println(res)
# # println("chiral part")
# # res=to_chiral_part(res,p)
# res=eliminate_epsilon(out)
# println(res)

# # t11=get_weighted_terms(get_weighted_terms(res,"θ",2),"Θ",1)

# # SUPERFIELDpp=(Tensor("j",[p =>lower, p => lower], x)
# # +Tensor("θ",[p=>upper],θ)*Tensor("S",[p => lower, p =>lower, p => lower], x)
# # +Tensor("Θ",[p=>upper],θbar)*Tensor("̄s",[p => lower, p =>lower, p => lower], x)
# # +Tensor("θ",[p=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("T",[p => lower,p => lower, p =>lower, p => lower], x))
# # SUPERFIELDll=(Tensor("j",[p =>upper, p => upper], x)
# # +Tensor("θ",[p=>lower],θ)*Tensor("S",[p => upper, p =>upper, p => upper], x)
# # +Tensor("Θ",[p=>lower],θbar)*Tensor("̄s",[p => upper, p =>upper, p => upper], x)
# # +Tensor("θ",[p=>lower],θ)*Tensor("Θ",[p=>lower],θbar)*Tensor("T",[p => upper,p => upper, p =>upper, p => upper], x))

# res2=Scalar×Tens


# final=get_weighted_terms(get_weighted_terms(res2,"θ",4),"Θ",0)
# # final=to_lc_expression(sort(final),Dict{Index, Index}())