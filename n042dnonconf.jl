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

spacetime_indices = IndexSet("spacetimelc", 1,lower)
flavour_indices = IndexSet("flavour_indices", 2,upper)
flavour2_indices = IndexSet("flavour2_indices", 2,upper)
lightcone_indices = IndexSet("lightcone", 2,lower)
± = Index("±", lightcone_indices, 1)
∓ = Index("∓", lightcone_indices, 1)
p = Index("+", lightcone_indices, 1)
m = Index("-", lightcone_indices, 1)

μ=Index("μ", spacetime_indices, 1)
ν=Index("ν", spacetime_indices, 1)
ρ=Index("ρ", spacetime_indices, 1)
σ=Index("σ", spacetime_indices, 1)
τ=Index("τ", spacetime_indices, 1)
λ=Index("λ", spacetime_indices, 1)


a = Index("a", flavour2_indices,1)
b = Index("b", flavour2_indices,1)
c = Index("c", flavour2_indices,1)
d = Index("d", flavour2_indices,1)
e = Index("e", flavour2_indices,1)
f = Index("f", flavour2_indices,1)
g = Index("g", flavour2_indices,1)
h = Index("h", flavour2_indices,1)
X = Index("x", flavour2_indices,1)
y = Index("y", flavour2_indices,1)
z = Index("z", flavour2_indices,1)
# h = Index("h", flavour2_indices,1)

i = Index("i", flavour_indices,0)
j = Index("j", flavour_indices,0)
k = Index("k", flavour_indices,0)
l = Index("l", flavour_indices,0)
n = Index("n", flavour_indices,0)
q = Index("q", flavour_indices,0)
r = Index("r", flavour_indices,0)
s = Index("s", flavour_indices,0)
t = Index("t", flavour_indices,0)
u = Index("u", flavour_indices,0)
v = Index("v", flavour_indices,0)
A = Index("A", flavour_indices,0)
B = Index("B", flavour_indices,0)
C = Index("C", flavour_indices,0)

x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices, flavour_indices, flavour2_indices])
# θbar = Coordinate("Θ", [spacetime_indices])

grassman_deriv=Derivative("𝛛",zero, [θ], [ i => upper, a =>upper])
# grassman_deriv_bar=Derivative("∂",zero, [θbar], [a => lower, i=>lower])
spatial_deriv=Derivative("∇",zero, [x])

# add_derivative_action(grassman_deriv_bar,Tensor("Θ",[c=>upper, j=>upper],θbar),KroneckerDelta([a =>lower, c=>upper])*KroneckerDelta([i =>lower, j=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[ j=>lower, b=>lower],θ),KroneckerDelta([j =>lower, i=>upper])*KroneckerDelta([b =>lower, a=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[ j=>upper, b=>lower],θ),EpsilonTensor([j =>upper, i=>upper])*KroneckerDelta([b =>lower, a=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[ j=>lower, b=>upper],θ),KroneckerDelta([j =>lower, i=>upper])*EpsilonTensor([b =>upper, a=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[ j=>upper, b=>upper],θ),EpsilonTensor([j =>upper, i=>upper])*EpsilonTensor([b =>upper, a=>upper]))


# covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[b=>upper, i=>lower],θ)*spatial_deriv
covar_unbar_expansion=1*grassman_deriv-im/2*Tensor("θ",[ i=>upper, a=>upper],θ)*spatial_deriv
covar=Derivative("D",zero, [x,θ], [ i => upper, a=>upper])
# covar_bar=Derivative("D",zero, [x,θbar], [a => lower, i=>lower])

# θgrass=Tensor("θ",[a=>upper],θ)
# θbargrass=Tensor("Θ",[b=>upper],θbar)

custom_sort_order(["θ","Θ","ϵ"])
# θci=Tensor("θ",[c=>upper, i=> lower],θ)
# barθci=Tensor("Θ",[c=>upper, i=> lower],θbar)

θia=Tensor("θ",[  i=> lower, a=>lower],θ)
θab=Tensor("θ",[  i=> lower, a=>lower],θ)*Tensor("θ",[ i=> upper, b=>lower],θ)
θij=Tensor("θ",[  i=> lower, a=>lower],θ)*Tensor("θ",[ j=> lower, a=>upper],θ)
θ3=Tensor("θ",[  j=> lower, a=>upper],θ)*Tensor("θ",[  i=> lower, a=>lower],θ)*Tensor("θ",[ i=> upper, b=>lower],θ)
θ4=Tensor("θ",[  j=> upper, b=>upper],θ)*Tensor("θ",[  j=> lower, a=>upper],θ)*Tensor("θ",[  i=> lower, a=>lower],θ)*Tensor("θ",[ i=> upper, b=>lower],θ)


#last two indices are the spacetime ones
Scalar=(Tensor("j", x)
+θia*Tensor("ψ",[   i=>upper, a=>upper],x)
+θab*Tensor("F",[   a=>upper, b=>upper],x)
+θij*Tensor("K",[   i=>upper, j=>upper],x)
+θ3*Tensor("S",[  b=>upper, j=>upper],x)
+θ4*Tensor("Tp",x)
)

Tens=(Tensor("Tm", x)
+θia*Tensor("A",[   i=>upper, a=>upper],x)
+θab*SymmetricTensor("B",[   a=>upper, b=>upper],x)
+θij*SymmetricTensor("C",[   i=>upper, j=>upper],x)
+θ3*Tensor("G",[b=>upper, j=>upper],x)
+θ4*Tensor("H",x)
)

W=-2*im*(Tensor("dT", x)
+θia*Tensor("dA",[   i=>upper, a=>upper],x)
# +θab*SymmetricTensor("B",[   a=>upper, b=>upper],x)
# +θij*SymmetricTensor("C",[   i=>upper, j=>upper],x)
# +θ3*Tensor("G",[b=>upper, j=>upper],x)
# +θ4*Tensor("H",x)
)

# Kc=(Tensor("J", [   k=>upper, c=>upper],x)
# +θia*Tensor("M",[   k=>upper, c=>upper, i =>upper, a=>upper],x)
# +θab*Tensor("N",[   k=>upper, c=>upper,   a=>upper, b=>upper],x)
# +θij*Tensor("P",[   k=>upper, c=>upper,   i=>upper, j=>upper],x)
# +θ3*Tensor("Q",[   k=>upper, c=>upper,b=>upper, j=>upper],x)
# +θ4*Tensor("R",[   k=>upper, c=>upper],x)
# )
# Kd=(Tensor("J", [   k=>upper, d=>upper],x)
# +θia*Tensor("M",[   k=>upper, d=>upper, i =>upper, a=>upper],x)
# +θab*Tensor("N",[   k=>upper, d=>upper,   a=>upper, b=>upper],x)
# +θij*Tensor("P",[   k=>upper, d=>upper,   i=>upper, j=>upper],x)
# +θ3*Tensor("Q",[   k=>upper, d=>upper,b=>upper, j=>upper],x)
# +θ4*Tensor("R",[   k=>upper, d=>upper],x)
# )


# dScalar=Derivative("∇",Scalar, [x], [τ => lower, λ=>upper])




# # D_SUPERFIELD=Derivative("D",SUPERFIELDpp, [x,θ], [g => lower, n=>upper])
# # println(covar)
# # println(covar_unbar_expansion)
# # println(covar_bar)
# # println(covar_expansion)
# # println(D_SUPERFIELD)
# # out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# # o=out
# # println(o)
# out=apply_derivative(dScalar)
# println("Found deriv: ", out)
# println("Elim delta: ")
# out=eliminate_kronecker(out)
# println(out)

# # res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
# res=to_lc_expression(out,Dict(τ=>m, λ=>m))
# println(res)
# println("chiral part")
# res=to_chiral_part(res,p)
# println(res)


#########################
#### tensor constraints

# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res1=eliminate_epsilon(out)

# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ r=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ s=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res2=eliminate_epsilon(out)

# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ s=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ r=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res3=eliminate_epsilon(out)

# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res4=eliminate_epsilon(out)


# res=1/2*(res1+res2+res3+res4)
# res=simplify(res)
# println(res)
#algreba^


#w_--

D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ r=>upper, f=>upper])
out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
o=out
out=apply_derivative(out)
out=eliminate_kronecker(out)
res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res1=eliminate_epsilon(out)

# out=EpsilonTensor([g=>lower, f=>lower])*(EpsilonTensor([s=>lower, r=>lower])*res1)
# res=eliminate_epsilon(out)


######### K constraints

# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res1=eliminate_epsilon(out)

# out=EpsilonTensor([g=>lower, f=>lower])*(EpsilonTensor([s=>lower, r=>lower])*res1)
# res=eliminate_epsilon(out)


# # D_SUPERFIELD=Derivative("D",Kd, [x,θ], [ r=>upper, c=>upper])
# # out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# # o=out
# # out=apply_derivative(out)
# # out=eliminate_kronecker(out)
# # res2=eliminate_epsilon(out)

# # res=1/2*(res1+res2)
# res=simplify(res)


#######SC2

# res=Scalar×Tens-W×W
# res2=simplify(res)
# println(res2)
# final=get_weighted_terms(get_weighted_terms(res2,"θ",4),"Θ",0)
# final=to_lc_expression(sort(final),Dict{Index, Index}())





# D_SUPERFIELD=Derivative("D",SUPERFIELDpp, [x,θ], [g => lower, n=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# # res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
# res=to_lc_expression(out,Dict(g=>p))
# res=to_chiral_part(res,p)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [h => lower, q=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# # res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
# res=to_lc_expression(out,Dict(h=>p))
# res=to_chiral_part(res,p)

# # t11=get_weighted_terms(get_weighted_terms(res,"θ",2),"Θ",1)

# # SUPERFIELDpp=(Tensor("j",[p =>lower, p => lower], x)
# # +Tensor("θ",[p=>upper],θ)*Tensor("S",[p => lower, p =>lower, p => lower], x)
# # +Tensor("Θ",[p=>upper],θbar)*Tensor("̄s",[p => lower, p =>lower, p => lower], x)
# # +Tensor("θ",[p=>upper],θ)*Tensor("Θ",[p=>upper],θbar)*Tensor("T",[p => lower,p => lower, p =>lower, p => lower], x))
# # SUPERFIELDll=(Tensor("j",[p =>upper, p => upper], x)
# # +Tensor("θ",[p=>lower],θ)*Tensor("S",[p => upper, p =>upper, p => upper], x)
# # +Tensor("Θ",[p=>lower],θbar)*Tensor("̄s",[p => upper, p =>upper, p => upper], x)
# # +Tensor("θ",[p=>lower],θ)*Tensor("Θ",[p=>lower],θbar)*Tensor("T",[p => upper,p => upper, p =>upper, p => upper], x))

# # res2=SUPERFIELDpp*SUPERFIELDpp


# # final=get_weighted_terms(get_weighted_terms(res2,"θ",1),"Θ",1)
# # final=to_lc_expression(sort(final),Dict{Index, Index}())