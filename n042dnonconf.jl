include("alakazam.jl")
include("superspace_old.jl")
#define index sets
#thesis specific code to get constraints
function strip_thetas(TEx::TensorExpression)
    used_inds=union(get_dummies(TEx),get_free_indices_nopos(TEx))
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
H = Index("H", flavour2_indices,1)
X2 = Index("X", flavour2_indices,1)
Y = Index("Y", flavour2_indices,1)
Z = Index("Z", flavour2_indices,1)
# # h = Index("h", flavour2_indices,1)

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
# D = Index("D", flavour_indices,0)
# E = Index("E", flavour_indices,0)
# F = Index("F", flavour_indices,0)
# G = Index("G", flavour_indices,0)

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
covar_unbar_expansion=1*grassman_deriv+im/2*Tensor("θ",[ i=>upper, a=>upper],θ)*spatial_deriv
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
+θab*SymmetricTensor("F",[   a=>upper, b=>upper],x)
+θij*SymmetricTensor("K",[   i=>upper, j=>upper],x)
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

Lmm=(Tensor("Λ", x)
+θia*Tensor("Ω",[   i=>upper, a=>upper],x)
+θab*SymmetricTensor("Γ",[   a=>upper, b=>upper],x)
+θij*SymmetricTensor("χ",[   i=>upper, j=>upper],x)
+θ3*Tensor("Ψ",[b=>upper, j=>upper],x)
+θ4*Tensor("Φ",x)
)

# Lmm=apply_derivative(Derivative("delmm",Scalar,[x]))




Wmmm=(Tensor("w",[   k=>upper, c=>upper], x)
+θia*Tensor("X",[  k=>upper, c=>upper,  i=>upper, a=>upper],x)
+θab*Tensor("Y",[ k=>upper, c=>upper,     a=>upper, b=>upper],x)
+θij*Tensor("Z",[ k=>upper, c=>upper,     i=>upper, j=>upper],x)
+θ3*Tensor("U",[ k=>upper, c=>upper,  b=>upper, j=>upper],x)
+θ4*Tensor("V",[ k=>upper, c=>upper],x)
)


Wm=(Tensor("J",[   l=>upper, d=>upper], x)
+θia*Tensor("M",[ l=>upper, d=>upper,  i=>upper, a=>upper],x)
+θab*Tensor("N",[ l=>upper, d=>upper,     a=>upper, b=>upper],x)
+θij*Tensor("P",[ l=>upper, d=>upper,     i=>upper, j=>upper],x)
+θ3*Tensor("Q",[ l=>upper, d=>upper,  b=>upper, j=>upper],x)
+θ4*Tensor("R",[ l=>upper, d=>upper],x)
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




##### d tensor
# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ l=>upper, d=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res=merge_sum(simplify(eliminate_epsilon(out)))






##############################
##### d tensor
# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ l=>upper, d=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# res=merge_sum(simplify(eliminate_epsilon(out)))
# # res=merge_sum(simplify(apply_derivative(Derivative("𝜵",res, [x]))))
# res2=merge_sum(simplify(apply_derivative(Derivative("delmm",Wm, [x]))))


# D_SUPERFIELD=Derivative("D",Scalar, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=drop_weight_above(out,"θ",4)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=drop_weight_above(out,"θ",4)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ t=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=drop_weight_above(out,"θ",4)

# # out=EpsilonTensor([g=>lower, f=>lower])*(EpsilonTensor([X=>lower, h=>lower])*(EpsilonTensor([u=>lower, s=>lower])*(EpsilonTensor([t=>lower, r=>lower])*out)))
# out=((EpsilonTensor([t=>lower, s=>lower])*(EpsilonTensor([g=>lower, f=>lower])*out)))
# # println(get_free_indices(out))
# # println(out)
# out=simplify(eliminate_epsilon(out))/9
# # println(out)
# # res2=get_weighted_terms(out,"θ",0)
# res2=merge_sum(simplify(apply_derivative(Derivative("𝛻",out, [x]))))
# # # println(res2)
# # # println("4 deriv")
# # # println(res2)
# # # res=res2-res















#########################
#### tensor constraints

# ob=Tens


# D_SUPERFIELD=Derivative("D",ob, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res1=eliminate_epsilon(out)

# D_SUPERFIELD=Derivative("D",ob, [x,θ], [ r=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ s=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res2=eliminate_epsilon(out)

# D_SUPERFIELD=Derivative("D",ob, [x,θ], [ s=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ r=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res3=eliminate_epsilon(out)

# D_SUPERFIELD=Derivative("D",ob, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
# D_SUPERFIELD2=Derivative("D",res, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res4=eliminate_epsilon(out)


# res=1/2*(res1+res2+res3+res4)
# res=simplify(res)
# println(res)

# #an extra:

# D_SUPERFIELDf=Derivative("D",res, [x,θ], [ v=>upper, z=>upper])
# out=expand_derivative(D_SUPERFIELDf,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=drop_weight_above(out,"θ",4)
# out=eliminate_kronecker(out)
# res=eliminate_epsilon(out)
#algreba^




################## W---


# D_SUPERFIELD=Derivative("D",Wmmm, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=EpsilonTensor([r=>lower, k=>lower])*(EpsilonTensor([f=>lower, c=>lower])*out)
# res=simplify(eliminate_epsilon(out))






##########w_-

#free inds l,d


# Wm2=(Tensor("J",[   r=>upper, f=>upper], x)
# +θia*Tensor("M",[ r=>upper, f=>upper,  i=>upper, a=>upper],x)
# +θab*Tensor("N",[ r=>upper, f=>upper,     a=>upper, b=>upper],x)
# +θij*Tensor("P",[ r=>upper, f=>upper,     i=>upper, j=>upper],x)
# +θ3*Tensor("Q",[ r=>upper, f=>upper,  b=>upper, j=>upper],x)
# +θ4*Tensor("R",[ r=>upper, f=>upper],x)
# )

# Wm3=(Tensor("J",[   r=>upper, d=>upper], x)
# +θia*Tensor("M",[ r=>upper, d=>upper,  i=>upper, a=>upper],x)
# +θab*Tensor("N",[ r=>upper, d=>upper,     a=>upper, b=>upper],x)
# +θij*Tensor("P",[ r=>upper, d=>upper,     i=>upper, j=>upper],x)
# +θ3*Tensor("Q",[ r=>upper, d=>upper,  b=>upper, j=>upper],x)
# +θ4*Tensor("R",[ r=>upper, d=>upper],x)
# )


# Wm4=(Tensor("J",[   l=>upper, f=>upper], x)
# +θia*Tensor("M",[ l=>upper, f=>upper,  i=>upper, a=>upper],x)
# +θab*Tensor("N",[ l=>upper, f=>upper,     a=>upper, b=>upper],x)
# +θij*Tensor("P",[ l=>upper, f=>upper,     i=>upper, j=>upper],x)
# +θ3*Tensor("Q",[ l=>upper, f=>upper,  b=>upper, j=>upper],x)
# +θ4*Tensor("R",[ l=>upper, f=>upper],x)
# )

# D_SUPERFIELD=Derivative("D",Wm, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)


# D_SUPERFIELD2=Derivative("D",Wm2, [x,θ], [ l=>upper, d=>upper])
# out2=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out2
# out2=apply_derivative(out2)
# out2=eliminate_kronecker(out2)


# D_SUPERFIELD3=Derivative("D",Wm3, [x,θ], [ l=>upper, f=>upper])
# out3=expand_derivative(D_SUPERFIELD3,covar,covar_unbar_expansion)
# out3=apply_derivative(out3)
# out3=eliminate_kronecker(out3)

# D_SUPERFIELD4=Derivative("D",Wm4, [x,θ], [ r=>upper, d=>upper])
# out4=expand_derivative(D_SUPERFIELD4,covar,covar_unbar_expansion)
# out4=apply_derivative(out4)
# out4=eliminate_kronecker(out4)



# res=simplify(out+out2+out3+out4)/4




# # res=simplify(out2+out3)/2



##### resitriction


# # Wm2=(Tensor("J",[   r=>upper, f=>upper], x)
# # +θia*Tensor("M",[ r=>upper, f=>upper,  i=>upper, a=>upper],x)
# # +θab*Tensor("N",[ r=>upper, f=>upper,     a=>upper, b=>upper],x)
# # +θij*Tensor("P",[ r=>upper, f=>upper,     i=>upper, j=>upper],x)
# # +θ3*Tensor("Q",[ r=>upper, f=>upper,  b=>upper, j=>upper],x)
# # +θ4*Tensor("R",[ r=>upper, f=>upper],x)
# # )

# Wm3=(Tensor("J",[   r=>upper, d=>upper], x)
# +θia*Tensor("M",[ r=>upper, d=>upper,  i=>upper, a=>upper],x)
# +θab*Tensor("N",[ r=>upper, d=>upper,     a=>upper, b=>upper],x)
# +θij*Tensor("P",[ r=>upper, d=>upper,     i=>upper, j=>upper],x)
# +θ3*Tensor("Q",[ r=>upper, d=>upper,  b=>upper, j=>upper],x)
# +θ4*Tensor("R",[ r=>upper, d=>upper],x)
# )


# Wm4=(Tensor("J",[   l=>upper, f=>upper], x)
# +θia*Tensor("M",[ l=>upper, f=>upper,  i=>upper, a=>upper],x)
# +θab*Tensor("N",[ l=>upper, f=>upper,     a=>upper, b=>upper],x)
# +θij*Tensor("P",[ l=>upper, f=>upper,     i=>upper, j=>upper],x)
# +θ3*Tensor("Q",[ l=>upper, f=>upper,  b=>upper, j=>upper],x)
# +θ4*Tensor("R",[ l=>upper, f=>upper],x)
# )



# D_SUPERFIELD=Derivative("D",Wm, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)




# # D_SUPERFIELD2=Derivative("D",Wm4, [x,θ], [ r=>upper, d=>upper])
# # out2=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# # o=out2
# # out2=apply_derivative(out2)
# # out2=eliminate_kronecker(out2)


# D_SUPERFIELD3=Derivative("D",Wm3, [x,θ], [ l=>upper, f=>upper])
# out3=expand_derivative(D_SUPERFIELD3,covar,covar_unbar_expansion)
# out3=apply_derivative(out3)
# out3=eliminate_kronecker(out3)

# D_SUPERFIELD4=Derivative("D",Wm4, [x,θ], [ r=>upper, d=>upper])
# out4=expand_derivative(D_SUPERFIELD4,covar,covar_unbar_expansion)
# out4=apply_derivative(out4)
# out4=eliminate_kronecker(out4)



# res=simplify(out+out2)/2

# res=simplify(out+out3)/2


# res=simplify(out2+out3)/2

# res=simplify(out+out4)/2


##################
##### D4W---=-W-


# D_SUPERFIELD=Derivative("D",Wmmm, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ t=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ u=>upper, X=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)


# out=EpsilonTensor([g=>lower, f=>lower])*(EpsilonTensor([X=>lower, h=>lower])*(EpsilonTensor([u=>lower, s=>lower])*(EpsilonTensor([t=>lower, r=>lower])*out)))
# res2=simplify(eliminate_epsilon(out))/144
# println("4 deriv")
# println(res2)
# # res=res2-res




##### D (DW)---=-W-


###Wmmm=(Tensor("w",[   k=>upper, c=>upper], x)

####Wm=(Tensor("J",[   l=>upper, d=>upper], x)

Wmmm2=(Tensor("w",[   l=>upper, d=>upper], x)
+θia*Tensor("X",[  l=>upper,d=>upper,  i=>upper, a=>upper],x)
+θab*Tensor("Y",[ l=>upper,d=>upper,     a=>upper, b=>upper],x)
+θij*Tensor("Z",[ l=>upper,d=>upper,     i=>upper, j=>upper],x)
+θ3*Tensor("U",[ l=>upper,d=>upper,  b=>upper, j=>upper],x)
+θ4*Tensor("V",[ l=>upper,d=>upper],x)
)

D_SUPERFIELD=Derivative("D",Wmmm, [x,θ], [ r=>upper, f=>upper])
out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
out=apply_derivative(out)
out=eliminate_kronecker(out)
D_SUPERFIELD=Derivative("D",out, [x,θ], [ l=>upper, d=>upper])
out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
out=apply_derivative(out)
out=eliminate_kronecker(out)


D_SUPERFIELD=Derivative("D",Wmmm, [x,θ], [ l=>upper, f=>upper])
out2=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
out2=apply_derivative(out2)
out2=eliminate_kronecker(out2)
D_SUPERFIELD=Derivative("D",out2, [x,θ], [ r=>upper, d=>upper])
out2=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
out2=apply_derivative(out2)
out2=eliminate_kronecker(out2)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ t=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ u=>upper, X=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# temp=copy(out)
out=simplify(eliminate_epsilon(EpsilonTensor([r=>lower, k=>lower])*(EpsilonTensor([f=>lower, c=>lower])*out)))
out2=simplify(eliminate_epsilon(EpsilonTensor([r=>lower, k=>lower])*(EpsilonTensor([f=>lower, c=>lower])*out2)))
res=Tensor("β")*out+Tensor("γ")*out2-Tensor("α")*apply_derivative(Derivative("delpp",Wmmm2,[x]))
res2=Wm

# out=EpsilonTensor([g=>lower, f=>lower])*(EpsilonTensor([X=>lower, h=>lower])*(EpsilonTensor([u=>lower, s=>lower])*(EpsilonTensor([t=>lower, r=>lower])*out)))
# res2=simplify(eliminate_epsilon(out))/144
# println("4 deriv")
# println(res2)
# res=res2-res


#DW---

# D_SUPERFIELD=Derivative("D",Wmmm, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=simplify(eliminate_epsilon(EpsilonTensor([r=>lower, k=>lower])*(EpsilonTensor([f=>lower, c=>lower])*out)))
# res=out

# res2=apply_derivative(Derivative("delmm",Scalar,[x]))


################# D2 D T


# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ r=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ t=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)

# ##no use!!
# # # D_SUPERFIELD=Derivative("D",out, [x,θ], [ u=>upper, X=>upper])
# # # out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# # # out=apply_derivative(out)
# # # out=eliminate_kronecker(out)
# #####

# out=EpsilonTensor([f=>lower, g=>lower])*(EpsilonTensor([t=>lower, s=>lower])*out)
# res2=im*simplify(eliminate_epsilon(out))/2

# # # #scalar con
# # #Δ_(3+) S
# D_SUPERFIELD=Derivative("D",Scalar, [x,θ], [ l=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ t=>upper, d=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)

# # #### not used
# # # # # D_SUPERFIELD=Derivative("D",out, [x,θ], [ u=>upper, X=>upper])
# # # # # out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# # # # # out=apply_derivative(out)
# # # # # out=eliminate_kronecker(out)
# # ######

# out=-1/9*(EpsilonTensor([t=>lower, s=>lower])*(EpsilonTensor([g=>lower, f=>lower])*out))
# res=Derivative("delmm",simplify(eliminate_epsilon(out)), [x])
# res=apply_derivative(res)


# #different contraction

# D_SUPERFIELD=Derivative("D",Scalar, [x,θ], [ t=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ l=>upper, d=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)

# β=Tensor("β")
# λ=Tensor("λ")

# out=-1/2*(EpsilonTensor([s=>lower, t=>lower])*(EpsilonTensor([g=>lower, f=>lower])*out))
# res2=Derivative("delmm",simplify(eliminate_epsilon(out)), [x])
# res2=-1*β*apply_derivative(res2)-λ*apply_derivative(Derivative("delpp",Wm, [x]))







# D_SUPERFIELD=Derivative("D",Wm, [x,θ], [ r=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
#  out=EpsilonTensor([r=>lower, l=>lower])*(EpsilonTensor([f=>lower, d=>lower])*out)
# res2=4*im*Derivative("delpp",simplify(eliminate_epsilon(out)), [x])
# res=apply_derivative(res)


#######################################################################################
############# w- different = D4Lmm
#l d free inds

# D_SUPERFIELD=Derivative("D",Wm, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
#  out=EpsilonTensor([r=>lower, l=>lower])*(EpsilonTensor([f=>lower, d=>lower])*out)
# res=4*im*Derivative("𝛻",simplify(eliminate_epsilon(out)), [x])
# res=apply_derivative(res)


# D_SUPERFIELD=Derivative("D",Lmm, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ t=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ u=>upper, X=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)


# out=EpsilonTensor([g=>lower, f=>lower])*(EpsilonTensor([X=>lower, h=>lower])*(EpsilonTensor([u=>lower, s=>lower])*(EpsilonTensor([t=>lower, r=>lower])*out)))
# res2=simplify(eliminate_epsilon(out))/144
# println("4 deriv")
# println(res2)
# res=res2-res


# a = Index("a", flavour2_indices,1)
# b = Index("b", flavour2_indices,1)
# c = Index("c", flavour2_indices,1)
# d = Index("d", flavour2_indices,1)
# e = Index("e", flavour2_indices,1)
# f = Index("f", flavour2_indices,1)
# g = Index("g", flavour2_indices,1)
# h = Index("h", flavour2_indices,1)
# X = Index("x", flavour2_indices,1)
# y = Index("y", flavour2_indices,1)
# z = Index("z", flavour2_indices,1)
# # H = Index("H", flavour2_indices,1)
# # X2 = Index("X", flavour2_indices,1)
# # Y = Index("Y", flavour2_indices,1)
# # Z = Index("Z", flavour2_indices,1)
# # # h = Index("h", flavour2_indices,1)

# i = Index("i", flavour_indices,0)
# j = Index("j", flavour_indices,0)
# k = Index("k", flavour_indices,0)
# l = Index("l", flavour_indices,0)
# n = Index("n", flavour_indices,0)
# q = Index("q", flavour_indices,0)
# r = Index("r", flavour_indices,0)
# s = Index("s", flavour_indices,0)
# t = Index("t", flavour_indices,0)
# u = Index("u", flavour_indices,0)
# v = Index("v", flavour_indices,0)






#######################################################################################
############# D4S
#l d free inds

# D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ l=>upper, d=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
#  out=im*(EpsilonTensor([r=>lower, l=>lower])*(EpsilonTensor([f=>lower, d=>lower])*out)/2)
# # res2=4*im*Derivative("𝛻",simplify(eliminate_epsilon(out)), [x])
# res=simplify(out)

# κ=Tensor("κ")
# #try D2 T +D2D2 T


# D_SUPERFIELD=Derivative("D",Scalar, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ t=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ u=>upper, X=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)


# out=EpsilonTensor([g=>lower, f=>lower])*(EpsilonTensor([X=>lower, h=>lower])*(EpsilonTensor([s=>lower, r=>lower])*(EpsilonTensor([u=>lower, t=>lower])*out)))
# res2=-1*simplify(eliminate_epsilon(out))/4
# println("4 deriv")
# println(res2)
# res2=apply_derivative(Derivative("delmm",res2,[x]))
# res=res-κ*res2




# D_SUPERFIELD=Derivative("D",Scalar, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=drop_weight_above(out,"θ",4)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=drop_weight_above(out,"θ",4)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ t=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=drop_weight_above(out,"θ",4)
# D_SUPERFIELD=Derivative("D",out, [x,θ], [ u=>upper, X=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# out=drop_weight_above(out,"θ",4)

# # # out=EpsilonTensor([g=>lower, f=>lower])*(EpsilonTensor([X=>lower, h=>lower])*(EpsilonTensor([u=>lower, s=>lower])*(EpsilonTensor([t=>lower, r=>lower])*out)))
#  out=EpsilonTensor([s=>lower, r=>lower])*(EpsilonTensor([f=>lower, X=>lower])*(EpsilonTensor([h=>lower, g=>lower])*(EpsilonTensor([u=>lower, t=>lower])*out)))
# # println(get_free_indices(out))
# # println(out)
# out=simplify(eliminate_epsilon(out))/144
# # println(out)
# # # res2=get_weighted_terms(out,"θ",0)
# res2=merge_sum(simplify(apply_derivative(Derivative("delmm",out, [x]))))
# # # println(res2)
# # # println("4 deriv")
# # # println(res2)
# # # res=res2-res



# D_SUPERFIELD=Derivative("D",Wm, [x,θ], [ l=>upper, d=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
# D_SUPERFIELD=Derivative("D",Wm, [x,θ], [ r=>upper, f=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# out=apply_derivative(out)
# out=eliminate_kronecker(out)
#  out=im*(EpsilonTensor([r=>lower, l=>lower])*(EpsilonTensor([f=>lower, d=>lower])*out)/2)
# # res=apply_derivative(Derivative("delpp",simplify(eliminate_epsilon(out)), [x]))
# res=simplify(out)







#######################################################################################

# out=EpsilonTensor([r=>lower, k=>lower])*(EpsilonTensor([f=>lower, c=>lower])*out)
# res=simplify(eliminate_epsilon(out))





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

# res=Tens×Scalar+EpsilonTensor([l=>lower, k=>lower])*(EpsilonTensor([d=>lower, c=>lower])*(Wm×Wmmm))
# # res2=simplify(res)
# # println(res2)
# # final=get_weighted_terms(get_weighted_terms(res2,"θ",4),"Θ",0)
# # final=to_lc_expression(sort(final),Dict{Index, Index}())
# project_component(res,θ,[i=>upper, a=>upper],4,0)

#  D_SUPERFIELD=Derivative("D",Tens, [x,θ], [ r=>upper, f=>upper])
#  out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# # #Wm=(Tensor("J",[   l=>upper, d=>upper], x)

####

# dL=Derivative("D",Wm, [x,θ], [r => upper, f=>upper])
# dL=expand_derivative(dL,covar,covar_unbar_expansion)
# dL=apply_derivative(dL)
# red=(EpsilonTensor([r=>lower, l=>lower])*(EpsilonTensor([f=>lower, d=>lower])*dL))
# println(red)
# red=simplify(red)
# rr=(red×Lmm)
# rr=simplify(rr)
# res=Tens×Scalar+rr
# # res2=simplify(res)
# # println(res2)
# # final=get_weighted_terms(get_weighted_terms(res2,"θ",4),"Θ",0)
# # final=to_lc_expression(sort(final),Dict{Index, Index}())
# res=simplify(res)
# println(res)
# res=project_component(res,θ,[j=>upper, b=>upper],4,0)
# println(to_latex(res))



# #alt:
# dL=Derivative("D",Lmm, [x,θ], [r => upper, f=>upper])
# dL=expand_derivative(dL,covar,covar_unbar_expansion)
# dL=apply_derivative(dL)
# dL=simplify(dL)
# # red=(EpsilonTensor([l=>lower, r=>lower])*(EpsilonTensor([d=>lower, f=>lower])*Wm))
# # println(red)
# # red=simplify(red)
# rr=(dL×Wm)
# # rr=(dL×apply_derivative(Derivative("delmm",Wm,[x])))
# println(get_free_indices(rr))
# rr=eliminate_epsilon(EpsilonTensor([l=>lower, r=>lower])*(EpsilonTensor([d=>lower, f=>lower])*rr))
# res=Tens×Scalar+rr
# # res2=simplify(res)
# # println(res2)
# # final=get_weighted_terms(get_weighted_terms(res2,"θ",4),"Θ",0)
# # final=to_lc_expression(sort(final),Dict{Index, Index}())
# res=simplify(res)
# println(res)
# res=project_component(res,θ,[j=>upper, b=>upper],4,0)
# println(to_latex(res))

# ###


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


############## another SC^2


# Wm2=(Tensor("J",[   r=>upper, f=>upper], x)
# +θia*Tensor("M",[ r=>upper, f=>upper,  i=>upper, a=>upper],x)
# +θab*Tensor("N",[ r=>upper, f=>upper,     a=>upper, b=>upper],x)
# +θij*Tensor("P",[ r=>upper, f=>upper,     i=>upper, j=>upper],x)
# +θ3*Tensor("Q",[ r=>upper, f=>upper,  b=>upper, j=>upper],x)
# +θ4*Tensor("R",[ r=>upper, f=>upper],x)
# )



# D_SUPERFIELD=Derivative("D",Wm×Wm2, [x,θ], [ t=>upper, h=>upper])
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# out=apply_derivative(out)
# out=eliminate_kronecker(out)


# D_SUPERFIELD2=Derivative("D",out, [x,θ], [ s=>upper, g=>upper])
# out2=expand_derivative(D_SUPERFIELD2,covar,covar_unbar_expansion)
# o=out2
# out2=apply_derivative(out2)
# out2=eliminate_kronecker(out2)

# # exp= out2×Wm2
# exp=out2

# println(exp)
# println(get_free_indices(exp))

# exp=Tensor("κ")*(EpsilonTensor([f=>lower,h=>lower])*(EpsilonTensor([d=>lower,g=>lower])*(EpsilonTensor([l=>lower,r=>lower])*(EpsilonTensor([s=>lower,t=>lower])*exp))))

# # exp=get_weighted_terms(exp,"θ",0)


# res=Tens×Scalar-exp
# # # res2=simplify(res)
# println(res)
# res=get_weighted_terms(get_weighted_terms(res,"θ",4),"Θ",0)
# # # final=to_lc_expression(sort(final),Dict{Index, Index}())
# println(res)
# println("starting simplify")
# res=simplify(res)
# println(res)
# println("starting projection")
# res=project_component(res,θ,[j=>upper, b=>upper],4,0)
# println(to_latex(res))

# D_SUPERFIELD3=Derivative("D",Wm3, [x,θ], [ l=>upper, f=>upper])
# out3=expand_derivative(D_SUPERFIELD3,covar,covar_unbar_expansion)
# out3=apply_derivative(out3)
# out3=eliminate_kronecker(out3)

# D_SUPERFIELD4=Derivative("D",Wm4, [x,θ], [ r=>upper, d=>upper])
# out4=expand_derivative(D_SUPERFIELD4,covar,covar_unbar_expansion)
# out4=apply_derivative(out4)
# out4=eliminate_kronecker(out4)



# res=simplify(out+out2)/2
