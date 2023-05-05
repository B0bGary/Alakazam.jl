include("alakazam.jl")

#define index sets

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


a = Index("a", flavour2_indices,0)
b = Index("b", flavour2_indices,0)
c = Index("c", flavour2_indices,0)
d = Index("d", flavour2_indices,0)
e = Index("e", flavour2_indices,0)
f = Index("f", flavour2_indices,0)

i = Index("i", flavour_indices,0)
j = Index("j", flavour_indices,0)
k = Index("k", flavour_indices,0)
l = Index("l", flavour_indices,0)
n = Index("n", flavour_indices,0)
q = Index("q", flavour_indices,0)
r = Index("r", flavour_indices,0)

x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices, flavour_indices, flavour2_indices])
# θbar = Coordinate("Θ", [spacetime_indices])

grassman_deriv=Derivative("𝛛",zero, [θ], [μ =>lower, i => upper, a =>upper])
# grassman_deriv_bar=Derivative("∂",zero, [θbar], [a => lower, i=>lower])
spatial_deriv=Derivative("∇",zero, [x], [μ => lower, ν => lower])

# add_derivative_action(grassman_deriv_bar,Tensor("Θ",[c=>upper, j=>upper],θbar),KroneckerDelta([a =>lower, c=>upper])*KroneckerDelta([i =>lower, j=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[ν=>upper, j=>lower, b=>lower],θ),KroneckerDelta([μ =>lower, ν=>upper])*KroneckerDelta([j =>lower, i=>upper])*KroneckerDelta([b =>lower, a=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[ν=>upper, j=>upper, b=>lower],θ),KroneckerDelta([μ =>lower, ν=>upper])*EpsilonTensor([i =>upper, j=>upper])*KroneckerDelta([b =>lower, a=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[ν=>upper, j=>lower, b=>upper],θ),KroneckerDelta([μ =>lower, ν=>upper])*KroneckerDelta([j =>lower, i=>upper])*EpsilonTensor([a =>upper, b=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[ν=>upper, j=>upper, b=>upper],θ),KroneckerDelta([μ =>lower, ν=>upper])*EpsilonTensor([i =>upper, j=>upper])*EpsilonTensor([a =>upper, b=>upper]))


# covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[b=>upper, i=>lower],θ)*spatial_deriv
covar_unbar_expansion=grassman_deriv-im/2*Tensor("θ",[ν=>upper, i=>lower, a=>lower],θ)*spatial_deriv
covar=Derivative("D",zero, [x,θ], [a => lower, i=>upper])
# covar_bar=Derivative("D",zero, [x,θbar], [a => lower, i=>lower])

# θgrass=Tensor("θ",[a=>upper],θ)
# θbargrass=Tensor("Θ",[b=>upper],θbar)

custom_sort_order(["θ","Θ","ϵ"])
# θci=Tensor("θ",[c=>upper, i=> lower],θ)
# barθci=Tensor("Θ",[c=>upper, i=> lower],θbar)

θia=Tensor("θ",[μ=>upper, i=> lower, a=>lower],θ)
θab=Tensor("θ",[μ=>upper, i=> lower, a=>lower],θ)*Tensor("θ",[ν=>upper, i=> upper, b=>lower],θ)
θij=Tensor("θ",[μ=>upper, i=> lower, a=>lower],θ)*Tensor("θ",[ν=>upper, j=> lower, a=>upper],θ)
θ3=Tensor("θ",[ρ=>upper, j=> lower, a=>upper],θ)*Tensor("θ",[μ=>upper, i=> lower, a=>lower],θ)*Tensor("θ",[ν=>upper, i=> upper, b=>lower],θ)
θ4=Tensor("θ",[σ=>upper, j=> upper, b=>upper],θ)*Tensor("θ",[ρ=>upper, j=> lower, a=>upper],θ)*Tensor("θ",[μ=>upper, i=> lower, a=>lower],θ)*Tensor("θ",[ν=>upper, i=> upper, b=>lower],θ)


#last two indices are the spacetime ones
Scalar=(Tensor("j", x)
+θia*Tensor("ψ",[μ=>lower, i=>upper, a=>upper],x)
+θab*Tensor("F",[μ=>lower,ν=>lower, a=>upper, b=>upper],x)
+θij*Tensor("K",[μ=>lower,ν=>lower, i=>upper, j=>upper],x)
+θ3*Tensor("S",[μ=>lower,ν=>lower,ρ=>lower, b=>upper, j=>lower],x)
+θ4*Tensor("T",[μ=>lower,ν=>lower,ρ=>lower, σ=>lower],x)
)

Tens=(Tensor("T", x)
+θia*Tensor("A",[μ=>lower, i=>upper, a=>upper],x)
+θab*Tensor("B",[μ=>lower,ν=>lower, a=>upper, b=>upper],x)
+θij*Tensor("C",[μ=>lower,ν=>lower, i=>upper, j=>upper],x)
+θ3*Tensor("G",[μ=>lower,ν=>lower,ρ=>lower, b=>upper, j=>lower],x)
+θ4*Tensor("H",[μ=>lower,ν=>lower,ρ=>lower, σ=>lower],x)
)


dScalar=Derivative("∇",Scalar, [x], [τ => lower, λ=>upper])




# D_SUPERFIELD=Derivative("D",SUPERFIELDpp, [x,θ], [g => lower, n=>upper])
# println(covar)
# println(covar_unbar_expansion)
# println(covar_bar)
# println(covar_expansion)
# println(D_SUPERFIELD)
# out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
# o=out
# println(o)
out=apply_derivative(dScalar)
println("Found deriv: ", out)
println("Elim delta: ")
out=eliminate_kronecker(out)
println(out)

# res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
res=to_lc_expression(out,Dict(τ=>m, λ=>m))
println(res)
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