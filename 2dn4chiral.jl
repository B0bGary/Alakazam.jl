include("alakazam.jl")

#define index sets

spacetime_indices = IndexSet("spacetimelc", 2,lower)
flavour_indices = IndexSet("flavour_indices", 2,upper)
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
n = Index("n", flavour_indices,0)
q = Index("q", flavour_indices,0)
r = Index("r", flavour_indices,0)

x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices])
θbar = Coordinate("Θ", [spacetime_indices])

grassman_deriv=Derivative("𝛛",zero, [θ], [a => lower, i =>upper])
grassman_deriv_bar=Derivative("∂",zero, [θbar], [a => lower, i=>lower])
spatial_deriv=Derivative("∇",zero, [x], [a => lower, b => lower])

add_derivative_action(grassman_deriv_bar,Tensor("Θ",[c=>upper, j=>upper],θbar),KroneckerDelta([a =>lower, c=>upper])*KroneckerDelta([i =>lower, j=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[c=>upper, j=>lower],θ),KroneckerDelta([a =>lower, c=>upper])*KroneckerDelta([j =>lower, i=>upper]))


covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[b=>upper, i=>lower],θ)*spatial_deriv
covar_unbar_expansion=grassman_deriv-im/2*Tensor("Θ",[b=>upper, i=>upper],θbar)*spatial_deriv
covar=Derivative("D",zero, [x,θ], [a => lower, i=>upper])
covar_bar=Derivative("D",zero, [x,θbar], [a => lower, i=>lower])

θgrass=Tensor("θ",[a=>upper],θ)
θbargrass=Tensor("Θ",[b=>upper],θbar)

custom_sort_order(["θ","Θ","ϵ"])
θci=Tensor("θ",[c=>upper, i=> lower],θ)
barθci=Tensor("Θ",[c=>upper, i=> lower],θbar)
#last two indices are the spacetime ones
SUPERFIELDpp=(Tensor("j", x)
+Tensor("θ",[c=>upper, i=> lower],θ)*Tensor("ψ",[c=>lower, i=>upper],x)
+Tensor("Θ",[c=>upper, i=> upper],θbar)*Tensor("̄ψ",[c=>lower, i=>lower],x)
+Tensor("θ",[c=>upper, i=> lower],θ)*Tensor("θ",[a=>upper, j=> lower],θ)*Tensor("F",[c=>lower, a=>lower, i=>upper,j=>upper],x)
+Tensor("Θ",[c=>upper, i=> upper],θbar)*Tensor("Θ",[a=>upper, j=> upper],θbar)*Tensor("̄f",[c=>lower, a=>lower, i=>lower,j=>lower],x)
+Tensor("θ",[c=>upper, i=> lower],θ)*Tensor("Θ",[a=>upper, j=> upper],θbar)*Tensor("K",[c=>lower, a=>lower, i=>upper,j=>lower],x)
+Tensor("θ",[c=>upper, i=> lower],θ)*Tensor("θ",[a=>upper, j=> lower],θ)*Tensor("Θ",[b=>upper, k=> upper],θbar)*Tensor("S",[c=>lower, a=>lower, b=>lower, i=>upper,j=>upper,k=>lower],x)
+Tensor("Θ",[c=>upper, i=> upper],θbar)*Tensor("Θ",[a=>upper, j=> upper],θbar)*Tensor("θ",[b=>upper, k=> lower],θ)*Tensor("̄s",[c=>lower, a=>lower, b=>lower, i=>lower,j=>lower,k=>upper],x)
+Tensor("θ",[c=>upper, i=> lower],θ)*Tensor("θ",[a=>upper, j=> lower],θ)*Tensor("Θ",[b=>upper, k=> upper],θbar)*Tensor("Θ",[d=>upper, l=> upper],θbar)*Tensor("T",[c=>lower, a=>lower, b=>lower, d=>lower, i=>upper,j=>upper,k=>lower,l=>lower],x)

# +Tensor("θ",[c=>upper],θ)*Tensor("S",[c => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("̄S",[c => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("Θ",[d=>upper],θbar)*Tensor("T",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("K",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("̄K",[c => lower,d => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
# +Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("θ",[e=>upper],θ)*Tensor("̄H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
# +Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("Θ",[f=>upper],θbar)*Tensor("F",[c => lower,d => lower,e => lower, f=>lower, a =>lower, b => lower], x)
)


D_SUPERFIELD=Derivative("D",SUPERFIELDpp, [x,θ], [g => lower, n=>upper])
println(covar)
println(covar_unbar_expansion)
println(covar_bar)
println(covar_expansion)
println(D_SUPERFIELD)
out=expand_derivative(D_SUPERFIELD,covar,covar_unbar_expansion)
o=out
println(o)
out=apply_derivative(out)
println("Found deriv: ", out)
println("Elim delta: ")
out=eliminate_kronecker(out)
println(out)

# res=to_lc_expression(out,Dict(a=>p,b=>p,i=>m))
res=to_lc_expression(out,Dict(g=>p))
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