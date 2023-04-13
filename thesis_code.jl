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
                            new_term*=KroneckerDelta([term.indices[indexes], unused_flav[2*tcounter-1]=>lower])
                        else
                            new_term*=KroneckerDelta([term.indices[indexes], unused_spin[tcounter]=>lower])
                        end
                    else
                        if term.indices[indexes][1].index_set==flavour_indices
                            new_term*=KroneckerDelta([term.indices[indexes], unused_flav[2*tcounter-1]=>upper])
                        else
                            new_term*=KroneckerDelta([term.indices[indexes], unused_spin[tcounter]=>upper])
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
                            new_term*=KroneckerDelta([term.indices[indexes], unused_flav[2*tbcounter]=>lower])
                        else
                            new_term*=KroneckerDelta([term.indices[indexes], unused_conj[tbcounter]=>lower])
                        end
                    else
                        if term.indices[indexes][1].index_set==flavour_indices
                            new_term*=KroneckerDelta([term.indices[indexes], unused_flav[2*tbcounter]=>upper])
                        else
                            new_term*=KroneckerDelta([term.indices[indexes], unused_conj[tbcounter]=>upper])
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
flavour_indices = IndexSet("flavour", 2,lower)
spinor_indices = IndexSet("spinor", 2)
spacetime_indices = IndexSet("spacetime", 4,lower)
spinor_conj_indices = IndexSet("spinor_conj", 2,lower)
flavour_conj_indices = IndexSet("flavour_conj", 2)


#define indices
i = Index("i", flavour_indices, 0)
j = Index("j", flavour_indices)
k = Index("k", flavour_indices, 0)
l = Index("l", flavour_indices)
m = Index("m", flavour_indices)
n = Index("n", flavour_indices)
p = Index("p", flavour_indices)
q = Index("q", flavour_indices)
r = Index("r", flavour_indices)
s = Index("s", flavour_indices)
t = Index("t", flavour_indices)
u = Index("u", flavour_indices)
v = Index("v", flavour_indices)
x = Index("x", flavour_indices)
y = Index("y", flavour_indices)
z = Index("z", flavour_indices)
A = Index("A", flavour_indices)
B = Index("B", flavour_indices)
C = Index("C", flavour_indices)

# ic = Index("̄i", flavour_conj_indices, 0)
# jc = Index("̄j", flavour_conj_indices)
# kc = Index("̄k", flavour_conj_indices, 0)
# lc = Index("̄l", flavour_conj_indices)
# mc = Index("̄m", flavour_conj_indices)
# nc = Index("̄n", flavour_conj_indices)
# pc = Index("̄p", flavour_conj_indices)
# qc = Index("̄q", flavour_conj_indices)

α = Index("α", spinor_indices, 1)
β = Index("β", spinor_indices, 1)
γ = Index("γ", spinor_indices, 1)
μ = Index("μ", spinor_indices, 1)
ν = Index("ν", spinor_indices, 1)
ρ = Index("ρ", spinor_indices, 1)
τ = Index("τ", spinor_indices, 1)
σ = Index("σ", spinor_indices, 1)
ω = Index("ω", spinor_indices, 1)
χ = Index("χ", spinor_indices, 1)
π = Index("π", spinor_indices, 1)
η = Index("η", spinor_indices, 1)

αdot = Index("̇α", spinor_conj_indices, 1)
βdot = Index("̇β", spinor_conj_indices, 1)
γdot = Index("̇γ", spinor_conj_indices, 1)
ρdot = Index("̇ρ", spinor_conj_indices, 1)
μdot = Index("̇μ", spinor_conj_indices, 1)
νdot = Index("̇ν", spinor_conj_indices, 1)
τdot = Index("̇τ", spinor_conj_indices, 1)
σdot = Index("̇σ", spinor_conj_indices, 1)
ωdot = Index("̇ω", spinor_conj_indices, 1)
πdot = Index("̇π", spinor_conj_indices, 1)
ηdot = Index("̇η", spinor_conj_indices, 1)

a = Index("a", spacetime_indices, 0)
b = Index("b", spacetime_indices)
c = Index("c", spacetime_indices)
d = Index("d", spacetime_indices)
e = Index("e", spacetime_indices)
f = Index("f", spacetime_indices)
g = Index("g", spacetime_indices)


#define coordinates, tensors will be functions of these coordinates (for use in derivatives)
x = Coordinate("x", spacetime_indices)
θ = Coordinate("θ", [flavour_indices, spinor_indices])
θbar = Coordinate("Θ", [flavour_conj_indices, spinor_conj_indices])

# define our grassman basis ̄θ 
θiα=Tensor("θ",[i=>lower, α=>upper],θ)
barθαi=Tensor("Θ",[i=>upper,αdot=>lower],θbar)
barθij=Tensor("Θ",[i=>lower, νdot=>upper],θbar)*Tensor("Θ",[j=>lower, νdot=>lower],θbar)
barθupij=Tensor("Θ",[i=>upper, μdot=>upper],θbar)*Tensor("Θ",[j=>upper, μdot=>lower],θbar)
θij=Tensor("θ",[i=>lower, ν=>upper],θ)*Tensor("θ",[j=>lower, ν=>lower],θ)
θupij=Tensor("θ",[i=>upper, μ=>upper],θ)*Tensor("θ",[j=>upper, μ=>lower],θ)
θjα=Tensor("θ",[j=>lower, α=>upper],θ)
θ4=Tensor("θ",[m=>lower, ν=>upper],θ)*Tensor("θ",[n=>lower, ν=>lower],θ)*Tensor("θ",[m=>upper, μ=>upper],θ)*Tensor("θ",[n=>upper, μ=>lower],θ)
θbar4=Tensor("Θ",[p=>lower, νdot=>upper],θbar)*Tensor("Θ",[q=>lower, νdot=>lower],θbar)*Tensor("Θ",[p=>upper, μdot=>upper],θbar)*Tensor("Θ",[q=>upper, μdot=>lower],θbar)
θᵅᵝ=Tensor("θ",[y=>upper, α=>upper],θ)*Tensor("θ",[y=>lower, β=>upper],θ)
θᵦᵧ =Tensor("θ",[y=>upper, β=>lower],θ)*Tensor("θ",[y=>lower, γ=>lower],θ)
barθᵅᵝ=Tensor("Θ",[z=>upper, αdot=>upper],θbar)*Tensor("Θ",[z=>lower, βdot=>upper],θbar)
barθᵦᵧ =Tensor("Θ",[z=>upper, βdot=>lower],θbar)*Tensor("Θ",[z=>lower, γdot=>lower],θbar)
θiβ=Tensor("θ",[i=>lower, β=>upper],θ)
# define a basic proof of concept derivative
grassman_deriv_normalrep=Derivative("∂",zero, [θ], [j => upper, β => lower])
grassman_deriv_conjrep=Derivative("𝝏",zero, [θbar], [j => lower, βdot => upper])
spatial_deriv=Derivative("∇",zero, [x], [b => lower])
#note that the derivative isn't defined on literally θiα, but on tensors with the same index patterns
add_derivative_action(grassman_deriv_normalrep,θiα,KroneckerDelta([i=>lower,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>lower, α=>lower],θ),KroneckerDelta([i=>lower,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>upper, α=>upper],θ),EpsilonTensor([i=>upper,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>upper, α=>lower],θ),EpsilonTensor([i=>upper,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))

add_derivative_action(grassman_deriv_conjrep,barθαi,KroneckerDelta([i=>upper,j=>lower])*KroneckerDelta([βdot=>upper,αdot=>lower]))


#define our superfield expansions tensors
# exp=(θjα*θupij)*Tensor("G",[i=>lower, α=>lower],x)
# println(exp)
# SUPERFIELD=(Tensor("j", x)+θiα*Tensor("ψ",[i=>upper, α=>lower],x)+θij*SymmetricTensor("F",[i=>upper, j=>upper],x)
# +(θjα*θupij)*Tensor("G",[i=>lower, α=>lower],x)+θij*θupij*Tensor("H",x)
# +θjα*barθαi*Tensor("K",[j=>upper, α=>lower,i => lower, αdot =>upper],x)
# +Tensor("θ",[k=>lower, α=>upper],θ)*barθupij*Tensor("L",[k=>upper, α =>lower, i=>lower,j => lower],x)
# )
SUPERFIELD=(Tensor("j", x)
+θiα*Tensor("ψ",[i=>upper, α=>lower],x)+θij*SymmetricTensor("F",[i=>upper, j=>upper],x)
+(θjα*θupij)*Tensor("G",[i=>lower, α=>lower],x)+θij*θupij*Tensor("H",x)
+θjα*barθαi*Tensor("K",[j=>upper, α=>lower,i => lower, αdot =>upper],x)
+Tensor("θ",[k=>lower, α=>upper],θ)*barθupij*Tensor("L",[k=>upper, α =>lower, i=>lower,j => lower],x)
+Tensor("θ",[k=>lower, α=>upper],θ)*barθαi*barθij*Tensor("M",[k=>upper, α =>lower, j=>lower, αdot=> upper],x)
+θiα*θbar4*Tensor("N",[i=>upper, α=>lower],x)
+θij*Tensor("Θ",[m=>lower, νdot=>upper],θbar)*Tensor("Θ",[n=>lower, νdot=>lower],θbar)*Tensor("P",[i=>upper,j=>upper,m=>upper,n=>upper],x)
+θij*Tensor("Θ",[k=>upper, αdot=>lower],θbar)*Tensor("Θ",[k=>lower, νdot=>upper],θbar)*Tensor("Θ",[l=>lower, νdot=>lower],θbar)*Tensor("Q",[i=>upper,j=>upper,l=>upper,αdot=>upper],x)
+θij*θbar4*SymmetricTensor("R",[i=>upper, j=>upper],x)
+θupij*Tensor("θ",[j=>lower, α=>upper],θ)*Tensor("Θ",[k=>upper, αdot=>lower],θbar)*Tensor("Θ",[k=>lower, νdot=>upper],θbar)*Tensor("Θ",[l=>lower, νdot=>lower],θbar)*Tensor("S",[i=>lower,α=>lower,l=>upper,αdot=>upper],x)
+Tensor("θ",[j=>lower, α=>upper],θ)*θupij*θbar4*Tensor("T",[i=>lower, α=>lower],x)
+θ4*θbar4*Tensor("V",x)
+barθαi*Tensor("ϕ",[i=>lower, αdot=>upper],x)
+barθij*Tensor("A",[i=>upper, j=>upper],x)
+barθαi*barθij*Tensor("B",[j=>upper, αdot=>upper],x)
+θbar4*Tensor("C",x)
+barθαi*Tensor("θ",[j=>upper, ν=>upper],θ)*Tensor("θ",[k=>upper, ν=>lower],θ)*Tensor("U",[j=>lower,k=>lower, i=>lower, αdot=>upper],x)
+barθαi*Tensor("θ",[k=>lower, β=>upper],θ)*Tensor("θ",[j=>upper, ν=>upper],θ)*Tensor("θ",[k=>upper, ν=>lower],θ)*Tensor("W",[j=>lower,β=>lower, i=>lower, αdot=>upper],x)
+barθαi*θ4*Tensor("X",[i=>lower, αdot=>upper],x)
+barθupij*Tensor("θ",[k=>upper, μ=>upper],θ)*Tensor("θ",[l=>upper, μ=>lower],θ)*Tensor("θ",[l=>lower, α=>upper],θ)*Tensor("Y",[k=>lower,α =>lower, i=>lower, j=>lower],x)
+barθupij*θ4*SymmetricTensor("Z",[i=>lower, j=>lower],x)
+Tensor("Θ",[j=>upper, αdot=>lower],θbar)*barθij*θ4*Tensor("E",[i=>lower, αdot=>upper],x)
+θᵦᵧ*SymmetricTensor("Ω",[β=>upper,γ=>upper],x)
+θiβ*θᵦᵧ*Tensor("Σ",[i=>upper,γ=>upper],x)
+barθᵦᵧ*Tensor("Λ",[βdot=>upper,γdot=>upper],x)
)

custom_sort_order(["θ","Θ","ϵ"])

# println(θupij*Tensor("θ",[j=>lower, α=>upper],θ)*Tensor("̄θ",[k=>upper, αdot=>lower],θbar)*Tensor("̄θ",[k=>lower, νdot=>upper],θbar)*Tensor("θ",[l=>lower, νdot=>lower],θbar)*Tensor("S",[i=>upper,α=>lower,l=>upper,αdot=>upper],x))
# break

covar_expansion=grassman_deriv_normalrep-im*Tensor("σ",[b =>upper, γdot=>upper, β =>lower])*Tensor("Θ",[j=>upper,γdot=>lower,],θbar)*spatial_deriv
covar=Derivative("D",zero, [x,θ], [j => upper, β => lower])
println("our superfield expansions is: ", SUPERFIELD)

#apply our derivative
println("Computing Symmetrised deriv: ")
D_SUPERFIELD=Derivative("D",SUPERFIELD, [x,θ], [s => upper, τ => lower])
out=expand_derivative(D_SUPERFIELD,covar,covar_expansion)
out=apply_derivative(out)
out=simplify(out)
DD_SUPERFIELD=Derivative("D",out,[x,θ], [r => upper, ρ => lower])

out=expand_derivative(DD_SUPERFIELD,covar,covar_expansion)
out=apply_derivative(out)
out=EpsilonTensor([τ=>upper,ρ=>upper])*out

out=eliminate_epsilon(out)
out=copy(out)
D_SUPERFIELD2=Derivative("D",SUPERFIELD, [x,θ], [r => upper, τ => lower])
out2=expand_derivative(D_SUPERFIELD2,covar,covar_expansion)
out2=apply_derivative(out2)
out2=simplify(out2)
DD_SUPERFIELD2=Derivative("D",out2,[x,θ], [s => upper, ρ => lower])

out2=expand_derivative(DD_SUPERFIELD2,covar,covar_expansion)
out2=apply_derivative(out2)
out2=EpsilonTensor([τ=>upper,ρ=>upper])*out2

out2=eliminate_epsilon(out2)
out2=copy(out2)
println("Result:")
println(out+out2)
res=simplify(out+out2)
println(res)
res=copy(res)

t11=get_weighted_terms(get_weighted_terms(res,"θ",1),"Θ",2)
println("t11 weight:")
println(merge_sum(t11))

t11=copy(t11)
println("rearrange")
t11=merge_sum(sort(rename_dummies(eliminate_epsilon(rename_dummies(sort(pop_metric(t11)))))))
println(t11)

println("elim delta")
t11=merge_sum(eliminate_kronecker(t11))
# for summand in t11.summands
#     println(summand)
#     println(get_free_indices(summand[2]))
# end
println(t11)
out=merge_sum(sort(pop_metric(t11)))
println(out)


println("stripping thetas")
out=strip_thetas(out)
println(out)
# for summand in out.summands
#         println(get_free_indices(summand[2]))
# end
println("constraint:")
println(copy(eliminate_epsilon(eliminate_kronecker(out))))








# println("We want to compute: ", DD_SUPERFIELD)
# println("Appyling the derivative:")
# out=apply_derivative(DD_SUPERFIELD)
# println(out)
# println("Getting rid of epsilon tensors")
# out=contract_metrics(out,EpsilonTensor)
# println(out)
# println("Eliminating deltas")
# out=eliminate_kronecker(out)
# println(out)
# println("Renaming dummies")
# rename_dummies(out)
# println(out)
# println("Simplifying")
# out=simplify(out)
# println(out)