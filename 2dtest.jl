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
                tens.indices=[ haskey(reps,ind[1]) ? reps[ind[1]]=>ind[2] : ind[1]=>ind[2]  for ind in tens.indices]
            end
            if !(typeof(term)<:SymorNum) && ((typeof(term)<:TensorTerm && tens in term.terms) || (typeof(term)<:TensorExpression && tens in term.summands[1][2].terms) )
                term=0
                break
            else
                term=term*tens
            end
        end
        new_exp+=term
    end


    #delete terms with repeated thetas



    return new_exp
end

function to_lc_expression(TEx::TensorExpression, reps::Dict{Index,Index})
    #replace non-dummies
    exp=specify_indices(TEx, reps)
    println(exp)
    #cookery for dummies....
    new_exp=zero

    for summand in exp.summands
        dums=get_dummies(summand[2])
        tempterm=TensorExpression(summand[2])
        while !isempty(dums)
            dum=pop!(dums)
            tempterm=specify_indices(tempterm, Dict(dum => ±))+specify_indices(tempterm, Dict(dum => ∓))
        end
        new_exp+=summand[1]*tempterm
    end
    # new_exp=zero
    # while !isempty(processing_stack)
    #     temp=pop!(processing_stack)

    # end
    return new_exp
end



#define index sets

spacetime_indices = IndexSet("spacetimelc", 2,lower)

lightcone_indices = IndexSet("lightcone", 2,lower)
± = Index("±", lightcone_indices, 1)
∓ = Index("∓", lightcone_indices, 1)


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
m = Index("m", spacetime_indices,1)
n = Index("n", spacetime_indices,1)


x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices])
θbar = Coordinate("Θ", [spacetime_indices])


grassman_deriv_bar=Derivative("∂",zero, [θbar], [a => lower])
spatial_deriv=Derivative("∇",zero, [x], [a => lower, b => lower])

add_derivative_action(grassman_deriv_bar,Tensor("Θ",[c=>upper],θbar),KroneckerDelta([a =>lower, c=>upper]))

covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[b=>upper],θ)*spatial_deriv
covar=Derivative("̄D",zero, [x,θbar], [a => lower])


θgrass=Tensor("θ",[a=>upper],θ)
θbargrass=Tensor("Θ",[b=>upper],θbar)

custom_sort_order(["θ","Θ","ϵ"])

#last two indices are the spacetime ones
SUPERFIELDpp=(Tensor("j",[a =>lower, b => lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("S",[c => lower, a =>lower, b => lower], x)
+Tensor("Θ",[c=>upper],θbar)*Tensor("̄S",[c => lower, a =>lower, b => lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("Θ",[d=>upper],θbar)*Tensor("T",[c => lower,d => lower, a =>lower, b => lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("K",[c => lower,d => lower, a =>lower, b => lower], x)
+Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("̄K",[c => lower,d => lower, a =>lower, b => lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
+Tensor("Θ",[c=>upper],θbar)*Tensor("Θ",[d=>upper],θbar)*Tensor("θ",[e=>upper],θ)*Tensor("̄H",[c => lower,d => lower,e => lower, a =>lower, b => lower], x)
+Tensor("θ",[c=>upper],θ)*Tensor("θ",[d=>upper],θ)*Tensor("Θ",[e=>upper],θbar)*Tensor("Θ",[f=>upper],θbar)*Tensor("F",[c => lower,d => lower,e => lower, f=>lower, a =>lower, b => lower], x)
)


D_SUPERFIELD=Derivative("̄D",SUPERFIELDpp, [x,θbar], [i => lower])

out=expand_derivative(D_SUPERFIELD,covar,covar_expansion)

out=apply_derivative(out)
println("Found deriv: ", out)
println("Elim delta: ")
out=eliminate_kronecker(out)
println(out)

to_lc_expression(out,Dict(a=>±,b=>±,i=>∓))






# +θbarp*Tensor("χ",[p =>lower],x)+θbarm*Tensor("χ",[m =>lower],x)
# +θp*θm*Tensor("A",[p =>lower,m=>lower],x)+θbarp*θbarm*Tensor("B",[p =>lower,m=>lower],x)
# +θbarp*θp*θm*Tensor("S",[p =>lower,p=>lower,m=>lower],x)+θbarm*θp*θm*Tensor("S",[m =>lower,p=>lower,m=>lower],x)
# +θp*θbarp*θbarm*Tensor("̄S",[p =>lower,p =>lower,m=>lower],x)+θp*θbarp*θbarm*Tensor("̄S",[p =>lower,p =>lower,m=>lower],x)
# +θbarp*θp*θm*Tensor("T",[p =>lower,p=>lower,m=>lower],x)
# +θp*θbarp*θbarm*Tensor("̄T",[p =>lower,p =>lower,m=>lower],x)+θp*θbarp*θbarm*Tensor("̄T",[p =>lower,p =>lower,m=>lower],x)
# )







# # define our grassman basis ̄θ 
# θiα=Tensor("θ",[i=>lower, α=>upper],θ)
# barθαi=Tensor("Θ",[i=>upper,αdot=>lower],θbar)
# barθij=Tensor("Θ",[i=>lower, νdot=>upper],θbar)*Tensor("Θ",[j=>lower, νdot=>lower],θbar)
# barθupij=Tensor("Θ",[i=>upper, μdot=>upper],θbar)*Tensor("Θ",[j=>upper, μdot=>lower],θbar)
# θij=Tensor("θ",[i=>lower, ν=>upper],θ)*Tensor("θ",[j=>lower, ν=>lower],θ)
# θupij=Tensor("θ",[i=>upper, μ=>upper],θ)*Tensor("θ",[j=>upper, μ=>lower],θ)
# θjα=Tensor("θ",[j=>lower, α=>upper],θ)
# θ4=Tensor("θ",[m=>lower, ν=>upper],θ)*Tensor("θ",[n=>lower, ν=>lower],θ)*Tensor("θ",[m=>upper, μ=>upper],θ)*Tensor("θ",[n=>upper, μ=>lower],θ)
# θbar4=Tensor("Θ",[p=>lower, νdot=>upper],θbar)*Tensor("Θ",[q=>lower, νdot=>lower],θbar)*Tensor("Θ",[p=>upper, μdot=>upper],θbar)*Tensor("Θ",[q=>upper, μdot=>lower],θbar)
# θᵅᵝ=Tensor("θ",[y=>upper, α=>upper],θ)*Tensor("θ",[y=>lower, β=>upper],θ)
# θᵦᵧ =Tensor("θ",[y=>upper, β=>lower],θ)*Tensor("θ",[y=>lower, γ=>lower],θ)
# barθᵅᵝ=Tensor("Θ",[z=>upper, αdot=>upper],θbar)*Tensor("Θ",[z=>lower, βdot=>upper],θbar)
# barθᵦᵧ =Tensor("Θ",[z=>upper, βdot=>lower],θbar)*Tensor("Θ",[z=>lower, γdot=>lower],θbar)
# θiβ=Tensor("θ",[i=>lower, β=>upper],θ)
# barθjα=Tensor("Θ",[j=>upper,αdot=>lower],θbar)
# barθjk=Tensor("Θ",[j=>lower, νdot=>upper],θbar)*Tensor("Θ",[k=>lower, νdot=>lower],θbar)
# θjk=Tensor("θ",[j=>lower, ν=>upper],θ)*Tensor("θ",[k=>lower, ν=>lower],θ)
# θʲᵏ=Tensor("θ",[j=>upper, ν=>upper],θ)*Tensor("θ",[k=>upper, ν=>lower],θ)
# # define a basic proof of concept derivative
# grassman_deriv_normalrep=Derivative("∂",zero, [θ], [j => upper, β => lower])
# grassman_deriv_conjrep=Derivative("𝝏",zero, [θbar], [j => lower, βdot => upper])
# spatial_deriv=Derivative("∇",zero, [x], [b => lower])
# #note that the derivative isn't defined on literally θiα, but on tensors with the same index patterns
# add_derivative_action(grassman_deriv_normalrep,θiα,KroneckerDelta([i=>lower,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
# add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>lower, α=>lower],θ),KroneckerDelta([i=>lower,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))
# add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>upper, α=>upper],θ),EpsilonTensor([i=>upper,j=>upper])*KroneckerDelta([β=>lower,α=>upper]))
# add_derivative_action(grassman_deriv_normalrep,Tensor("θ",[i=>upper, α=>lower],θ),EpsilonTensor([i=>upper,j=>upper])*EpsilonTensor([α=>lower,β=>lower]))

# add_derivative_action(grassman_deriv_conjrep,barθαi,KroneckerDelta([i=>upper,j=>lower])*KroneckerDelta([βdot=>upper,αdot=>lower]))


#define our superfield expansions tensors
# exp=(θjα*θupij)*Tensor("G",[i=>lower, α=>lower],x)
# println(exp)
# SUPERFIELD=(Tensor("j", x)+θiα*Tensor("ψ",[i=>upper, α=>lower],x)+θij*SymmetricTensor("F",[i=>upper, j=>upper],x)
# +(θjα*θupij)*Tensor("G",[i=>lower, α=>lower],x)+θij*θupij*Tensor("H",x)
# +θjα*barθαi*Tensor("K",[j=>upper, α=>lower,i => lower, αdot =>upper],x)
# +Tensor("θ",[k=>lower, α=>upper],θ)*barθupij*Tensor("L",[k=>upper, α =>lower, i=>lower,j => lower],x)
# # )
# SUPERFIELD=(Tensor("j", x)
# +θiα*Tensor("ψ",[i=>upper, α=>lower],x)+θij*SymmetricTensor("F",[i=>upper, j=>upper],x)
# +(θjα*θupij)*Tensor("G",[i=>lower, α=>lower],x)+θij*θupij*Tensor("H",x)
# +θjα*barθαi*Tensor("K",[j=>upper, α=>lower,i => lower, αdot =>upper],x)
# +Tensor("θ",[k=>lower, α=>upper],θ)*barθupij*Tensor("L",[k=>upper, α =>lower, i=>lower,j => lower],x)
# +Tensor("θ",[k=>lower, α=>upper],θ)*barθαi*barθij*Tensor("M",[k=>upper, α =>lower, j=>lower, αdot=> upper],x)
# +θiα*θbar4*Tensor("N",[i=>upper, α=>lower],x)
# +θij*Tensor("Θ",[m=>lower, νdot=>upper],θbar)*Tensor("Θ",[n=>lower, νdot=>lower],θbar)*Tensor("P",[i=>upper,j=>upper,m=>upper,n=>upper],x)
# +θij*Tensor("Θ",[k=>upper, αdot=>lower],θbar)*Tensor("Θ",[k=>lower, νdot=>upper],θbar)*Tensor("Θ",[l=>lower, νdot=>lower],θbar)*Tensor("Q",[i=>upper,j=>upper,l=>upper,αdot=>upper],x)
# +θij*θbar4*SymmetricTensor("R",[i=>upper, j=>upper],x)
# +θupij*Tensor("θ",[j=>lower, α=>upper],θ)*Tensor("Θ",[k=>upper, αdot=>lower],θbar)*Tensor("Θ",[k=>lower, νdot=>upper],θbar)*Tensor("Θ",[l=>lower, νdot=>lower],θbar)*Tensor("S",[i=>lower,α=>lower,l=>upper,αdot=>upper],x)
# +Tensor("θ",[j=>lower, α=>upper],θ)*θupij*θbar4*Tensor("T",[i=>lower, α=>lower],x)
# +θ4*θbar4*Tensor("V",x)
# +barθαi*Tensor("ϕ",[i=>lower, αdot=>upper],x)
# +barθij*Tensor("A",[i=>upper, j=>upper],x)
# +barθαi*barθij*Tensor("B",[j=>upper, αdot=>upper],x)
# +θbar4*Tensor("C",x)
# +barθαi*Tensor("θ",[j=>upper, ν=>upper],θ)*Tensor("θ",[k=>upper, ν=>lower],θ)*Tensor("U",[j=>lower,k=>lower, i=>lower, αdot=>upper],x)
# +barθαi*Tensor("θ",[k=>lower, β=>upper],θ)*Tensor("θ",[j=>upper, ν=>upper],θ)*Tensor("θ",[k=>upper, ν=>lower],θ)*Tensor("W",[j=>lower,β=>lower, i=>lower, αdot=>upper],x)
# +barθαi*θ4*Tensor("X",[i=>lower, αdot=>upper],x)
# +barθupij*Tensor("θ",[k=>upper, μ=>upper],θ)*Tensor("θ",[l=>upper, μ=>lower],θ)*Tensor("θ",[l=>lower, α=>upper],θ)*Tensor("Y",[k=>lower,α =>lower, i=>lower, j=>lower],x)
# +barθupij*θ4*SymmetricTensor("Z",[i=>lower, j=>lower],x)
# +Tensor("Θ",[j=>upper, αdot=>lower],θbar)*barθij*θ4*Tensor("E",[i=>upper, αdot=>upper],x)
# +θᵦᵧ*SymmetricTensor("Ω",[β=>upper,γ=>upper],x)
# # +θiβ*θᵦᵧ*Tensor("Σ",[i=>upper,γ=>upper],x)
# +barθᵦᵧ*Tensor("Λ",[βdot=>upper,γdot=>upper],x)
# # +barθαi*barθᵅᵝ*Tensor("Ξ",[i=>lower,βdot=>lower],x)
# +θᵦᵧ*barθαi*Tensor("𝒜",[β=>upper, γ=>upper,i=>lower,αdot=>upper],x)
# +θᵦᵧ*barθij*Tensor("ℬ",[β=>upper, γ=>upper,i=>upper,j=>upper],x)
# +θᵦᵧ*barθαi*barθij*Tensor("𝒢",[β=>upper, γ=>upper,j=>upper,αdot=>upper],x)
# +θᵦᵧ*θbar4*SymmetricTensor("𝒞",[β=>upper, γ=>upper],x)
# +θᵦᵧ*barθᵦᵧ*Tensor("ℋ",[β=>upper, γ=>upper,βdot=>upper, γdot=>upper],x)
# # +θᵦᵧ*barθαi*barθᵅᵝ*Tensor("ℱ",[β=>upper, γ=>upper,i=>lower,βdot=>lower],x)
# # +θiβ*θᵦᵧ*barθjα*Tensor("𝒲",[i=>upper, γ=>upper,j=>lower,αdot=>upper],x)
# # +θiβ*θᵦᵧ*barθjk*Tensor("𝒵",[i=>upper, γ=>upper,j=>upper,k=>upper],x)
# # +θiβ*θᵦᵧ*barθjα*barθjk*Tensor("ℰ",[i=>upper, γ=>upper,αdot=>upper,k=>upper],x)
# # +θiβ*θᵦᵧ*θbar4*Tensor("ℛ",[i=>upper, γ=>upper],x)
# # +θiβ*θᵦᵧ*barθᵦᵧ*Tensor("𝒴",[i=>upper, γ=>upper,βdot=>upper,γdot=>upper],x)
# # +θiβ*θᵦᵧ*barθjα*barθᵅᵝ*Tensor("𝒦",[i=>upper, γ=>upper,j=>lower,βdot=>lower],x)
# +barθᵦᵧ*θiα*Tensor("ℒ",[i=>upper, α =>lower,βdot=>upper,γdot=>upper],x)
# +barθᵦᵧ*θij*Tensor("ℳ",[i=>upper, j =>upper,βdot=>upper,γdot=>upper],x)
# +barθᵦᵧ*θjα*θupij*Tensor("𝒩",[i=>lower, α =>lower,βdot=>upper,γdot=>upper],x)
# +barθᵦᵧ*θ4*SymmetricTensor("𝒮",[βdot=>upper,γdot=>upper],x)
# # +barθαi*barθᵅᵝ*θjα*Tensor("𝒫",[j=>upper,α =>lower,i=>lower,βdot=>lower],x)
# # +barθαi*barθᵅᵝ*θjk*Tensor("𝒬",[j=>upper,k =>upper,i=>lower,βdot=>lower],x)
# # +barθαi*barθᵅᵝ*θjα*θʲᵏ*Tensor("𝒳",[k =>lower,α =>lower,i=>lower,βdot=>lower],x)
# # +barθαi*barθᵅᵝ*θ4*Tensor("𝒱",[i=>lower,βdot=>lower],x)
# )



# println(θupij*Tensor("θ",[j=>lower, α=>upper],θ)*Tensor("̄θ",[k=>upper, αdot=>lower],θbar)*Tensor("̄θ",[k=>lower, νdot=>upper],θbar)*Tensor("θ",[l=>lower, νdot=>lower],θbar)*Tensor("S",[i=>upper,α=>lower,l=>upper,αdot=>upper],x))
# break

# covar_expansion=grassman_deriv_normalrep-im*Tensor("σ",[b =>upper, γdot=>upper, β =>lower])*Tensor("Θ",[j=>upper,γdot=>lower,],θbar)*spatial_deriv
# covar=Derivative("D",zero, [x,θ], [j => upper, β => lower])
# println("our superfield expansions is: ", SUPERFIELD)

# #apply our derivative
# println("Computing Symmetrised deriv: ")
# D_SUPERFIELD=Derivative("D",SUPERFIELD, [x,θ], [s => upper, τ => lower])
# out=expand_derivative(D_SUPERFIELD,covar,covar_expansion)
# out=apply_derivative(out)
# out=simplify(out)
# DD_SUPERFIELD=Derivative("D",out,[x,θ], [r => upper, ρ => lower])

# out=expand_derivative(DD_SUPERFIELD,covar,covar_expansion)
# out=apply_derivative(out)
# out=EpsilonTensor([τ=>upper,ρ=>upper])*out

# out=eliminate_epsilon(out)
# out=copy(out)
# D_SUPERFIELD2=Derivative("D",SUPERFIELD, [x,θ], [r => upper, τ => lower])
# out2=expand_derivative(D_SUPERFIELD2,covar,covar_expansion)
# out2=apply_derivative(out2)
# out2=simplify(out2)
# DD_SUPERFIELD2=Derivative("D",out2,[x,θ], [s => upper, ρ => lower])

# out2=expand_derivative(DD_SUPERFIELD2,covar,covar_expansion)
# out2=apply_derivative(out2)
# out2=EpsilonTensor([τ=>upper,ρ=>upper])*out2

# out2=eliminate_epsilon(out2)
# out2=copy(out2)
# println("Result:")
# println(out+out2)
# res=simplify(out+out2)
# println(res)

# t11=get_weighted_terms(get_weighted_terms(res,"θ",1),"Θ",2)
# println("t11 weight:")
# println(merge_sum(t11))

# t11=copy(t11)
# println("rearrange")
# t11=merge_sum(sort(rename_dummies(eliminate_epsilon(rename_dummies(sort(pop_metric(t11)))))))
# println(t11)

# println("elim delta")
# t11=merge_sum(eliminate_kronecker(t11))
# # for summand in t11.summands
# #     println(summand)
# #     println(get_free_indices(summand[2]))
# # end
# println(t11)
# out=merge_sum(sort(pop_metric(t11)))
# println(out)


# println("stripping thetas")
# out=strip_thetas(out)
# println(out)
# # for summand in out.summands
# #         println(get_free_indices(summand[2]))
# # end
# println("constraint:")
# println(copy(eliminate_epsilon(eliminate_kronecker(out))))








# # println("We want to compute: ", DD_SUPERFIELD)
# # println("Appyling the derivative:")
# # out=apply_derivative(DD_SUPERFIELD)
# # println(out)
# # println("Getting rid of epsilon tensors")
# # out=contract_metrics(out,EpsilonTensor)
# # println(out)
# # println("Eliminating deltas")
# # out=eliminate_kronecker(out)
# # println(out)
# # println("Renaming dummies")
# # rename_dummies(out)
# # println(out)
# # println("Simplifying")
# # out=simplify(out)
# # println(out)