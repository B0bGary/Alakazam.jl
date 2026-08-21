using Alakazam, Alakazam.Lightcone
################### Setup ###################
# coordinates
x = Coordinate("x", [lightcone_indices,lightcone_indices])
θ,Θ=Coordinate("θ", lightcone_indices),Coordinate("̄θ", lightcone_indices)
θ⁺,Θ⁺=Tensor("θ",ℓ⁺=>↑,θ),Tensor("Θ",ℓ⁺=>↑,Θ)
# basic partial derivatives
∂₊₊,∂₋₋= PD("∇",x,ℓ⁺=>↓,ℓ⁺=>↓),PD("∇",x,ℓ⁻=>↓,ℓ⁻=>↓)
∂θ,∂Θ=PD("∂",θ,ℓ⁺=>↓),PD("𝛛",Θ,ℓ⁺=>↓)
@action ∂θ(θ⁺) = 1; @action ∂Θ(Θ⁺) = 1
# covariant derivatives
D,Dbar=Derivative("D",[x, θ],ℓ⁺=>↓),Derivative("̄D",[x, Θ],ℓ⁺=>↓)
@covariant D=∂θ-im//2*Θ⁺*∂₊₊;@covariant Dbar=-∂Θ+im//2*θ⁺*∂₊₊
superspace = Superspace([θ⁺, Θ⁺], [D, Dbar])
# generate the superfields
sf(n,inds)=generate_superfield_naive(n,inds,superspace,x)
S,T=sf("𝒮", [ℓ⁺=>↓,ℓ⁺=>↓]),sf("𝒯", [(ℓ⁻ => ↓) for _ in 1:4])
W,Wb = sf("𝒲", [ℓ⁻=>↓]),sf("̄𝒲", [ℓ⁻=>↓])
# coefficient of θᵃΘᵇ in the expansion of the superfield F
sfcoef(F,a,b)=isTensor(grassmann_coefficient(F,superspace,[a,b]))
# substitute target => value rules via our replacement pattern
apply_rules(e, rules) = foldl((a, r) -> replace_pattern(a,
    TensorTerm(r[1]), r[2]),rules; init = TensorExpression(e))

############### Defining the Problem ###############
# The constraints, each = 0, on the set F = {𝒮, 𝒯, 𝒲, ̄𝒲}
constraints = [(F)->Dbar(F[3]),(F)->D(F[4]),(F)->D(F[2])-1//2*∂₋₋(F[4]),
    (F)->Dbar(F[2])-1//2*∂₋₋(F[3]),(F)->∂₋₋(F[1])-D(F[3])+Dbar(F[4]),]
# Helper for ith constraint for each superfield at order θʲΘᵏ
expand(i,j,k)=simplify(derive(apply_rules(grassmann_coefficient(
    derive(constraints[i]([apply_rules(F, sols) 
    for F in (S, T, W, Wb)]),superspace),superspace, j, k),diffs)))

################### Solving ###################
# Brute force solve the constraints wherever a component appears
sols, diffs = Pair[], Pair[]
for F∈(S,T,W,Wb),i∈1:length(constraints),j∈0:1,k∈0:1,
    ∂∈(nothing, ∂₊₊, ∂₋₋),a ∈ 0:1,b∈0:1
        target = ∂===nothing ? sfcoef(F,a,b) : ∂(sfcoef(F,a,b))
        (comp,val)=solve_linear(expand(i, j, k), target)
        comp===nothing || push!(∂===nothing ? sols : diffs,comp=>val)
end
# simplify constraints, drop redundant ones
diffs=[r for r∈(t=>simplify(derive(apply_rules(v,sols))) 
    for (t,v)∈diffs) if !isZeroTens(simplify(r[2]-r[1]))]
# apply the constraints, then integrate over superspace
Ss, Ts, Ws, Wbs = (apply_rules(F, sols) for F ∈ (S, T, W, Wb))
𝒪 = ∫(Ts * Ss - Wbs * Ws, superspace)
O = reduce_bp_with(𝒪, diffs)   # by parts, drop total derivs

################### Display ###################
println("solved components (", length(sols), ")")
foreach(r -> println("   ", r[1], " = ", r[2]), sols)
println("relations (", length(diffs), ")")
foreach(r -> println("   ", r[1], " = ", r[2]), diffs)
println("\nconstrained superfields")
foreach((nm, F) -> println("  ", rpad(nm, 6), " = ", F),
    ("𝒮₊₊", "𝒯₋₋₋₋", "𝒲₋", "̄𝒲₋"), (Ss, Ts, Ws, Wbs))
@show 𝒪
@show O