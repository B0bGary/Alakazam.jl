using Alakazam, Symbolics
st = IndexSet("spacetime", 4, upper)
@indices μ:ν ∈ st ; @indices a:d ∈ st ; @indices f:k ∈ st
@syms m e s1 s2   #s1 = k₁·k₃, s2 = k₂·k₃ (scalar invars)
γ(i, p=upper) = Gamma(i => p)
slash(n,i) = Tensor("k",i=>upper, "$n" => lower)*γ(i,lower)# k̸ₙ
Id = IdMatrix(spinor_indexset(st))# the spinor space tr(𝟙) runs over, not st itself
# s and u channel vertex structs, denoms (k₃+k₁)²-m² = 2s1 and (k₃-k₂)²-m² = -2s2
Γ=-im*e^2*(γ(ν)*(slash(3, a)+slash(1, b) + m*Id)*γ(μ)*(1//2) / s1 -
              γ(μ)*(slash(3, c)-slash(2, d) + m*Id)*γ(ν)*(1//2) / s2)
Γbar=im*e^2*(γ(μ,lower)*(slash(3, f)+slash(1, g)+m*Id)*γ(ν,lower)*(1//2) / s1 -
              γ(ν,lower)*(slash(3, h)-slash(2, i)+m*Id)*γ(μ,lower)*(1//2) / s2)

kernel = Trace((slash(4, j) + m * Id) * Γ * (slash(3, k) + m * Id) * Γbar)
traced = simplify(eliminate_trace(kernel))

# k₄ = k₁ + k₃ - k₂, in both index positions
onshell = replace_pattern(traced,  T"k_4^?w", T"k_1^?w" + T"k_3^?w" - T"k_2^?w")
onshell = replace_pattern(onshell, T"k_4_?w", T"k_1_?w" + T"k_3_?w" - T"k_2_?w")

# replacement rules for the scalar invariants # kₙ, wildcard index helper
kx(n, p) = Tensor("k", [DummyPatternIndex("x", st) => p], "$n" => lower)
prods = ((1,1,0), (2,2,0), (3,3,m^2),# k₁² = k₂² = 0,  k₃² = m²
       (1,3,s1), (2,3,s2), (1,2,s1-s2))# k₁·k₃ = s1, k₂·k₃ = s2, k₁·k₂ = s1-s2
rules = [(kx(x,p)*kx(y,op_pos(p)), r) for (x,y,r)∈prods for p∈(lower,upper)]
# replace momenta by the scalar invariants
for  (pat, rep) in rules
    global onshell = replace_pattern(onshell, pat, rep)
end
length(summands(onshell))==1 || error("Simplify failed.") ; 
M2=simplify(onshell)
println("¼Σ|M|²      = ", simplify(expand(first_coef(1/4*M2))))
kn = 2e^4*(s1/s2 + s2/s1 + 2m^2*(1/s1-1/s2) + m^4*(1/s1-1/s2)^2) #from PS (5.87) 
println("Our form- Klein-Nishina = ", simplify(expand(first_coef(1/4*M2) - kn)))