module Lightcone

using ..Alakazam

const lightcone_indices = IndexSet("lightcone",2, lower, fermionic)
const ℓ⁺⁻, ℓ⁻⁺, ℓ⁺, ℓ⁻ = Index(["±", "∓", "+", "-"], lightcone_indices)

export ℓ⁺⁻, ℓ⁻⁺, ℓ⁺, ℓ⁻, lightcone_indices,expand_lc


"""
Expands ± and ∓ into + and -
"""
function expand_lc(T::Tensorial)
    out1=copy(T)
    out2=copy(T)
    set_indices_deep!(out1, Dict((ℓ⁺⁻ => upper)=>(ℓ⁺=>upper),(ℓ⁺⁻ => lower)=>(ℓ⁺=>lower),(ℓ⁻⁺ => upper)=>(ℓ⁻=>upper),(ℓ⁻⁺ => lower)=>(ℓ⁻=>lower)))
    set_indices_deep!(out2, Dict((ℓ⁺⁻ => upper)=>(ℓ⁻=>upper),(ℓ⁺⁻ => lower)=>(ℓ⁻=>lower),(ℓ⁻⁺ => upper)=>(ℓ⁺=>upper),(ℓ⁻⁺ => lower)=>(ℓ⁺=>lower)))
    return out1+out2
end




end

