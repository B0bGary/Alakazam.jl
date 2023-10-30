#add the package
include("alakazam.jl")



#create some index sets
flavour_indices = IndexSet("flavour", 2)
spinor_indices = IndexSet("spinor", 2)
spacetime_indices = IndexSet("spacetime", 4)


# create some index objects
i = Index("i", flavour_indices, 0)
j = Index("j", flavour_indices)
k = Index("k", flavour_indices, 0)
l = Index("l", flavour_indices)

α = Index("α", spinor_indices, 1)
β = Index("β", spinor_indices, 1)
γ = Index("γ", spinor_indices, 1)

μ = Index("μ", spacetime_indices)
ν = Index("ν", spacetime_indices)
ρ = Index("ρ", spacetime_indices)
τ = Index("τ", spacetime_indices)


# create some coordinates (for taking deriviatves w/r to, etc)
x = Coordinate("x", spacetime_indices)
θ = Coordinate("θ", [flavour_indices, spinor_indices])



J = Tensor("J", [i=>upper, j=>upper])
K = Tensor("K", [j=>upper,k=>lower])

# exp=J+K

η = Metric("η", [k=>upper, i=>upper])

exp = J+ η*K

exp=contract_metrics(exp)


