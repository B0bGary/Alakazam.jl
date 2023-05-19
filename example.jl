include("alakazam.jl")
spacetime_indices = IndexSet("spacetime", 4,lower)
# lightcone_indices = IndexSet("lightcone", 2,lower)
#     ± = Index("±", lightcone_indices, 1)
#     ∓ = Index("∓", lightcone_indices, 1)
#     p = Index("+", lightcone_indices, 1)
#     m = Index("-", lightcone_indices, 1)

#     a = Index("a", spacetime_indices, 1)
#     b = Index("b", spacetime_indices,1)
#     c = Index("c", spacetime_indices,1)
    
#     x = Coordinate("x", [spacetime_indices, spacetime_indices])
#     θ = Coordinate("θ", [spacetime_indices])

    
#     field=Tensor("θ",[c=>upper],θ)*Tensor("S",[c => lower, a =>lower, b => lower], x)
#     Dfield=Derivative("∇",field, [x], [i => lower, j => lower])
#     out=apply_derivative(Dfield)
#     println(out)
#     res=to_lc_expression(out,Dict(a=>m,b=>m,i=>p, j=> p))


    x = Coordinate("x", [spacetime_indices])
    X=Tensor("X",x)
    exp= 2*X*X+8*X+9*X*X*X
    get_weighted_terms(exp,"x",3)
