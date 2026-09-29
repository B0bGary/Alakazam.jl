# Alakazam

<div align="center">
  <img src="https://gitlab.com/B0bGary/alakazam.jl/-/raw/main/newlogoalakazamtest.png" width="400" height="400" alt="Alakazam!">
</div>

Symbolic tensor algebra for field theory, in Julia.

Alakazam is designed to perform symbolic algebraic manipulations of tensors and tensor equations, with field theories and quantum mechanics in mind. It offers notable speed improvements over analogous operations in other packages like Cadabra and xAct.
In particular, Alakazam features the ability to raise and lower tensor indices with a metric, simplify expressions using symmetries of tensors, and perform covariant differentiation.
For those working on supersymmetric theories, Alakazam supports grading natively, and can handle covariant derivatives with multiple indices, perform superspace integration, and more. 


## Getting started

To get started import the package after installing by running
```julia
using Alakazam
```
The simplest thing you can do is create a tensor,
```julia
T=Tensor("T")
```
This creates a tensor named T. This however, isn't particularly useful. To create a tensor with indices, we can firstly create an IndexSet,
```julia
spacetime_indices = IndexSet("spacetime", 4)
```
Then, we can create some index objects
```julia
@indices μ:τ∈spacetime_indices
```
Now we can create more interesting tensor objects,
```julia
T = Tensor("T", μ => upper,  ν=> lower, τ=>upper)
    # T^{μ∘τ}_{∘ν∘}
```
Alakazam can also parse tensors in a variety of ways, including from LaTeX, using macros
```julia
T = @Tensor T_μν^τ
T = T"T_μν^τ"
T = LT"T_\mu\nu^\tau"
    # T^{μ∘τ}_{∘ν∘}
```

Arbitrary symmetries are supported by Young tableaux, and are recognised by simplification algorithms
```julia
riemann_tableau=YoungTableau([2,2],[1,3,2,4])
@show riemann_tableau
    # ┌───┬───┐
    # │ 1 │ 3 │
    # ├───┼───┤
    # │ 2 │ 4 │
    # └───┴───┘
R(i,j,k,l) = Tensor("R", i=>lower, j=>lower, k=>lower, l=>lower,riemann_tableau)
bianchi = R(μ,ν,σ,τ) + R(μ,σ,τ,ν) + R(μ,τ,ν,σ)
@show hard_simplify(bianchi; multiterm=true)
    #   0
```

Non-commuting matrix objects are also supported, with commutation relations
```julia
X=Matrix("X",μ=>upper)
P=Matrix("P",ν=>upper)
@commutator [X, P] = im * Metric(μ=>upper, ν=>upper)
simplify(X*P - P*X)                 # iη^{μ,ν}
```


There is also pattern matching with wildcard indices, Lie algebras, traces of gamma matrices, superspace, a lightcone module, and serialisation...

... and much more!


## Going further

`demo.jl` has a short suite of features:

```
julia --project=. demo.jl
```

The `demo/` folder has complete, non-trivial calculations:

- **`demo/compton_spinor.jl`** — Compton scattering end to end, with matrix reduction rules and substitutions
- **`demo/ttbar2DN02.jl`** — an 𝒩=(0,2) superspace constraint analysis, with covariant derivatives, superfield expansions and Grassmann coefficients.

For a more comprehensive description of these demos, see the paper [_"Alakazam: Symbolic Tensor Algebra for Field Theory in Julia"_](https://arxiv.org/abs/2608.20452).

## Support

If you'd like to support me and future development of this project, please donate [here](https://www.paypal.com/donate/?business=L23N9VSVR9KB4&no_recurring=0&item_name=Thanks+for+supporting+the+development+of+great+tools%21&currency_code=USD).

## License
Alakazam is available under the conditions of the GNU General Public License. 
If you use this software in part of your work, please cite: [_"Alakazam: Symbolic Tensor Algebra for Field Theory in Julia"_](https://arxiv.org/abs/2608.20452), and [_"Decomposing Grassmann Monomials for Superfield Expansions"_](https://arxiv.org/abs/2608.20455) if using the `generate_superfield` function.

## Contributing
If you'd like to contribute, or have any questions in general, please get in touch at <jesse.woods@unibe.ch> :)

## Authors and acknowledgment
Written by me.
Many thanks to Cian Luke Martin for the suggestion of macros during the early stages of development. 
Many thanks to Alejandro Soto Franco for several useful merge requests. 
Credits to Kasper Peeters, author of Cadabra, whose project inspired me to create Alakazam.
