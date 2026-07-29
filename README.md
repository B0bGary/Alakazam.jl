# Alakazam

<div align="center">
  <img src="newlogo3.png" width="400" height="400">
</div>

Alakazam is a Julia-based solution to performing symbolic algebraic manipulations of tensors and tensor equations.
It was designed primarily for studying field theories in physics.
In particular, it includes the ability to raise and lower tensor indices with a metric, simplify expressions using symmetries of tensors, and perform covariant differentiation.

For those working on supersymmetric theories, Alakazam can handle covariant derivatives with multiple indices, perform superspace integration, and more.

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
spacetime_indices = IndexSet("spacetime", 4,lower)
```
Then, we can create some index objects
```julia
μ,ν,ρ,τ = Index(["μ","ν","ρ","τ"], spacetime_indices)
```
Now we can create more interesting tensor objects,
```julia
T = Tensor("T", [μ => upper,  ν=> lower, τ=>upper])
    # T^{μ∘τ}_{∘ν∘}
```
We can now begin to perform some basic operations, for example raising and lowering indices with a metric,
```julia
η_νρ = Metric([ν=>upper, ρ=>upper])
expression = η_νρ*T
contract_metrics(expression)
    # T^{μρτ}
```
More important functionality includes taking covariant derivatves, pattern matching and replacements in expressions, and simplifications using symmetries of tensors.
There is also inbuilt support for Lie algebras, matrices, supersymmetric/fermionic objects, and some aspects of superspace.

This section will be updated when an upcoming paper on this package is published.

For more info, see demo.jl


## Support

If you'd like to support me and future development of this project, please donate [here](https://www.paypal.com/donate/?business=L23N9VSVR9KB4&no_recurring=0&item_name=Thanks+for+supporting+the+development+of+great+tools%21&currency_code=USD).

## License
Alakazam is available under the conditions of the GNU General Public License. 
If you use this software in part of your work, please cite: (publication coming soon)

## Contributing
If you'd like to contribute, or have any questions in general, please get in touch at <jesse.woods@unibe.ch> :)

## Authors and acknowledgment
Written by me.
Credits to Kasper Peeters, author of Cadabra, whose project inspired me to create Alakazam.

