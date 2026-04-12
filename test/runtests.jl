using Alakazam
using Test
using Symbolics
zero=Alakazam.zero
# To run the tests in a terminal, type:
# julia --project=. test/runtests.jl

# include("test_new_macros.jl")



include("test_arithmetic.jl")
include("test_tableaux.jl")


#include("test_tensor_term.jl")
# include("test_tterm_equivs.jl")

# include("test_opsimps.jl")
# include("test_thesis.jl")
# ...

return