# using Alakazam
# using Test
# using Symbolics

# # To run the tests in a terminal, type:
# # julia --project=. test/runtests.jl
# println("initialised.")
# @testset "Alakazam Tests" begin    
#     #basic arithmetic tests, to make sure basic manipulation is working
#     include("test_arithmetic.jl")
#     #test tensor properties, metrics, eps, deltas etc
#     include("test_basic.jl")
#     #test replacement patterns
#     include("test_replacements.jl")
#     #test DG, appension etc
#     include("test_arithmetic_adv.jl")
#     # derivatives and calculus
#     include("test_derivatives.jl")
#     # properties involving YT
#     include("test_tableaux.jl")
#     # tests of tensor term simplifications
#     include("test_tterm_equivs.jl")
#     # tests of nested derivs, terms, and contractions with susy objects etc.
#     include("test_commuting.jl")
#     # test of eps identities
#     include("test_epsdelta.jl")
#     # tests matrix and commutator relations and lie algebra
#     include("test_matrices.jl")
#     # lie bracket explicit
#     include("test_lie.jl")

#     include("test_serialise.jl")

#     include("test_serialise_load.jl")

#     include("test_algebras.jl")
#     include("test_variational.jl")

#     include("test_new_macros.jl")
#     include("test_superspace.jl")
#         include("test_superspace_adv.jl")
# include("test_lightcone.jl")

# end;




# # include("test_superspace.jl")
# # include("test_graph.jl")
# #include("test_tensor_term.jl")
# # include("test_tterm_equivs.jl")
# # include("test_opsimps.jl")
# # include("test_thesis.jl")
# # ...
# return


# using Symbolics,Test,Alakazam
# using SnoopCompileCore
# invs = @snoop_invalidations (using Symbolics, Test, Alakazam);
# using SnoopCompile
# trees = invalidation_trees(invs)
# println("total: ", length(uinvalidated(invs)))
# for t in sort(trees, by=SnoopCompile.countchildren)[max(1, end-10):end]
#     println(SnoopCompile.countchildren(t), "  ", t.method.module, ".", t.method.name)
# end
using #=Symbolics,=# Test, Alakazam, Alakazam.Lightcone
if true
    const TEST_FILES = [
        "test_arithmetic.jl",
        "test_basic.jl",
        "test_replacements.jl",
        "test_arithmetic_adv.jl",
        "test_derivatives.jl",
        "test_tableaux.jl",
        "test_tterm_equivs.jl",
        "test_commuting.jl",
        "test_epsdelta.jl",
        "test_matrices.jl",
        "test_lie.jl",
        "test_serialise.jl",
        "test_serialise_load.jl",
        "test_algebras.jl",
        "test_variational.jl",
        "test_new_macros.jl",
        "test_superspace.jl",
        "test_superspace_adv.jl",
        "test_lightcone.jl",
    ]

    @testset verbose=true "Alakazam Tests" begin
        for f in TEST_FILES
            print(rpad(f, 32))
            @time @testset "$f" begin
                include(f)
            end
        end
    end
    # using Symbolics, Test, Alakazam
    # for f in TEST_FILES
    #     print(rpad(f, 32))
    #     @time include(f)          # no @testset wrapping
    # end
end

# using Profile
# @profview include("test_superspace_adv.jl")        # where time goes
# @profview_allocs include("test_superspace_adv.jl") # where allocations come from