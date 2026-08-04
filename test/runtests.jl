
using Symbolics, Test, Alakazam, Alakazam.Lightcone
if true
    const TEST_FILES = [
        "test_arithmetic.jl",
        "test_basic.jl",
        "test_indices.jl",
        "test_replacements.jl",
        "test_arithmetic_adv.jl",
        "test_derivatives.jl",
        "test_parts.jl",
        "test_tableaux.jl",
        "test_tterm_equivs.jl",
        "test_commuting.jl",
        "test_epsdelta.jl",
        "test_identities.jl",
        "test_matrices.jl",
        "test_lie.jl",
        "test_trace.jl",
        "test_serialise.jl",
        "test_serialise_load.jl",
        "test_algebras.jl",
        "test_variational.jl",
        "test_new_macros.jl",
        "test_superspace.jl",
        "test_superspace_adv.jl",
        "test_lightcone.jl",
        "test_n02.jl",
        "test_forms.jl",
        "test_exterior.jl",
        "test_compton.jl",
    ]
    _TESTDIR = mktempdir(cleanup=false)
    _TESTFILE1 = joinpath(_TESTDIR, "test_expression")
    _TESTFILE2 = joinpath(_TESTDIR, "test_expression_headers")
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