using InteractiveUtils: subtypes

function _concretes(T)
    out = Type[]
    for S in subtypes(T)
        isconcretetype(S) && push!(out, S)
        isabstracttype(S) && append!(out, _concretes(S))
    end
    unique(out)
end


module _PC
using Symbolics
using Alakazam
using Alakazam.Lightcone

macro test(ex)
    ;
    esc(ex);
end
macro test_throws(T, ex)
    ;
    esc(:(
        try
            $ex
        catch
        end
    ));
end
macro testset(args...)
    ;
    esc(args[end]);
end
macro test_nowarn(ex)
    ;
    esc(ex);
end
macro test_warn(msg, ex)
    ;
    esc(ex);
end
macro test_logs(args...)
    ;
    esc(args[end]);
end
macro test_broken(ex)
    ;
    esc(:(
        try
            $ex
        catch
        end
    ));
end
macro test_skip(ex)
    ;
    nothing;
end
macro inferred(ex)
    ;
    esc(ex);
end

const _DIR = joinpath(pkgdir(Alakazam), "test")

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

_TESTDIR = mktempdir(cleanup=false)
_TESTFILE1 = joinpath(_TESTDIR, "test_expression")
_TESTFILE2 = joinpath(_TESTDIR, "test_expression_headers")
for f in TEST_FILES
    include(joinpath(_DIR, f))
end
end


let
    precompile(Base.deepcopy_internal, (IndexSet, IdDict))
    fi = IndexSet("fi", 2, upper);
    @indices a:c ∈ fi;
    x = Coordinate("x", fi)
    T = Tensor("T", [a=>lower, b=>lower, c=>lower], x)
    S2 = Tensor("S", [a=>upper, b=>upper, c=>upper], x)
    T * S2
    for A in [_concretes(TensorSuperType);]
        precompile(TensorTerm, (Vector{A},))
    end

    precompile(Tensor, (String, Coordinate))

    TENS = [_concretes(TensorSuperType); TensorTerm; TensorExpression]
    for A in TENS
        precompile(copy, (A,))
        precompile(deepcopy, (A,))
        precompile(TensorTerm, (Vector{A},))
        precompile(-, (A,))
        precompile(/, (A, SymorNum))
        precompile(//, (A, SymorNum))
        precompile(+, (SymorNum, A))
        precompile(+, (A, SymorNum))
        precompile(simplify, (A,))
        for B in TENS
            precompile(*, (A, B))
            precompile(+, (A, B))
            precompile(-, (A, B))
        end
    end
    precompile(TensorTerm, (Vector{TensorSuperType},))
end

let
    for T in (Pair{IndexSuperType,IndexPosition}, Pair{String,IndexPosition},
        Pair{Index,IndexPosition}, Pair{DummyPatternIndex,IndexPosition},
        Coordinate, IndexSet, YoungTableau, TensorSuperType, TensorTerm)
        precompile(Memory{T}, (UndefInitializer, Int))
        precompile(Vector{T}, (UndefInitializer, Int))
    end
    for I in (Index, DummyPatternIndex)
        precompile(Pair, (I, IndexPosition))
    end
    precompile(Pair, (String, IndexPosition))
end
let
    for I in (Index, DummyPatternIndex)
        precompile(Pair, (I, IndexPosition))
        precompile(Memory{Pair{I,IndexPosition}}, (UndefInitializer, Int))
        precompile(Vector{Pair{I,IndexPosition}}, (UndefInitializer, Int))
    end
    for T in (Pair{IndexSuperType,IndexPosition}, Pair{String,IndexPosition},
        Coordinate, IndexSet, YoungTableau, TensorSuperType)
        precompile(Memory{T}, (UndefInitializer, Int))
        precompile(Vector{T}, (UndefInitializer, Int))
    end
    precompile(Pair, (String, IndexPosition))
end

let
    IDX = (Index, DummyPatternIndex)
    # element types both appear in your inference dump — cover both
    PR = (Pair{Index,IndexPosition}, Pair{DummyPatternIndex,IndexPosition},
        Pair{IndexSuperType,IndexPosition})
    IOS = (IOBuffer, IOContext{IOBuffer}, Base.TTY, IOContext{Base.TTY})

    # ── construction ──────────────────────────────────
    precompile(IndexSet, (String,))
    precompile(IndexSet, (String, Int))
    precompile(IndexSet, (String, Int, IndexPosition))
    precompile(IndexSet, (String, Int, IndexGrading))
    precompile(IndexSet, (String, Int, IndexPosition, IndexGrading))
    precompile(IndexSet, (String, Tuple{Int64,Int64}))
    precompile(IndexSet, (String, Int, Type{Metric}))
    precompile(IndexSet, (String, Int, IndexPosition, Type{Metric}))

    precompile(Index, (String,))
    precompile(Index, (String, IndexSet))
    precompile(Index, (Char, IndexSet))
    precompile(Index, (Vector{Char}, IndexSet))      # ← the @indices path
    precompile(Index, (Vector{String}, IndexSet))
    precompile(valid_name, (String,))

    precompile(DummyPatternIndex, (String,))
    precompile(DummyPatternIndex, (String, IndexSet))
    precompile(DummyPatternIndex, (Vector{String}, IndexSet))

    precompile(Coordinate, (String, IndexSet))
    precompile(Coordinate, (String, Vector{IndexSet}))

    # ── scalar accessors ──────────────────────────────
    for I in IDX
        precompile(grading, (I,))
        precompile(det, (I,))
        precompile(name, (I,))
        precompile(dimension, (I,))
        precompile(isLightcone, (I,))
        precompile(==, (I, I))
        precompile(Base.hash, (I, UInt))
        precompile(Base.isless, (I, I))
        precompile(string, (I,))
    end

    precompile(==, (IndexSet, IndexSet))
    precompile(==, (IndexPosition, IndexPosition))
    precompile(==, (Coordinate, Coordinate))
    precompile(Base.hash, (IndexSet, UInt))
    precompile(Base.hash, (IndexPosition, UInt))
    precompile(Base.hash, (Coordinate, UInt))
    precompile(Base.isless, (IndexPosition, IndexPosition))

    # ── pair / position machinery ─────────────────────
    precompile(op_pos, (IndexPosition,))
    for P in PR
        precompile(op_pos, (P,))
        precompile(op_pos, (Vector{P},))
        precompile(op_pos, (OrderedSet{P},))
        precompile(in_default_pos, (P,))
        precompile(dimension, (P,))
        precompile(isLightcone, (P,))
        precompile(Base.isless, (P, P))
        precompile(sort!, (Vector{P},))
        precompile(copy, (Vector{P},))
    end

    # ── index-set operations ──────────────────────────
    precompile(get_index, (String, IndexSet))
    precompile(get_index, (String, String))
    precompile(generate_indices, (IndexSet,))
    precompile(safe_pop!, (IndexSet, Int))
    for I in IDX
        precompile(safe_pop!, (OrderedSet{I}, IndexSet))
    end
    precompile(safe_pop!, (OrderedSet{IndexSuperType}, IndexSet))
    precompile(clear_coords!, ())
    precompile(clear_indexsets!, ())
    precompile(Base.deepcopy_internal, (IndexSet, IdDict{Any,Any}))

    # ── printing ──────────────────────────────────────
    for IO_ in IOS, I in IDX
        precompile(Base.show, (IO_, I))
    end
    for IO_ in IOS
        precompile(Base.show, (IO_, IndexSet))
        precompile(Base.show, (IO_, IndexPosition))
        precompile(Base.show, (IO_, Coordinate))
        precompile(Base.show, (IO_, Vector{Coordinate}))
    end

end


let
    TENS = (Tensor, Fermion, Boson, DummyPatternTensor)
    IDX = (Index, DummyPatternIndex)
    IPV = (Vector{IndexPair}, Vector{Pair{Index,IndexPosition}},
        Vector{Pair{DummyPatternIndex,IndexPosition}})
    NUMS = (Int, Rational{Int}, Rational{BigInt}, Float64, Complex{Int})

    #TensorData ────────────────────────────────────
    precompile(TensorData, (String,))
    for V in IPV
        precompile(TensorData, (String, V))
        precompile(TensorData, (String, V, Coordinate))
        precompile(TensorData, (String, V, Vector{Any}, Vector{Coordinate}))  # @Tensor macro shape
        precompile(TensorData, (String, V, Vector{YoungTableau}, Vector{Coordinate}))
        # precompile(_absorb!, (V, Vector{IndexPair}, OrderedSet{Coordinate},
        #                       Dict{String,Number}, OrderedSet{YoungTableau}, Vector{LabelPair}))
    end

    # leaf tensors ──────────────────────────────────
    for T in TENS
        precompile(T, (TensorData,))
        precompile(T, (String,))
        for V in IPV
            precompile(T, (String, V))
            precompile(T, (String, V, Coordinate))
        end
    end
    for V in IPV
        precompile(SymmetricTensor, (String, V))
        precompile(AntisymmetricTensor, (String, V))
        precompile(RiemannTensor, (String, V))
    end

    # ── terms ─────────────────────────────────────────
    precompile(TensorTerm, (Vector{TensorSuperType},))
    for T in TENS
        precompile(TensorTerm, (Vector{T},))
        precompile(TensorTerm, (T,))
    end
    precompile(TensorTerm, (TensorTerm,))

    # ── expressions ───────────────────────────────────
    precompile(TensorExpression, ())
    precompile(TensorExpression, (Vector{CoeffPair},))
    precompile(TensorExpression, (TensorTerm,))
    precompile(TensorExpression, (TensorExpression,))
    precompile(TensorExpression, (Vector{TensorTerm},))
    precompile(TensorExpression, (TensorExpression, TensorExpression))
    for T in TENS
        precompile(TensorExpression, (T,))
        precompile(TensorExpression, (Vector{T},))
        precompile(TensorExpression, (TensorExpression, T))
        for N in NUMS
            precompile(TensorExpression, (N, T))
            precompile(TensorExpression, (Pair{N,T},))
            precompile(TensorExpression, (Vector{Pair{N,T}},))
        end
    end
    for N in NUMS
        precompile(TensorExpression, (N,))
        precompile(TensorExpression, (Pair{N,TensorTerm},))
        precompile(TensorExpression, (Vector{Pair{N,TensorTerm}},))
    end
    precompile(merge_sum, (TensorExpression,))

    # ── tableaux ──────────────────────────────────────
    precompile(generate_symmetric_tableau, (Vector{IndexPair},))
    precompile(generate_antisymmetric_tableau, (Vector{IndexPair},))
    precompile(generate_symmetric_tableau, (Int,))
    precompile(generate_antisymmetric_tableau, (Int,))
    precompile(generate_riemann_tableau, ())
    precompile(reduce_tableaux, (Set{YoungTableau},))
    precompile(Base.push!, (Set{YoungTableau}, YoungTableau))
end
clear_global!()
