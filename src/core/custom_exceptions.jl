

struct FreeIndexMismatch <: Exception
    A::Any
    B::Any
    FreeIndexMismatch() = new(nothing, nothing)
    FreeIndexMismatch(A, B) = new(A, B)
end
Base.showerror(io::IO, e::FreeIndexMismatch) = begin
    print(io, "FreeIndexMismatch: terms in a sum must carry the same free indices")
    if e.A !== nothing
        print(io, "\n  first term:  ", e.A, "   free: ", free_indices(e.A))
        print(io, "\n  second term: ", e.B, "   free: ", free_indices(e.B))
    end
end

struct OverusedIndices <: Exception
    T1::Vector{Tensorial}
    OverusedIndices(T1::Tensorial,T2::Tensorial)=new([T1,T2])
    OverusedIndices(T1::Tensorial)=new([T1])
end

Base.showerror(io::IO, e::OverusedIndices) = begin
    print(io, "OverusedIndices: an index appears more than twice, so it cannot be a valid dummy pair")
    for t in e.T1
        print(io, "\n  ", t)
    end
end

struct InvalidKroneckerIndices <: Exception
    inds::Any
end
Base.showerror(io::IO, e::InvalidKroneckerIndices) = begin
    println(io,"Invalid index configuration:",e.inds)
end

struct InvalidEpsilonIndices <: Exception
    inds::Any
    reason::String
    InvalidEpsilonIndices() = new(nothing, "")
    InvalidEpsilonIndices(inds, reason::String) = new(inds, reason)
end

Base.showerror(io::IO, e::InvalidEpsilonIndices) = begin
    print(io, "InvalidEpsilonIndices")
    isempty(e.reason) || print(io, ": ", e.reason)
    e.inds === nothing || print(io, "\n  given indices: ", e.inds)
end

struct OperatorPrecedenceError <: Exception
end

struct OutOfIndices <: Exception
    I::IndexSet
end
Base.showerror(io::IO, e::OutOfIndices) = begin
    print(io, "OutOfIndices: index set '", e.I.name, "' has no unused indices left.",
              "\n  Auto-generation of new indices failed (out of text characters). Index set size: ",length(e.I.indices),".")
end

struct InvalidReplacementPattern <: Exception
    reason::String
    detail::Any
    InvalidReplacementPattern() = new("", nothing)
    InvalidReplacementPattern(reason::String) = new(reason, nothing)
    InvalidReplacementPattern(reason::String, detail) = new(reason, detail)
end

Base.showerror(io::IO, e::InvalidReplacementPattern) = begin
    print(io, "InvalidReplacementPattern")
    isempty(e.reason) || print(io, ": ", e.reason)
    e.detail === nothing || print(io, "\n  ", e.detail)
end

struct MixedWeightOperandTerms <: Exception
end

struct SameLevelDummyError <: Exception
    T1::Union{Tensorial,Vector{IndexPair}}
    T2::Union{Tensorial,Nothing}
    SameLevelDummyError(T1)=new(T1,nothing)
    SameLevelDummyError(T1,T2)=new(T1,T2)
end
Base.showerror(io::IO, e::SameLevelDummyError) = begin
    print(io, "SameLevelDummyError: a repeated index appears at the same position (both raised or both lowered);",
              " a dummy pair must have one raised and one lowered index to obey Einstein summation convention.")
    print(io, "\n  ", e.T1)
    e.T2 === nothing || print(io, "\n  ", e.T2)
end

struct InvalidDerivativeBasis <: Exception
    D::Any
    exp::Any
    InvalidDerivativeBasis() = new(nothing, nothing)
    InvalidDerivativeBasis(D, exp) = new(D, exp)
end

Base.showerror(io::IO, e::InvalidDerivativeBasis) = begin
    print(io, "InvalidDerivativeBasis: a derivative action must be declared on a single tensor, ",
              "not on a product or sum.")
    if e.exp !== nothing
        print(io, "\n  derivative: ", e.D)
        print(io, "\n  given:      ", e.exp)
    end
end

struct InvalidDifferentialForm <: Exception
    T::String
end
Base.showerror(io::IO, e::InvalidDifferentialForm) = begin
    println(io, e.T)
end

struct InvalidCommutatorDefinition <: Exception
    A::Any
    B::Any
    C::Any
end

Base.showerror(io::IO, e::InvalidCommutatorDefinition) = begin
    print(io, "InvalidCommutatorDefinition: the result of a commutator must carry the same ",
              "free indices as the product of its operands, unless it vanishes.")
    print(io, "\n  [", name(e.A), ", ", name(e.B), "] was given as: ", e.C)
    print(io, "\n  free indices of ", name(e.A), "*", name(e.B), ": ", free_indices(e.A * e.B))
    print(io, "\n  free indices of result:            ", free_indices(e.C))
end