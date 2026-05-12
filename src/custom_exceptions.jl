

struct FreeIndexMismatch <: Exception
    #println("Free indices in terms are different.")
end


struct OverusedIndices <: Exception
    #println("One or more free indices in the expression is already used as a dummy index pair.")
end


struct InvalidKroneckerIndices <: Exception
    inds::Any
end


struct InvalidEpsilonIndices <: Exception
end

struct OperatorPrecedenceError <: Exception
end

struct OutOfIndices <: Exception
    I::IndexSet
end

struct InvalidReplacementPattern <: Exception
end

struct MixedWeightOperandTerms <: Exception
end

struct SameLevelDummieError <: Exception
end

struct InvalidDerivativeBasis <: Exception
end

struct InvalidDifferentialForm <: Exception
    T::String
end

Base.showerror(io::IO, e::InvalidDifferentialForm) = begin
    println(io, e.T)
    # println(io, e.T)
    # println(io, " should be totally antisymmetric.")
end
Base.showerror(io::IO, e::OutOfIndices) = begin
    println("Ran out of indices in index set '"*e.I.name*"'.")
end
Base.showerror(io::IO, e::InvalidKroneckerIndices) = begin
    println(io,"Invalid index configuration:",e.inds)
end