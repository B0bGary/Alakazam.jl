

struct FreeIndexMismatch <: Exception
end


struct OverusedIndices <: Exception
    T1::Vector{Tensorial}
    OverusedIndices(T1::Tensorial,T2::Tensorial)=new([T1,T2])
    OverusedIndices(T1::Tensorial)=new([T1])
end

Base.showerror(io::IO, e::OverusedIndices) = begin
    println(io, "Indices are used too many times in: ", t for t in e.T1)
end


struct InvalidKroneckerIndices <: Exception
    inds::Any
end
Base.showerror(io::IO, e::InvalidKroneckerIndices) = begin
    println(io,"Invalid index configuration:",e.inds)
end

struct InvalidEpsilonIndices <: Exception
end

struct OperatorPrecedenceError <: Exception
end

struct OutOfIndices <: Exception
    I::IndexSet
end
Base.showerror(io::IO, e::OutOfIndices) = begin
    println(io, "Ran out of indices in IndexSet ",e.I.name,".")
end


struct InvalidReplacementPattern <: Exception
end

struct MixedWeightOperandTerms <: Exception
end

struct SameLevelDummieError <: Exception
    T1::Union{Tensorial,Vector{IndexPair}}
    T2::Union{Tensorial,Nothing}
    SameLevelDummieError(T1)=new(T1,nothing)
    SameLevelDummieError(T1,T2)=new(T1,T2)
end
Base.showerror(io::IO, e::SameLevelDummieError) = begin
    println(io, "Repeated indices at the same level: ",e.T1, (e.T2===nothing ? "" : e.T2))
end


struct InvalidDerivativeBasis <: Exception
end

struct InvalidDifferentialForm <: Exception
    T::String
end
Base.showerror(io::IO, e::InvalidDifferentialForm) = begin
    println(io, e.T)
end

struct InvalidCommutatorDefinition <: Exception
    # T::String
end
