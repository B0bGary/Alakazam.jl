

struct FreeIndexMismatch <: Exception
    #println("Free indices in terms are different.")
end


struct OverusedIndices <: Exception
    #println("One or more free indices in the expression is already used as a dummy index pair.")
end


struct InvalidKroneckerIndices <: Exception
end


struct InvalidEpsilonIndices <: Exception
end

struct OperatorPrecedenceError <: Exception
end

struct OutOfIndices <: Exception
end

struct InvalidReplacementPattern <: Exception
end

struct MixedWeightOperandTerms <: Exception
end

struct SameLevelDummieError <: Exception
end