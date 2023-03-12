

struct FreeIndexMismatch <: Exception
    #println("Free indices in terms are different.")
end


struct OverusedIndices <: Exception
    #println("One or more free indices in the expression is already used as a dummy index pair.")
end