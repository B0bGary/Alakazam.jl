###########################################
############ commutator behaviour



"""
Might make a version that automatically detects bosonic vs fermionic
"""
mutable struct Commutator
    grading::Int #±1 for commutator or anticommutator
    A::TensorSuperType 
    B::TensorSuperType #these probs should be matrices.. problem is i want to have operators later
    C::TensorExpression
end

#keys should be a symmetric pair eg (σ,σ), (M,P) for gens etc (they're sorted into lexo-order)
#stores sign for com/acom, indexed pair and their com rel. These should be matrix type
const commutators = Dict{Tuple{String,String},Commutator}()


# ==(C1::Commutator,C2::Commutator)= f1==f2 for fields in C1,C2
==(a::Commutator, b::Commutator) = begin
    a.grading == b.grading &&
    a.A == b.A &&
    a.B == b.B &&
    a.C == b.C
end

function clear_commutators!()
    empty!(commutators)
end

function remove_commutator!(A::MatrixSuperType,B::MatrixSuperType)
    key=minmax(name(A),name(B))
    if haskey(commutators,key)
        delete!(commutators,key)
    end
end

function remove_commutator!(L::LieBracketSuperType)
    A=first(L)
    B=second(L)
    key=minmax(name(isTensor(A)),name(isTensor(B)))
    if haskey(commutators,key)
        delete!(commutators,key)
    end
end

# first(c::Commutator)=c.A
# second(c::Commutator)=c.B

"""
Non-commuting tensor must be defined as matrixes for detection purposes.
"""
function add_commutator!(grading::Integer, A::MatrixSuperType,B::MatrixSuperType,C::TensorExpression )
    temp_term=A*B
    if free_indices(temp_term)!=free_indices(C)
        if !isZeroTens(C)
            throw(InvalidCommutatorDefinition)#they gotta have same index pattern, unless commutes
        end
    end
    name_key=minmax(name(A),name(B))
    #add the commutator
    commutators[name_key]=Commutator(grading,A,B,C)
end
add_commutator!(grading::Integer, A::MatrixSuperType,B::MatrixSuperType,C::Int64 )=add_commutator!(grading, A,B,TensorExpression(C) )
add_commutator!(k::Tuple{String,String},va::Commutator)=commutators[k]=va
function add_commutator!(va::Commutator)
    k=minmax(name(va.A),name(va.B))
    commutators[k]=va
end


"""Gets the result of (anti)commuting A with B"""
function get_commutator(A::MatrixSuperType,B::MatrixSuperType)::Tuple{Integer,TensorExpression}
    key=minmax(name(A),name(B))
    if haskey(commutators,key)
        comrel=commutators[key]
        dictout=Dict{Pair{IndexTypes,IndexPosition}, Pair{IndexTypes,IndexPosition}}()
        sign=1
        #check which order the rel corresponds to
        if name(comrel.A)==name(A)
            merge!(dictout,get_pattern(comrel.A,A))
            merge!(dictout,get_pattern(comrel.B,B))
        else
            sign*=(-comrel.grading)#fermions: {B,A}={A,B}, bosons [A,B]=-[B,A]
            merge!(dictout,get_pattern(comrel.A,B))
            merge!(dictout,get_pattern(comrel.B,A))
        end
        
        out=sign*set_from_pattern(comrel.C,dictout)
        return (comrel.grading,out)
        
    else #no commutation rels defined, trivial swap
        return (1, zero())
    end
end

"""
Based on the actual tensor, figure out the index structure result from commutator
"""
function get_pattern(Asym::MatrixSuperType,Aact::MatrixSuperType)
    dictout=Dict{Pair{IndexSuperType,IndexPosition}, Pair{IndexSuperType,IndexPosition}}()
    for (i,j) in zip(indices(Asym),indices(Aact))
        dictout[i]=j
    end
    return dictout
end
function set_from_pattern(CSym::TensorExpression,indreps::Dict{Pair{M,IndexPosition}, Pair{M,IndexPosition}}) where {M<:IndexSuperType}
    out=copy(CSym)
    set_indices_deep!(out,indreps)
    return out
end

function Base.show(io::IO,com::Commutator)
    stringout=""
    if com.grading==1
        stringout*="["
    else com.grading==-1
        stringout*="{"
    end
    stringout*=name(com.A)
    stringout*=","
    stringout*=name(com.B)
    if com.grading==1
        stringout*="]"
    else com.grading==-1
        stringout*="}"
    end
    stringout*="="
    stringout*=string(com.C)
    print(io, stringout)
end


"""
Adds a commutation relation [A,B]=C.
Also accepts {A,B}=C for convenience.
"""
macro commutator(expr)
    @assert expr.head == :(=)
    lhs = expr.args[1]
    rhs = expr.args[2]
    @assert (lhs.head == :vect || lhs.head == :braces) && length(lhs.args) == 2
    A = lhs.args[1]
    B = lhs.args[2]
    grading = lhs.head == :vect ? 1 : -1
    return esc(quote
        add_commutator!($grading, $A, $B, $rhs)
    end)
end


"""
Adds an anticommutation relation {A,B}=C.
"""
macro anticommutator(expr)
    @assert expr.head == :(=)
    lhs = expr.args[1]
    rhs = expr.args[2]
    @assert lhs.head == :braces && length(lhs.args) == 2
    A = lhs.args[1]
    B = lhs.args[2]
    return esc(quote
        add_commutator!(-1, $A, $B, $rhs)
    end)
end



