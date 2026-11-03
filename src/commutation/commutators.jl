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
"""
Key under which a commutation relation is stored: name, concrete type, and index signature.
"""
const CommutatorKey = Tuple{String,Symbol,String}

"""
Sorted multiset of the KINDS of index an object carries: `v` spacetime, `s` spinor, `i` internal.
Infers the type from metrics and grading etc. Hard because we can't tie to the like indexset names,
as users can change this to random things
"""
function _ckey_sig(inds)
    tags = Vector{String}(undef, length(inds))
    @inbounds for k in eachindex(inds)
        idx = inds[k][1]
        S = idx.index_set
        kind = _is_spinor_set(S) ? "s" : (_is_spacetime_set(S) ? "v" : "i")
        tags[k] = kind * _accent_of(idx.name)
    end
    sort!(tags)
    return join(tags, ",")
end

"""
The combining mark decorating an index name, or `""` when it carries none.
"""
function _accent_of(nm::AbstractString)
    for (_, mark) in INDEX_SUFFIXES
        occursin(mark, nm) && return mark
    end
    return ""
end

_ckey(A)::CommutatorKey = (name(A), nameof(typeof(A)), _ckey_sig(indices(A)))

"""
Operators key on their OWN indices, not on `indices`.
"""
_ckey(A::OperatorSuperType)::CommutatorKey =
    (name(A), nameof(typeof(A)), _ckey_sig(op_indices(A)))

const commutators = Dict{Tuple{CommutatorKey,CommutatorKey},Commutator}()


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
    key=minmax(_ckey(A),_ckey(B))
    if haskey(commutators,key)
        delete!(commutators,key)
    end
end

function remove_commutator!(L::LieBracketSuperType)
    A=first(L)
    B=second(L)
    key=minmax(_ckey(isTensor(A)),_ckey(isTensor(B)))
    if haskey(commutators,key)
        delete!(commutators,key)
    end
end

# first(c::Commutator)=c.A
# second(c::Commutator)=c.B

"""
Non-commuting tensor must be defined as matrixes for detection purposes.
"""
function add_commutator!(grading::Integer, A::MatrixSuperType, B::MatrixSuperType, C::TensorExpression)
    if free_indices(A * B) != free_indices(C) && !isZeroTens(C)
        throw(InvalidCommutatorDefinition(A, B, C))
    end
    commutators[minmax(_ckey(A), _ckey(B))] = Commutator(grading, A, B, C)
end

add_commutator!(grading::Integer, A::MatrixSuperType,B::MatrixSuperType,C::Int64 )=add_commutator!(grading, A,B,TensorExpression(C) )
add_commutator!(k::Tuple{CommutatorKey,CommutatorKey},va::Commutator)=commutators[k]=va
function add_commutator!(va::Commutator)
    k=minmax(_ckey(va.A),_ckey(va.B))
    commutators[k]=va
end


"""Gets the result of (anti)commuting A with B"""
function get_commutator(A::MatrixSuperType,B::MatrixSuperType)::Tuple{Integer,TensorExpression}
    key=minmax(_ckey(A),_ckey(B))
    if haskey(commutators,key)
        comrel=commutators[key]
        dictout=Dict{Pair{IndexTypes,IndexPosition}, Pair{IndexTypes,IndexPosition}}()
        sign=1
        if _ckey(comrel.A)==_ckey(A)
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
    stringout*=_display_name(name(com.A))
    stringout*=","
    stringout*=_display_name(name(com.B))
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



