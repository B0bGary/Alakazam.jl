
##############################################################################################################################
##############################################################################################################################
### Mathematical ops, arithmetic 
##############################################################################################################################
##############################################################################################################################


########################
#addition
@inline function _free_masks(T)
    up = UInt64(0);
    dn = UInt64(0);
    nl = UInt64(0)
    inds = indices(T)
    @inbounds for k in eachindex(inds)
        idx, pos = inds[k]
        b = _bit(idx)
        m = -UInt64(pos.position)
        dn |= b & m
        up |= b & ~m
        isLightcone(idx) || (nl |= b)
    end
    # free = appears at exactly one position; non-lightcone only
    return ((up | dn) & ~(up & dn)) & nl
end

# function _addable_nonlightcone(A, B)
#     _free_masks(A) == _free_masks(B) || return false     # cheap reject
#     return _addable_exact(A, B)                          # allocation-free confirm
# end
_addable_nonlightcone(A, B) = _addable_exact(A, B)
@inline function _is_dummy_of(idx, inds)
    hi = false;
    lo = false
    @inbounds for p in inds
        p[1].uid == idx.uid || continue
        p[2].position ? (lo = true) : (hi = true)
        hi && lo && return true
    end
    return false
end

function _addable_exact(A, B)
    ia = indices(A);
    ib = indices(B)
    nA = 0;
    nB = 0
    @inbounds for p in ia
        idx = p[1]
        isLightcone(idx) && continue
        _is_dummy_of(idx, ia) && continue
        nA += 1
        # must appear in B at the same position, and be free there
        found = false
        for q in ib
            (q[1].uid == p[1].uid && q[2] === p[2]) || continue          # was: q === p
            _is_dummy_of(idx, ib) || (found = true)
            break
        end
        found || return false
    end
    @inbounds for p in ib
        idx = p[1]
        isLightcone(idx) && continue
        _is_dummy_of(idx, ib) && continue
        nB += 1
    end
    return nA == nB
end

+(A::TensorExpression, B::SymorNum) = +(B, A)
+(A::TensorTerm, B::SymorNum) = +(B, TensorExpression(A))
+(A::SymorNum, B::TensorTerm) = +(A, TensorExpression(B))
+(A::TensorSuperType, B::SymorNum) = +(B, TensorExpression(A))
+(A::SymorNum, B::TensorSuperType) = +(A, TensorExpression(B))

function +(A::TensorSuperType, B::TensorSuperType)
    _addable_nonlightcone(A, B) || throw(FreeIndexMismatch())
    return TensorExpression(A, B)
end

+(A::TensorExpression, B::TensorSuperType) = +(B, A)
function +(A::TensorSuperType, B::TensorExpression)
    isZeroTens(B) && return TensorExpression(A)
    if !_addable_nonlightcone(A, B)
        @show A, B
        throw(FreeIndexMismatch())
    end
    return TensorExpression(A, B)
end
function +(A::SymorNum, B::TensorExpression)
    isZero(A) && return B
    inds = indices(B)
    @inbounds for p in inds
        idx = p[1]
        isLightcone(idx) && continue
        _is_dummy_of(idx, inds) && continue
        throw(FreeIndexMismatch())      # a non-lightcone free index → can't add a scalar
    end
    return TensorExpression(CoeffPair[A => TensorTerm(TensorSuperType[]); B.summands])
end
+(B::TensorSuperType, A::TensorTerm) = +(A, B)
function +(A::TensorTerm, B::TensorSuperType)
    _addable_nonlightcone(A, B) || throw(FreeIndexMismatch())
    return TensorExpression(A, B)
end



function +(A::TensorTerm, B::TensorTerm)
    _addable_nonlightcone(A, B) || throw(FreeIndexMismatch())
    return TensorExpression(A, B)
end


+(B::TensorExpression, A::TensorTerm) = +(A, B)
function +(A::TensorTerm, B::TensorExpression)
    isZeroTens(B) && return TensorExpression(A)
    if _addable_nonlightcone(A, B)
        return TensorExpression(A, B)
    else
        @show A, B
        throw(FreeIndexMismatch())
    end
end

function +(A::TensorExpression, B::TensorExpression)
    isZeroTens(A) && return B
    isZeroTens(B) && return A
    if _addable_nonlightcone(A, B)
       return TensorExpression(A, B)
    else
        @show A, B
        throw(FreeIndexMismatch())
    end
end

"""Unchecked addition, only use if certain that these things can be added"""
⨥(A::TensorExpression, B::TensorExpression) = TensorExpression(A, B)
⨥(A::TensorExpression, B::TensorTerm) = TensorExpression(A, B)
⨥(A::TensorExpression, B::TensorSuperType) = TensorExpression(A, B)
⨥(A::TensorTerm, B::TensorExpression) = TensorExpression(A, B)
⨥(A::TensorSuperType, B::TensorExpression) = TensorExpression(A, B)

⨥(A::TensorTerm, B::TensorTerm)             = TensorExpression(A, B)
⨥(A::TensorTerm, B::TensorSuperType)        = TensorExpression(A, B)
⨥(A::TensorSuperType, B::TensorTerm)        = TensorExpression(B, A)   
⨥(A::TensorSuperType, B::TensorSuperType)   = TensorExpression(A, B)

⨥(A::SymorNum, B::TensorExpression)         = isZero(A) ? B : A * consttens() ⨥ B
⨥(A::TensorExpression, B::SymorNum)         = isZero(B) ? A : B * consttens() ⨥ A
⨥(A::SymorNum, B::TensorTerm)               = A ⨥ TensorExpression(B)
⨥(A::TensorTerm, B::SymorNum)               = TensorExpression(A) ⨥ B
⨥(A::SymorNum, B::TensorSuperType)          = A ⨥ TensorExpression(B)
⨥(A::TensorSuperType, B::SymorNum)          = TensorExpression(A) ⨥ B



########################
#multiplication
##

#basic mult
*(N::SymorNum, B::TensorSuperType) = TensorExpression(N, copy(B))
*(N::SymorNum, B::TensorTerm) = TensorExpression(N => B)
*(B::TensorSuperType, N::SymorNum) = TensorExpression(N, copy(B)) #copy so one can reuse in a different expression
*(B::TensorTerm, N::SymorNum) = TensorExpression(N => B)

"""
I think one should change this to return a new TEXP
Implicitly there are grading swaps... one should ensure N is a c-num
"""
# function *(B::TensorExpression, N::SymorNum)

#     if isZeroTens(B) || isZero(N)
#         return zero()
#     end
#     c = convert(SymorNum, N)
#     out = TensorExpression(_RAW,Pair{SymorNum,TensorTerm}[c * summand[1] => summand[2] for summand in summands(B)])
#     return out
# end
function *(B::TensorExpression, N::SymorNum)
    (isZeroTens(B) || isZero(N)) && return zero()
    c = convert(SymorNum, N)
    s = summands(B)
    out = Vector{CoeffPair}(undef, length(s))
    @inbounds for i in eachindex(s)
        out[i] = (c * s[i][1]) => s[i][2]
    end
    return TensorExpression(_RAW,out)
end
*(N::SymorNum, B::TensorExpression) = *(B, convert(SymorNum, N)) #fix for grading....?

function *(T::TensorSuperType, TEx::TensorExpression)
    new_exp = zero()
    for sum in summands(TEx)
        new_exp += sum[1] * (T * sum[2])#should stop overflow?
    end
    return (new_exp)
end
function *(TEx::TensorExpression, T::TensorSuperType)
    new_exp = zero()
    for sum in summands(TEx)
        new_exp += sum[1] * (sum[2] * T)#should stop overflow?
    end
    return (new_exp)
    # end
end

function *(T::TensorTerm, TEx::TensorExpression)
    return (TensorExpression(Pair{SymorNum,TensorTerm}[sum[1] => T * sum[2] for sum in summands(TEx)]))
end
function *(TEx::TensorExpression, T::TensorTerm)
    return (TensorExpression(Pair{SymorNum,TensorTerm}[sum[1] => sum[2] * T for sum in summands(TEx)]))
end
function *(T1::TensorSuperType, T2::TensorSuperType)
    if valid_index_placement(T1, T2)
        return (TensorTerm(copy(T1), copy(T2)))
    end
end
function *(TEx::TensorTerm, T::TensorSuperType)
    if valid_index_placement(TEx, T)
        return (TensorTerm(vcat(TEx.terms, copy(T))))
    end
end
function *(T::TensorSuperType, TEx::TensorTerm)
    if valid_index_placement(T, TEx)
        return (TensorTerm(vcat(copy(T), TEx.terms)))
    end
end
function *(TEx1::TensorExpression, TEx2::TensorExpression)
    new_tex = zero()
    if TEx1 == zero() || TEx2 == zero()
        return zero()
    end
    for T in summands(TEx1)
        for T2 in summands(TEx2)
            new_tex += TensorExpression(T[1] * T2[1] => T[2] * T2[2])
        end
    end
    return (new_tex)
end

function *(T::TensorTerm, T2::TensorTerm)
    if valid_index_placement(T, T2)
        return (TensorTerm(vcat(terms(T), terms(T2))))
    end
end

########################
#subtraction (defined via mult and add)
##

-(A::TensorSuperType, B::TensorSuperType) = TensorExpression(-1 * B + A)
-(A::TensorSuperType, B::TensorExpression) = TensorExpression(A + (-1 * B))
-(A::TensorExpression, B::TensorSuperType) = TensorExpression(A + (-1 * B))
-(A::TensorExpression, B::TensorExpression) = TensorExpression(A + (-1 * B))
-(A::TensorTerm, B::TensorExpression) = 1 * A - B
-(A::TensorTerm, B::TensorTerm) = 1 * A - B
-(A::TensorExpression, B::TensorTerm) = A - (1 * B)
-(A::TensorTerm, B::TensorSuperType) = (1 * A) - (1 * B)
-(A::TensorSuperType, B::TensorTerm) = (1 * A) - (1 * B)

#negation
-(A::TensorSuperType) = -1 * A
-(A::TensorTerm) = -1 * A
-(A::TensorExpression) = -1 * A


########################
#division (by numerical factors)
##

/(B::TensorSuperType, N::SymorNum) = *(1 / N, B)
/(B::TensorExpression, N::SymorNum) = *(1 / N, B)
/(B::TensorTerm, N::SymorNum) = *(1 / N, B)

(//)(B::TensorSuperType, N::SymorNum) = *(1 // N, B)
(//)(B::TensorExpression, N::SymorNum) = *(1 // N, B)
(//)(B::TensorTerm, N::SymorNum) = *(1 // N, B)



##########################
#### powers

Base.:^(a::TensorSuperType, p::Integer) = begin
    p < 0 && error("Negative tensor powers not implemented")
    p == 0 && return 1
    p == 1 && return a
    pair = a * op_pos(a)
    npairs = p ÷ 2
    result = pair
    for _ in 2:npairs
        result = result × pair
    end
    if isodd(p)
        result = a × result
    end
    result
end