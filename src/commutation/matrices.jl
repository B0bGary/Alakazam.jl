
"""
Basic Matrix object. Generically, is non-commuting: two matrices are left in the order written unless a relation between them has been registered. 
Declaring the commutator to vanish, with `@commutator [A, B] = zero()`, is what makes a pair freely exchangeable; declaring a non-zero one
lets them be exchanged at the cost of the remainder.
    """
mutable struct Mat <: MatrixSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Mat(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Mat(name::String, args...)
        return new(build_tensor(name, ℂ, args...)...)
    end
end

#explicit fermionic "matrix"
mutable struct Spinor <: MatrixSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Spinor(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Spinor(name::String, args...)
        return new(build_tensor(name, ℂ, args...)...)
    end
end
grading(::Spinor)=1

"""
Back-compatible alias. Renamed to mat to avoid conflict with base
"""
const Matrix = Mat

"""
The identity matrix. 
Commutes with everything, and absorbs into any adjacent matrix under `contract_adjacent()`.
"""
mutable struct IdMatrix <: MatrixSuperType
    data::TensorData
    indices::Vector{IndexPair}
    index_set::Union{IndexSet,Nothing}# The space it is the identity on
    IdMatrix(d::TensorData, i::Vector{IndexPair}) = new(d, i, nothing)
    IdMatrix(d::TensorData, i::Vector{IndexPair}, st::Union{IndexSet,Nothing}) = new(d, i, st)
    IdMatrix(st::IndexSet, args...) = new(build_tensor("𝟙", args...)..., st)
    IdMatrix(args...) = new(build_tensor("𝟙", args...)..., nothing)
end

# specialised copy with indexset
copy(m::IdMatrix) = IdMatrix(m.data, m.indices, m.index_set)

# Two identities on different spaces are different objects, even where the dims are the same
function ==(a::IdMatrix, b::IdMatrix)
    a.index_set == b.index_set || return false
    return indices(a) == indices(b) && labels(a) == labels(b)
end
Base.hash(m::IdMatrix, h::UInt) =
    hash(m.index_set === nothing ? "" : m.index_set.name, hash(IdMatrix, h))
_extra_identity(m::IdMatrix) = m.index_set === nothing ? nothing : m.index_set.name


"""Ket used in QM"""
mutable struct Ket <: KetSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Ket(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Ket(name::String, args...)
        return new(build_tensor(name, ℂ, args...)...)
    end
    Ket(d::TensorData) = new(d, IndexPair[])
end
"""Bra used in QM"""
mutable struct Bra <: BraSuperType
    data::TensorData
    indices::Vector{IndexPair}
    Bra(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Bra(name::String, args...)
        return new(build_tensor(name, ℂ, args...)...)
    end
    Bra(d::TensorData) = new(d, IndexPair[])
end

function Base.show(io::IO, T::BraSuperType)
    print(io, "⟨", tens2str(T), "|")
end
function Base.show(io::IO, T::KetSuperType)
    print(io, "|", tens2str(T), "⟩")
end

"""Custom printing for braket notation."""
function show_factor(io::IO, T::KetSuperType, pos::Int64, TT::TensorTerm)
    after_bra = pos > 1 && terms(TT)[pos - 1] isa BraSuperType
    print(io, after_bra ? "" : "|", tens2str(T), "⟩")
end

mutable struct PauliM <: PauliSuperType
    data::TensorData
    indices::Vector{IndexPair}
    PauliM(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function PauliM(name::String, args...)
        return new(build_tensor(name, ℂ, args...)...)
    end
    function PauliM(args...)
        return new(build_tensor("σ", ℂ, args...)...)
    end

end

"""
Rank of a Clifford basis element read straight off an index list: the number of spacetime vector
indices.
"""
_clifford_rank_inds(i) =
    count(p -> !_is_spinor_set(p[1]) && _is_spacetime_set(p[1].index_set), i)

mutable struct GammaM <: GammaSuperType
    data::TensorData
    indices::Vector{IndexPair}
    GammaM(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function GammaM(name::String, args...)
        return new(build_tensor(name, ℂ, args...)...)
    end
    function GammaM(args...)
        return new(build_tensor("γ", ℂ, args...)...)
    end
end

# Component forms, these commute.
mutable struct GammaT <: GammaTensorSuperType
    data::TensorData
    indices::Vector{IndexPair}
    GammaT(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function GammaT(name::String, args...)
        return new(build_tensor(name, ℂ, args...)...)
    end
    function GammaT(args...)
        return new(build_tensor("γ", ℂ, args...)...)
    end
end

mutable struct PauliT <: PauliTensorSuperType
    data::TensorData
    indices::Vector{IndexPair}
    PauliT(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function PauliT(name::String, args...)
        return new(build_tensor(name, ℂ, args...)...)
    end
    function PauliT(args...)
        return new(build_tensor("σ", ℂ, args...)...)
    end
end

"""
True when the indices span more than one index set.

This is what separates the two forms. The matrix form carries only the vector indices it is
named for, with the spinor indices implicit in the matrix product. The component form carries
spinor indices explicitly, and those necessarily belong to a different set.
"""
function _is_component_form(t)
    inds = indices(t)
    length(inds) < 2 && return false
    s = inds[1][1].index_set
    return any(p -> p[1].index_set !== s, inds)
end

"""
True when an index set carries spinors rather than spacetime vectors.
"""
function _is_spinor_set(S::IndexSet)
    S.metric_type !== nothing && S.metric_type <: AntisymMetricSuperType && return true
    return S.grading == fermionic.grading  
end
_is_spinor_set(idx) = _is_spinor_set(idx.index_set)

"""
    clifford_rank(t)

Rank of a Clifford basis element: the number of spacetime vector indices it carries, so 1 for
γ^{a}, 2 for γ^{ab} or Σ_{ab}, and d for the top form.
"""
clifford_rank(t) = _clifford_rank_inds(indices(t))

"""
A set named "spacetime" is taken at its word, that being what the usual set-up calls it.

Failing the name, spacetime sets are those carrying a symmetric metric, which is what
distinguishes them from internal ones such as an SU(2) R-symmetry set. A set with no
declared metric is taken to be spacetime, since that is the common case for a plain vector
index
"""
function _is_spacetime_set(S::IndexSet)
    S.name == "spacetime" && return true
    S.metric_type === nothing && return true
    return S.metric_type <: SymMetricSuperType
end

"""
    Gamma(args...)

Gamma matrices. Returns the matrix form `GammaM` when the indices lie in a single set, and the
component form `GammaT` when spinor indices are given explicitly.
"""
Gamma(d::TensorData, i::Vector{IndexPair}) = GammaM(d, i)
function Gamma(name::String, args...)
    t = GammaM(name, args...)
    return _is_component_form(t) ? GammaT(data(t), indices(t)) : t
end
function Gamma(args...)
    t = GammaM(args...)
    return _is_component_form(t) ? GammaT(data(t), indices(t)) : t
end

"""
    Pauli(args...)

As `Gamma`, for the Pauli matrices. The component form (σ^{μ})_{α α̇} carries one undotted and
one dotted spinor index, which live in different sets again, so the same test selects it.
"""
Pauli(d::TensorData, i::Vector{IndexPair}) = PauliM(d, i)
function Pauli(name::String, args...)
    t = PauliM(name, args...)
    return _is_component_form(t) ? PauliT(data(t), indices(t)) : t
end
function Pauli(args...)
    t = PauliM(args...)
    return _is_component_form(t) ? PauliT(data(t), indices(t)) : t
end

mutable struct Gamma5 <: Gamma5SuperType
    data::TensorData
    indices::Vector{IndexPair}
    Gamma5(d::TensorData, i::Vector{IndexPair}) = new(d, i)
    function Gamma5(name::String, args...)
        return new(build_tensor(name, ℂ, "5" => upper, args...)...)
    end
    function Gamma5(args...)
        return new(build_tensor("γ5", ℂ, "5" => upper, args...)...)
    end

end

"""
    @matrix MatrixType

Macro to define a  new MatrixType<:MatrixSuperType, with basic constructors.
For example,

    @matrix GMMatrix

    A = GMMatrix("A", [μ => lower])
"""
macro matrix(name)
    name isa Symbol || error("@matrix expects a type name symbol, got $name")
    esc(quote
        mutable struct $name <: MatrixSuperType
            data::$(TensorData)
            indices::Vector{$(IndexPair)}
            $name(args...) = new($(build_tensor)(args...)...)
            $name(d::$(TensorData), i::Vector{$(IndexPair)}) = new(d, i)
        end
    end)
end

####
#reduce special objects
#####
function contract_sigma(T::Tensorial)
    result = contract_adjacent(T, PauliSuperType)
    return eliminate_kronecker(contract_metrics(eliminate_epsilon(result)))
end

"""Eliminate matrices from an expression using matrix product rules."""
eliminate_matrices(T::Tensorial,mat_type::Type=MatrixSuperType)=contract_adjacent(T,mat_type)

"""
Generic reduce matrix function, can be overloaded to define new identities.
"""
reduce_matrix(m::MatrixSuperType, m2::MatrixSuperType, tterm::TensorTerm) = (false,  zero())

# the identity absorbs either side; the Id-Id method is needed to break the ambiguity between them
#removed because may belong to differnt spaces. fix at somepoint when relevant
reduce_matrix(::IdMatrix, m::MatrixSuperType, ::TensorTerm) = (true, m) #careful, need to decide best way to handle
reduce_matrix(m::MatrixSuperType, ::IdMatrix, ::TensorTerm) = (true, m)# careful, need to decide best way to handle
reduce_matrix(i::IdMatrix, ::IdMatrix, ::TensorTerm) = (false, zero()) #makes sure bad prods dont vanish. really, need to check that theyd id on same space

"""
Largest gamma matrix product to define a rule for. 
We can stil cancel higher products, but there just won't be explicit rules, have to recurse down
"""
const MAX_GAMMA_CONTRACTION = 8

"""
Recursive helper for contracting gamma products `γ^μ γ^{ν1}...γ^{νn} γ_μ` in `d` dimensions.

The two- and three-gamma interiors are given in closed form,

    γ^μ γ^ν γ^ρ γ_μ      =  4η^{νρ} + (d-4) γ^ν γ^ρ
    γ^μ γ^ν γ^ρ γ^σ γ_μ  = -2γ^σ γ^ρ γ^ν + (4-d) γ^ν γ^ρ γ^σ

rather than left to the recursion. The recursion only ever anticommutes the last gamma past `γ_μ` -
the three-gamma rule applied repeatedly - so on its own it lands on equivalent but unreduced forms
such as `2γ^ργ^ν + (d-2)γ^νγ^ρ`. Longer interiors still recurse, but onto these, and the four-gamma
interior then comes out as the standard `2(γ^λγ^νγ^ργ^σ + γ^σγ^ργ^νγ^λ) + (d-4)γ^νγ^ργ^σγ^λ` with
no rule of its own.
"""
function _gamma_sandwich(mid::Vector{TensorSuperType}, d)
    n = length(mid)
    n == 0 && return TensorExpression(d)
    if n == 2
        (i1, p1), (i2, p2) = indices(mid[1])[1], indices(mid[2])[1]
        out = TensorExpression(4 * _gamma_pair(i1, p1, i2, p2))
        return d == 4 ? out : out + (d - 4) * TensorTerm(mid)
    elseif n == 3
        out = TensorExpression(-2 * TensorTerm(TensorSuperType[mid[3], mid[2], mid[1]]))
        return d == 4 ? out : out + (4 - d) * TensorTerm(mid)
    end
    lastg = mid[end]
    head = TensorSuperType[mid[i] for i in 1:(length(mid)-1)]
    lead = isempty(head) ? TensorExpression(2 * lastg) : 2 * lastg * TensorTerm(head)
    return lead - _gamma_sandwich(head, d) * lastg
end

"""
Reduce a run of gamma matrices whose outermost two are contracted with one another.

The dimension is read from the contracted index at reduction time rather than baked in when a rule
is registered.
"""
function _reduce_gamma_contraction(gs::Vector{TensorSuperType}, ::TensorTerm)
    length(gs) >= 2 || return (false, zero())
    all(length(indices(g)) == 1 for g in gs) || return (false, zero())

    (i1, p1) = indices(gs[1])[1]
    (i2, p2) = indices(gs[end])[1]
    (i1 == i2 && p1 != p2) || return (false, zero())   # not a contracted pair

    mid = TensorSuperType[gs[k] for k in 2:(length(gs)-1)]
    return (true, _gamma_sandwich(mid, i1.index_set.range))
end
for _n in 2:MAX_GAMMA_CONTRACTION
    _args = [:($(Symbol(:g, i))::GammaSuperType) for i in 1:_n]
    _names = [Symbol(:g, i) for i in 1:_n]
    @eval reduce_matrix($(_args...), tterm::TensorTerm) =
        _reduce_gamma_contraction(TensorSuperType[$(_names...)], tterm)
end

# γ5 γ5 = 1, true in any dimension 
reduce_matrix(::Gamma5SuperType, ::Gamma5SuperType, ::TensorTerm) = (true, 1)

function reduce_matrix(σ1::PauliSuperType, σ2::PauliSuperType, tterm::TensorTerm)
    if length(indices(σ2)) == length(indices(σ1))#same #inds for contraction.
        if length(indices(σ1)) == 1#basic pauli
            refindex = indices(σ1)[1]
            if refindex[2] != indices(σ2)[1][2] #different pos
                return (false,  zero())
            end
            spare = get_unused_indices(tterm, refindex[1].index_set)
            newind = safe_pop!(spare,refindex[1].index_set)
            mt = generate_metric(indices(σ1)[1], indices(σ2)[1])
            if mt === nothing
                error("Reduction of index set not supported.")
            end
            metind=safe_pop!(spare,refindex[1].index_set)
            res = mt + im * EpsilonTensor(vcat(indices(σ1), indices(σ2), newind => refindex[2])) * generate_metric(newind => op_pos(refindex[2]),metind=>op_pos(refindex[2]))*Pauli(metind => (refindex[2]))
            return (true, res)
        elseif length(indices(σ1)) == 1#explicit spinor inds 

        end
    else
        return (false,  zero())
    end
end

"""
Accept the entry points `Gamma` and `Pauli` where a type is expected.
"""
contract_adjacent(To::Tensorial, ::typeof(Gamma)) = contract_adjacent(To, GammaSuperType)
contract_adjacent(To::Tensorial, ::typeof(Pauli)) = contract_adjacent(To, PauliSuperType)

"""
Moves matrices of the given type to be adjacent, then contracts them when possible
using simplification rules. When no identity applies to the full set of target
matrices, progressively smaller subsets are dragged adjacent and tried instead.
"""
function contract_adjacent(To::Tensorial, target_type::Type=MatrixSuperType)
    T = copy(1 * To)

    resultout =  zero()
    any_contraction = false
    changed = true

    while changed
        changed = false
        new_resultout =  zero()

        for k in eachindex(summands(T))
            summand = summands(T)[k]
            all_terms = summand[2].terms
            target_indices = findall(t -> t isa target_type, all_terms)

            if length(target_indices) < 2
                new_resultout += summand[1] * TensorTerm(all_terms)
                continue
            end
            summand_contracted = false
            # Largest subsets first
            for subset in subsets_decreasing(target_indices)
                anchor = subset[1]
                m = length(subset)

                dragged = drag_subset_adjacent(summand[1], all_terms, subset, target_type)

                candidate =  zero()
                subset_hit = false

                for (c, trms) in summands(dragged)
                    t = terms(trms)
                    run_ok = (anchor + m - 1 <= length(t)) &&
                             all(t[anchor + j - 1] isa target_type for j in 1:m)

                    reduced_here = false
                    if run_ok
                        try
                            (elim, result) = reduce_matrix(t[anchor:(anchor + m - 1)]..., trms)
                            if elim
                                candidate += c *
                                             TensorTerm(t[1:(anchor - 1)]) *
                                             result *
                                             TensorTerm(t[(anchor + m):end])
                                reduced_here = true
                                subset_hit = true
                            end
                        catch e
                            e isa MethodError || rethrow(e)
                        end
                    end

                    # commutator side-terms
                    if !reduced_here
                        candidate += c * trms
                    end
                end

                if subset_hit
                    # commit the dragged+reduced form only on success
                    new_resultout += candidate
                    summand_contracted = true
                    any_contraction = true
                    changed = true
                    break
                end
            end

            if !summand_contracted
                new_resultout += summand[1] * TensorTerm(all_terms)
            end
        end
        T = new_resultout
        resultout = new_resultout
    end
    any_contraction || return To
    return eliminate_kronecker(
        contract_metrics(
            eliminate_epsilon(resultout)
        )
    )
end

"""
Generalization of `drag_to_adjacent_n` to an arbitrary (sorted, ascending) subset
of target positions. Drags each subset element, in order, to sit directly after
the first subset element (the anchor), accumulating commutator terms.
"""
function drag_subset_adjacent(summand_coeff, all_terms, subset::Vector{Int}, target_type::Type)
    anchor = subset[1]
    n = length(subset)
    all_occ = findall(t -> t isa target_type, all_terms)
    ranks = [findfirst(==(i), all_occ) for i in subset]
    term_list = summand_coeff * TensorTerm(all_terms)
    n == 1 && return term_list
    for k in 2:n
        target_pos = anchor + k - 1
        new_term_list =  zero()
        for (c, tterm) in summands(term_list)
            trms = terms(tterm)
            occurrences = findall(t -> t isa target_type, trms)
            # commutator side-terms may have consumed some targets: skip those
            if length(occurrences) < ranks[k]
                new_term_list += c * tterm
                continue
            end
            current_pos = occurrences[ranks[k]]
            if current_pos == target_pos
                new_term_list += c * tterm
                continue
            end
            dragged = drag_to_adjacent(c, trms, target_pos - 1, current_pos)
            # a drag blocked by a non-commuting pair leaves its term as it stood
            new_term_list += dragged === nothing ? c * tterm : dragged
        end
        term_list = new_term_list
    end

    return term_list
end

"""
All subsets of `v` of size >= 2, largest first. Within a size class, subsets are
ordered by span (fewest interleaved non-subset elements first), so cheaper drags
that generate fewer commutator terms are attempted before long-range ones.
"""
function subsets_decreasing(v::Vector{Int})
    n = length(v)
    out = Vector{Vector{Int}}()
    for r in n:-1:2
        group = Vector{Vector{Int}}()
        _combinations!(group, v, r)
        sort!(group, by = s -> s[end] - s[1])   # prefer tight subsets
        append!(out, group)
    end
    return out
end

function _combinations!(out::Vector{Vector{Int}}, v::Vector{Int}, r::Int)
    n = length(v)
    r > n && return
    idx = collect(1:r)
    while true
        push!(out, v[idx])
        i = r
        while i > 0 && idx[i] == n - r + i
            i -= 1
        end
        i == 0 && return
        idx[i] += 1
        for j in (i + 1):r
            idx[j] = idx[j - 1] + 1
        end
    end
end
#######
"""
In a vector of tensorterms, move j to i+1, keeping track of grading and accumulated tensors from commutation relations

Returns `nothing` when the move would have to cross a pair of matrices that neither commute nor
carry a registered relation. Such a pair cannot be exchanged at any price, so the factors can never
be brought together and the caller should leave its term alone.
"""
function drag_to_adjacent(summand_coeff, all_terms, i, j)
    target = all_terms[j]
    current_terms = copy(all_terms)
    # drag target from position j to i+1 one step at a time
    pos = j
    coeff = summand_coeff
    term_list = coeff * TensorTerm(current_terms) #return tensor expression
    if j == i + 1
        return term_list
    end
    while pos > i + 1
        new_term_list =  zero()
        for (c, tterm) in summands(term_list)
            trms = terms(tterm)
            # Remainders generated on the way have lost the target to a commutator, so there is nothing left to drag in them.
            # they are carried through untouched and left for the caller's next pass. only summands still holding the target are moved.
            if pos > length(trms) || trms[pos] !== target
                new_term_list += c * tterm
                continue
            end
            left = trms[pos-1]
            right = trms[pos]  # this is target
            # check if they have non-trivial commutation
            if can_swap(left, right)
                grade = 0
                extra = nothing
            elseif !has_commutation_rel(left, right)
                return nothing        # non-commuting and unrelated: the drag cannot be made
            else
                (grade, extra) = get_commutator(left, right)  # returns (sign*swap, extra commutator terms)
            end
            # simple swap term
            new_terms = copy(trms)
            new_terms[pos-1] = right
            new_terms[pos] = left
            swaps = grading(trms[pos-1]) * grading(trms[pos]) + (grade < 0)
            new_term_list += (c * (-1)^(swaps) * TensorTerm(new_terms))
            # extra commutator terms if any
            if extra !== nothing
                extra_terms = TensorTerm(trms[1:(pos-2)]) * (extra) * TensorTerm(trms[(pos+1):end])
                new_term_list += (c * extra_terms)
            end
        end
        term_list = new_term_list
        pos -= 1
    end
    return term_list
end

function drag_to_adjacent_n(summand_coeff, all_terms, indices::Vector{Int})
    anchor = indices[1]
    n = length(indices)
    # start with a single summand
    term_list = summand_coeff * TensorTerm(all_terms)
    for k in 2:n
        target_pos = anchor + k - 1
        new_term_list =  zero()
        for (c, tterm) in summands(term_list)
            trms = terms(tterm)
            # find current position of k-th target by searching for target_type
            # it should be the k-th occurrence of target_type in trms
            occurrences = findall(t -> t isa typeof(all_terms[indices[k]]), trms)
            if length(occurrences) < k || occurrences[k] == target_pos
                new_term_list += c * tterm
                continue
            end
            current_pos = occurrences[k]
            dragged = drag_to_adjacent(c, trms, target_pos - 1, current_pos)
            new_term_list += dragged === nothing ? c * tterm : dragged
        end
        term_list = new_term_list
    end
    # @show term_list,"outer"
    return term_list
end

"""
If there is a simplification rule found, manipulate by raising and lowering indicies until the layout matches
"""
function match_and_apply_rule(pattern_terms, actual_terms, rhs_template, full_term)
    if length(free_indices_ordered(TensorTerm(pattern_terms))) != length(free_indices_ordered(TensorTerm(actual_terms)))
        return (false,  zero())
    end
    substitution = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()

    corrected_terms = TensorSuperType[]
    used_dummies = Set{IndexSuperType}()

    for (pat, act) in zip(pattern_terms, actual_terms)
        index_map = Dict{Pair{IndexSuperType,IndexPosition},Pair{IndexSuperType,IndexPosition}}()  # a_idx => dummy at p_pos, for replacing on act
        for (p_idx, a_idx) in zip(indices(pat), indices(act))
            p_name, p_pos = p_idx
            a_name, a_pos = a_idx

            if p_pos != a_pos
                dummies = get_unused_indices(full_term, a_name.index_set)
                filter!(x -> !(x in used_dummies), dummies)
                dummy = safe_pop!(dummies, a_name.index_set)
                push!(used_dummies, dummy)
                substitution[p_idx] = dummy => p_pos
                
                push!(corrected_terms, Metric([a_name => a_pos, dummy => a_pos]))
             
                index_map[a_idx] = (dummy => p_pos)  # act gets dummy at pattern's level
            else
                substitution[p_idx] = a_idx
                index_map[a_idx] = a_idx  # no change
            end
        end
        # replace indices on act with corrected versions
        push!(corrected_terms, set_indices_deep!(act, index_map))
    end
    # @show corrected_terms
    #handle dummies that clash in the pattern
    rhs_dummies = dummies(rhs_template)
    actual_free = Set(i[1] for i in union(indices(full_term), indices(TensorTerm(corrected_terms))))
    clashing = filter(d -> any(x -> x== d, actual_free), rhs_dummies)
    if !isempty(clashing)
        rhs_template = copy(rhs_template)
        for clash_idx in clashing
            fresh = get_unused_indices(full_term, clash_idx.index_set)
            filter!(x -> !(x in actual_free), fresh)
            replacement = safe_pop!(fresh,clash_idx.index_set)
            push!(used_dummies, replacement)
            substitution[clash_idx=>upper] = replacement => upper
            substitution[clash_idx=>lower] = replacement => lower
        end
    end
    rhs = set_from_pattern(rhs_template, substitution)
    target_indices = findall(x -> x isa typeof(actual_terms[1]), corrected_terms)
    dragged = drag_to_adjacent_n(1, corrected_terms, target_indices)
    # replace adjacent block with rhs, keeping surrounding metrics
    result =  zero()
    n = length(actual_terms)
    for (c, tterm) in summands(dragged)
        t = terms(tterm)
        new_anchor = findfirst(x -> x isa typeof(actual_terms[1]), t)
        result += c * TensorTerm(t[1:(new_anchor-1)]) * rhs * TensorTerm(t[(new_anchor+n):end])
    end
    return (true, result)
end

"""
The dimension a rule actually matched: the range of the nearest spacetime-like index
the matched terms carry. This is what `:dim` in an `@mrule` right-hand side resolves to.

A set named "spacetime" settles it outright, so every slot is searched for one of those
before `_is_spacetime_set` is consulted. Sweeping the whole term first is what makes the
difference: that test also passes an internal set which never declared a metric, and on a
term carrying both such a set would otherwise win merely by sitting first.
"""
function _matched_dimension(actual_terms)
    for pick in (S -> S.name == "spacetime", _is_spacetime_set)
        for t in actual_terms
            for (idx, _) in indices(t)
                pick(idx.index_set) && return idx.index_set.range
            end
        end
    end
    #Nothing spacetime-like: fall back to whatever slot there is
    for t in actual_terms
        inds = indices(t)
        isempty(inds) || return inds[1][1].index_set.range
    end
    return nothing
end

"""
Macro to generate a function reduce_matrix that reduces a product of matrices.
This also defines rules for reducing traces.

The right-hand side may use `:dim` for anything dimension-dependent. It is a placeholder,
not a variable. At each application it becomes the dimension of the terms the rule
actually matched, the range of the nearest spacetime-like index they carry.

    @mrule g1*g3 = :dim                 # d, whichever set the gammas came from
    @mrule g1*g2*g3 = (2 - :dim)*g2
"""
macro mrule(expr)
    @assert expr.head == :(=)

    lhs = expr.args[1]
    rhs = expr.args[2]
    is_trace = lhs isa Expr && lhs.head == :call && lhs.args[1] === :Trace
    is_trace && (lhs = length(lhs.args) == 2 ? lhs.args[2] :
                 error("@mrule: a Trace on the left must wrap a single product of matrices"))

    function extract_terms(ex)
        if ex isa Expr && ex.head == :call && ex.args[1] == :*
            out = Any[]
            for arg in ex.args[2:end]
                append!(out, extract_terms(arg))
            end
            return out
        else
            return [ex]
        end
    end
    lhs_terms = extract_terms(lhs)
    n = length(lhs_terms)

    #`:dim` on the right stands for a dimension not known until the rule is applied.
    dimvar = gensym(:dim)
    function sub_dim(ex)
        ex isa QuoteNode && ex.value === :dim && return dimvar
        ex isa Expr && return Expr(ex.head, map(sub_dim, ex.args)...)
        return ex
    end
    rhs = sub_dim(rhs)
    #Only pay for the lookup when a rule actually asks for it.
    dim_assign = rhs == expr.args[2] ? Expr(:block) :
                 :($(esc(dimvar)) = Alakazam._matched_dimension(actual_terms))

    arg_names = [gensym() for _ in 1:n]
    type_exprs = [:(typeof($(esc(t)))) for t in lhs_terms]

    args_with_types = [
        :($(arg_names[i])::$(type_exprs[i]))
        for i in 1:n
    ]

    fname = is_trace ? :(Alakazam.reduce_trace) : :(Alakazam.reduce_matrix)
    return quote
        function $fname($(args_with_types...), tterm::Alakazam.TensorTerm)
            pattern_terms = TensorSuperType[$(map(esc, lhs_terms)...)]
            actual_terms = TensorSuperType[$(arg_names...)]
            $dim_assign
            rhs_template = TensorExpression($(esc(rhs)))
            return match_and_apply_rule(
                pattern_terms,
                actual_terms,
                rhs_template,
                tterm
            )
        end
        # Force this specialisation into the pkgimage.
        precompile($fname,
                   ($(type_exprs...), Alakazam.TensorTerm))
    end
end


function to_latex(T::KetSuperType)
    out="\\ket{"
    out*=_to_latex(T)
    out*="}"
    return out
end

function to_latex(T::BraSuperType)
    out="\\bra{"
    out*=_to_latex(T)
    out*="}"
    return out
end

#A bra immediately followed by a ket is one inner product, and the braket package
#writes it as a single \braket{ψ|ϕ}
_ket_follows_bra(pos::Int64, TT::TensorTerm) =
    pos < length(terms(TT)) && terms(TT)[pos + 1] isa KetSuperType

_bra_precedes_ket(pos::Int64, TT::TensorTerm) =
    pos > 1 && terms(TT)[pos - 1] isa BraSuperType

function to_latex_factor(T::BraSuperType, pos::Int64, TT::TensorTerm)
    _ket_follows_bra(pos, TT) || return to_latex(T)
    return "\\braket{" * _to_latex(T) * "|"
end

function to_latex_factor(T::KetSuperType, pos::Int64, TT::TensorTerm)
    _bra_precedes_ket(pos, TT) || return to_latex(T)
    return _to_latex(T) * "}"
end
