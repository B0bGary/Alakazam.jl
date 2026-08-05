##############################################################################################
# Default identities, registered into a session by load_defaults!.
# @commutator/@anticommutator write into the global commutators table, which is emptied by clear_global!
# @mrule defines methods on reduce_matrix/reduce_trace, which clear_global! cannot reach.
##############################################################################################
"""
    load_gamma_defaults!()

Register the parts of the gamma algebra that hold in every dimension and for every index set.
These are just patterns, so the indexset is a holder
"""
function load_gamma_defaults!()
    st = IndexSet("_alakazam_clifford", 4, upper, Metric)
    a, b = Index("_γa", "_γb", st)
    ga = Gamma(a => upper)
    gb = Gamma(b => upper)
    @anticommutator {ga, gb} = 2 * Metric([a => upper, b => upper])
    @anticommutator {Gamma5(), ga} = zero()
    delete!(SET_REGISTRY, st.name) # the temp indexset not needed after pattern is made
    return nothing
end

"""
Load default package commutators and other properties to restore to a clean session
"""
function load_defaults!()
   load_gamma_defaults!()
    return nothing
end


##############################################################################################
# Identities
##############################################################################################


############################################
# Gamma-matrix traces.
"""
    spinor_dimension(st) -> Int

The number of components a Dirac spinor has over the spacetime `st`, namely `2^⌊d/2⌋`
"""
spinor_dimension(st::IndexSet) = 2^(st.range ÷ 2)

"""
    spinor_indexset(st) -> IndexSet

A helper indexset. The spinor space belonging to the spacetime `st`: a real index set whose `range` is the number of
components a Dirac spinor carries there. Built once per spacetime and reused, so that identities raised over the 
same spacetime compare equal while those over genuinely different spaces stay apart.
"""
function spinor_indexset(st::IndexSet)
    name = "spinor_" * st.name
    haskey(SET_REGISTRY, name) && return SET_REGISTRY[name]
    return IndexSet(name, spinor_dimension(st), upper)
end

"""The size of the space a gamma matrix acts on, read from the spacetime it is built over."""
_spinor_dimension(g::MatrixSuperType) = spinor_indexset(indices(g)[1][1].index_set).range

const MAX_GAMMA_TRACE = 8

"""
The metric pairing two gamma indices. `η` wants both at the same level, so one raised against one
lowered is a Kronecker delta instead.
"""
_gamma_pair(i1, p1, i2, p2) =
    p1 == p2 ? Metric([i1 => p1, i2 => p2]) : KroneckerDelta([i1 => p1, i2 => p2])

"""
Recursive helper for `tr(γ^{μ1}⋯γ^{μn})` in a space of spinor dimension `dim`.
"""
function _gamma_trace(gs::Vector{TensorSuperType}, dim)
    n = length(gs)
    isodd(n) && return zero()
    n == 0 && return TensorExpression(dim)
    out = zero()
    (i1, p1) = indices(gs[1])[1]
    for k in 2:n
        (ik, pk) = indices(gs[k])[1]
        rest = TensorSuperType[gs[j] for j in 2:n if j != k]
        out += (iseven(k) ? 1 : -1) * _gamma_pair(i1, p1, ik, pk) * _gamma_trace(rest, dim)
    end
    return out
end

"""Trace of a product of gammas, each carrying exactly one index."""
function _reduce_gamma_trace(gs::Vector{TensorSuperType}, ::TensorTerm)
    all(length(indices(g)) == 1 for g in gs) || return (false, zero())
    return (true, _gamma_trace(gs, _spinor_dimension(gs[1])))
end

for _n in 2:MAX_GAMMA_TRACE
    _args = [:($(Symbol(:g, i))::GammaSuperType) for i in 1:_n]
    _names = [Symbol(:g, i) for i in 1:_n]
    @eval reduce_trace($(_args...), tterm::TensorTerm) =
        _reduce_gamma_trace(TensorSuperType[$(_names...)], tterm)
end

"""
`tr(𝟙)` is the range of the space the identity was built on. An `IdMatrix` carrying no
space has no size to report, so the trace is not reduced.
"""
function reduce_trace(m::IdMatrix, ::TensorTerm)
    m.index_set === nothing && return (false, zero())
    return (true, TensorExpression(m.index_set.range))
end

# a single gamma, below the recursion's range
reduce_trace(::GammaSuperType, ::TensorTerm) = (true, zero())
# γ5 alone, or against two gammas
reduce_trace(::Gamma5SuperType, ::TensorTerm) = (true, zero())
reduce_trace(::Gamma5SuperType, ::GammaSuperType, ::GammaSuperType, ::TensorTerm) = (true, zero())

