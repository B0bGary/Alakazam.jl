###########################################################################

# Native Julia replacement for the PyCall + networkx graph isomorphisms we used to use.
#replacements:
#   ISO.is_isomorphic()        ->  is_isomorphic(ISO)
#   ISO.isomorphisms_iter()    ->  isomorphisms_iter(ISO)   (yields Dict{Int,Int})
#   nx.vf2pp_isomorphism(a,b)  ->  first_isomorphism(a, b)  (Dict or nothing)
###########################################################################

import Graphs
import Graphs: add_edge!            # method for LabeledMultiGraph
const _GE = Graphs.Experimental

# ------------------------------------------------------------------
# new graph struct

"""
Undirected multigraph with arbitrary integer node ids (e.g. `-2^i * 3^j`),
string node names, weighted parallel edges, and weighted self-loops.
"""
mutable struct LabeledMultiGraph
    ids::Vector{Int}                          # original ids, insertion order
    id2v::Dict{Int,Int}                       # original id -> internal 1..n
    names::Vector{String}                     # node "name" attribute
    loops::Vector{Vector{Float64}}            # sorted self-loop weights
    adj::Vector{Vector{Int}}                  # internal adjacency (no loops)
    edges::Dict{Tuple{Int,Int},Vector{Float64}} # (min,max) -> sorted weights
    skel::Union{Nothing,Graphs.SimpleGraph{Int}} # unweighted skeleton handed to vf2, built on demand and dropped whenever the graph changes. should save some time
end

"""
    LabeledMultiGraph([nnodes, nedges])

Empty graph, optionally sized for what is about to be put in it.
"""
function LabeledMultiGraph(nnodes::Int=0, nedges::Int=0)
    ids   = Int[];                        nnodes > 0 && sizehint!(ids, nnodes)
    names = String[];                     nnodes > 0 && sizehint!(names, nnodes)
    loops = Vector{Vector{Float64}}();    nnodes > 0 && sizehint!(loops, nnodes)
    adj   = Vector{Vector{Int}}();        nnodes > 0 && sizehint!(adj, nnodes)
    id2v  = Dict{Int,Int}();              nnodes > 0 && sizehint!(id2v, nnodes)
    edges = Dict{Tuple{Int,Int},Vector{Float64}}()
    nedges > 0 && sizehint!(edges, nedges)
    return LabeledMultiGraph(ids, id2v, names, loops, adj, edges, nothing)
end

Base.copy(g::LabeledMultiGraph) = LabeledMultiGraph(
    copy(g.ids), copy(g.id2v), copy(g.names),
    [copy(l) for l in g.loops], [copy(a) for a in g.adj],
    Dict(k => copy(v) for (k, v) in g.edges), nothing)

_insert_sorted!(v::Vector{Float64}, x::Float64) =
    insert!(v, searchsortedfirst(v, x), x)

"""
Stand-in for "this node has no self-loops", shared by every such node instead of giving each its
own empty vector. Saves on allocations
"""
const _NO_LOOPS = Float64[]

"""
    add_node!(g, id; name="")

Mirrors `nx.add_node`: creates the node if new, otherwise updates its name.
"""
function add_node!(g::LabeledMultiGraph, id::Integer; name="")
    g.skel = nothing
    v = get(g.id2v, id, 0)
    if v == 0
        push!(g.ids, Int(id))
        push!(g.names, string(name))
        push!(g.loops, _NO_LOOPS)          # replaced with a private vector if a loop arrives
        push!(g.adj, Int[])
        g.id2v[Int(id)] = length(g.ids)
    else
        g.names[v] = string(name)
    end
    return g
end

function add_nodes_from!(g::LabeledMultiGraph, itr)
    for id in itr
        add_node!(g, id)
    end
    return g
end

"""
    add_edge!(g, u, v; weight)

Mirrors `nx.add_edge` on a MultiGraph: auto-creates missing nodes and
accumulates parallel edges.  Self-loops go into the node colour.
"""
function add_edge!(g::LabeledMultiGraph, u::Integer, v::Integer; weight::Real)
    g.skel = nothing
    haskey(g.id2v, u) || add_node!(g, u)
    haskey(g.id2v, v) || add_node!(g, v)
    a, b = g.id2v[u], g.id2v[v]
    if a == b
        # never mutate the shared empty
        g.loops[a] === _NO_LOOPS && (g.loops[a] = Float64[])
        _insert_sorted!(g.loops[a], Float64(weight))
    else
        key = a < b ? (a, b) : (b, a)
        ws = get!(g.edges, key) do
            # first edge between this pair -> register adjacency once
            push!(g.adj[a], b)
            push!(g.adj[b], a)
            Float64[]
        end
        _insert_sorted!(ws, Float64(weight))
    end
    return g
end

"Text dump replacement for the old matplotlib `draw_ind_graph` (debug only)."
function draw_ind_graph(g::LabeledMultiGraph)
    println("nodes:")
    for (v, id) in enumerate(g.ids)
        println("  ", id, "  name=\"", g.names[v], "\"",
            isempty(g.loops[v]) ? "" : "  loops=$(g.loops[v])")
    end
    println("edges:")
    for ((a, b), ws) in sort!(collect(g.edges))
        println("  ", g.ids[a], " -- ", g.ids[b], "  weights=", ws)
    end
end

# ------------------------------------------------------- 
# isomorphism check

function _skeleton(g::LabeledMultiGraph)
    g.skel === nothing || return g.skel
    sg = Graphs.SimpleGraph(length(g.ids))
    for (a, b) in keys(g.edges)
        Graphs.add_edge!(sg, a, b)
    end
    g.skel = sg
    return sg
end

"""
Holds the two graphs of an isomorphism problem.  Mappings run
`gsource -> gtarget`, matching the old `nx.GraphMatcher(g2, g1)` direction.
"""
struct TensorGraphIso
    gsource::LabeledMultiGraph
    gtarget::LabeledMultiGraph
    ssource::Graphs.SimpleGraph{Int}
    starget::Graphs.SimpleGraph{Int}
end

"Same call signature as before: mappings go from `g2`'s ids to `g1`'s ids."
find_tensor_graph_iso(g1::LabeledMultiGraph, g2::LabeledMultiGraph) =
    TensorGraphIso(g2, g1, _skeleton(g2), _skeleton(g1))

@inline _ekey(e) = Graphs.src(e) < Graphs.dst(e) ?
                   (Graphs.src(e), Graphs.dst(e)) : (Graphs.dst(e), Graphs.src(e))

# vertex colour = (name, self-loop weight multiset); edge colour = weight multiset
_vrel(I::TensorGraphIso) = (u, v) ->
    (I.gsource.names[u] == I.gtarget.names[v] &&
     I.gsource.loops[u] == I.gtarget.loops[v])

_erel(I::TensorGraphIso) = (e1, e2) ->
    I.gsource.edges[_ekey(e1)] == I.gtarget.edges[_ekey(e2)]

# pre-filter before VF2 algo - automatically reject things that cant be iso
function _quick_reject(gs::LabeledMultiGraph, gt::LabeledMultiGraph)
    length(gs.ids) == length(gt.ids) || return true
    length(gs.edges) == length(gt.edges) || return true
    sum(length, values(gs.edges); init=0) ==
    sum(length, values(gt.edges); init=0) || return true
    sort(gs.names) == sort(gt.names) || return true
    sort(gs.loops) == sort(gt.loops) || return true
    sort!(map(length, gs.adj)) ==
    sort!(map(length, gt.adj)) || return true
    return false
end

function is_isomorphic(I::TensorGraphIso)
    _quick_reject(I.gsource, I.gtarget) && return false
    return _GE.has_isomorph(I.ssource, I.starget;
        vertex_relation=_vrel(I), edge_relation=_erel(I))
end

"""
    isomorphisms_iter(I) -> iterator of Dict{Int,Int}

Apply iso. Each dict maps original node ids of g2 to original node
ids of g1, exactly like networkx's `isomorphisms_iter`.
"""
function isomorphisms_iter(I::TensorGraphIso)
    _quick_reject(I.gsource, I.gtarget) && return (Dict{Int,Int}() for _ in 1:0)
    ch = _GE.all_isomorph(I.ssource, I.starget;
        vertex_relation=_vrel(I), edge_relation=_erel(I))
    return (Dict{Int,Int}(I.gsource.ids[a] => I.gtarget.ids[b] for (a, b) in m)
            for m in ch)
end

"""
    first_isomorphism(gsource, gtarget) -> Dict{Int,Int} or nothing

Replacement for nx.vf2pp_isomorphism(g2, g1).  
NOTE: unlike vf2pp, this DOES match edge weights.
"""
function first_isomorphism(gsource::LabeledMultiGraph, gtarget::LabeledMultiGraph)
    I = TensorGraphIso(gsource, gtarget, _skeleton(gsource), _skeleton(gtarget))
    for m in isomorphisms_iter(I)
        return m
    end
    return nothing
end

# ------------------------------------------------------------- 
# util

"""
Extracts the powers (p2,p3) from a number 2^{p2}3^{p3}
"""
function factor_2_3(n::Int)
    n < 1 && error("n must be positive")
    n2 = 0
    n3 = 0
    while iseven(n)
        n >>>= 1
        n2 += 1
    end
    while n % 3 == 0
        n ÷= 3
        n3 += 1
    end
    n == 1 || error("has primes other than 2 and 3")
    return n2, n3
end

"""
    edges_with_weight(g, vertex_id, target_weight)

Original ids of neighbours connected to `vertex_id` by at least one edge
carrying `target_weight` (self-loops excluded).
"""
function edges_with_weight(g::LabeledMultiGraph, vertex_id::Integer, target_weight::Real)
    out = Int[]
    v = get(g.id2v, Int(vertex_id), 0)
    v == 0 && return out
    w = Float64(target_weight)
    for u in g.adj[v]
        key = u < v ? (u, v) : (v, u)
        (w in g.edges[key]) && push!(out, g.ids[u])
    end
    return out
end

@inline function _has_weight(g::LabeledMultiGraph, a::Int, b::Int, w::Float64)
    key = a < b ? (a, b) : (b, a)
    ws = get(g.edges, key, nothing)
    return ws !== nothing && (w in ws)
end

"""
    wl_colours(g) -> Vector{UInt}

Weisfeiler-Leman colour refinement. vf2 can only exclude a candidate pairing that
its vertex relation rejects, and the raw labelling is coarse. Iterates. Cells that survive refinement 
are genuinely interchangeable slots of a symmetric tensor.

Computed once per graph and reused across every comparison that graph takes part in.
"""
function wl_colours(g::LabeledMultiGraph)
    n = length(g.ids)
    col = Vector{UInt}(undef, n)
    @inbounds for v in 1:n
        col[v] = hash(g.names[v], hash(g.loops[v]))
    end
    prev = length(unique(col))
    prev == n && return col            
    buf = Vector{UInt}(undef, n)
    for _ in 1:n
        @inbounds for v in 1:n
            acc = UInt(0)
            for u in g.adj[v]
                acc += hash(col[u], hash(g.edges[v < u ? (v, u) : (u, v)]))
            end
            buf[v] = hash(col[v], acc)
        end
        copyto!(col, buf)
        c = length(unique(col))
        (c == prev || c == n) && break  # stable, or discrete
        prev = c
    end
    return col
end

"""
Invariant hash from an already-computed refinement, unequal ones prove non-isomorphism.
"""
invariant_from(col::Vector{UInt}) = hash(sort(col))

"""
    graph_invariant(g) -> UInt

A hash that isomorphic graphs always agree on, used to reject pairs before paying for the full
isomorphism search. 
"""
function graph_invariant(g::LabeledMultiGraph)
    return invariant_from(wl_colours(g))
end

"""
    find_weighted_square(g) -> Bool

True iff the graph contains a 4-cycle whose consecutive edges carry
weights `-1, 2, 1, 2` (up to rotation/reflection). This corresponds to terms that are equal to  zero(), since symmetric contracted against antisymmetric.
 Direct search this time instead of the old DFS.
"""
function find_weighted_square(g::LabeledMultiGraph)
    pat = (-1.0, 2.0, 1.0, 2.0)
    n = length(g.ids)
    for a in 1:n
        for b in g.adj[a]
            _has_weight(g, a, b, pat[1]) || continue
            for c in g.adj[b]
                c == a && continue
                _has_weight(g, b, c, pat[2]) || continue
                for d in g.adj[c]
                    (d == a || d == b) && continue
                    _has_weight(g, c, d, pat[3]) || continue
                    _has_weight(g, d, a, pat[4]) && return true
                end
            end
        end
    end
    return false
end

"""
True when a permutation of a term's factors only exchanges factors that are free to commute.
Only matrices can fail `can_swap`, so this costs nothing for ordinary tensor expressions.
"""
function perm_respects_commutation(perm_dict::Dict{Int,Int}, tt::TensorTerm)
    ts = terms(tt)
    n = length(ts)
    @inbounds for i in 1:(n-1)
        ts[i] isa MatrixSuperType || continue
        for j in (i+1):n
            ts[j] isa MatrixSuperType || continue
            # the pair is exchanged exactly when the permutation inverts their relative order
            if perm_dict[i] > perm_dict[j] && !can_swap(ts[i], ts[j])
                return false
            end
        end
    end
    return true
end

"""
Finds the overall sign change when permuting graded tensors with a particular perm.
Returns 0 if even number fermionic exchanges, 1 if odd number
"""
function perm_grading(perm_dict::Dict{Int,Int}, tt::TensorTerm)
    n = length(tt.terms)
    # Initialize inversion count
    inv_count = 0
    # Loop over each pair of indices in the original permutation
    for i in 1:(n-1)
        for j in (i+1):n
            # Check if the elements are in reverse order in the permutation (i.e., an inversion)
            # println(perm_dict)
            if perm_dict[i] > perm_dict[j]
                # Check if the grading of elements at i and j are different
                if grading(tt.terms[i]) * grading(tt.terms[j]) == 1
                    inv_count += 1
                end
            end
        end
    end
    # If the number of inversions is even, return true (even permutation), else false (odd permutation)
    return inv_count % 2 # to then do (-1)^perm_grading

end
# =====================================================================
#  Equality checkers
# =====================================================================

"""
Potentially faster do_swaps that doesn't require explicitly computing all symmetries of the tableau.
"""
function do_swaps(T1::TensorSuperType, patt::OrderedSet, movers::Dict)
    if length(patt) <= 1
        return 1 - length(patt)
    end

    orig_inds = [first(p) for p in indices(T1)]
    target_inds = copy(orig_inds)
    for (pos, idx) in movers
        try
            target_inds[pos] = idx
        catch e
            @show T1
            @show patt
            @show movers
            rethrow(e)
        end
    end
    rho = derive_rho(orig_inds, target_inds)
    rho === nothing && return 0
    for tab in tableaux(T1)
        kind = permutation_kind(tab, rho)
        kind === :row && return 1
        kind === :col && return -1
        kind === :identity && return 1
    end
    return 0
end

# the "simple algo" in question
function dumless_same_tensor_mod_syms(T1::TensorSuperType, T2::TensorSuperType)
    #new version to use the tableau swapping algo, which should be more general and work for any symmetries
    # (sign, new_inds) = find_canonical_rep(T)(T1) #just to make sure the tableaux are in the correct order, since the swapping algo relies on this. also should be a no-op
    (sign1, new_inds1) = find_canonical_rep(T1)
    (sign2, new_inds2) = find_canonical_rep(T2)
    return new_inds1 == new_inds2 ? sign1 * sign2 : 0
end

# the "simple algo" in question (again)
function to_dumless_same_tensor_mod_syms(T1::TensorSuperType, T2::TensorSuperType)
    #name,indpos
    inds1 = get_index_index(T1)
    dums1 = dummies(T1)
    dums2 = dummies(T2)
    inds2 = get_index_index(T2)
    inds1 = OrderedSet([pair for pair in (inds1) if !(pair[1] in dums1)])
    inds2 = OrderedSet([pair for pair in (inds2) if !(pair[1] in dums2)])
    needperm1 = setdiff(inds1, inds2) #indices of T1 that don't match T2
    needperm2 = setdiff(inds2, inds1) #indices of T2 that don't match T1
    inds_to_move = index_index_to_dict(needperm2)
    out = do_swaps(T1, needperm1, inds_to_move)
    return out
end

"""Special op version."""
function same_tensor_mod_syms(O1::OperatorSuperType, O2::OperatorSuperType)
    #names don't match
    if name(O1) != name(O2) || function_of(O1) != function_of(O2)
        return (false, 0)
    end
    return signed_equals(O1, O2)
end

function same_tensor_mod_syms(T1::TensorSuperType, T2::TensorSuperType)
    if !(typeof(T1) <: OperatorSuperType && typeof(T2) <: OperatorSuperType) && (typeof(T1) != typeof(T2) || name(T1) != name(T2) || function_of(T1) != function_of(T2) || tableaux(T1) != tableaux(T2) || labels(T1) != labels(T2) || _extra_identity(T1) != _extra_identity(T2))
        return (false, 0)
    end
    dums1 = dummies(T1)
    dums2 = dummies(T2)
    #if we dont need to worry about dummies, do simple algo
    if isempty(dums1) && isempty(dums2)
        # println("no dummies")
        out = dumless_same_tensor_mod_syms(T1, T2)
        if out != 0
            return (true, out)
        end
    else
        out = same_tensor_mod_dummies(T1, T2)
        if out != 0
            return (true, out)
        end
    end
    return (false, 0)
end

"""
Vertex colour for a tensor node in the isomorphism graph (encodes distinguishing properties)
"""
function _graph_vertex_name(t)
    base = string(nameof(typeof(t)), ":", name(t))
    ex = _extra_identity(t)
    ex === nothing || (base = string(base, "|", ex))
    ls = labels(t)
    isempty(ls) && return base
    return string(base, "|", join((string(l[1], l[2] === upper ? "^" : "_") for l in ls), ","))
end

"""
Explicitly construct a labeled graph from a tensor term.
"""
function construct_tensorterm_graph(T1::TensorTerm, debug=false)
    # Size the graph before filling it. One node per factor and per index slot; one edge per slot
    # (slot to its factor), one per symmetry pair a tableau contributes, and one per contracted
    # pair. Operators add a few more, so this is a lower bound and the containers may still grow
    _nt = length(terms(T1))
    _ni = 0
    _ne = 0
    for tens in terms(T1)
        r = length(indices(tens))
        _ni += r
        _ne += r
        for tab in tableaux(tens)
            d = size(tab.adj_matrix, 1)
            _ne += (d * (d - 1)) ÷ 2
        end
    end
    g1 = LabeledMultiGraph(_nt + _ni, _ne + _ni ÷ 2)
    node_counter = 0
    dums = dummies(T1)   # don't name dummies in the graph

    for t1pl in 1:length(terms(T1))
        tens1 = T1.terms[t1pl]

        if typeof(tens1) <: OperatorSuperType
            chn = get_operator_chain(tens1)
            for ch in 1:length(chn)
                opnm = name(chn[ch])
                opp = chn[ch]
                # prime-power decomposition to assign unique ids to op nodes
                add_node!(g1, -2^t1pl * 3^ch, name=opnm)
                for indexx in 1:length(op_indices(opp))
                    if !(op_indices(opp)[indexx][1] in dums)
                        add_node!(g1, 2 * (node_counter + indexx) + 1,
                            name=name(op_indices(opp)[indexx][1]))
                    else
                        add_node!(g1, 2 * (node_counter + indexx) + 1, name="")
                    end
                    add_edge!(g1, -2^t1pl * 3^ch, 2 * (node_counter + indexx) + 1, weight=0)
                end
                if ch > 1
                    for bk in 1:(ch-1)
                        add_edge!(g1, -2^t1pl * 3^bk, -2^t1pl * 3^ch, weight=3)
                    end
                end
                node_counter += length(op_indices(opp))

                oprnd = operand(chn[ch])
                singl = isTensor(oprnd)
                if !(singl === nothing)
                    opr_name = _graph_vertex_name(singl)
                    if ch == length(chn)   # last operator in chain: add operand node
                        add_node!(g1, 2 * t1pl, name=opr_name)
                        for bk in 1:ch
                            add_edge!(g1, -2^t1pl * 3^bk, 2 * t1pl, weight=3)
                        end
                        for indexx in 1:length(indices(singl))
                            if !(indices(singl)[indexx][1] in dums)
                                add_node!(g1, 2 * (node_counter + indexx) + 1,
                                    name=name(indices(singl)[indexx][1]))
                            else
                                add_node!(g1, 2 * (node_counter + indexx) + 1, name="")
                            end
                            add_edge!(g1, 2 * t1pl, 2 * (node_counter + indexx) + 1, weight=0)
                        end
                        for tab in tableaux(singl)
                            adm = tab.adj_matrix
                            dims = size(adm)[1]
                            for i in 1:dims, j in (i+1):dims
                                add_edge!(g1, 2 * (node_counter + i) + 1,
                                    2 * (node_counter + j) + 1, weight=adm[i, j])
                            end
                        end
                        node_counter += length(indices(singl))
                    end
                end
            end
        else
            add_node!(g1, 2 * t1pl, name=_graph_vertex_name(tens1))
        end

        if tens1 isa OperatorSuperType
            continue
        else
            for indexx in 1:length(indices(tens1))
                if !(indices(tens1)[indexx][1] in dums)
                    add_node!(g1, 2 * (node_counter + indexx) + 1,
                        name=indices(tens1)[indexx][1].name)
                else
                    add_node!(g1, 2 * (node_counter + indexx) + 1, name="")
                end
                add_edge!(g1, 2 * t1pl, 2 * (node_counter + indexx) + 1, weight=0)
            end
            for tab in tableaux(tens1)
                adm = tab.adj_matrix
                dims = size(adm)[1]
                for i in 1:dims, j in (i+1):dims
                    add_edge!(g1, 2 * (node_counter + i) + 1,
                        2 * (node_counter + j) + 1, weight=adm[i, j])
                end
            end
            node_counter += length(indices(tens1))
        end
    end
    # contractions between dummies
    dums11 = dummies_with_pos_index_pos_ul(T1)
    contractions1 = [(x, y) for x in dums11 for y in dums11
                                                if x[2][1] == y[2][1] && x != y]
    if debug
        println("dummies")
        println(T1)
        println(indices(T1))
        println(dummies(T1))
        println(dums11)
        println(contractions1)
    end
    for pair in contractions1
        add_edge!(g1, 2 * pair[1][1] + 1, 2 * pair[2][1] + 1, weight=2)
        if pair[1][2][2] == upper
            add_edge!(g1, 2 * pair[1][1] + 1, 2 * pair[1][1] + 1, weight=-2)
        else
            add_edge!(g1, 2 * pair[2][1] + 1, 2 * pair[2][1] + 1, weight=-2)
        end
    end

    return g1
end

function same_tensorterm_mod_syms(T1::TensorTerm, T2::TensorTerm, debug=false)
    eq, gr = same_tensorterm_mod_syms_grade(T1, T2, debug)
    return eq && gr == 1
end

"""
Grading is the sign flip to turn T2 into T1.
"""
function same_tensorterm_mod_syms_grade(T1::TensorTerm, T2::TensorTerm, debug=false;
                                        g1=nothing, g2=nothing, ws1=nothing, ws2=nothing)
    length(terms(T1)) == length(terms(T2)) || return (false, 0)
    # update_indices!(T1)
    # update_indices!(T2)
    if free_indices(T1) != free_indices(T2)
        return (false, 0)
    end
    # A caller comparing one term against many can build each graph once and pass it in.
    g1 === nothing && (g1 = construct_tensorterm_graph(T1))
    g2 === nothing && (g2 = construct_tensorterm_graph(T2))
    if debug
        draw_ind_graph(g1)
        draw_ind_graph(g2)
    end
    ISO = find_tensor_graph_iso(g1, g2)
    # Only the weighted-square shortcut needs isomorphism established up front
    if (ws2 === nothing ? find_weighted_square(g2) : ws2) &&
       (ws1 === nothing ? find_weighted_square(g1) : ws1)
        return is_isomorphic(ISO) ? (true, 1) : (false, 0)
    end
    #usually just 1, sometimes 2. biggest ive seen is 12.
    for isom_map in isomorphisms_iter(ISO)
        # decode the graph id scheme back into term/index/operator positions
        filtered_even = Dict(k ÷ 2 => v ÷ 2 for (k, v) in isom_map
                                                if k > 0 && v > 0 && iseven(k) && iseven(v))
        # filtered_odd = Dict((k - 1) ÷ 2 => (v - 1) ÷ 2 for (k, v) in isom_map
        #                                                    if k > 0 && v > 0 && isodd(k) && isodd(v))
        # filtered_odd_undiv = Dict(k => v for (k, v) in isom_map
        #                                      if k > 0 && v > 0 && isodd(k) && isodd(v))
        
        #faster deriving 
        filtered_odd_undiv = Dict(k => v for (k, v) in isom_map if k > 0 && v > 0 && isodd(k) && isodd(v))
        filtered_odd = Dict((k-1)÷2 => (v-1)÷2 for (k, v) in filtered_odd_undiv)

        filtered_nega = Dict(-k => -v for (k, v) in isom_map if k < 0 && v < 0)
        filtered_ops = Dict(factor_2_3(k) => factor_2_3(v) for (k, v) in filtered_nega)

        TT2 = copy_owned(T2) ## safe to hard-reassign on indices
        sign = 1
        flag = false
        treeswaps = Dict{Int,Int}()

        for perm_pair in filtered_ops
            if haskey(treeswaps, perm_pair[1][1])
                if treeswaps[perm_pair[1][1]] != perm_pair[2][1]
                    println("Nested operators are switching trees!")
                    flag = true
                    break
                end
            else
                treeswaps[perm_pair[1][1]] = perm_pair[2][1]
            end

            # only a chain of partial derivatives can be reordered by operator index: a trace has
            # no operator indices, and a Lie bracket does not commute with a derivative
            chain = get_operator_chain(TT2.terms[perm_pair[1][1]], PartialSuperType)
            if !(chain[1] isa PartialSuperType)
                flag = true
                break
            end
            rt = root_oprnd(chain[1])
            swp = graded_swap(chain, [perm_pair[1][2], perm_pair[2][2]])
            rb = rebuild_nested(swp[1], rt)
            sign *= swp[2]
            TT2.terms[perm_pair[1][1]] = rb
            update_indices!(TT2)
        end
        if flag
            continue
        end
        # reassign dummies so they match T1
        counter = 0
        for termid in 1:length(TT2.terms)
            ttensor = TT2.terms[termid]
            for indid in 1:length(indices(ttensor))
                wgtd = edges_with_weight(g2, 2 * (indid + counter) + 1, 2)
                if !isempty(wgtd)
                    if typeof(TT2.terms[termid]) <: OperatorSuperType
                        ind_id_counter = indid
                        nested = isTensor(TT2.terms[termid])
                        while (isTensor(nested) isa OperatorSuperType) &&
                            length(op_indices(isTensor(nested))) < ind_id_counter
                            ind_id_counter -= length(op_indices(isTensor(nested)))
                            nested = operand(isTensor(nested))
                        end
                        if (isTensor(nested) isa OperatorSuperType)
                            op_indices(isTensor(nested))[ind_id_counter] =
                                indices(T1)[Integer((filtered_odd_undiv[2*(indid+counter)+1] - 1) // 2)]
                        else
                            indices(isTensor(nested))[ind_id_counter] =
                                indices(T1)[Integer((filtered_odd_undiv[2*(indid+counter)+1] - 1) // 2)]
                        end
                    else
                        indices(TT2.terms[termid])[indid] =
                            indices(T1)[Integer((filtered_odd_undiv[2*(indid+counter)+1] - 1) // 2)]
                    end
                else
                    continue
                end
            end
            counter += length(indices(ttensor))
        end
        update_indices!(TT2)
        filtered_even = merge!(filtered_even, treeswaps)
        # the graph does not record commutation relations, reject invalid swaps
        perm_respects_commutation(filtered_even, TT2) || continue
        tens_grad = perm_grading(filtered_even, TT2)
        newts = Vector{TensorSuperType}(undef, length(TT2.terms))
        for perm_pair in filtered_even
            newts[perm_pair[2]] = TT2.terms[perm_pair[1]]
        end
        newtenter = TensorTerm(newts)
        newinds = Vector{Pair{IndexTypes,IndexPosition}}(undef, length(indices(T2)))
        for perm_pair in filtered_odd
            #removing try catch to allow for compiler optimising
            # try
                newinds[perm_pair[2]] = indices(TT2)[perm_pair[1]]
            # catch e
            #     println("Error occurred while checking isomorphisms: ", e)
            #     g1 = construct_tensorterm_graph(T1, true)
            #     draw_ind_graph(g1)
            # end
        end
        counter = 1
        totgrweight = 1
        flag = false
        for tensid in 1:length(newtenter.terms)
            ntensinds = indices(newtenter.terms[tensid])
            rearrinds = newinds[counter:(counter+length(ntensinds)-1)]
            curtens = terms(newtenter)[tensid]
            comptens = copy(curtens)
            replace_guide = Dict(r => n for (r, n) in zip(ntensinds, rearrinds))
            set_indices_deep!(comptens, replace_guide)

            (same, grweight) = same_tensor_mod_syms(comptens, curtens)
            if !same
                # this isomorphism relies on a symmetry the tensor doesn't have
                flag = true
                break
            end
            totgrweight *= grweight
            counter += length(ntensinds)
        end
        if flag
            continue   # try a different isomorphism
        end
        oversign = totgrweight * ((-1)^tens_grad) * sign
        return (true, oversign)
    end
    return (false, 0)
end


function same_tensor_mod_dummies(T1::TensorSuperType, T2::TensorSuperType)
    typeof(T1) != typeof(T2) && return 0
    name(T1) != name(T2) && return 0   # early exit
    length(indices(T1)) != length(indices(T2)) && return 0
    labels(T1) != labels(T2) && return 0
    dums11 = dummies_with_pos_index_pos_ul(T1)
    contractions1 = [(x, y) for x in dums11 for y in dums11
                                                if x[2][1] == y[2][1] && x != y]
    dums22 = dummies_with_pos_index_pos_ul(T2)
    contractions2 = [(x, y) for x in dums22 for y in dums22
                                                if x[2][1] == y[2][1] && x != y]
    # index graph: tableau adjacency weights + weight-2 contraction edges +
    # weight -2 self-loop marking the upper index of each contraction
    g1 = LabeledMultiGraph()
    add_nodes_from!(g1, 1:length(indices(T1)))   # NB: nx range(1,n) dropped node n; fixed here
    for tab in tableaux(T1)
        adm = tab.adj_matrix
        dims = size(adm)[1]
        for i in 1:dims, j in (i+1):dims
            add_edge!(g1, i, j, weight=adm[i, j])
        end
    end
    g2 = copy(g1)

    for pair in contractions1
        add_edge!(g1, pair[1][1], pair[2][1], weight=2)
        if pair[1][2][2] == upper
            add_edge!(g1, pair[1][1], pair[1][1], weight=-2)
        else
            add_edge!(g1, pair[2][1], pair[2][1], weight=-2)
        end
    end
    for pair in contractions2
        add_edge!(g2, pair[1][1], pair[2][1], weight=2)
        if pair[1][2][2] == upper
            add_edge!(g2, pair[1][1], pair[1][1], weight=-2)
        else
            add_edge!(g2, pair[2][1], pair[2][1], weight=-2)
        end
    end
    ISO = first_isomorphism(g2, g1)   # Dict g2-node => g1-node, or nothing
    ISO === nothing && return 0
    T3 = copy_owned(T2)
    for (oldpos, newpos) in pairs(ISO)
        indices(T3)[newpos] = indices(T2)[oldpos]
    end
    out = bubble_sort_index_with_grading!(copy_owned(T2), ISO)
    grad = out[1]
    grad *= to_dumless_same_tensor_mod_syms(T1, T3)
    return grad
end