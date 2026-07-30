"""
Young Tableau diagram that stores information about symmetries and irreducible representations.
"""
mutable struct YoungTableau
    row_lengths::Vector{Int64}
    index_assignment::Vector{Int64}
    tableau::Vector{Vector{Int64}}
    adj_matrix::Base.Matrix{Int64}
    # symmetries::Vector{String}
    # antisymmetries::Vector{String}
end

function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64})
    YoungTableau(row_lengths, index_assignment, maximum(index_assignment))
end


function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64}, tot_inds::Int64)
    #build tableau
    tableau = build_rows(index_assignment, row_lengths) #changed to this form, idk if old form was used anywhere else
    #new adjm construction without using the full symmetry sets
    adjm2 = zeros(Int, tot_inds, tot_inds)
    for i in 1:tot_inds
        adjm2[i, i] = 1
    end
    for ((i, j), s) in tableau_single_swaps(tableau)
        adjm2[i, j] = s
        adjm2[j, i] = s
    end
    return YoungTableau(row_lengths, index_assignment, tableau, adjm2)#, symmetries, antisymmetries
end


function ==(a::YoungTableau,b::YoungTableau)
    return a.row_lengths==b.row_lengths && a.index_assignment == b.index_assignment
end
function Base.hash(a::YoungTableau, h::UInt)
    h = hash(a.row_lengths, h)
    h = hash(a.index_assignment, h)
    return h
end




function build_rows(index_assignment::Vector, row_lengths::Vector)
    rows = Vector{Vector{eltype(index_assignment)}}(undef, length(row_lengths))

    idx = 1
    for r in 1:length(row_lengths)
        len = row_lengths[r]
        rows[r] = index_assignment[idx:(idx+len-1)]
        idx += len
    end
    return rows
end

function tableau_single_swaps(tableau)

    n = maximum(vcat(tableau...))

    gens = Tuple{Vector{Int},Int}[]

    for g in row_generators(tableau)
        push!(gens, (dict_to_array(g, n), 1))
    end

    for g in column_generators(tableau)
        push!(gens, (dict_to_array(g, n), -1))
    end


    # permutation generators without signs
    pgens = [g for (g, s) in gens]

    swaps = Dict{Tuple{Int,Int},Int}()

    # find generator transpositions
    for (g, sgn) in gens

        moved = findall(i -> g[i] != i, 1:n)

        # only true single swaps
        if length(moved)==2
            a, b = moved

            if g[a]==b && g[b]==a

                # conjugate through the group
                orbit = permutation_orbit_pairs(pgens, a, b)

                for (i, j) in orbit
                    swaps[(min(i, j), max(i, j))] = sgn
                end
            end
        end
    end

    swaps
end

function row_generators(rows)
    gens = []

    # build columns explicitly
    cols = Dict{Int,Vector{Int}}()

    for r in 1:length(rows)
        for c in 1:length(rows[r])
            push!(get!(cols, c, Int[]), rows[r][c])
        end
    end

    # NO row symmetry exists unless two entries share a row AND column position
    # in Young diagrams this only happens in rectangular row segments

    minlen = minimum(length.(rows))

    for c in 1:(minlen-1)
        perm = Dict{Int,Int}()

        for r in 1:length(rows)
            perm[rows[r][c]] = rows[r][c+1]
            perm[rows[r][c+1]] = rows[r][c]
        end

        push!(gens, perm)
    end

    gens
end
function column_generators(rows)
    cols = Dict{Int,Vector{Int}}()

    for r in 1:length(rows)
        for c in 1:length(rows[r])
            cols[c] = get(cols, c, Int[])
            push!(cols[c], rows[r][c])
        end
    end

    gens = []

    for (_, col) in cols
        for i in 1:(length(col)-1)
            for j in (i+1):length(col)
                push!(gens, Dict(col[i] => col[j], col[j] => col[i]))
            end
        end
    end

    gens
end

function dict_to_array(p, n)
    a = collect(1:n)
    for (k, v) in p
        a[k] = v
    end
    a
end

function compose(p, q)
    [p[q[i]] for i in 1:length(p)]
end

function generate_group_with_paths(rows)
    n = maximum(vcat(rows...))

    sym_gens = [dict_to_array(g, n) for g in row_generators(rows)]
    antisym_gens = [dict_to_array(g, n) for g in column_generators(rows)]

    gens = vcat(sym_gens, antisym_gens)
    gen_types = vcat(fill(0, length(sym_gens)), fill(1, length(antisym_gens)))

    id = collect(1:n)

    group = Dict{Vector{Int},Set{Int}}()
    group[id] = Set([0])

    changed = true
    while changed
        changed = false
        for (g, signs) in collect(group)
            for (h, t) in zip(gens, gen_types)
                new = compose(g, h)
                new_signs = Set{Int}()
                for s in signs
                    push!(new_signs, (s + t) % 2)
                end

                if !haskey(group, new)
                    group[new] = new_signs
                    changed = true
                else
                    old = group[new]
                    unioned = union(old, new_signs)
                    if unioned != old
                        group[new] = unioned
                        changed = true
                    end
                end
            end
        end
    end

    group
end

function dict_to_cycles(p::Dict{Int,Int})
    visited = Set{Int}()
    cycles = []
    for start in sort(collect(keys(p)))
        if start ∉ visited
            cycle = []
            x = start
            while x ∉ visited
                push!(cycle, x)
                push!(visited, x)
                x = p[x]
            end
            if length(cycle) > 1
                push!(cycles, cycle)
            end
        end
    end
    cycles
end

function cycles_to_string(cycles)
    isempty(cycles) && return "()"
    join(["(" * join(c, " ") * ")" for c in cycles], "")
end

"""
Currently unused, from old version which stored this data but was too slow for large tableaux. 
Generates all permutations permitted by a tableau. Exponentially grows.
"""
function young_symmetry_sets(rows)
    n = maximum(vcat(rows...))
    group = generate_group_with_paths(rows)

    sym = String[]
    antisym = String[]

    for (g, signs) in group
        p = Dict(i => g[i] for i in 1:n)
        str = cycles_to_string(dict_to_cycles(p))

        if 0 in signs
            push!(sym, str)
        end
        if 1 in signs
            push!(antisym, str)
        end
    end

    sort!(unique(sym))
    sort!(unique(antisym))

    sym, antisym
end


"""Is one tableau embedded in the other?"""
function is_subsumed(t1::YoungTableau, t2::YoungTableau)
    is_row(t) = length(t.row_lengths) == 1
    is_col(t) = all(r -> r == 1, t.row_lengths)

    if is_row(t1) && is_row(t2)
        return t1.row_lengths[1] <= t2.row_lengths[1] &&
               issubset(t1.index_assignment, t2.index_assignment)

    elseif is_col(t1) && is_col(t2)
        n1 = length(t1.row_lengths)
        n2 = length(t2.row_lengths)
        return n1 <= n2 &&
               issubset(t1.index_assignment, t2.index_assignment)
    end

    return false
end


#custom push to stop redundant YT
function Base.push!(tabs::Set{YoungTableau}, t::YoungTableau)
    tabs.dict[t] = nothing
    reduced = reduce_tableaux(tabs)
    reduced === tabs && return tabs
    empty!(tabs.dict)
    for r in reduced
        tabs.dict[r] = nothing
    end
    return tabs
end

function add_tableau!(T::TensorDataSuperType, tab::YoungTableau)
    tabs = tableaux(T)
    if tabs === notableaux()
        tabs = Set{YoungTableau}()
        T.tableaux = tabs
        tabs.dict[tab] = nothing
        return tabs
    end
    push!(tabs, tab)
    return T
end

"""Merge tableaux, if they are compatible"""
function reduce_tableaux(tableaux::Set{YoungTableau})
    n = length(tableaux)
    n <= 1 && return tableaux
    

    ts = collect(keys(tableaux.dict))

    keep = trues(n)
    @inbounds for i in 1:n
        keep[i] || continue
        for j in 1:n
            (i == j || !keep[j]) && continue
            if is_subsumed(ts[i], ts[j])
                keep[i] = false
                break
            end
        end
    end

    result = Set{YoungTableau}()
    sizehint!(result.dict, count(keep))
    @inbounds for i in 1:n
        keep[i] && (result.dict[ts[i]] = nothing)
    end
    return result
end

"""
Shift the labels of a tableau by an offset.
"""
function offset_tableau(tab::YoungTableau, offset::Int)

    new_assignment = [
        x + offset for x in tab.index_assignment
    ]

    return YoungTableau(
        copy(tab.row_lengths),
        new_assignment,
        maximum(new_assignment)
    )
end


"""
Helper function to create a vector of young tableaux, by appending single
tableau and relabeling. `base` shifts all assignments by a constant (the
number of leading extra indices on ψ that carry no symmetry).
"""
function append_offset!(tabs, tab; base::Int=0)
    tab === nothing && return length(tabs)
    off = base + (isempty(tabs) ? 0 : sum(length(t.index_assignment) for t in tabs))
    push!(tabs, offset_tableau(tab, off))
end





function row_col_lookup(tab::YoungTableau)
    lookup = Dict{Int,Tuple{Int,Int}}()
    for (r, row) in enumerate(tab.tableau), (c, pos) in enumerate(row)
        lookup[pos] = (r, c)   # tab.tableau entries are the actual position numbers
    end
    lookup
end

# rho[p] = original position whose value ends up at position p
function derive_rho(orig_inds, target_inds)
    val_to_pos = Dict(orig_inds[k] => k for k in eachindex(orig_inds))
    n = length(orig_inds)
    rho = Vector{Int}(undef, n)
    for p in 1:n
        v = target_inds[p]
        haskey(val_to_pos, v) || return nothing
        rho[p] = val_to_pos[v]
    end
    rho
end

function permutation_kind(tab::YoungTableau, rho::Vector{Int})
    lookup = row_col_lookup(tab)
    is_row, is_col, touched = true, true, false
    for p in tab.index_assignment
        q = rho[p]
        q == p && continue
        touched = true
        (!haskey(lookup, q)) && return :outside
        rp, cp = lookup[p]
        rq, cq = lookup[q]
        is_row &= (rp == rq)
        is_col &= (cp == cq)
    end
    !touched && return :identity
    is_row && return :row
    is_col && return :col
    return :neither
end


function inverse(p)
    ip = similar(p)
    for i in eachindex(p)
        ;
        ip[p[i]] = i;
    end
    ip
end

# --- orbit + transversal (Schreier tree) ---
function orbit_transversal(base_pt::Int, gens::Vector{Vector{Int}}, N::Int)
    id = collect(1:N)
    transversal = Dict{Int,Vector{Int}}(base_pt => id)
    orbit = [base_pt];
    queue = [base_pt]
    while !isempty(queue)
        pt = popfirst!(queue)
        u = transversal[pt]
        for g in gens
            img = g[pt]
            if !haskey(transversal, img)
                transversal[img] = compose(g, u)
                push!(orbit, img);
                push!(queue, img)
            end
        end
    end
    orbit, transversal
end

# --- Schreier-Sims: build a base + strong generating set ---
struct BSGSNode
    base_pt::Int
    orbit::Vector{Int}
    transversal::Dict{Int,Vector{Int}}
end

function build_bsgs(gens::Vector{Vector{Int}}, N::Int)
    chain = BSGSNode[]
    current_gens = gens
    remaining = collect(1:N)
    id = collect(1:N)

    while !isempty(current_gens) && any(g -> g != id, current_gens)
        base_pt = findfirst(pt -> any(g -> g[pt] != pt, current_gens), remaining)
        base_pt === nothing && break
        bp = remaining[base_pt]

        orbit, transversal = orbit_transversal(bp, current_gens, N)
        push!(chain, BSGSNode(bp, orbit, transversal))

        # Schreier generators of the stabilizer of bp
        seen = Set{Vector{Int}}()
        new_gens = Vector{Int}[]
        for y in orbit
            u_y = transversal[y]
            for g in current_gens
                yg = g[y]
                u_yg = transversal[yg]
                sg = compose(inverse(u_yg), compose(g, u_y))   # fixes bp
                if sg != id && !(sg in seen)
                    push!(seen, sg);
                    push!(new_gens, sg)
                end
            end
        end
        current_gens = new_gens
        remaining = filter(!=(bp), remaining)
    end
    chain
end

function permutation_orbit_pairs(gens, a, b)

    seen = Set{Tuple{Int,Int}}()
    queue = [(a, b)]

    while !isempty(queue)

        x, y = popfirst!(queue)

        key = (min(x, y), max(x, y))
        key in seen && continue
        push!(seen, key)

        for g in gens
            push!(queue, (g[x], g[y]))
        end
    end

    seen
end

function symcopy(s::Set{YoungTableau})
    out = Set{YoungTableau}()
    for t in s
        out.dict[t] = nothing        # direct insert — skip push!/reduce_tableaux
    end
    return out
end
