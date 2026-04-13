mutable struct YoungTableau
    row_lengths::Vector{Int64}
    index_assignment::Vector{Int64}
    tableau::Vector{Vector{Int64}}
    adj_matrix::Matrix{Int64}
    symmetries::Vector{String}
    antisymmetries::Vector{String}
end

function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64})
    YoungTableau(row_lengths, index_assignment, maximum(index_assignment))
end


function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64}, tot_inds::Int64)

    nrows = length(row_lengths)
    ncols = maximum(row_lengths)

    # ---- build tableau ----
    tableau = build_rows(index_assignment, row_lengths) #changed to this form, idk if old form was used anywhere else
    # tableau=fill(-1, nrows, ncols)
    # cum = 1
    # for r in 1:nrows
    #     for c in 1:row_lengths[r]
    #         tableau[r, c] = index_assignment[cum]
    #         cum += 1
    #     end
    # end

    # -------------------------------------------------
    symmetries, antisymmetries = young_symmetry_sets(tableau)

    # println("Calculating symmetries and antisymmetries...")
    # println(tableau)
    # println(collect(eachrow(tableau)))
    # symmetries, antisymmetries = young_symmetry_sets(collect(eachrow(tableau)))




    # ---- adjacency matrix (default: all zeros) ----
    # adjm = zeros(Int, tot_inds, tot_inds)
    # adjm =Matrix{Int}(I, tot_inds, tot_inds)
    # adjm = zeros(Int, tot_inds, tot_inds)
    # for i in 1:tot_inds
    #     adjm[i, i] = 1
    # end
    # # # -------------------------------------------------
    # # # detect shape type
    # # # -------------------------------------------------
    # # is_single_row = all(r -> r == 0, row_lengths[2:end])  # only first row non-empty
    # # is_single_col = all(r -> r <= 1 for r in row_lengths)

    # # # NOTE: stronger and safer detection:
    # # is_single_row = nrows == 1 || all(x -> x == row_lengths[1], row_lengths)
    # # is_single_col = all(x -> x <= 1 for x in row_lengths)

    # # BUT better: check structure directly
    # if nrows == 1
    #     is_single_row = true
    #     is_single_col = false

    # elseif all(r -> r == 1, row_lengths)
    #     is_single_row = false
    #     is_single_col = true

    # else
    #     is_single_row = false
    #     is_single_col = false
    # end

    # # -------------------------------------------------
    # # fill adjacency ONLY if pure symmetry type
    # # -------------------------------------------------
    # if is_single_row
    #     # totally symmetric
    #     for r in 1:nrows
    #         for c1 in 1:row_lengths[r]
    #             i = tableau[r, c1]
    #             for c2 in c1+1:row_lengths[r]
    #                 j = tableau[r, c2]
    #                 adjm[i, j] = 1
    #                 adjm[j, i] = 1
    #             end
    #         end
    #     end

    # elseif is_single_col
    #     # totally antisymmetric
    #     for c in 1:ncols
    #         col_entries = Int[]

    #         for r in 1:nrows
    #             if c <= row_lengths[r]
    #                 push!(col_entries, tableau[r, c])
    #             end
    #         end

    #         for a in 1:length(col_entries)
    #             i = col_entries[a]
    #             for b in a+1:length(col_entries)
    #                 j = col_entries[b]
    #                 adjm[i, j] = -1
    #                 adjm[j, i] = -1
    #             end
    #         end
    #     end
    # end
    adjm = zeros(Int, tot_inds, tot_inds)
    for i in 1:tot_inds
        adjm[i, i] = 1
    end

    # -------------------------------------------------
    # detect shape type
    # -------------------------------------------------
    if nrows == 1
        is_single_row = true
        is_single_col = false

    elseif all(r -> r == 1, row_lengths)
        is_single_row = false
        is_single_col = true

    else
        is_single_row = false
        is_single_col = false
    end

    # -------------------------------------------------
    # fill adjacency ONLY if pure symmetry type
    # -------------------------------------------------
    if is_single_row
        # totally symmetric (row-wise)
        for r in 1:length(tableau)
            row = tableau[r]
            for c1 in 1:length(row)
                i = row[c1]
                for c2 in c1+1:length(row)
                    j = row[c2]
                    adjm[i, j] = 1
                    adjm[j, i] = 1
                end
            end
        end

    elseif is_single_col
        # totally antisymmetric (column-wise)
        ncols = maximum(length.(tableau))

        for c in 1:ncols
            col_entries = Int[]

            for r in 1:length(tableau)
                if c <= length(tableau[r])
                    push!(col_entries, tableau[r][c])
                end
            end

            for a in 1:length(col_entries)
                i = col_entries[a]
                for b in a+1:length(col_entries)
                    j = col_entries[b]
                    adjm[i, j] = -1
                    adjm[j, i] = -1
                end
            end
        end
    end
    return YoungTableau(row_lengths, index_assignment, tableau, adjm, symmetries, antisymmetries)
end



function build_rows(index_assignment::Vector, row_lengths::Vector)
    rows = Vector{Vector{eltype(index_assignment)}}(undef, length(row_lengths))

    idx = 1
    for r in 1:length(row_lengths)
        len = row_lengths[r]
        rows[r] = index_assignment[idx:idx+len-1]
        idx += len
    end

    return rows
end

function column_permutation(rows, σ)
    perm = Dict{Int,Int}()
    for r in 1:length(rows)
        for c in 1:length(rows[r])
            perm[rows[r][c]] = rows[r][σ[c]]
        end
    end
    perm
end

function row_swap(rows, i, j)
    perm = Dict{Int,Int}()
    maxlen = maximum(length.(rows))
    for c in 1:maxlen
        if c <= length(rows[i]) && c <= length(rows[j])
            perm[rows[i][c]] = rows[j][c]
            perm[rows[j][c]] = rows[i][c]
        end
    end
    perm
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

    for c in 1:minlen-1
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
        for i in 1:length(col)-1
            for j in i+1:length(col)
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
    n = sum(length.(rows))

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

function young_symmetry_sets(rows)
    n = sum(length.(rows))
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




































# mutable struct YoungTableau
#     row_lengths::Vector{Integer}
#     index_assignment::Vector{Integer}
#     tableau::Matrix{Integer}
#     adj_matrix::Matrix{Integer}
# end
#     function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64})
#         YoungTableau(row_lengths,index_assignment,maximum(index_assignment))
#     end


#     function YoungTableau(row_lengths::Vector{Int64}, index_assignment::Vector{Int64},tot_inds::Int64)
#         #create tabelau in matri form
#         matrix =[-1 for _ in 1:size(row_lengths)[1], _ in 1:row_lengths[1]]
#         cum=1
#         for row in 1:size(row_lengths)[1]
#             for i in 1:row_lengths[row]
#                 matrix[row,i]=index_assignment[cum]
#                 cum+=1

#             end
#         end
#         # println(matrix)
#         # println("here")
#         #add symmetric perms
#         adjm =[0 for _ in 1:tot_inds, _ in 1:tot_inds]
#         init=-1
#         for row in 1:size(matrix)[1]
#             init*=-1
#             cur=matrix[row,1]


#             if cur!=-1
#             for col in 1:size(matrix)[2]
#                 nex=matrix[row,col]
#                 if nex==-1
#                     break
#                 else
#                     last_col=col
#                 end



#                 while last_col>0
#                     cr=matrix[row,last_col]
#                     for cc in 1:last_col
#                         nex=matrix[row,cc]
#                         adjm[cr,nex]=1
#                         adjm[nex,cr]=1
#                     end
#                     last_col-=1
#                 end
#             end
#             end


#         for col in 1:size(matrix)[2]
#             entry=matrix[row,col]
#             if entry==-1
#                 break
#             else
#                 for row2 in row+1:size(matrix)[1]
#                     for col2 in 1:size(matrix)[2] 
#                         nex=matrix[row2,col2]
#                         if nex==-1
#                             break
#                         else
#                             adjm[entry,nex]=-1
#                             adjm[nex,entry]=-1
#                         end
#                     end
#                 end
#             end
#         end
#     end
#         # println(matrix)
#         # println(adjm)
#         YoungTableau(row_lengths, index_assignment,matrix,adjm)
#     end