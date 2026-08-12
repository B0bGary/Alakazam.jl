############################################################
# Author: Jesse Woods
# Please cite ... if you use this script in your work.


############################################################
# This script splits a product of n Grassmann numbers θ^i_α (spinor index α, flavour
# index i) into independent group-invariant structures, essentially so that a superfield
# expansion can be expressed in terms of coefficient fields with a maximally reduced number
# of indices.
#############################################################
# Usage:
#
# One entry point is via the function 'independent_structures(n::Int, d_S::Int, d_F::Int, S_type::Symbol, F_type::Symbol)'
# Here:
#       'n' is the monomial order (length of the product of Grassmanns),
#       'd_S' the dimension of the spin group
#       'd_F' the dimension of the flavour group
#       'S_type' a symbol :none, :GL, :U, :SL, :SU, :O, :SO, :Sp determining the spin group
#       'F_type' a symbol :none, :GL, :U, :SL, :SU, :O, :SO, :Sp determining the flavour group
# 
# This function returns a vector of custom `GrassmannStructure` objects containing all relevant information, 
# which can be passed to another function if needed. These objects can also be passed to `explain_structure` to
# obtain an explicit printout.
#
# Alternatively, one could use 'render_compact(n::Int, d_S::Int, d_F::Int, S_type::Symbol, F_type::Symbol)'
# to directly print to the console the decomposition in a somewhat readable form.
#
###############################################################################
#
# The procedure follows as in the paper.
#
# Starting from the space of antisymmetric grassmanns Λ^n(V⊗W), we split into Young-diagram 
# blocks using the Cauchy identity. For :Sp/:O/:SO each block is further split into that group's
# actual irreducible pieces. Those branchings are computed by comparing  characters in exact 
# rational arithmetic and are re-checked on every call (the result must reproduce the character exactly, 
# be whole numbers, and have dimensions that add up); on any failure the code stops with an error.
# Explicit example witness contractions (ε, ω, η) are attached in the printout wherever a simple one exists;
# if none can be found, it is reported as "composite".
#
# Worth knowing:
# * Odd SO(N): the candidate list spans both parities of the trace count, because contracting with the
#   odd-rank ε changes a diagram's size by an odd amount (e.g. Λ²C³ ≅ C³).
# * :O is computed on the identity component, so it really gives SO  content; 
#   ε is excluded for :O (it flips sign under reflections).
# * Totally antisymmetric contractions serve :SL/:SU, :Sp and :SO. For :SL/:SO they are always the full-rank ε.
#   For :Sp the invariant on 2k antisymmetrised slots is ω^∧k, which is ε only at full valence 2k = dim (ε ∝ ω∧...∧ω);
#   at lower valence it prints as ω or ω^∧k, not ε.
#
############################################################



############## partitions & dimension formulas ###################
function integer_partitions(n::Int, max_part::Int=n)
    n == 0 && return [Int[]]
    out = Vector{Vector{Int}}()
    for k in min(n, max_part):-1:1
        for rest in integer_partitions(n - k, k)
            push!(out, vcat([k], rest))
        end
    end
    return out
end
function conjugate(p::Vector{Int})
    isempty(p) && return Int[]
    return [count(x -> x > j, p) for j in 0:(p[1]-1)]
end
"Dimension of the GL(m) irrep with Young diagram `lam`."
function gl_dim(lam::Vector{Int}, m::Int)
    isempty(lam) && return 1
    length(lam) > m && return 0
    conj = conjugate(lam)
    num = 1
    den = 1
    for i in 1:length(lam)
        for j in 1:lam[i]
            num *= (m + j - i)
            den *= lam[i] - j + conj[j] - i + 1
        end
    end
    return num ÷ den
end
# Shared core of the Weyl dimension formulas:
# Π_{i<j} (l_i²-l_j²)/(m_i²-m_j²), optionally times Π_i l_i/m_i.
function weyl_dim_ratio(l, m, with_singles::Bool)
    num = 1 // 1
    den = 1 // 1
    for i in 1:length(l), j in (i+1):length(l)
        num *= (l[i] - l[j]) * (l[i] + l[j])
        den *= (m[i] - m[j]) * (m[i] + m[j])
    end
    if with_singles
        for i in 1:length(l)
            num *= l[i]
            den *= m[i]
        end
    end
    return num / den
end
"Dimension of the Sp(2r) irrep labelled by `mu`."
function sp_dim(mu::Vector{Int}, r::Int)
    mupad = vcat(mu, zeros(Int, r - length(mu)))
    l = [mupad[i] + r - i + 1 for i in 1:r]
    m = [r - i + 1 for i in 1:r]
    return Int(weyl_dim_ratio(l, m, true))
end

"""
    o_dim(mu, dim)

Dimension of the O(dim)/SO(dim) irrep labelled by `mu`. For even dim, a
label with dim/2 rows describes a self-dual/anti-self-dual pair and counts
both (hence the factor 2).
"""
function o_dim(mu::Vector{Int}, dim::Int)
    r = dim ÷ 2
    mupad = vcat(mu, zeros(Int, r - length(mu)))
    if dim % 2 == 0
        l = [mupad[i] + r - i for i in 1:r]
        m = [r - i for i in 1:r]
        val = weyl_dim_ratio(l, m, false)
        (length(mu) == r && !isempty(mu) && mu[end] != 0) && (val *= 2)
        return Int(val)
    else
        l = [mupad[i] + r - i + 1 // 2 for i in 1:r]
        m = [r - i + 1 // 2 for i in 1:r]
        return Int(weyl_dim_ratio(l, m, true))
    end
end

fits_in_box(p::Vector{Int}, rows::Int, cols::Int) =
    length(p) <= rows && (isempty(p) || p[1] <= cols)

# Cauchy-block tableaux

"""
    dual_tableaux(lam)

A matched pair of tableaux (spinor shape `lam` and flavour shape `lam'`)
sharing the labels 1..n, one label per θ copy. `(nothing, nothing)` at n=0.
"""
function dual_tableaux(lam::Vector{Int})
    n = sum(lam; init=0)
    n == 0 && return nothing, nothing
    spinor_tab = YoungTableau(copy(lam), collect(1:n), n)
    lamp = conjugate(lam)
    flavour_labels = Int[]
    for c in 1:(isempty(lam) ? 0 : lam[1])
        for row in spinor_tab.tableau
            length(row) >= c && push!(flavour_labels, row[c])
        end
    end
    return spinor_tab, YoungTableau(lamp, flavour_labels, n)
end

row_groups_rows(rows::Vector{Vector{Int}}) = [row for row in rows if length(row) > 1]
function column_groups_rows(rows::Vector{Vector{Int}})
    cols = Vector{Int}[]
    for c in 1:(isempty(rows) ? 0 : maximum(length.(rows); init=0))
        col = Int[row[c] for row in rows if length(row) >= c]
        length(col) > 1 && push!(cols, col)
    end
    return cols
end
column_groups(tab::YoungTableau) = column_groups_rows(tab.tableau)

############## branching by exact character comparison ####################
#
# Everything below runs in exact rational arithmetic at deterministic prime
# sample points. Numerical error is impossible by construction.

"First `count` primes, via a sieve."
function primes_upto_count(count::Int)
    count <= 6 && return Int[2, 3, 5, 7, 11, 13][1:max(count, 0)]
    limit = ceil(Int, count * (log(count) + log(log(count)))) + 10
    while true
        sieve = trues(limit + 1)
        sieve[1] = false
        for i in 2:isqrt(limit)
            if sieve[i]
                for j in (i*i):i:limit
                    sieve[j] = false
                end
            end
        end
        primes = findall(sieve)
        length(primes) >= count && return primes[1:count]
        limit *= 2
    end
end

"`num` sample points in Rational{BigInt}^r built from consecutive primes.
`offset` shifts which primes are used so separate calls never overlap."
function prime_sample_points(num::Int, r::Int, offset::Int=0)
    primes = primes_upto_count(num * r + offset + 5)
    points = Vector{Vector{Rational{BigInt}}}()
    idx = offset
    for _ in 1:num
        pt = Rational{BigInt}[]
        for _ in 1:r
            idx += 1
            push!(pt, Rational{BigInt}(primes[idx]))
        end
        push!(points, pt)
    end
    return points
end

"Exact determinant (Bareiss fraction-free elimination, O(n³))."
function exact_det(M::Base.Matrix{Rational{BigInt}})
    n = size(M, 1)
    n == 0 && return Rational{BigInt}(1)
    A = copy(M)
    sign = 1
    prev_pivot = Rational{BigInt}(1)
    for k in 1:(n-1)
        if A[k, k] == 0
            swap_row = 0
            for i in (k+1):n
                if A[i, k] != 0
                    swap_row = i
                    break
                end
            end
            swap_row == 0 && return Rational{BigInt}(0)
            A[[k, swap_row], :] = A[[swap_row, k], :]
            sign = -sign
        end
        for i in (k+1):n, j in (k+1):n
            A[i, j] = (A[i, j] * A[k, k] - A[i, k] * A[k, j]) / prev_pivot
        end
        prev_pivot = A[k, k]
    end
    return sign * A[n, n]
end

"Schur polynomial of `lam` at the point `xs`, as a ratio of determinants."
function schur_poly_exact(lam::Vector{Int}, xs::Vector{Rational{BigInt}})
    m = length(xs)
    lampad = vcat(lam, zeros(Int, m - length(lam)))
    A = Rational{BigInt}[xs[i]^(lampad[j] + m - j) for i in 1:m, j in 1:m]
    V = Rational{BigInt}[xs[i]^(m - j) for i in 1:m, j in 1:m]
    return exact_det(A) / exact_det(V)
end

"Sp(2r) Weyl character. Throws `DivideError` at a degenerate sample point;
callers catch it and retry with fresh points."
function sp_character_exact(mu::Vector{Int}, xs::Vector{Rational{BigInt}})
    r = length(xs)
    mupad = vcat(mu, zeros(Int, r - length(mu)))
    entry(x, k) = x^k - x^(-k)
    A = Rational{BigInt}[entry(xs[i], mupad[j] + r - j + 1) for i in 1:r, j in 1:r]
    D = Rational{BigInt}[entry(xs[i], r - j + 1) for i in 1:r, j in 1:r]
    return exact_det(A) / exact_det(D)
end

"""
    o_character_exact(mu, xs, dim)

O(dim)/SO(dim) Weyl character. For odd dim, `xs` plays the role of
y = sqrt(x): the half-integer exponents are tracked as odd integers 2k+1,
so no fractional arithmetic is ever needed. Same DivideError contract as
`sp_character_exact`.
"""
function o_character_exact(mu::Vector{Int}, xs::Vector{Rational{BigInt}}, dim::Int)
    r = length(xs)
    mupad = vcat(mu, zeros(Int, r - length(mu)))
    if dim % 2 == 0
        entry_even = (x, k) -> k == 0 ? Rational{BigInt}(2) : x^k + x^(-k)
        A = Rational{BigInt}[entry_even(xs[i], mupad[j] + r - j) for i in 1:r, j in 1:r]
        D = Rational{BigInt}[entry_even(xs[i], r - j) for i in 1:r, j in 1:r]
        val = exact_det(A) / exact_det(D)
        (length(mu) == r && mu[end] != 0) && (val *= 2)
        return val
    else
        entry_odd = (y, k2) -> y^k2 - y^(-k2)   # k2 = 2k+1
        A = Rational{BigInt}[entry_odd(xs[i], 2 * (mupad[j] + r - j) + 1) for i in 1:r, j in 1:r]
        D = Rational{BigInt}[entry_odd(xs[i], 2 * (r - j) + 1) for i in 1:r, j in 1:r]
        return exact_det(A) / exact_det(D)
    end
end

"Exact Gauss-Jordan solve of a square rational system."
function exact_solve(A::Base.Matrix{Rational{BigInt}}, b::Vector{Rational{BigInt}})
    n = size(A, 1)
    M = hcat(A, b)
    for col in 1:n
        piv = 0
        for row in col:n
            if M[row, col] != 0
                piv = row
                break
            end
        end
        piv == 0 && error("exact_solve: singular system")
        M[[col, piv], :] = M[[piv, col], :]
        M[col, :] ./= M[col, col]
        for row in 1:n
            if row != col && M[row, col] != 0
                M[row, :] .-= M[row, col] .* M[col, :]
            end
        end
    end
    return M[:, end]
end

"Exact least squares via the normal equations A'A x = A'b. The system is
over-determined on purpose, so one unlucky sample point cannot make it
singular."
function exact_normal_equations_solve(A::Base.Matrix{Rational{BigInt}}, b::Vector{Rational{BigInt}})
    At = permutedims(A)
    return exact_solve(At * A, At * b)
end

"Candidate labels: partitions of n, n-kstep, n-2·kstep, ... with at most r rows."
branch_candidates(n::Int, r::Int, kstep::Int) =
    Vector{Int}[mu for k in 0:kstep:n for mu in integer_partitions(n - k) if length(mu) <= r]

"""
    branch_via_characters(name, lam, dim, r, candidates, charfn, full_eigs, dimfn; base_offset)

Shared engine for the two branchings. Fits the Schur function of `lam`
(evaluated on `full_eigs(xs)`) against the candidate characters at exact
sample points, then verifies everything: the fit must be exact at every
point used AND at one independent extra point, the multiplicities must be
whole numbers, and their dimensions must sum to gl_dim. Errors on any failure.
"""
function branch_via_characters(name::String, lam::Vector{Int}, dim::Int, r::Int,
                               candidates::Vector{Vector{Int}}, charfn, full_eigs, dimfn;
                               base_offset::Int=0)
    isempty(candidates) && return Dict{Vector{Int},Int}()
    num_samples = length(candidates) + 8
    A = Base.Matrix{Rational{BigInt}}(undef, num_samples, length(candidates))
    b = Vector{Rational{BigInt}}(undef, num_samples)
    samples = Vector{Vector{Rational{BigInt}}}()
    attempt = 0
    while true
        samples = prime_sample_points(num_samples, r, base_offset + attempt * num_samples * r)
        ok = true
        try
            for (si, xs) in enumerate(samples)
                b[si] = schur_poly_exact(lam, full_eigs(xs))
                for (ci, mu) in enumerate(candidates)
                    A[si, ci] = charfn(mu, xs)
                end
            end
        catch e
            e isa DivideError || rethrow()
            ok = false
        end
        ok && break
        attempt += 1
        attempt > 50 && error("$name: could not find non-degenerate sample points after 50 attempts for lam=$lam, dim=$dim - this should not happen mathematically; please report.")
    end
    coeffs = exact_normal_equations_solve(A, b)
    for (si, xs) in enumerate(samples)
        rhs = sum(coeffs[ci] * charfn(mu, xs) for (ci, mu) in enumerate(candidates))
        rhs != b[si] && error("$name: exact fit failed at a sample point for lam=$lam, dim=$dim: refusing to trust result")
    end
    result = Dict{Vector{Int},Int}()
    for (mu, c) in zip(candidates, coeffs)
        c.den != 1 && error("$name: non-integer coefficient $c for mu=$mu, lam=$lam, dim=$dim: refusing to trust result")
        rc = Int(c.num)
        rc != 0 && (result[mu] = rc)
    end
    check_attempt = 0
    lhs = Rational{BigInt}(0)
    rhs = Rational{BigInt}(0)
    while true
        xs_check = prime_sample_points(1, r, base_offset + (attempt + 1) * num_samples * r + check_attempt * r)[1]
        ok = true
        try
            lhs = schur_poly_exact(lam, full_eigs(xs_check))
            rhs = sum(mult * charfn(mu, xs_check) for (mu, mult) in result; init=Rational{BigInt}(0))
        catch e
            e isa DivideError || rethrow()
            ok = false
        end
        ok && break
        check_attempt += 1
        check_attempt > 50 && error("$name: could not find a non-degenerate verification point for lam=$lam, dim=$dim: please report.")
    end
    lhs != rhs && error("$name: independent verification point failed for lam=$lam, dim=$dim: refusing to trust result")
    total_dim = sum(mult * dimfn(mu) for (mu, mult) in result; init=0)
    true_dim = gl_dim(lam, dim)
    total_dim != true_dim && error("$name: DIMENSION CHECK FAILED for lam=$lam, dim=$dim: branch dims sum to $total_dim, expected $true_dim. Refusing to trust result.")
    return result
end

"GL(dim) -> Sp(dim) branching multiplicities, self-verified."
function branch_gl_to_sp(lam::Vector{Int}, dim::Int)
    dim % 2 == 0 || error("branch_gl_to_sp: dim must be even (Sp(dim) needs dim=2r), got dim=$dim")
    r = dim ÷ 2
    n = sum(lam; init=0)
    return branch_via_characters("branch_gl_to_sp", lam, dim, r, branch_candidates(n, r, 2),
        (mu, xs) -> sp_character_exact(mu, xs),
        xs -> vcat(xs, [1 / x for x in xs]),
        mu -> sp_dim(mu, r))
end

"""
    branch_gl_to_o(lam, dim)

GL(dim) -> O(dim) branching multiplicities, self-verified. At odd dim the
candidate list spans both parities of the trace count, because ε-folding
of long labels changes a diagram's size by an odd amount (Λ²C³ ≅ C³ shows
up this way). Note: characters are sampled on the identity component, so
the result is really SO content.
"""
function branch_gl_to_o(lam::Vector{Int}, dim::Int)
    r = dim ÷ 2
    n = sum(lam; init=0)
    full_eigs = iseven(dim) ?
        (xs -> vcat(xs, [1 / x for x in xs])) :
        (xs -> begin
            xs2 = [y^2 for y in xs]
            vcat(xs2, [1 / x for x in xs2], [Rational{BigInt}(1)])
        end)
    return branch_via_characters("branch_gl_to_o", lam, dim, r,
        branch_candidates(n, r, iseven(dim) ? 2 : 1),
        (mu, xs) -> o_character_exact(mu, xs, dim), full_eigs,
        mu -> o_dim(mu, dim); base_offset=100_000)
end

################# explicit contraction witnesses ####################
function skew_cells(lam::Vector{Int}, mu::Vector{Int})
    cells = Tuple{Int,Int}[]
    for i in 0:(length(lam)-1)
        lo = i < length(mu) ? mu[i+1] : 0
        for j in lo:(lam[i+1]-1)
            push!(cells, (i, j))
        end
    end
    return cells
end
"""
    lr_witnesses(lam, mu, delta, max_count)

Up to `max_count` Littlewood-Richardson fillings of the skew shape
`lam/mu` with content `delta`. Only used to produce EXAMPLE contractions
for branches whose existence and multiplicity are already certified by the
character fit. A failed search here can never change any count.
"""
function lr_witnesses(lam::Vector{Int}, mu::Vector{Int}, delta::Vector{Int}, max_count::Int)
    if length(mu) > length(lam) || any(mu[i] > lam[i] for i in 1:length(mu))
        return Dict{Tuple{Int,Int},Int}[]
    end
    cells = skew_cells(lam, mu)
    length(cells) != sum(delta; init=0) && return Dict{Tuple{Int,Int},Int}[]
    K = length(delta)
    filled = Dict{Tuple{Int,Int},Int}()
    counts_used = zeros(Int, K + 1)
    found = Dict{Tuple{Int,Int},Int}[]
    function lattice_ok_full()
        word = Int[]
        for i in 0:(length(lam)-1)
            lo = i < length(mu) ? mu[i+1] : 0
            for j in (lam[i+1]-1):-1:lo
                push!(word, filled[(i, j)])
            end
        end
        counts = zeros(Int, K + 1)
        for v in word
            counts[v+1] += 1
            for k in 1:(v-1)
                counts[k+1] < counts[k+2] && return false
            end
        end
        return true
    end
    function backtrack(idx::Int)
        length(found) >= max_count && return
        if idx > length(cells)
            lattice_ok_full() && push!(found, copy(filled))
            return
        end
        i, j = cells[idx]
        lo_val = 1
        left_bound = i < length(mu) ? mu[i+1] : 0
        j - 1 >= left_bound && (lo_val = max(lo_val, filled[(i, j - 1)]))
        above_lo = 1
        if i - 1 >= 0
            mu_im1 = (i - 1) < length(mu) ? mu[i] : 0
            mu_im1 <= j < lam[i] && (above_lo = filled[(i - 1, j)] + 1)
        end
        for v in max(lo_val, above_lo):K
            length(found) >= max_count && return
            counts_used[v+1] >= delta[v] && continue
            filled[(i, j)] = v
            counts_used[v+1] += 1
            backtrack(idx + 1)
            counts_used[v+1] -= 1
            delete!(filled, (i, j))
        end
    end
    backtrack(1)
    return found
end

"Reads a witness filling with content (j,j) as ω-pairs: it is a balanced
bracket word, so the pairing is unambiguous."
function decode_elementary_antisymmetric_pairs(filled::Dict{Tuple{Int,Int},Int},
        lam::Vector{Int}, mu::Vector{Int}, tab_rows::Vector{Vector{Int}}, j::Int)
    word_cells = Tuple{Int,Int}[]
    for i in 0:(length(lam)-1)
        lo = i < length(mu) ? mu[i+1] : 0
        for c in (lam[i+1]-1):-1:lo
            push!(word_cells, (i, c))
        end
    end
    stack = Tuple{Int,Int}[]
    pairs = Vector{Vector{Int}}()
    for cell in word_cells
        haskey(filled, cell) || continue
        v = filled[cell]
        if v == 1
            push!(stack, cell)
        elseif v == 2
            i1, jx1 = pop!(stack)
            i2, jx2 = cell
            push!(pairs, [tab_rows[i1+1][jx1+1], tab_rows[i2+1][jx2+1]])
        end
    end
    return pairs
end

"""
Reads a witness filling as η-pairs: cells sharing a value are contracted.
"""
function decode_elementary_symmetric_pairs(filled::Dict{Tuple{Int,Int},Int},
        tab_rows::Vector{Vector{Int}}, j::Int)
    pairs = Vector{Vector{Int}}()
    isempty(filled) && return pairs
    for v in 1:maximum(values(filled))                    
        cells = sort([c for (c, val) in filled if val == v])  # (row, col) order
        isempty(cells) && continue
        for t in 1:2:(length(cells) - 1)                   
            i1, c1 = cells[t]
            i2, c2 = cells[t+1]
            push!(pairs, sort([tab_rows[i1+1][c1+1], tab_rows[i2+1][c2+1]]))
        end
    end
    return pairs
end

"""true iff every label of `filled` occupies a single row.
A δ-row of length 2a confined to one row of λ/μ sits on mutually
symmetrised boxes, so the symmetrised η^{⊗a} collapses to a bare product
of a metrics and the filling is a legitimate witness."""
function filling_row_confined(filled::Dict{Tuple{Int,Int},Int})
    seen = Dict{Int,Int}()
    for ((i, _), v) in filled
        haskey(seen, v) && seen[v] != i && return false
        seen[v] = i
    end
    return true
end
"""the δ partitions the η search may use. Littlewood's rule
sums over all δ with even rows; the pairwise (2,...,2) is only one of them
and is the only one that encodes a pairing directly, so it is tried first.
The others are admitted only through the row-confined test above."""
function symmetric_delta_candidates(j::Int)
    j == 0 && return Vector{Int}[Int[]]
    pairwise = fill(2, j)
    out = Vector{Int}[pairwise]
    for d in integer_partitions(2 * j)
        (all(iseven, d) && d != pairwise) || continue
        push!(out, d)
    end
    return out
end

"Every way to remove `f` boxes from `mu`, at most one per row. The places
where an ε's leftover free slots can sit inside the coefficient."
function vertical_strips_removed(mu::Vector{Int}, f::Int)
    f == 0 && return [copy(mu)]
    results = Vector{Vector{Int}}()
    L = length(mu)
    f > L && return results
    choice = zeros(Int, L)
    function rec(i::Int, left::Int)
        if i > L
            left == 0 || return
            sigma = [mu[t] - choice[t] for t in 1:L]
            for t in 1:L
                sigma[t] < 0 && return
                t < L && sigma[t] < sigma[t+1] && return
            end
            push!(results, [p for p in sigma if p > 0])
            return
        end
        rec(i + 1, left)
        if left > 0 && mu[i] >= 1
            choice[i] = 1
            rec(i + 1, left - 1)
            choice[i] = 0
        end
        return
    end
    rec(1, f)
    return results
end

########### Grassmann structures ###################

"""
    SideData(dim, kind, tableau, residual_rows, residual_shape, pair_groups, is_composite[, eps_groups])

One explicit term for one side (spinor or flavour) of a Cauchy block:
which θ copies its ε's and pairs consume, and which labels the coefficient
keeps (`residual_rows`, shape `residual_shape`). `residual_rows === nothing`
marks a composite (the piece is certified to exist with the given shape and
dimension, but no simple explicit contraction was found for it). NEGATIVE
labels are ε free slots carried by the coefficient (rendered β on the
spinor side, m on the flavour side; each pairs up with the matching ε).
The 7-argument constructor leaves `eps_groups` empty.
"""
struct SideData
    dim::Int
    kind::Symbol
    tableau::Union{Any,Nothing}
    residual_rows::Union{Vector{Vector{Int}},Nothing}
    residual_shape::Vector{Int}
    pair_groups::Vector{Vector{Int}}
    is_composite::Bool
    eps_groups::Vector{Vector{Int}}
end

SideData(dim, kind, tableau, residual_rows, residual_shape, pair_groups, is_composite) =
    SideData(dim, kind, tableau, residual_rows, residual_shape, pair_groups, is_composite, Vector{Int}[])

"The θ-label groups an ε consumes: full-height columns for :SL, the
recorded `eps_groups` for :Sp/:SO, nothing for :O (where ε is not
invariant) or other kinds."
function side_epsilon_columns(side::SideData)
    side.tableau === nothing && return Vector{Int}[]
    side.kind == :SL && return [c for c in column_groups(side.tableau) if length(c) == side.dim]
    side.kind in (:Sp, :SO) && return side.eps_groups
    return Vector{Int}[]
end

struct GrassmannStructure
    spinor::Union{SideData,Nothing}
    flavour::Union{SideData,Nothing}
end

"""
    sp_elementary_terms(lam, mu, tab, dim, mult)

:Sp witness search for one certified branch. Sp has one invariant ω^k for
every even valence 2k <= dim, and each such factor must act inside ONE
antisymmetrised column (that antisymmetry is what makes the contraction
nonzero). A strategy therefore assigns each original column c an even
bottom-strip size x_c <= min(height_c, dim), fed to one ω^(x_c/2), subject
to the leftover heights being weakly decreasing, exactly the condition
for the leftover cells to form a genuine sub-tableau whose columns are the
top parts of original columns, so the remaining ω-pair search stays
faithful to the symmetry structure. This covers full column strips
(x_c = height_c), odd columns (which have no full strip), and the pair
enumeration's one blind spot: several ω-pairs inside a single column are
nonzero for Sp yet unreachable as tableau fillings.
"""
function sp_elementary_terms(lam::Vector{Int}, mu::Vector{Int}, tab, dim::Int, mult::Int)
    n = sum(lam; init=0)
    rows0 = tab.tableau
    terms = SideData[]
    isempty(lam) && return terms
    ncols = lam[1]
    heights = conjugate(lam)
    budget = n - sum(mu; init=0)
    # enumerate strip vectors: x_c even, <= min(height, dim, budget),
    # leftover heights weakly decreasing
    xs_list = Vector{Vector{Int}}()
    x = zeros(Int, ncols)
    function rec(c::Int, prev_leftover::Int, total::Int)
        if c > ncols
            push!(xs_list, copy(x))
            return
        end
        h = heights[c]
        for xc in 0:2:min(h, dim, budget - total)
            h - xc > prev_leftover && continue
            x[c] = xc
            rec(c + 1, h - xc, total + xc)
            x[c] = 0
        end
    end
    rec(1, typemax(Int), 0)
    # pure pairs first, then fewer/larger strips
    sort!(xs_list, by = s -> (count(!=(0), s), [-v for v in s]))
    for xv in xs_list
        length(terms) >= mult && break
        two_j = n - sum(xv) - sum(mu; init=0)
        (two_j < 0 || isodd(two_j)) && continue
        j = two_j ÷ 2
        eps_groups = Vector{Int}[]
        for c in 1:ncols
            xv[c] == 0 && continue
            h = heights[c]
            push!(eps_groups, Int[rows0[i][c] for i in (h-xv[c]+1):h])
        end
        lam_a = Int[]
        rows_a = Vector{Vector{Int}}()
        for i in 1:length(lam)
            len = count(c -> heights[c] - xv[c] >= i, 1:lam[i])
            len > 0 || continue
            push!(lam_a, len)
            push!(rows_a, Int[rows0[i][c] for c in 1:len])
        end
        need = mult - length(terms)
        delta = j == 0 ? Int[] : [j, j]
        for filled in lr_witnesses(lam_a, mu, delta, need)
            pairs = j == 0 ? Vector{Int}[] :
                decode_elementary_antisymmetric_pairs(filled, lam_a, mu, rows_a, j)
            res = Vector{Int}[copy(rows_a[i][1:mu[i]]) for i in 1:length(mu)]
            push!(terms, SideData(dim, :Sp, tab, res, mu, pairs, false, eps_groups))
        end
    end
    return terms
end

"""
    elementary_invariant_terms(lam, mu, tab, dim, kind, mult)

Explicit example contractions for one certified branch `lam -> mu`, never
more than `mult` of them (so this search can never change any count).
:Sp delegates to `sp_elementary_terms` (per-column ω^k strips plus ω-pair
fillings). :O/:SO use η-pair fillings; :SO additionally uses ε, which has
one FIXED valence, dim: a full column of height dim strips entirely, and a
shorter column may borrow ε partially, leaving the unused slots as extra
antisymmetric indices on the coefficient. Whatever the search cannot
resolve falls through to a composite term.

The η search also tries every even-row δ of Littlewood's rule, not
only the pairwise (2,...,2). Non-pairwise δ are accepted only when
row-confined (see `filling_row_confined`), which is exactly when the
symmetrised metrics collapse to a bare product.
"""
function elementary_invariant_terms(lam::Vector{Int}, mu::Vector{Int}, tab,
                                    dim::Int, kind::Symbol, mult::Int)
    kind == :Sp && return sp_elementary_terms(lam, mu, tab, dim, mult)
    n = sum(lam; init=0)
    rows0 = tab.tableau
    terms = SideData[]
    strategies = Tuple{Int,Int}[(0, 0)]   # (full ε strips, partial-ε flag)
    if kind == :SO
        max_full = (!isempty(lam) && length(lam) == dim) ? lam[end] : 0
        for a in 0:max_full, x in 0:1
            (a, x) == (0, 0) && continue
            push!(strategies, (a, x))
        end
        sort!(strategies, by = t -> (t[1] + t[2], t[2]))
    end
    for (a, x) in strategies
        length(terms) >= mult && break

        lam_a = [p - a for p in lam if p > a]
        rows_a = [row[(a+1):end] for row in rows0 if length(row) > a]
        eps_groups = Vector{Int}[Int[rows0[i][c] for i in 1:length(lam) if lam[i] >= c] for c in 1:a]
        consumed = a * length(lam)
        f = 0
        partial_labels = Int[]
        if x == 1
            # partial ε: consume the whole next column (height h < dim),
            # leaving f = dim - h free slots for the coefficient
            h = length(lam_a)
            (h == 0 || h >= dim) && continue
            f = dim - h
            partial_labels = Int[rows_a[i][1] for i in 1:h]
            consumed += h
            lam_a = [p - 1 for p in lam_a if p > 1]
            rows_a = [row[2:end] for row in rows_a if length(row) > 1]
        end

        two_j = n - consumed + f - sum(mu; init=0)
        (two_j < 0 || isodd(two_j)) && continue
        j = two_j ÷ 2

        for sigma in (f > 0 ? vertical_strips_removed(mu, f) : [copy(mu)])
            need = mult - length(terms)
            need <= 0 && break
            # loop over every admissible δ, not just the pairwise (2,...,2)
            for delta in symmetric_delta_candidates(j)
                need = mult - length(terms)
                need <= 0 && break
                pairwise = (delta == fill(2, j))
                for filled in lr_witnesses(lam_a, sigma, delta, need)
                    # a non-pairwise δ only yields a bare product of metrics when its labels are row-confined, otherwise would give a sum of invariants
                    (pairwise || filling_row_confined(filled)) || continue
                    pairs = j == 0 ? Vector{Int}[] :
                        decode_elementary_symmetric_pairs(filled, rows_a, j)
                    sigpad = vcat(sigma, zeros(Int, max(length(mu) - length(sigma), 0)))
                    res = Vector{Int}[]
                    sent = 0
                    for i in 1:length(mu)
                        row = i <= length(rows_a) ? copy(rows_a[i][1:sigpad[i]]) : Int[]
                        if f > 0 && mu[i] == sigpad[i] + 1
                            sent -= 1
                            push!(row, sent)
                        end
                        push!(res, row)
                    end
                    alleps = copy(eps_groups)
                    f > 0 && push!(alleps, vcat(partial_labels, collect(-1:-1:-f)))
                    push!(terms, SideData(dim, kind, tab, res, mu, pairs, false, alleps))
                    length(terms) >= mult && break
                end
            end
        end
    end
    return terms
end

"""
    side_terms(lam, tab, dim, kind)

One or more `SideData` for shape `lam` under `kind`:
* `:none`/`:GL`/`:U` -> one term, no contraction (`:U` behaves as `:GL`:
  its central phase forbids any invariant here).
* `:SL`/`:SU` -> one term, ε on every full-height column (`:SU` behaves as
  `:SL`: same complexified group).
* `:Sp`/`:O`/`:SO` -> one term per certified branch of the group's
  branching, with an explicit ω/η/ε contraction where one exists, else a
  composite. `:Sp` at dim 2 delegates to `:SL`.
Assumes θ carries a single holomorphic index (no barred partner).
"""
function side_terms(lam::Vector{Int}, tab, dim::Int, kind::Symbol)
    kind = kind == :U ? :GL : kind == :SU ? :SL : kind
    kind = (kind == :Sp && dim == 2) ? :SL : kind
    if tab === nothing
        return [SideData(dim, kind, nothing, Vector{Int}[], Int[], Vector{Int}[], false)]
    end
    if kind in (:none, :GL)
        return [SideData(dim, kind, tab, tab.tableau, tab.row_lengths, Vector{Int}[], false)]
    end
    if kind == :SL
        eps_cols = [c for c in column_groups(tab) if length(c) == dim]
        consumed = Set{Int}()
        for c in eps_cols
            union!(consumed, c)
        end
        residual = Vector{Int}[]
        for row in tab.tableau
            free = [k for k in row if !(k in consumed)]
            !isempty(free) && push!(residual, free)
        end
        return [SideData(dim, kind, tab, residual, length.(residual), Vector{Int}[], false)]
    end
    if kind in (:Sp, :O, :SO)
        branch = kind == :Sp ? branch_gl_to_sp(lam, dim) : branch_gl_to_o(lam, dim)
        terms = SideData[]
        for (mu, mult) in branch
            found = elementary_invariant_terms(lam, mu, tab, dim, kind, mult)
            append!(terms, found)
            remaining = mult - length(found)
            if remaining > 0 && isempty(mu)
                #  zero() free indices is unambiguous whatever invariant achieves it
                all_labels = sort(vcat(tab.tableau...))
                for _ in 1:remaining
                    push!(terms, SideData(dim, kind, tab, Vector{Int}[], mu, [all_labels], true))
                end
            elseif remaining > 0
                for _ in 1:remaining
                    push!(terms, SideData(dim, kind, tab, nothing, mu, Vector{Int}[], true))
                end
            end
        end
        return terms
    end
    error("side_terms: :$(kind) is not a supported group type.")
end

"""
    independent_structures(n, d_S, d_F, S_type, F_type)

All independent Grassmann structures at degree n: for every Cauchy-block
shape fitting the d_S*d_F box, every combination of a spinor term and a
flavour term. Types per side: :none, :GL, :U, :SL, :SU, :O, :SO, :Sp.
"""
function independent_structures(n::Int, d_S::Int, d_F::Int, S_type::Symbol, F_type::Symbol)
    out = GrassmannStructure[]
    if d_F == 0
        lam = fill(1, n)
        fits_in_box(lam, d_S, typemax(Int)) || return out
        spinor_tab, _ = dual_tableaux(lam)
        for sside in side_terms(lam, spinor_tab, d_S, S_type)
            push!(out, GrassmannStructure(sside, nothing))
        end
        return out
    end
    for lam in integer_partitions(n)
        fits_in_box(lam, d_S, d_F) || continue
        spinor_tab, flavour_tab = dual_tableaux(lam)
        sterms = side_terms(lam, spinor_tab, d_S, S_type)
        fterms = side_terms(conjugate(lam), flavour_tab, d_F, F_type)
        for sside in sterms, fside in fterms
            push!(out, GrassmannStructure(sside, fside))
        end
    end
    return out
end
independent_structures(n::Int, d_S::Int, S_type::Symbol) = independent_structures(n, d_S, 0, S_type, :none)

"""
    side_dimension(side)

True dimension of this term's coefficient. For :Sp/:O/:SO the coefficient
is the group's own traceless tensor of the residual shape, so its
dimension comes from `sp_dim`/`o_dim` (generally smaller than the plain
`gl_dim` of the same shape).
"""
function side_dimension(side::Union{SideData,Nothing})
    side === nothing && return 1
    side.kind == :Sp && return sp_dim(side.residual_shape, side.dim ÷ 2)
    side.kind in (:O, :SO) && return o_dim(side.residual_shape, side.dim)
    return gl_dim(side.residual_shape, side.dim)
end

######################## rendering #####################

"""
Positive labels are θ-copy indices (`prefix * k`). Negative labels are ε
free slots on the coefficient, drawn with a side-specific dummy letter so
the two sides never collide: spinor (α) slots become β1, β2, ...; flavour
(i) slots become m1, m2, ... Each appears twice, on C and inside its ε.
"""
const SLOT_PREFIX = Dict("α" => "β", "i" => "m")
idx_label(prefix::AbstractString, k::Int) =
    k > 0 ? string(prefix, k) : string(get(SLOT_PREFIX, prefix, "m" * prefix), -k)

"""
    antisym_invariant_name(len, kind, dim)

Name of the totally antisymmetric invariant consuming `len` antisymmetrised
slots. For :Sp that invariant is ω^∧(len/2). It is totally antisymmetric at
any valence, but it equals the volume form ε ONLY at full valence
len == dim, where ε ∝ ω^∧(dim/2); at len < dim it is ω (len=2) or a lower
wedge power, and naming it ε would misstate its rank. :SL and :SO groups are
always full length dim (an :SO partial ε carries its unused slots as
negative labels), so those are genuine ε.
"""
function antisym_invariant_name(len::Int, kind::Symbol, dim::Int)
    (kind == :Sp && len != dim) || return "ε"
    return len == 2 ? "ω" : "ω^{∧$(len ÷ 2)}"
end

function epsilon_symbol(indices, prefix, upper=true, kind::Symbol=:SL, dim::Int=length(indices))
    body = join([idx_label(prefix, k) for k in indices], "")
    sym = antisym_invariant_name(length(indices), kind, dim)
    occursin('^', sym) && (sym = "($sym)")   # avoid a double superscript
    return upper ? "$(sym)^{$body}" : "$(sym)_{$body}"
end

"ω (:Sp) or η (:O/:SO) pair; Ω for a fully-contracted composite, where
`pair` lists every consumed label."
function pair_symbol(pair::Vector{Int}, prefix::String, kind::Symbol, composite::Bool=false)
    body = join([idx_label(prefix, k) for k in pair], "")
    composite && return "Ω_{$body}"
    return (kind == :Sp ? "ω" : "η") * "_{" * body * "}"
end

"Labels the coefficient still carries (negatives are ε slots); `nothing`
only for a composite whose explicit labels are unknown."
function free_indices_structure(side::Union{SideData,Nothing}, n::Int)
    side === nothing && return Int[]
    side.tableau === nothing && return collect(1:n)
    side.residual_rows === nothing && return nothing
    return sort(vcat(Int[], side.residual_rows...))
end

"Symmetrization notation for the coefficient: one row -> (sym), one column
-> [antisym], mixed -> each row in its own parentheses."
function symmetrized_str(side::Union{SideData,Nothing}, prefix::String)
    (side === nothing || side.residual_rows === nothing || isempty(side.residual_rows)) && return ""
    rows = side.residual_rows
    length(rows) == 1 &&
        return "(" * join([idx_label(prefix, k) for k in rows[1]], "") * ")"
    all(r -> length(r) == 1, rows) &&
        return "[" * join([idx_label(prefix, k) for k in vcat(rows...)], "") * "]"
    return join(["(" * join([idx_label(prefix, k) for k in row], "") * ")" for row in rows], "")
end

"All invariant-tensor factors on one side: ε groups, then ω/η/Ω."
function invariant_factors(side::Union{SideData,Nothing}, prefix::String, upper::Bool)
    (side === nothing || side.tableau === nothing) && return String[]
    out = String[epsilon_symbol(col, prefix, upper, side.kind, side.dim)
                 for col in side_epsilon_columns(side)]
    for pair in side.pair_groups
        push!(out, pair_symbol(pair, prefix, side.kind, side.is_composite))
    end
    return out
end

"The coefficient field: C with its symmetrized free indices, empty when
fully contracted, or a composite marker."
function coefficient_str(s::GrassmannStructure)
    n = sum(s.spinor.tableau.row_lengths; init=0)
    free_alpha = free_indices_structure(s.spinor, n)
    free_i = s.flavour === nothing ? Int[] : free_indices_structure(s.flavour, n)
    (free_alpha === nothing || free_i === nothing) && return "C(composite, unresolved)"
    (isempty(free_alpha) && isempty(free_i)) && return ""
    up = isempty(free_alpha) ? "" : "^{" * symmetrized_str(s.spinor, "α") * "}"
    dn = isempty(free_i) ? "" : "_{" * symmetrized_str(s.flavour, "i") * "}"
    return "C" * up * dn
end

"Coefficient tableau and true dimension for one side; \"traceless\" flags a
dimension below the plain gl_dim of the shown shape. Empty for a scalar."
function coeff_note(side::Union{SideData,Nothing}, prefix::String)
    side === nothing && return ""
    d = side_dimension(side)
    side.residual_rows === nothing &&
        return "$prefix: shape $(side.residual_shape) dim $d (unresolved)"
    isempty(side.residual_shape) && return ""
    tless = (side.kind in (:Sp, :O, :SO) && gl_dim(side.residual_shape, side.dim) != d) ?
        ", traceless" : ""
    return "$prefix: $(side.residual_rows) dim $d$tless"
end

"True iff at least one side of `s` carries an actual invariant contraction."
function has_forced_invariant(s::GrassmannStructure)
    hit(side) = side !== nothing && side.tableau !== nothing &&
        (!isempty(side_epsilon_columns(side)) || !isempty(side.pair_groups) || side.is_composite)
    return hit(s.spinor) || hit(s.flavour)
end

"""
    render_expression(s)

Fully explicit rendering of one structure: shape labels, coefficient,
invariant tensors, all θ factors, and per-side coefficient annotations.
"""
function render_expression(s::GrassmannStructure)
    (s.spinor === nothing || s.spinor.tableau === nothing) && return "1"
    n = sum(s.spinor.tableau.row_lengths; init=0)
    tensors = vcat(invariant_factors(s.spinor, "α", true),
                   invariant_factors(s.flavour, "i", false))
    thetas = ["θ^{i$k}_{α$k}" for k in 1:n]
    coeff = coefficient_str(s)
    notes = filter(!isempty, [coeff_note(s.spinor, "α"), coeff_note(s.flavour, "i")])
    ann = isempty(notes) ? "" : "   [coeff " * join(notes, "; ") * "]"
    lbl(side, name) = (side === nothing || side.tableau === nothing) ? "" :
        " " * name * string(side.tableau.row_lengths) *
        (side.residual_shape == side.tableau.row_lengths ? "" : "->" * string(side.residual_shape))
    shape = strip(lbl(s.spinor, "spinor") * lbl(s.flavour, "flavour"))
    body = join(vcat(coeff == "" ? String[] : [coeff], tensors, thetas), " ")
    return "[$shape]  " * body * ann
end

"""
    render_compact(n, d_S, d_F, S_type, F_type) -> Vector{String}

Compact expansion at degree n: a header with the expected dimension
binomial(d_S·N, n), then one line per term showing the invariant tensors,
the coefficient, and its tableau with the true dimension. If any block has
no reduction at all, the whole degree collapses to the single generic
term. A MISMATCH line is appended if the dimensions fail to add up.
Print with `foreach(println, render_compact(...))`.
"""
function render_compact(n::Int, d_S::Int, d_F::Int, S_type::Symbol, F_type::Symbol)
    expected = binomial(d_S * d_F, n)
    lines = String["n=$n  (dim $expected)"]
    n == 0 && return push!(lines, "  1")
    structures = independent_structures(n, d_S, d_F, S_type, F_type)
    isempty(structures) && return lines
    if !all(has_forced_invariant, structures)
        push!(lines, "  C^{α1..α$n}_{i1..i$n}  (generic Λ^$n(V⊗W) term, no reduction at this degree)")
        return lines
    end
    total = 0
    for s in structures
        d = side_dimension(s.spinor) * side_dimension(s.flavour)
        total += d
        parts = vcat(let c = coefficient_str(s); c == "" ? String[] : [c] end,
                     invariant_factors(s.spinor, "α", true),
                     invariant_factors(s.flavour, "i", false))
        notes = filter(!isempty, [coeff_note(s.spinor, "α"), coeff_note(s.flavour, "i")])
        note = isempty(notes) ? "scalar" : join(notes, " × ")
        push!(lines, "  " * (isempty(parts) ? "1" : join(parts, " ")) * "   [" * note * "]  dim $d")
    end
    total != expected &&
        push!(lines, "  !! dimension MISMATCH: total $total vs expected $expected: real bug, please report")
    return lines
end

#################### explanation ###############################

"Prints a plain-language account of one side's contractions."
function explain_side(side::SideData, label::String)
    if side.tableau === nothing
        println("  ($label: none)")
        return
    end
    println("  $label tableau  ", side.tableau.tableau, "  [", side.kind, "]")
    for col in side_epsilon_columns(side)
        sym = antisym_invariant_name(length(col), side.kind, side.dim)
        if any(x -> x < 0, col)
            println("    $label $sym contracts θ-copies (negative entries = free coefficient slots): ", col)
        else
            println("    $label $sym contracts θ-copies: ", col)
        end
    end
    for pair in side.pair_groups
        println("    $label ", side.kind == :Sp ? "ω" : "η", " contracts θ-copies: ", pair)
    end
    side.is_composite &&
        println("    $label composite pairing (shape $(side.residual_shape), dim $(side_dimension(side)), unresolved)")
    if side.kind == :SL && side.residual_rows !== nothing
        for col in column_groups_rows(side.residual_rows)
            println("    $label antisymmetrizes θ-copies (free indices): ", col)
        end
        for row in row_groups_rows(side.residual_rows)
            println("    $label symmetrizes θ-copies (free indices): ", row)
        end
    end
end

function explain_structure(s::GrassmannStructure)
    println("cauchy block")
    if s.spinor === nothing || s.spinor.tableau === nothing
        println("  (degenerate n=0 case: scalar, no tableau)")
    else
        explain_side(s.spinor, "spinor")
    end
    s.flavour === nothing ? println("  (no flavour)") : explain_side(s.flavour, "flavour")
end