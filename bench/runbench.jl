# coverage + timing for the simplification API (not tests, nothing fails here)
#
#   julia --project=. bench/runbench.jl
#
# coverage: each case is an expression we already know the answer to. every
# entry point gets applied and the result recorded. what comes out is a matrix:
# a few cases only reduce under one entry point and it's a different one each
# time
#
# timing: canonicalise against the number of contracted dummy pairs, plus how
# many dummies a set hands out before it runs dry
#
# expected values are what xAct/xTensor gives for the same input, written in as
# consts (doesn't shell out to Mathematica)

using Alakazam
using Printf

const ENTRIES = ["simplify" => Alakazam.simplify,
                 "hard_simplify" => Alakazam.hard_simplify,
                 "canonicalise" => Alakazam.canonicalise]

is_zero(v) = try
    isZeroTens(v) || string(v) == "0"
catch
    string(v) == "0"
end

"""One case: a name, a thunk that builds the expression, and what it should reduce to."""
struct Case
    name::String
    build::Function
    expect::String          # "0", a literal like "2", or "" when only timing matters
end

function run_case(c::Case)
    row = String[]
    for (_, f) in ENTRIES
        got = try
            v = f(c.build())
            if c.expect == "0"
                is_zero(v) ? "reduced" : "no"
            else
                string(v) == c.expect ? "reduced" : "no"
            end
        catch
            "ERR"
        end
        push!(row, got)
    end
    row
end

# --------------------------------------------------------------------- cases
surf = IndexSet("surface", 4, lower)
a, b, c, d = Index(["a", "b", "c", "d"], surf)

S(nm, is) = SymmetricTensor(nm, is)
A(nm, is) = AntisymmetricTensor(nm, is)
R(i, j, k, l) = RiemannTensor("R", [i => lower, j => lower, k => lower, l => lower])

CASES = Case[
    Case("symmetric slot swap: h_ab - h_ba",
         () -> S("h", [a => lower, b => lower]) - S("h", [b => lower, a => lower]), "0"),
    Case("antisymmetric slot swap: F_ab + F_ba",
         () -> A("F", [a => lower, b => lower]) + A("F", [b => lower, a => lower]), "0"),
    Case("antisymmetric against symmetric: F^ab S_ab",
         () -> A("F", [a => upper, b => upper]) * S("S", [a => lower, b => lower]), "0"),
    Case("antisymmetric against the metric: F^ab g_ab",
         () -> A("F", [a => upper, b => upper]) * Metric([a => lower, b => lower]), "0"),
    Case("dummy relabelling: A^ab B_ab - A^cd B_cd",
         () -> S("A", [a => upper, b => upper]) * S("B", [a => lower, b => lower]) -
               S("A", [c => upper, d => upper]) * S("B", [c => lower, d => lower]), "0"),
    Case("Riemann pair exchange: R_abcd - R_cdab",
         () -> R(a, b, c, d) - R(c, d, a, b), "0"),
    Case("Riemann first-pair antisymmetry: R_abcd + R_bacd",
         () -> R(a, b, c, d) + R(b, a, c, d), "0"),
    Case("Ricci symmetry: R^a_bad - R^a_dab",
         () -> contract_metrics(
             RiemannTensor("R", [a => upper, b => lower, a => lower, d => lower]) -
             RiemannTensor("R", [a => upper, d => lower, a => lower, b => lower])), "0"),
    Case("metric trace: g^ab g_ab",
         () -> Metric([a => upper, b => upper]) * Metric([a => lower, b => lower]), "4"),
    Case("delta pair trace: d^a_b d^b_a",
         () -> KroneckerDelta([a => upper, b => lower]) *
               KroneckerDelta([b => upper, a => lower]), "4"),
    Case("delta self trace: d^a_a",
         () -> KroneckerDelta([a => upper, a => lower]), "4"),
    Case("first Bianchi: R_abcd + R_acdb + R_adbc",
         () -> R(a, b, c, d) + R(a, c, d, b) + R(a, d, b, c), "0"),
]

println("Alakazam simplification coverage")
println("dimension 4; expected values are those xAct/xTensor gives for the same input\n")
@printf("  %-46s %-13s %-14s %-13s\n", "case", "simplify", "hard_simplify", "canonicalise")
println("  ", "-"^88)

rows = Tuple{Case,Vector{String}}[]
for cs in CASES
    row = run_case(cs)
    push!(rows, (cs, row))
    @printf("  %-46s %-13s %-14s %-13s\n", first(cs.name, 46), row[1], row[2], row[3])
end

println("\n  summary")
none = [c.name for (c, r) in rows if !any(==("reduced"), r)]
one = [c.name for (c, r) in rows if count(==("reduced"), r) == 1]
errs = [(c.name, i) for (c, r) in rows for (i, v) in enumerate(r) if v == "ERR"]
println("  cases reduced by exactly one entry point: $(length(one))")
for n in one
    println("      $n")
end
println("  cases no entry point reduces: $(length(none))")
for n in none
    println("      $n")
end
if !isempty(errs)
    println("  entry points that raised: $(length(errs))")
end
println("""
  first Bianchi is a multiterm symmetry. xAct's ToCanonical doesn't impose it
  either, so "no" there matches xAct.
  a case only one entry point reduces means the caller has to know which one up
  front, and it isn't the same one across cases.""")

# ------------------------------------------------------------------- timing
println("\nCanonicalisation timing against the number of contracted dummy pairs")

tsurf = IndexSet("timing", 4, lower)
TIDX = Index(string.(collect("abcdefghijklmnopqrstuvwx")), tsurf)

function chain(n::Int; reversed::Bool=false)
    ex = 1
    rng = reversed ? (n:-1:1) : (1:n)
    for k in rng
        p, q = reversed ? (TIDX[2k], TIDX[2k-1]) : (TIDX[2k-1], TIDX[2k])
        ex = ex * S("A", [p => upper, q => upper]) * S("B", [p => lower, q => lower])
    end
    ex
end

#tens of microseconds each, so one time() reading is noise (an earlier version
#of this had n=6 coming out faster than n=2). best of REPS with time_ns instead
const REPS = 200

function best_of(f, reps::Int=REPS)
    best = Inf
    for _ in 1:reps
        t0 = time_ns()
        f()
        best = min(best, (time_ns() - t0) / 1e3)   # microseconds
    end
    best
end

best_of(() -> Alakazam.canonicalise(chain(1)), 5)   #warm up so row one isnt compilation
@printf("  %-7s %-14s %-16s %-8s\n", "pairs", "build (us)", "canonicalise (us)", "agree")
for n in 1:12
    tb = best_of(() -> chain(n))
    e1 = chain(n)
    e2 = chain(n; reversed=true)
    tc = best_of(() -> Alakazam.canonicalise(chain(n)))
    agree = Alakazam.canonicalise(e1) == Alakazam.canonicalise(e2)
    @printf("  %-7d %-14.1f %-16.1f %-8s\n", n, tb, tc, string(agree))
end
println("""
  canonicalise is timed on a fresh expression, so its column includes build.
  both orderings land on the same canonical form, so the work happened.""")

# ------------------------------------------------------- index pool headroom
println("\nIndices an index set can supply in total, declared plus generated")

#cap only bounds the runtime. print where the supply stops and pick the closing
#note from that number, since it moves with how generate_indices searches
const SUPPLY_CAP = 4000
const ONE_ALPHABET = length(first(Alakazam.ORDERED_POOLS))
const ALL_POOLS = length(unique(reduce(vcat, Alakazam.ORDERED_POOLS)))

supplies = Int[]
for (label, names) in ("4 declared" => "abcd", "12 declared" => "abcdefghklmn")
    iset = IndexSet("pool_$label", 4, lower)
    Index(string.(collect(names)), iset)
    n = 0
    msg = "hit this benchmark's cap of $SUPPLY_CAP, so at least that many"
    try
        pool = copy(iset.indices)
        while n < SUPPLY_CAP
            safe_pop!(pool, iset)
            n += 1
        end
    catch err
        msg = replace(first(sprint(showerror, err), 90), "\n" => " ")
    end
    push!(supplies, n)
    @printf("  %-14s supplied %5d in total, then: %s\n", label, n, msg)
end

#template name thats in no pool. the Index constructor takes one, so generation
#has to do something with it either way
odd = IndexSet("pool_odd", 4, lower)
Index(["i1", "i2"], odd)
oddmsg = try
    g = Alakazam.generate_indices(odd)
    g === nothing ? "returned nothing" : "generated \"$(g.name)\""
catch err
    replace(first(sprint(showerror, err), 90), "\n" => " ")
end
@printf("  %-14s template in no pool: %s\n", "\"i1\"", oddmsg)

@printf("\n  one alphabet is %d names; ORDERED_POOLS holds %d across %d pools\n",
        ONE_ALPHABET, ALL_POOLS, length(Alakazam.ORDERED_POOLS))
if all(<=(ONE_ALPHABET), supplies)
    println("""
  both stop at one alphabet, so declaring more indices buys no extra dummies.
  generate_indices takes the pool holding the FIRST declared index and draws
  from that one only, so a set whose names fill an alphabet cant hand out
  another dummy even with the other pools untouched.""")
else
    println("""
  supply runs past one alphabet, so generation falls through to the remaining
  pools and the ceiling is the whole of ORDERED_POOLS.""")
end
