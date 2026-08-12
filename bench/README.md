# bench

Coverage and timing for the simplification API.

```
julia --project=. bench/runbench.jl
```

It needs nothing beyond the package's own dependencies and always exits zero.
Nothing here passes or fails, so it does not belong in `runtests.jl`.

## What it reports

**A coverage matrix.** Each row is an expression whose answer is already known.
Every public entry point is applied to it and the outcome recorded, so the
output says which entry point reduces which kind of expression. The count in the
summary is the useful part: a few rows reduce under only one entry point, and it
is a different one each time, so the caller has to know which to reach for.

**Canonicalisation timing** against the number of contracted dummy pairs. Each
figure is the best of `REPS` runs taken with `time_ns`. A single call here is
tens of microseconds, and one `time()` reading of that is noise. Both orderings
of each expression are canonicalised and compared, so the timing covers work
that happened.

**Index supply.** How many indices a set hands out in total, declared plus
generated, and what it does when it runs dry.

## Reading the index supply

That section prints what it measures, since the answer depends on how
`generate_indices` searches `ORDERED_POOLS`. It reports **all 1425 names of
`ORDERED_POOLS`, across 55 pools**, because generation falls through to the
remaining pools once one is exhausted, and raises `OutOfIndices` when every
pool is used up. A template name belonging to no pool is served from the
start of the order, so `"i1"` yields `"a"`.

The script prints the other outcome, **one alphabet, 26**, whenever a set stops
at the pool holding its first declared index. Both notes follow the
measurement, so read the printed one.

## Where the expected values come from

They are what xAct/xTensor gives for the same expressions, written in as
constants. Nothing shells out to Mathematica.

One row is expected to stay unreduced and is not a shortfall: **first Bianchi**,
`R_abcd + R_acdb + R_adbc`. It is a multiterm symmetry, which xAct's
`ToCanonical` does not impose either. Monoterm canonicalisation leaves it alone
and reducing it takes a Young projection.

## Adding a row

Append a `Case` to `CASES` with a name, a thunk that builds the expression, and
the string it should reduce to. `"0"` is checked with `isZeroTens`, anything
else by string comparison against the printed form, so a dimension-valued row in
dimension four expects `"4"`.
