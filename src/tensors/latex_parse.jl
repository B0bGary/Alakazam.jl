###############################################################################
# LT"..." read TeX markup back into tensors.
#
# Accepts what a person would actually type, eg: `A^\mu_\nu`, `F_{\mu\nu}`, `\bar\psi_\alpha`.
###############################################################################

"""
Spellings a person may type that `latexify` never produces, so they cannot be recovered by
inverting its table. Only aliases belong here
"""
const _TEX_EXTRA = Dict{String,String}(
    "hbar" => "ħ",          # latexify writes ħ as \hslash, so \hbar is never produced
)

const _TEX_TO_UNI = let d = Dict{String,String}()
    if isdefined(Latexify, :unicodedict)
        pat = Regex("^" * "\\\\" * "([A-Za-z]+)\$")
        for (uni, tex) in Latexify.unicodedict
            tex isa AbstractString || continue
            m = match(pat, tex)
            m === nothing && continue
            get!(d, String(m.captures[1]), string(uni))
        end
    end
    for (k, v) in _TEX_EXTRA
        get!(d, k, v)
    end
    d
end

"""Accents TeX writes around a name, and the combining mark each stands for. `\\overline` and
`\\bar` are the same mark; so are `\\dot` and `\\Dot`, which `latexify` may emit."""
const _TEX_ACCENTS = Dict{String,Char}(
    "dot" => COMBINING_DOT, "Dot" => COMBINING_DOT,
    "overline" => COMBINING_MACRON, "bar" => COMBINING_MACRON,
    "hat" => COMBINING_HAT, "widehat" => COMBINING_HAT,
    "tilde" => COMBINING_TILDE, "widetilde" => COMBINING_TILDE,
)

"""Spacing and delimiter commands that carry no meaning here and are dropped."""
const _TEX_IGNORED = Set(["quad", "qquad", "left", "right", "!", ",", ";", ":", " ", "displaystyle", "cdot"])

const _TEX_MARK_TO_SUFFIX = Dict(COMBINING_DOT => "dot", COMBINING_MACRON => "bar",
                                 COMBINING_HAT => "hat", COMBINING_TILDE => "tilde")

#------------------------------------------------------------------ character-level reader

mutable struct _TeXReader
    s::String
    i::Int
end
_TeXReader(s::AbstractString) = _TeXReader(String(s), 1)

_eof(r::_TeXReader) = r.i > lastindex(r.s)
_peek(r::_TeXReader) = _eof(r) ? '\0' : r.s[r.i]
_bump!(r::_TeXReader) = (c = r.s[r.i]; r.i = nextind(r.s, r.i); c)

function _skipspace!(r::_TeXReader)
    while !_eof(r) && isspace(_peek(r))
        _bump!(r)
    end
end

"""Read a control word, the backslash having been consumed. `\\alpha` gives `"alpha"`; a control
symbol such as `\\,` gives that single character."""
function _read_command!(r::_TeXReader)
    _eof(r) && return ""
    isletter(_peek(r)) || return string(_bump!(r))
    start = r.i
    while !_eof(r) && isletter(_peek(r))
        _bump!(r)
    end
    return r.s[start:prevind(r.s, r.i)]
end

"""Read the text of a `{...}` group, the opening brace having been consumed, honouring nesting."""
function _read_group!(r::_TeXReader)
    depth = 1
    start = r.i
    while !_eof(r)
        c = _peek(r)
        if c == '{'
            depth += 1
        elseif c == '}'
            depth -= 1
            depth == 0 && break
        end
        _bump!(r)
    end
    txt = r.s[start:prevind(r.s, r.i)]
    _eof(r) || _bump!(r)          # closing brace
    return txt
end

"""The argument of an accent or a script: either a braced group or the single token after it."""
function _read_argument!(r::_TeXReader)
    _skipspace!(r)
    _eof(r) && return ""
    if _peek(r) == '{'
        _bump!(r)
        return _read_group!(r)
    elseif _peek(r) == '[' || _peek(r) == '('
        # an unbraced symmetry group, as the tensor strings write it: A_[ab] as well as A_{[ab]}
        open = _bump!(r)
        close = open == '[' ? ']' : ')'
        start = r.i
        while !_eof(r) && _peek(r) != close
            _bump!(r)
        end
        body = r.s[start:prevind(r.s, r.i)]
        _eof(r) || _bump!(r)
        return string(open) * body * string(close)
    elseif _peek(r) == '\\'
        _bump!(r)
        return "\\" * _read_command!(r)
    else
        return string(_bump!(r))
    end
end

#------------------------------------------------------------------ names

"""
    _tex_name(text) -> String

One TeX name as the unicode string Alakazam stores. Accents become a leading combining mark, which
is the same convention `uni_to_TeX` decodes, so `\\dot{\\alpha}` and `̇α` are the same name.
"""
function _tex_name(text::AbstractString)
    r = _TeXReader(strip(text))
    _skipspace!(r)
    while !_eof(r) && _peek(r) == '$'
        _bump!(r)
        _skipspace!(r)
    end
    _eof(r) && return ""
    if _peek(r) == '\\'
        _bump!(r)
        cmd = _read_command!(r)
        if haskey(_TEX_ACCENTS, cmd)
            return string(_TEX_ACCENTS[cmd]) * _tex_name(_read_argument!(r))
        end
        cmd in _TEX_IGNORED && return _tex_name(r.s[r.i:end])
        return get(_TEX_TO_UNI, cmd, cmd)
    end
    if _peek(r) == '{'                       # a redundant brace group around a name
        _bump!(r)
        return _tex_name(_read_group!(r))
    end
    return string(_bump!(r))
end

"""
    _identifier_spelling(name) -> Union{Symbol,Nothing}

How the index would be written as a Julia variable, so that the caller's scope can be consulted.
"""
function _identifier_spelling(name::AbstractString)
    isempty(name) && return nothing
    c = first(name)
    if haskey(_TEX_MARK_TO_SUFFIX, c)
        cand = chop(name, head=1, tail=0) * _TEX_MARK_TO_SUFFIX[c]
        return Base.isidentifier(cand) ? Symbol(cand) : nothing
    end
    return Base.isidentifier(name) ? Symbol(name) : nothing
end

#------------------------------------------------------------------ index groups

"""
    _parse_script(text, slot0) -> (slots, groups, nextslot)

Read the body of one script into its slots and any symmetry groups written over them.

`{\\mu\\nu}` is two slots, and so is `{ab}`: a run of plain letters is one index each, the rule the
tensor-string parser follows. `?` marks the slot after it as a wildcard.

Round brackets make the enclosed slots symmetric and square brackets antisymmetric, as they are
written in the literature, and a pair of vertical bars holds slots out of the group being
accumulated, so that `R_{a[bc|d|e]}` antisymmetrises over `b`, `c` and `e` while leaving `d` alone.
Slot numbers are counted from `slot0` across the whole tensor, since a tableau indexes the tensor's
slots and not the script's.
"""
function _parse_script(text::AbstractString, slot0::Int)
    slots = Tuple{String,Bool}[]
    groups = Tuple{Bool,Vector{Int}}[]
    r = _TeXReader(text)
    slotnum = slot0
    acc = Int[]
    open_sym = nothing          # true for a row (symmetric), false for a column
    paused = false

    function take!(nm, wild)
        slotnum += 1
        push!(slots, (nm, wild))
        (open_sym === nothing || paused) || push!(acc, slotnum)
    end

    while true
        _skipspace!(r)
        _eof(r) && break
        c = _peek(r)

        if c == '(' || c == '['
            _bump!(r)
            open_sym === nothing ||
                error("LT: a symmetry group is already open in \"$text\"; nested groups are not " *
                      "part of the notation, close the first before opening another.")
            open_sym = (c == '(')
            acc = Int[]
            paused = false
            continue
        elseif c == ')' || c == ']'
            _bump!(r)
            open_sym === nothing &&
                error("LT: `$c` closes a symmetry group that was never opened, in \"$text\".")
            (open_sym == (c == ')')) ||
                error("LT: mismatched symmetry brackets in \"$text\"; a group opened with one " *
                      "kind must close with the same kind.")
            push!(groups, (open_sym, copy(acc)))
            open_sym = nothing
            continue
        elseif c == '|'
            _bump!(r)
            paused = !paused
            continue
        elseif c == '$'
            _bump!(r)          # math delimiters carry no slot
            continue
        end

        wild = false
        if c == '?'
            _bump!(r)
            wild = true
            _skipspace!(r)
            _eof(r) && break
        end

        if _peek(r) == '\\'
            _bump!(r)
            cmd = _read_command!(r)
            if haskey(_TEX_ACCENTS, cmd)
                take!(string(_TEX_ACCENTS[cmd]) * _tex_name(_read_argument!(r)), wild)
            elseif cmd in _TEX_IGNORED
                continue
            else
                take!(get(_TEX_TO_UNI, cmd, cmd), wild)
            end
        elseif _peek(r) == '{'
            _bump!(r)
            inner = _read_group!(r)
            (sl, gr, slotnum) = _parse_script(inner, slotnum)
            for (nm, w) in sl
                push!(slots, (nm, w || wild))
            end
            isempty(gr) ||
                error("LT: a symmetry group inside a nested brace in \"$text\" is not read; " *
                      "write the brackets in the script itself, as in A_{[ab]}.")
        elseif isdigit(_peek(r))
            start = r.i                      # a label such as the 1 of k_1, kept whole
            while !_eof(r) && isdigit(_peek(r))
                _bump!(r)
            end
            take!(r.s[start:prevind(r.s, r.i)], wild)
        else
            take!(string(_bump!(r)), wild)
        end
    end

    open_sym === nothing ||
        error("LT: a symmetry group opened in \"$text\" is never closed. Antisymmetrisation " *
              "written across separate factors, as in \\partial_{[a}A_{b]}, is not a symmetry of " *
              "any one tensor and so cannot be read as one; expand it into terms instead.")

    return (slots, groups, slotnum)
end

#------------------------------------------------------------------ the parse plan
#
# The reader produces a plan rather than either tensors or code, so that the macro and the runtime
# function share every line of the parsing and differ only in how a name is resolved.

struct _LTSlot
    name::String
    wild::Bool
    pos::IndexPosition
end

struct _LTFactor
    name::String
    slots::Vector{_LTSlot}
    coords::Vector{String}
    groups::Vector{Tuple{Bool,Vector{Int}}}   # (symmetric?, the slots it covers)
end

"""
Read one factor: a name, then any number of scripts and an optional coordinate list. Returns
`nothing` when the reader is not looking at a factor.
"""
function _parse_factor!(r::_TeXReader)
    _skipspace!(r)
    _eof(r) && return nothing

    name = ""
    while true
        _skipspace!(r)
        _eof(r) && break
        c = _peek(r)
        if c == '{'
            save = r.i
            _bump!(r)
            body = _read_group!(r)
            if isempty(strip(body))
                r.i = save                   # `{ }` introduces a script, it is not a name
                break
            end
            name = _tex_name(body)
            break
        elseif c == '\\'
            _bump!(r)
            cmd = _read_command!(r)
            if cmd in _TEX_IGNORED
                continue
            elseif haskey(_TEX_ACCENTS, cmd)
                name = string(_TEX_ACCENTS[cmd]) * _tex_name(_read_argument!(r))
                break
            else
                name = get(_TEX_TO_UNI, cmd, cmd)
                break
            end
        elseif c == '(' || c == ')'
            _bump!(r)          # a grouping paren, from \left( ... \right) or plain
            continue
        elseif c == '?' || isletter(c) || isdigit(c)
            start = r.i
            _peek(r) == '?' && _bump!(r)
            _peek(r) == '?' && _bump!(r)     # ??A adopts the matched object's data
            _eof(r) || _bump!(r)
            name = r.s[start:prevind(r.s, r.i)]
            break
        else
            return nothing
        end
    end
    isempty(name) && return nothing

    slots = _LTSlot[]
    coords = String[]
    groups = Tuple{Bool,Vector{Int}}[]
    while true
        save = r.i
        _skipspace!(r)
        if _eof(r)
            r.i = save
            break
        end
        c = _peek(r)
        if c == '{'                                    # the `{ }` standing before a printed slot
            _bump!(r)
            body = _read_group!(r)
            isempty(strip(body)) && continue
            r.i = save
            break
        elseif c == '^' || c == '_'
            _bump!(r)
            pos = c == '^' ? upper : lower
            (sl, gr, _) = _parse_script(_read_argument!(r), length(slots))
            for (nm, wild) in sl
                push!(slots, _LTSlot(nm, wild, pos))
            end
            append!(groups, gr)
        elseif c == '('
            _bump!(r)
            for piece in split(_read_group_paren!(r), ',')
                nm = _tex_name(piece)
                isempty(nm) || push!(coords, nm)
            end
        else
            r.i = save
            break
        end
    end

    return _LTFactor(name, slots, coords, groups)
end

"""Read to the matching `)`, the opening parenthesis having been consumed."""
function _read_group_paren!(r::_TeXReader)
    depth = 1
    start = r.i
    while !_eof(r)
        c = _peek(r)
        if c == '('
            depth += 1
        elseif c == ')'
            depth -= 1
            depth == 0 && break
        end
        _bump!(r)
    end
    txt = r.s[start:prevind(r.s, r.i)]
    _eof(r) || _bump!(r)
    return txt
end

"""A coefficient, if one stands here"""
function _parse_coefficient!(r::_TeXReader)
    _skipspace!(r)
    _eof(r) && return nothing
    if _peek(r) == '\\'
        save = r.i
        _bump!(r)
        cmd = _read_command!(r)
        if cmd == "frac"
            num = _read_argument!(r)
            den = _read_argument!(r)
            n = tryparse(Int, strip(num))
            d = tryparse(Int, strip(den))
            (n === nothing || d === nothing) &&
                error("LT: only a numeric fraction is read as a coefficient, got " *
                      "a numerator `$num` and denominator `$den`.")
            return n // d
        end
        r.i = save
        return nothing
    end
    if isdigit(_peek(r))
        start = r.i
        while !_eof(r) && (isdigit(_peek(r)) || _peek(r) == '.')
            _bump!(r)
        end
        txt = r.s[start:prevind(r.s, r.i)]
        _skipspace!(r)
        if !_eof(r) && (_peek(r) == '^' || _peek(r) == '_')
            r.i = start          # a digit carrying a script is a tensor name, not a coefficient
            return nothing
        end
        v = tryparse(Int, txt)
        return v === nothing ? parse(Float64, txt) : v
    end
    return nothing
end

"""
    _lt_plan(expr) -> Vector{Tuple{Any,Vector{_LTFactor}}}

The whole string as a list of terms, each a coefficient together with the factors multiplied.
"""
function _lt_plan(expr::String)
    r = _TeXReader(expr)
    out = Tuple{Any,Vector{_LTFactor}}[]
    sgn = 1
    coeff = nothing
    factors = _LTFactor[]

    function flush!()
        isempty(factors) || push!(out, (sgn * (coeff === nothing ? 1 : coeff), copy(factors)))
        empty!(factors)
        coeff = nothing
        sgn = 1
    end

    while true
        _skipspace!(r)
        _eof(r) && break
        c = _peek(r)
        if c == '+'
            _bump!(r); flush!(); continue
        elseif c == '-'
            _bump!(r); flush!(); sgn = -1; continue
        elseif c == '$' || c == '&' || c == '(' || c == ')'
            # grouping parens, typically from \left( ... \right); a coordinate list is
            # consumed by the factor that owns it, so any paren reaching here is decoration
            _bump!(r); continue
        end
        if isempty(factors) && coeff === nothing
            v = _parse_coefficient!(r)
            if v !== nothing
                coeff = v
                continue
            end
        end
        f = _parse_factor!(r)
        if f === nothing
            _eof(r) && break        # trailing decoration, nothing left to read
            error("LT: cannot read `" * r.s[r.i:end] * "` in \"" * expr * "\".")
        end
        push!(factors, f)
    end
    flush!()
    return out
end

#------------------------------------------------------------------ resolution

"""
    _lt_slot(scoped, name, wild, pos) -> Pair

One parsed slot as the pair a tensor constructor expects. A wildcard is always a
`DummyPatternIndex`. Otherwise the caller's own binding wins, then any index of that name already
registered in a set, and failing both the slot is a label
"""
function _lt_slot(scoped, name::String, wild::Bool, pos::IndexPosition)
    wild && return Pair{IndexTypes,IndexPosition}(DummyPatternIndex(name), pos)
    scoped isa IndexSuperType && return Pair{IndexTypes,IndexPosition}(scoped, pos)
    idx = _find_registered_index(name)
    idx === nothing && return Pair{String,IndexPosition}(name, pos)
    return Pair{IndexTypes,IndexPosition}(idx, pos)
end

"""
The index of this name in any registered set, or `nothing`. This is what lets a printed expression
be read back without the indices being in scope under the same variable names.
"""
function _find_registered_index(name::AbstractString)
    for setname in sort!(collect(keys(SET_REGISTRY)))
        for idx in SET_REGISTRY[setname].indices
            idx.name == name && return idx
        end
    end
    return nothing
end

"""A coordinate named in a function-of list, from the caller's scope or the registry."""
function _lt_coord(scoped, name::String)
    scoped isa Coordinate && return scoped
    haskey(COORD_REGISTRY, name) && return COORD_REGISTRY[name]
    error("LT: no Coordinate named `" * name * "`. A parenthesised list after a tensor gives " *
          "the coordinates it depends on, so each entry must already be a declared Coordinate.")
end

"""A symmetry group as the tableau standing for it: a row when symmetric, a column when not."""
_lt_tableau(sym::Bool, acc::Vector{Int}) =
    sym ? YoungTableau([length(acc)], acc) : YoungTableau(fill(1, length(acc)), acc)

_lt_tableaux(f::_LTFactor) = [_lt_tableau(sym, acc) for (sym, acc) in f.groups]

_lt_build(f::_LTFactor) = Tensor(f.name,
    [_lt_slot(nothing, s.name, s.wild, s.pos) for s in f.slots]...,
    _lt_tableaux(f),
    [_lt_coord(nothing, c) for c in f.coords]...)

"""Combine the parsed terms into one object: a bare tensor when the markup holds a single factor
with unit coefficient, so that the simple case behaves as `T"..."` does, and a sum otherwise."""
function _lt_assemble(plan, build)
    isempty(plan) && return zero()
    if length(plan) == 1
        (c, fs) = plan[1]
        c == 1 && length(fs) == 1 && return build(fs[1])
        p = foldl(*, [build(f) for f in fs])
        return c == 1 ? TensorTerm(p) : c * TensorExpression(p)
    end
    total = nothing
    for (c, fs) in plan
        p = foldl(*, [build(f) for f in fs])
        piece = c == 1 ? TensorExpression(p) : c * TensorExpression(p)
        total = total === nothing ? piece : total + piece
    end
    return total
end

"""The counterpart of [`_lt_assemble`](@ref) for the macro, combining expressions not values."""
function _lt_assemble_expr(plan, build)
    isempty(plan) && return :($(zero)())
    if length(plan) == 1
        (c, fs) = plan[1]
        c == 1 && length(fs) == 1 && return build(fs[1])
        p = foldl((x, y) -> :($x * $y), [build(f) for f in fs])
        return c == 1 ? :($(TensorTerm)($p)) : :($c * $(TensorExpression)($p))
    end
    total = nothing
    for (c, fs) in plan
        p = foldl((x, y) -> :($x * $y), [build(f) for f in fs])
        piece = c == 1 ? :($(TensorExpression)($p)) : :($c * $(TensorExpression)($p))
        total = total === nothing ? piece : :($total + $piece)
    end
    return total
end

"""
    parse_latex(s) -> tensor

Read TeX markup held in a `String` at run time, as [`@LT_str`](@ref) does at parse time. Names
resolve against the registered index sets and coordinates only, there being no caller scope to
consult, which is what is wanted when reading back something [`to_latex`](@ref) wrote.
"""
parse_latex(s::AbstractString) = _lt_assemble(_lt_plan(String(s)), _lt_build)

"""
    LT"..."

Read genuine TeX markup as a tensor, inverting [`to_latex`](@ref).

What we print round-trips, so an expression can be copied out of a document and pasted back:

    LT"A { }_{\\mu}{ }^{\\nu}"

and so does what one would ordinarily type, where the `{ }` are absent and one script may carry
several indices:

    LT"F_{\\mu\\nu}"        LT"A^\\mu B_\\mu"        LT"\\bar\\psi_{\\dot\\alpha}"

Accents become the combining mark Alakazam stores, so that `\\dot{\\alpha}`, `\\dot\\alpha` and the
bare decorated character all name one index. A slot is an index if the caller has a variable of
that name, or if an index of that name is registered in some set; failing both it is a label, which
is the rule `T"..."` follows. `?` marks a wildcard, a parenthesised list gives the coordinates a
tensor depends on, and terms may be added or subtracted with an integer or fractional coefficient.

Symmetry is written as it is in the literature. Round brackets make the enclosed slots symmetric
and square brackets antisymmetric, and a pair of bars holds slots out of the group:

    LT"A_{(ab)}"      LT"F_{[\\mu\\nu]}"      LT"R_{a[bc|d|e]}"

These build the same tableaux the tensor strings do, and the braces around them are optional. Note
that `to_latex` does not write the brackets, so symmetry is the one thing that does not survive a
round trip through the printed form; it has to be read from markup that states it.
Antisymmetrisation spread across two factors, as in `\\partial_{[a}A_{b]}`, is not a symmetry of
either tensor and is refused rather than misread.

Only a numeric coefficient is read as one. A named factor is always a tensor, because in the
printed form a symbolic coefficient and a tensor carrying no indices are the same markup: both
`c * A_mu` and `Tensor("c") * A_mu` print as `c A { }_{\\mu}`, differing by a space. So

    parse_latex(to_latex(c * A)) == Tensor("c") * A

Use [`parse_latex`](@ref) for markup that is only known at run time.
"""
macro LT_str(expr::String)
    plan = _lt_plan(expr)
    build = function (f)
        slots = map(f.slots) do s
            sym = s.wild ? nothing : _identifier_spelling(s.name)
            scoped = sym === nothing ? :(nothing) :
                     :($(esc(:(@isdefined $sym))) ? $(esc(sym)) : nothing)
            :($(_lt_slot)($scoped, $(s.name), $(s.wild), $(s.pos)))
        end
        coords = map(f.coords) do c
            sym = _identifier_spelling(c)
            scoped = sym === nothing ? :(nothing) :
                     :($(esc(:(@isdefined $sym))) ? $(esc(sym)) : nothing)
            :($(_lt_coord)($scoped, $c))
        end
        tabs = _lt_tableaux(f)
        :($(Tensor)($(f.name), $(slots...), $tabs, $(coords...)))
    end
    return _lt_assemble_expr(plan, build)
end
