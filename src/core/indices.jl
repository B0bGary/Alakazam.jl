"""position is 0 for upper index, 1 for lower index"""
struct IndexPosition
   position::Bool
end

struct IndexGrading
   grading::Int64
end
const bosonic=IndexGrading(0)
const fermionic=IndexGrading(1)


mutable struct IndexSet
   name::String
   range::Int
   signature::Tuple{Int64,Int64}
   metric_type::Union{Type{<:MetricSuperType},Nothing}
   indices::Vector{IndexSuperType}
   default_position::IndexPosition
   grading::Int
   conjugate::Union{IndexSet,Nothing}
   algebra::Union{LieAlgebraSuperType,Nothing}
   islightcone::Bool
   function IndexSet(name::String, args...)
      range = 4
      metric_type = Metric
      signature = nothing
      default_position = upper
      grading=0
      conjugate=nothing
      algebra=nothing
      islightcone=(name=="lightcone" ? true : false)

      for arg in args
         if arg isa Int
            range = arg
         elseif arg isa Type{<:MetricSuperType}
            metric_type = arg
         elseif arg isa Tuple{Int64,Int64}
            signature = arg
         elseif arg isa Vector{Int64}
            if length(arg)!=2
               println("Signature "*string(arg)*" is invalid. Signatures should have two numbers.")
            else
               signature = (arg[1], arg[2])
            end
         elseif arg isa IndexPosition
            default_position = arg
         elseif arg isa IndexGrading
            grading=arg.grading
         elseif arg isa Nothing
            nothing
         elseif arg isa LieAlgebraSuperType
            metric_type=metric(arg)
            algebra=arg
            range=dimension(arg)
         elseif arg isa Bool
            islightcone=arg
         else
            error("Unexpected argument: $arg")
         end
      end
      metric_type == MetricSuperType && error("Explicit metric subtype of MetricSuperType must be specified")
      if signature isa Nothing
         if metric_type!==nothing && metric_type <: SymMetricSuperType
            signature = (1, range - 1) #default to lorentzian
         else
            signature = (0, range) #otherwise, euclidean
         end
      else
         range = signature[1] + signature[2]
      end
      res=new(name, range, signature, metric_type, Vector{IndexSuperType}(), default_position, grading, conjugate,algebra,islightcone)
      global SET_REGISTRY
      SET_REGISTRY[name]=res
      return res
   end

end

"""Names that are forbidden, reserved for internal pattern matching logic."""
const forbidden=Set(["*","?","_","{","}","[","]","^","(",")","∘"])

valid_name(name::String)= (name in forbidden) ? error("Index name \"$name\" is reserved for internal logic.") : name


mutable struct Index <: IndexSuperType
   name::String
   index_set::IndexSet
   uid::UInt64
   function Index(name::String, index_set::IndexSet)
      for idx in index_set.indices
         if idx.name == name #stop double creation
            return idx
         end
      end
      temp = new(valid_name(name), index_set, hash(name, hash(index_set.name)))
      push!(index_set.indices, temp)
      return temp
   end
   function Index(name::String)
      s = IndexSet("")
      return new(valid_name(name), s, hash(name, hash("")))
   end

   function Index(args::Union{String,IndexSet,}...)
      if args[end] isa IndexSet
         index_set = args[end]
         return Tuple(Index(name, index_set) for name in args[1:(end-1)])
      else
         return Tuple(Index(name) for name in args)
      end
   end

   function Index(names::Vector{String}, index_set::IndexSet)
      return [Index(name, index_set) for name in names]
   end
      function Index(names::Vector{Char}, index_set::IndexSet)
      return [Index(name, index_set) for name in names]
   end
   function Index(name::Char, index_set::IndexSet)
      return Index(string(name), index_set)
   end
   function Index(names::String...; index_set::IndexSet)
      if length(names) == 1
         return Index(first(names), index_set)  # fall through to original
      end
      last_arg = names[end]
      # check if last arg was meant to be used differently
      return [Index(name, index_set) for name in names]
   end
end



mutable struct DummyPatternIndex <: IndexSuperType
   name::String
   index_set::IndexSet
   uid::UInt64
   function DummyPatternIndex(name::String, index_set::IndexSet)
      for idx in index_set.indices
         if idx.name == name
            return idx
         end
      end
      temp = new(name, index_set, hash(name, hash(index_set.name)))
      push!(index_set.indices, temp)
      return temp
   end

   function DummyPatternIndex(name::String)
      return new(name, IndexSet(""), hash(name, hash("")))
   end
   function DummyPatternIndex(names::Vector{String}, index_set::IndexSet)
      return [DummyPatternIndex(name, index_set) for name in names]
   end
end

mutable struct Coordinate
   name::String
   index_set::Vector{IndexSet}
   function Coordinate(name::String, index_set2::IndexSet)
      c=new(name, [index_set2])
      COORD_REGISTRY[name]=c
      return c
   end
   function Coordinate(name::String, index_set::Vector{IndexSet})
      c=new(name, index_set)
      COORD_REGISTRY[name]=c
      return c
   end
end

#consts
const upper = IndexPosition(0)
const lower = IndexPosition(1)
const ↑ = upper
const ↓ = lower
const SET_REGISTRY=Dict{String,IndexSet}()
const COORD_REGISTRY=Dict{String,Coordinate}()
const IndexTypes = Union{Index,DummyPatternIndex}
const IndexPair = Pair{IndexTypes,IndexPosition}

#funcs

clear_coords!() = empty!(COORD_REGISTRY)
clear_indexsets!() = empty!(SET_REGISTRY)


grading(i::I) where I<:IndexSuperType = i.index_set.grading


function get_index(iname::String, set::IndexSet)
   for i in set.indices
      i.name == iname && return i
   end
end
function get_index(iname::String, setname::String)
   haskey(SET_REGISTRY, setname) || error("Unknown index set: $setname")

   for i in SET_REGISTRY[setname].indices
      i.name == iname && return i
   end

   error("No index named $iname in set $setname")
end

"""Returns the det of a metric associated with this index type."""
function det(I::IndexSuperType)
   (-1)^(I.index_set.signature[1])
end

function name(I::IndexSuperType)
   return I.name
end

function in_default_pos(i::Pair{I,IndexPosition}) where {I<:IndexSuperType}
   return i[2] == i[1].index_set.default_position

end

function op_pos(ip::IndexPosition)
   if ip == upper
      return lower
   else
      return upper
   end
end

function op_pos(i::Pair{I,IndexPosition}) where {I<:IndexSuperType}
   return i[1]=>op_pos(i[2])
end
function op_pos(v::Vector{Pair{I,IndexPosition}}) where {I<:IndexSuperType}
   return [i[1]=>op_pos(i[2]) for i in v]
end

function op_pos(v::T) where {T<:AbstractSet{<:Pair{I,IndexPosition}}} where {I<:IndexSuperType}
   return T(i[1] => op_pos(i[2]) for i in v)
end


function dimension(I::IndexSuperType)
   return I.index_set.range
end
function dimension(IP::Pair{I,IndexPosition}) where {I<:IndexSuperType}
   return IP[1].index_set.range
end

function ==(a::Coordinate, b::Coordinate)
   return a.name == b.name && a.index_set == b.index_set # index set ....  would cause issues since index set is a vecot. would have to use name but good for now
end

# Base.hash(a::Coordinate, h::UInt) = hash((a.name, a.index_set), h)


Base.hash(a::Coordinate, h::UInt) = hash(a.index_set, hash(a.name, h))

Base.hash(a::IndexSet, h::UInt) = hash(a.default_position, hash(a.range, hash(a.name, h)))
#uid compare should be faster
function ==(a::Index, b::Index)
    a.uid == b.uid && return true
    # wildcard: empty set name matches anything with the same name
    (isempty(a.index_set.name) || isempty(b.index_set.name)) && return a.name == b.name
    return false
end
# function ==(a::Index, b::Index)
#    return a.name == b.name && ((a.index_set.name=="" || b.index_set.name=="")||(a.index_set.name == b.index_set.name)) # index set ....  would cause issues since index set is a vecot. would have to use name but good for now
# end#can;t use index set name due to replacement pattern matching

# Base.hash(a::Index, h::UInt) = hash((a.name,a.index_set.name), h)
Base.hash(a::Index, h::UInt) = hash(a.uid, h)
Base.hash(a::DummyPatternIndex, h::UInt) = hash(a.uid, h)
function ==(a::DummyPatternIndex, b::DummyPatternIndex)
   return a.name == b.name #&& a.index_set.name  == b.index_set.name  # index set ....  would cause issues since index set is a vecot. would have to use name but good for now
end
# Base.hash(a::DummyPatternIndex, h::UInt) = hash((a.name,a.index_set.name), h)


function ==(a::IndexSet, b::IndexSet)
   return a.name == b.name && a.range == b.range && a.default_position == b.default_position
end
# Base.hash(a::IndexSet, h::UInt) = hash((a.name, a.range, a.default_position), h)


function ==(a::IndexPosition, b::IndexPosition)
   return a.position == b.position
end
Base.hash(a::IndexPosition, h::UInt) = hash((a.position), h)


isLightcone(I::IndexSuperType) = I.index_set.islightcone
function isLightcone(i::Pair{I,IndexPosition}) where {I<:IndexSuperType}
    return i[1].index_set.islightcone
end


function Base.show(io::IO, C::Coordinate)
   print(io, C.name)
   for ind in C.index_set
      print(io, "[")
      print(io, ind.indices[1].name)
      print(io, "]")
   end

end

function Base.show(io::IO, C::Vector{Coordinate})
   print(io, "{")
   it = Iterators.Stateful(C)
   for i in it
      print(io, i)
      if peek(it) != nothing
         print(io, ", ")
      end
   end
   print(io, "}")
end

function Base.show(io::IO, I::IndexSuperType)
   if typeof(I) == DummyPatternIndex
      print(io, "?")
   end
   print(io, I.name)
end

function Base.show(io::IO, I::IndexPosition)
   print(io, I.position == 0 ? "upper" : "lower")
end

function Base.show(io::IO, I::IndexSet)
   print(io, "[")
   it = Iterators.Stateful(I.indices)
   for i in it
      print(io, i.name)
      if peek(it) != nothing
         print(io, ", ")
      end
   end
   print(io, "]")
end

function Base.deepcopy_internal(x::IndexSet, dict::IdDict)
   return x
end

function Base.isless(a::IndexTypes, b::IndexTypes)
   a.name < b.name
end

function Base.isless(a::Pair{<:IndexSuperType,IndexPosition}, b::Pair{<:IndexSuperType,IndexPosition})
   if a[2] != b[2]
      return a[2] < b[2]
   else
      return a[1] < b[1]
   end
end

function Base.isless(a::IndexPosition, b::IndexPosition)
   a.position < b.position
end



# macro indices(ex)
#     @assert ex.head == :(=) "Usage: @indices a:z = expr"

#     lhs, rhs = ex.args

#     @assert lhs.head == :call && lhs.args[1] == :(:)"LHS must be a symbol range, e.g. a:z"

#     first = string(lhs.args[2])
#     last  = string(lhs.args[3])

#     @assert length(first) == 1 && length(last) == 1 "Range endpoints must be single-character identifiers"

#     vars = Symbol.(collect(first[1]:last[1]))

#     :( $(esc(Expr(:tuple, vars...))) = $(esc(rhs)) )
# end

# macro indices(ex)
#     @assert ex.head == :(=)

#     lhs, rhs = ex.args
#     @assert lhs.head == :call && lhs.args[1] == :(:)

#     first = only(String(lhs.args[2]))
#     last  = only(String(lhs.args[3]))

#     chars = collect(first:last)
#     vars  = map(Symbol, chars)
#     tuple = Expr(:tuple, vars...)

#     quote
#         $(esc(tuple)) = Index($chars, $(esc(rhs)))
#     end
# end
macro indices(ex)
    @assert ex.head == :call
    @assert ex.args[1] == :∈

    lhs, rhs = ex.args[2], ex.args[3]

    @assert lhs.head == :call
    @assert lhs.args[1] == :(:)

    first = only(String(lhs.args[2]))
    last  = only(String(lhs.args[3]))

    chars = collect(first:last)
    vars  = Symbol.(chars)
    tuple = Expr(:tuple, vars...)

    quote
        $(esc(tuple)) = Index($chars, $(esc(rhs)))
    end
end


"""
Gets an unused index from a set, or generate a new index if we've run out.
"""
function safe_pop!(unused::OrderedSet{I}, IS::IndexSet) where I<:IndexSuperType
   if isempty(unused)
      nid = generate_indices(IS)
      nid === nothing && throw(OutOfIndices(IS))
      push!(unused, nid)
   end
   return popfirst!(unused) #pop first or pop??
end
function safe_pop!(unused::Vector{I}, IS::IndexSet) where I<:IndexSuperType
    if isempty(unused)
        nid = generate_indices(IS)
        nid === nothing && throw(OutOfIndices(IS))
        push!(unused, nid)
    end
    return popfirst!(unused)   # popfirst! works on Vector
end
"""
Gets a a list of indices from the indexset of length n, generating more if needed.
"""
function safe_pop!(is::IndexSet,n::Int)
   indpool=copy(is.indices)
   return [safe_pop!(indpool,is)=>lower for i in 1:n]
end

const COMBINING_DOT = '\u0307'
const COMBINING_MACRON = '\u0304'
const COMBINING_HAT = '\u0302'
const COMBINING_TILDE = '\u0303'

dotted(c::Char) = string(COMBINING_DOT) * string(c)
barred(c::Char) = string(COMBINING_MACRON) * string(c)
hatted(c::Char) = string(COMBINING_HAT) * string(c)
tilded(c::Char) = string(COMBINING_TILDE) * string(c)


const LATIN_LOWER = [string(c) for c in 'a':'z']
const LATIN_UPPER = [string(c) for c in 'A':'Z']
const GREEK_LOWER = [string(Char(i)) for i in 0x3B1:0x3C9]
const GREEK_UPPER = [string(Char(i)) for i in 0x391:0x3A9]
const HEBREW = [string(Char(i)) for i in 0x05D0:0x05EA]

const LATIN_DOTTED_LOWER = [dotted(c) for c in 'a':'z']
const LATIN_DOTTED_UPPER = [dotted(c) for c in 'A':'Z']
const LATIN_BARRED_LOWER = [barred(c) for c in 'a':'z']
const LATIN_BARRED_UPPER = [barred(c) for c in 'A':'Z']
const LATIN_HATTED_LOWER = [hatted(c) for c in 'a':'z']
const LATIN_HATTED_UPPER = [hatted(c) for c in 'A':'Z']
const LATIN_TILDE_LOWER = [tilded(c) for c in 'a':'z']
const LATIN_TILDE_UPPER = [tilded(c) for c in 'A':'Z']

const GREEK_DOTTED_LOWER = [dotted(Char(i)) for i in 0x3B1:0x3C9]
const GREEK_DOTTED_UPPER = [dotted(Char(i)) for i in 0x391:0x3A9]
const GREEK_BARRED_LOWER = [barred(Char(i)) for i in 0x3B1:0x3C9]
const GREEK_BARRED_UPPER = [barred(Char(i)) for i in 0x391:0x3A9]
const GREEK_HATTED_LOWER = [hatted(Char(i)) for i in 0x3B1:0x3C9]
const GREEK_HATTED_UPPER = [hatted(Char(i)) for i in 0x391:0x3A9]
const GREEK_TILDE_LOWER = [tilded(Char(i)) for i in 0x3B1:0x3C9]
const GREEK_TILDE_UPPER = [tilded(Char(i)) for i in 0x391:0x3A9]

const HEBREW_DOTTED = [dotted(Char(i)) for i in 0x05D0:0x05EA]
const HEBREW_BARRED = [barred(Char(i)) for i in 0x05D0:0x05EA]
const HEBREW_HATTED = [hatted(Char(i)) for i in 0x05D0:0x05EA]
const HEBREW_TILDE = [tilded(Char(i)) for i in 0x05D0:0x05EA]

const FRAKTUR_LOWER = [string(Char(i)) for i in 0x1D51E:0x1D537]
# Fraktur uppercase exceptions: C=\u212D, H=\u210C, I=\u2111, R=\u211C, Z=\u2128
const FRAKTUR_UPPER = [
   "\U1D504", "\U1D505", "\u212D", "\U1D507", "\U1D508", "\U1D509",
   "\U1D50A", "\u210C", "\u2111", "\U1D50D", "\U1D50E", "\U1D50F",
   "\U1D510", "\U1D511", "\U1D512", "\U1D513", "\U1D514", "\u211C",
   "\U1D516", "\U1D517", "\U1D518", "\U1D519", "\U1D51A", "\U1D51B",
   "\U1D51C", "\u2128"
]
# Script uppercase has exceptions: B=\u212C, E=\u2130, F=\u2131, H=\u210B, I=\u2110, L=\u2112, M=\u2133, R=\u211B
const SCRIPT_UPPER = [
   "\U1D49C", "\u212C", "\U1D49E", "\U1D49F", "\u2130", "\u2131",
   "\U1D4A2", "\u210B", "\u2110", "\U1D4A5", "\U1D4A6", "\u2112",
   "\u2133", "\U1D4A9", "\U1D4AA", "\U1D4AB", "\U1D4AC", "\u211B",
   "\U1D4AE", "\U1D4AF", "\U1D4B0", "\U1D4B1", "\U1D4B2", "\U1D4B3",
   "\U1D4B4", "\U1D4B5"
]

# Script lowercase has exceptions: e=\u212F, g=\u210A, o=\u2134
const SCRIPT_LOWER = [
   "\U1D4B6", "\U1D4B7", "\U1D4B8", "\U1D4B9", "\u212F", "\U1D4BB",
   "\u210A", "\U1D4BD", "\U1D4BE", "\U1D4BF", "\U1D4C0", "\U1D4C1",
   "\U1D4C2", "\U1D4C3", "\u2134", "\U1D4C5", "\U1D4C6", "\U1D4C7",
   "\U1D4C8", "\U1D4C9", "\U1D4CA", "\U1D4CB", "\U1D4CC", "\U1D4CD",
   "\U1D4CE", "\U1D4CF"
]

# Double-struck uppercase has exceptions: C=\u2102, H=\u210D, N=\u2115, P=\u2119, Q=\u211A, R=\u211D, Z=\u2124
const DOUBLESTRUCK_UPPER = [
   "\U1D538", "\U1D539", "\u2102", "\U1D53B", "\U1D53C", "\U1D53D",
   "\U1D53E", "\u210D", "\U1D540", "\U1D541", "\U1D542", "\U1D543",
   "\U1D544", "\u2115", "\U1D546", "\u2119", "\u211A", "\u211D",
   "\U1D54A", "\U1D54B", "\U1D54C", "\U1D54D", "\U1D54E", "\U1D54F",
   "\U1D550", "\u2124"
]

# Double-struck lowercase is contiguous
const DOUBLESTRUCK_LOWER = [string(Char(i)) for i in 0x1D552:0x1D56B]

const FRAKTUR_DOTTED_LOWER = [dotted(Char(i)) for i in 0x1D51E:0x1D537]
const FRAKTUR_DOTTED_UPPER = [dotted(Char(i)) for i in 0x1D504:0x1D51D]
const FRAKTUR_BARRED_LOWER = [barred(Char(i)) for i in 0x1D51E:0x1D537]
const FRAKTUR_BARRED_UPPER = [barred(Char(i)) for i in 0x1D504:0x1D51D]
const FRAKTUR_HATTED_LOWER = [hatted(Char(i)) for i in 0x1D51E:0x1D537]
const FRAKTUR_HATTED_UPPER = [hatted(Char(i)) for i in 0x1D504:0x1D51D]
const FRAKTUR_TILDE_LOWER = [tilded(Char(i)) for i in 0x1D51E:0x1D537]
const FRAKTUR_TILDE_UPPER = [tilded(Char(i)) for i in 0x1D504:0x1D51D]

const SCRIPT_DOTTED_LOWER = [dotted(Char(i)) for i in 0x1D4B6:0x1D4CF]
const SCRIPT_DOTTED_UPPER = [dotted(Char(i)) for i in 0x1D49C:0x1D4B5]
const SCRIPT_BARRED_LOWER = [barred(Char(i)) for i in 0x1D4B6:0x1D4CF]
const SCRIPT_BARRED_UPPER = [barred(Char(i)) for i in 0x1D49C:0x1D4B5]
const SCRIPT_HATTED_LOWER = [hatted(Char(i)) for i in 0x1D4B6:0x1D4CF]
const SCRIPT_HATTED_UPPER = [hatted(Char(i)) for i in 0x1D49C:0x1D4B5]
const SCRIPT_TILDE_LOWER = [tilded(Char(i)) for i in 0x1D4B6:0x1D4CF]
const SCRIPT_TILDE_UPPER = [tilded(Char(i)) for i in 0x1D49C:0x1D4B5]

const DOUBLESTRUCK_DOTTED_LOWER = [dotted(Char(i)) for i in 0x1D552:0x1D56B]
const DOUBLESTRUCK_DOTTED_UPPER = [dotted(Char(i)) for i in 0x1D538:0x1D551]
const DOUBLESTRUCK_BARRED_LOWER = [barred(Char(i)) for i in 0x1D552:0x1D56B]
const DOUBLESTRUCK_BARRED_UPPER = [barred(Char(i)) for i in 0x1D538:0x1D551]
const DOUBLESTRUCK_HATTED_LOWER = [hatted(Char(i)) for i in 0x1D552:0x1D56B]
const DOUBLESTRUCK_HATTED_UPPER = [hatted(Char(i)) for i in 0x1D538:0x1D551]
const DOUBLESTRUCK_TILDE_LOWER = [tilded(Char(i)) for i in 0x1D552:0x1D56B]
const DOUBLESTRUCK_TILDE_UPPER = [tilded(Char(i)) for i in 0x1D538:0x1D551]

const ORDERED_POOLS = [
   LATIN_LOWER, LATIN_UPPER,
   GREEK_LOWER, GREEK_UPPER,
   HEBREW,
   FRAKTUR_LOWER, FRAKTUR_UPPER,
   SCRIPT_LOWER, SCRIPT_UPPER,
   DOUBLESTRUCK_LOWER, DOUBLESTRUCK_UPPER,
   LATIN_DOTTED_LOWER, LATIN_DOTTED_UPPER,
   LATIN_BARRED_LOWER, LATIN_BARRED_UPPER,
   LATIN_HATTED_LOWER, LATIN_HATTED_UPPER,
   LATIN_TILDE_LOWER, LATIN_TILDE_UPPER,
   GREEK_DOTTED_LOWER, GREEK_DOTTED_UPPER,
   GREEK_BARRED_LOWER, GREEK_BARRED_UPPER,
   GREEK_HATTED_LOWER, GREEK_HATTED_UPPER,
   GREEK_TILDE_LOWER, GREEK_TILDE_UPPER,
   HEBREW_DOTTED, HEBREW_BARRED,
   HEBREW_HATTED, HEBREW_TILDE,
   FRAKTUR_DOTTED_LOWER, FRAKTUR_DOTTED_UPPER,
   FRAKTUR_BARRED_LOWER, FRAKTUR_BARRED_UPPER,
   FRAKTUR_HATTED_LOWER, FRAKTUR_HATTED_UPPER,
   FRAKTUR_TILDE_LOWER, FRAKTUR_TILDE_UPPER,
   SCRIPT_DOTTED_LOWER, SCRIPT_DOTTED_UPPER,
   SCRIPT_BARRED_LOWER, SCRIPT_BARRED_UPPER,
   SCRIPT_HATTED_LOWER, SCRIPT_HATTED_UPPER,
   SCRIPT_TILDE_LOWER, SCRIPT_TILDE_UPPER,
   DOUBLESTRUCK_DOTTED_LOWER, DOUBLESTRUCK_DOTTED_UPPER,
   DOUBLESTRUCK_BARRED_LOWER, DOUBLESTRUCK_BARRED_UPPER,
   DOUBLESTRUCK_HATTED_LOWER, DOUBLESTRUCK_HATTED_UPPER,
   DOUBLESTRUCK_TILDE_LOWER, DOUBLESTRUCK_TILDE_UPPER,
]

"""First unused name in `pool` as a new index, or nothing if the pool is full."""
function _first_unused_in(pool, used::Set{String}, template, I::IndexSet)
   for name in pool
      if name ∉ used
         #dont need push because index constructor creates
         return typeof(template)(name, I)
      end
   end
   nothing
end

"""
Gets a new index for a set, or nothing when every pool is used up.
"""
function generate_indices(I::IndexSet)
   isempty(I.indices) && error("IndexSet has no indices to infer character set from")
   # @show objectid(I)
   template = I.indices[1]#get the type of character
   used = Set(idx.name for idx in I.indices)

   #which pool the template is from, if any. a freely chosen name like "i1" is in
   #none of them and the Index constructor allows it, so this can come back empty
   template_index = findfirst(pool -> template.name in pool, ORDERED_POOLS)

   #own pool first, so the names are unchanged when one alphabet is enough
   if template_index !== nothing
      idx = _first_unused_in(ORDERED_POOLS[template_index], used, template, I)
      idx === nothing || return idx
   end

   #then the rest, else a set caps at one alphabet with the other pools untouched
   for (p, pool) in enumerate(ORDERED_POOLS)
      p == template_index && continue
      idx = _first_unused_in(pool, used, template, I)
      idx === nothing || return idx
   end

   #every caller does `nid === nothing && throw(OutOfIndices(...))`
   nothing
end
