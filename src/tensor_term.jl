# include("indices.jl")
# include("custom_exceptions.jl")

SymorNum = Union{SymbolicUtils.Symbolic,Number}
abstract type TensorSuperType end


#pattern matching struct... don't touch
mutable struct DummyPatternTensor <: TensorSuperType
    name::String
    indices::Vector{Pair{I,IndexPosition}} where I<:IndexSuperType
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}

    function DummyPatternTensor(name::String, args...)
        indices = Pair{IndexSuperType,IndexPosition}[]  # Initialize indices as an empty array
        function_of = Coordinate[]  # Initialize function_of as an empty array
        weights = Dict{String,Number}()  # Initialize an empty weights dictionary

        for arg in args
            if arg isa Pair{I,IndexPosition} where I<:IndexSuperType
                push!(indices, arg)
            elseif arg isa Vector{Pair{I,IndexPosition}} where I<:IndexSuperType
                append!(indices, arg)
            elseif arg isa Coordinate
                push!(function_of, arg)
            elseif arg isa Vector{Coordinate}
                append!(function_of, arg)
            elseif arg isa Dict{String,Number}
                weights = arg
            end
        end

        # Ensure function_of is not overwritten if it's not passed
        if isempty(function_of)
            function_of = [Coordinate()]  # Use a default if necessary
        end

        # Ensure weights is not empty and infer if necessary
        if isempty(weights) && !isempty(function_of)
            for coord in function_of
                weights[coord.name] = 1  # Default weight
            end
        end

        new(name, indices, function_of, weights)
    end
end


#####################
####### tensor data struct. more universal.
#####################
"""
Mutual tensor data that all tensors should have in common.
"""
mutable struct TensorData
    name::String
    indices::Vector{Pair{M,IndexPosition}} where {M<:IndexSuperType}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}
    tableaux::Vector{YoungTableau}

    function TensorData(name::String, args...)
        indices = Vector{Pair{IndexSuperType,IndexPosition}}()
        function_of = Vector{Coordinate}()
        weights = Dict{String,Number}()
        tableaux = Vector{YoungTableau}()
        for arg in args
            if arg isa Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
                indices = vcat(indices, arg)
            elseif arg isa Pair{I,IndexPosition} where {I<:IndexSuperType}
                push!(indices, arg)
            elseif arg isa Coordinate
                function_of = [arg]
            elseif arg isa Vector{Coordinate}
                function_of = arg
            elseif arg isa Dict{String,Number}
                weights = arg
            elseif arg isa YoungTableau
                tableaux = [arg]
            elseif arg isa Vector{YoungTableau}
                tableaux = arg
            end
        end
        w = Dict{String,Number}()
        for coord in function_of
            w[coord.name] = 1
        end
        return new(name, indices, function_of,
            merge(w, weights), tableaux)
    end
end


#####################
####### primitive tensor types
#####################

mutable struct Tensor <: TensorSuperType
    name::String
    indices::Vector{Pair{M,IndexPosition}} where {M<:IndexSuperType}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}
    tableaux::Vector{YoungTableau}

    function Tensor(name::String, args...)
        indices = Vector{Pair{IndexSuperType,IndexPosition}}()
        function_of = Vector{Coordinate}()
        weights = Dict{String,Number}()
        tableaux = Vector{YoungTableau}()

        for arg in args
            if arg isa Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
                indices = vcat(indices, arg)
            elseif arg isa Pair{I,IndexPosition} where {I<:IndexSuperType}
                push!(indices, arg)
            elseif arg isa Coordinate
                function_of = [arg]
            elseif arg isa Vector{Coordinate}
                function_of = arg
            elseif arg isa Dict{String,Number}
                weights = arg
            elseif arg isa YoungTableau
                tableaux = [arg]
            elseif arg isa Vector{YoungTableau}
                tableaux = arg
            end
        end

        w = Dict{String,Number}()

        for coord in function_of
            w[coord.name] = 1
        end

        return new(name, indices, function_of,
            merge(w, weights), tableaux)
    end
end

function SymmetricTensor(name::String, args...)
    t = Tensor(name, args...)
    t.tableaux = generate_symmetric_tableau(t.indices)
    return t
end

function AntisymmetricTensor(name::String, args...)
    t = Tensor(name, args...)
    t.tableaux = generate_antisymmetric_tableau(t.indices)
    return t
end

function RiemannTensor(args...)
    return RiemannTensor("R", args...)
end

function RiemannTensor(name::String, args...)
    t = Tensor(name, args...)
    if length(t.indices) != 4
        throw(ArgumentError("Riemann tensor must have exactly 4 indices"))
    end
    t.tableaux = generate_riemann_tableau()
    return t
end

abstract type MatrixSuperType <: TensorSuperType end #supertype which may not commute
abstract type PauliSuperType <: MatrixSuperType end #new supertype for pauli matrices, since they have special properties that we want to be able to check for in simplification rules
abstract type GammaSuperType <: MatrixSuperType end #new supertype for pauli matrices, since they have special properties that we want to be able to check for in simplification rules


# mutable struct Pauli <: PauliSuperType
#     name::String
#     indices::Vector{Pair{M,IndexPosition}} where {M<:IndexSuperType}
#     function_of::Vector{Coordinate}
#     weights::Dict{String,Number}
#     tableaux::Vector{YoungTableau}


#     t = Tensor(name, args...)
#     if length(t.indices) != 3 || length(t.indices) != 1
#         throw(ArgumentError("Pauli matrices must have exactly 1 or 3 indices"))
#     end
#     t.tableaux = generate_symmetric_tableau(t.indices)
#     return t

#     function Tensor(name::String, args...)
#         indices = Vector{Pair{IndexSuperType,IndexPosition}}()
#         function_of = Vector{Coordinate}()
#         weights = Dict{String,Number}()
#         tableaux = Vector{YoungTableau}()

#         for arg in args
#             if arg isa Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
#                 indices = vcat(indices, arg)
#             elseif arg isa Pair{I,IndexPosition} where {I<:IndexSuperType}
#                 push!(indices, arg)
#             elseif arg isa Coordinate
#                 function_of = [arg]
#             elseif arg isa Vector{Coordinate}
#                 function_of = arg
#             elseif arg isa Dict{String,Number}
#                 weights = arg
#             elseif arg isa YoungTableau
#                 tableaux = [arg]
#             elseif arg isa Vector{YoungTableau}
#                 tableaux = arg
#             end
#         end

#         w = Dict{String,Number}()

#         for coord in function_of
#             w[coord.name] = 1
#         end

#         return new(name, indices, function_of,
#             merge(w, weights), tableaux)
#     end
# end


# """
# Pauli Matrices. 3 index version are numbers and can commute. The 1 index version are matrices and do not commute.
# """
# function Pauli(name::String="σ", args...)
#     <:PauliSuperType
#     t = Tensor(name, args...)
#     if length(t.indices) != 3 || length(t.indices) != 1
#         throw(ArgumentError("Pauli matrices must have exactly 1 or 3 indices"))
#     end
#     t.tableaux = generate_symmetric_tableau(t.indices)
#     return t
# end


# """
# Gamma Matrices. 3 index version are numbers and can commute. The 1 index version are matrices and do not commute.
# """
# function Gamma(name::String="γ", args...)
#     <:GammaSuperType
#     t = Tensor(name, args...)
#     if length(t.indices) != 3 || length(t.indices) != 1
#         throw(ArgumentError("Pauli matrices must have exactly 1 or 3 indices"))
#     end
#     t.tableaux = generate_symmetric_tableau(t.indices)
#     return t
# end

# function Gamma5(name::String="γ5", args...)
#     <:GammaSuperType
#     t = Tensor(name, args...)
#     if length(t.indices) != 3 || length(t.indices) != 1
#         throw(ArgumentError("Pauli matrices must have exactly 1 or 3 indices"))
#     end
#     t.tableaux = generate_symmetric_tableau(t.indices)
#     return t
# end


# keep these, since we need to check if metrics in some checks
abstract type MetricSuperType <: TensorSuperType end
abstract type SymMetricSuperType <: MetricSuperType end
abstract type AntisymMetricSuperType <: MetricSuperType end
abstract type EpsilonTensorSuperType <: AntisymMetricSuperType end
abstract type KroneckerSuperType <: TensorSuperType end

"""
    Creates a metric tensor with (-,+) signature. By default, the signature is (1,d-1) where d is the range of the IndexSet of the attached indices, i.e. we use mostly-plus Minkowski convention unless otherwise specified.
"""
mutable struct Metric <: SymMetricSuperType
    name::String
    indices::Vector{Pair{IndexSuperType,IndexPosition}}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}
    tableaux::Vector{YoungTableau}

    signature::Pair{Int64} #added metric signature

    Metric(indices::Vector{Pair{I,IndexPosition}}) where {I<:IndexSuperType} = Metric("η", indices)


    function Metric(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(), signature=nothing) where {I<:IndexSuperType}
        if length(indices) != 2
            throw(ArgumentError("Metric must have exactly 2 indices"))
        end
        if signature === nothing
            signature = (1 => indices[1][1].index_set.range - 1)
        end
        new(name, indices, Vector{Coordinate}(), Dict{String,Number}(), generate_symmetric_tableau(indices), signature)
    end
end


mutable struct EpsilonTensor <: EpsilonTensorSuperType
    name::String
    indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
    weights::Dict{String,Number}
    function_of::Vector{Coordinate}
    tableaux::Vector{YoungTableau}


    EpsilonTensor(indices::Vector{Pair{I,IndexPosition}}) where {I<:IndexSuperType} = EpsilonTensor("ϵ", indices)
    function EpsilonTensor(name::String, indices::Vector{Pair{I,IndexPosition}}) where {I<:IndexSuperType}#lets you have scalars with no indices
        if indices[1][2] != indices[2][2]
            throw(InvalidEpsilonIndices()) #should be both upper or both lower
        end
        if indices[1][1].index_set != indices[2][1].index_set
            println(indices[1][1])
            println(indices[2][1])
            throw(InvalidEpsilonIndices()) #should be from the same index set
        end
        new(name, indices, Dict{String,Number}(), [], generate_antisymmetric_tableau(indices))
    end
end


mutable struct KroneckerDelta <: KroneckerSuperType
    name::String
    indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
    weights::Dict{String,Number}
    function_of::Vector{Coordinate}
    tableaux::Vector{YoungTableau}

    KroneckerDelta(indices::Vector{Pair{I,IndexPosition}}) where {I<:IndexSuperType} = KroneckerDelta("δ", indices)

    function KroneckerDelta(name::String, indices::Vector{Pair{I,IndexPosition}}) where {I<:IndexSuperType}#lets you have scalars with no indices
        if indices[1][2] == indices[2][2] #should be one upper and one lower
            throw(InvalidKroneckerIndices(indices))
        end
        if indices[1][1].index_set != indices[2][1].index_set #should be from the same index set
            throw(InvalidKroneckerIndices(indices))
        end
        new(name, indices, Dict{String,Number}(), [], generate_symmetric_tableau(indices))
    end
end



#######
# common tableaux generation functions
######

function generate_symmetric_tableau(indices::Vector)
    n = length(indices)
    return tableaux = [YoungTableau([n], collect(1:n))]
end
function generate_antisymmetric_tableau(indices::Vector)
    n = length(indices)
    return tableaux = [YoungTableau(fill(1, n), collect(1:n))]
end
function generate_symmetric_tableau(n::Integer)
    return tableaux = [YoungTableau([n], collect(1:n))]
end
function generate_antisymmetric_tableau(n::Integer)
    return tableaux = [YoungTableau(fill(1, n), collect(1:n))]
end
function generate_riemann_tableau()
    return tableaux = [YoungTableau([2, 2], [1, 3, 2, 4])]
end



#####
# access properties of tensors
#####

#will implement this soon... but will require work to change across everything...

# name(t::TensorSuperType) = t.data.name
# indices(t::TensorSuperType) = t.data.indices
# function_of(t::TensorSuperType) = t.data.function_of
# weights(t::TensorSuperType) = t.data.weights
# tableaux(t::TensorSuperType) = t.data.tableaux

######
# compound tensor terms and expressions
######


mutable struct TensorTerm
    terms::Vector{T} where {T<:TensorSuperType}
    indices::Vector{Pair{IndexSuperType,IndexPosition}}
    function TensorTerm(terms::Vector{T} where {T<:TensorSuperType})
        inds = []
        for t in terms
            for i in t.indices
                push!(inds, i)

            end
        end
        new(terms, inds)
    end
end
function TensorTerm(term::TensorSuperType)
    TensorTerm([term])
end

mutable struct TensorExpression
    summands::Vector{Pair{SymorNum,TensorTerm}}
    function TensorExpression(summands::Vector{TensorTerm})
        merge_sum(new(Vector((1, summand) for summand in summands)))
    end
    function TensorExpression(summands::Vector{T} where {T<:TensorSuperType})
        merge_sum(new([1 => TensorTerm(summand) for summand in summands]))
    end
    function TensorExpression(exp::TensorExpression, T::TensorSuperType)
        merge_sum(new([exp.summands; [1 => TensorExpression(T)]]))
    end
    function TensorExpression(exp::TensorExpression, T::TensorExpression)
        merge_sum(new([exp.summands; T.summands]))
    end
    function TensorExpression(num::SymorNum, T::TensorSuperType)
        merge_sum(new([num => TensorTerm(T)]))
    end
    function TensorExpression(vec::Vector{<:Pair{<:SymorNum,<:TensorSuperType}})
        merge_sum(new([v[1] => TensorTerm(v[2]) for v in vec]))
    end
    function TensorExpression(T::TensorSuperType)
        merge_sum(new([1 => TensorTerm(T)]))
    end
    function TensorExpression(T::TensorTerm)
        merge_sum(new([1 => T]))
    end
    function TensorExpression(p::Pair{T,TensorTerm}) where {T<:SymorNum}
        merge_sum(new([p]))
    end
    function TensorExpression(v::Vector{Pair{T,TensorTerm}} where {T<:SymorNum})
        merge_sum(new(v))
    end
    function TensorExpression(p::Pair{N,T}) where {N<:SymorNum,T<:TensorSuperType}
        merge_sum(new([p[1] => TensorTerm(p[2])]))
    end
    function TensorExpression(tex::TensorExpression)
        tex
    end
    function TensorExpression()
        new([0 => TensorTerm(Tensor(""))])
    end
    function TensorExpression(n::SymorNum)
        new([n => TensorTerm(Tensor(""))])
    end
end

# function get_dummies_positional(T::TensorExpression)
#     dum::Set{Pair{Index, IndexPosition}}=Set()
#     for d in get_dummies(T.summands[1]) 
#         push!(dum, d=>upper)
#         push!(dum, d=>lower)
#     end
#     return dum
# end









#t.upper_indices for t in terms

#  function update_indices(T::TensorTerm)
#     for t in T.terms

#  end


#  *(A::Tensor, B::TensorSuperType)= TensorTerm([A,B])

#  function remove_element(vector::Vector{T}, element::T) where T
#     index = findfirst(vector .== element)
#     if index === nothing
#         return false  # element not found
#     else
#         splice!(vector, index)
#         return true  # element removed successfully
#     end
# end

#  function getFreeIndices(Ts::TensorTerm)
#     freeup::OrderedSet{Index}=OrderedSet([])
#     freelow::OrderedSet{Index}=OrderedSet([])
#     for term in Ts.terms
#         for ind in term.upper_indices
#         if !(ind in freeup)
#             push!(freeup, ind)
#         else
#             delete!(freeup, ind)
#         end


#         end

#         for ind in term.lower_indices
#             if !(ind in freeup) && !(ind in freelow)
#                 push!(freelow, ind)
#             elseif (ind in freeup)
#                 delete!(freeup, ind)
#             else 
#                 delete!(freelow, ind)
#             end
#             end
#     end
#     return [freeup, freelow]
#  end

#  function getDummieIndices(Ts::TensorTerm)
#     dumup::OrderedSet{Index}=OrderedSet([])
#     dumlow::OrderedSet{Index}=OrderedSet([])
#     tempup::OrderedSet{Index}=OrderedSet([])
#     templow::OrderedSet{Index}=OrderedSet([])

#     for term in Ts.terms
#         for ind in term.upper_indices
#         if !(tempup in freeup)
#             push!(freeup, ind)
#         else
#             delete!(freeup, ind)
#         end


#         end

#         for ind in term.lower_indices
#             if !(ind in freeup) && !(ind in freelow)
#                 push!(freelow, ind)
#             elseif (ind in freeup)
#                 delete!(freeup, ind)
#             else 
#                 delete!(freelow, ind)
#             end
#             end
#     end
#     return [freeup, freelow]
#  end




# #  function Base.show(io::IO, T::TensorTerm)
# #     for t in T
# #         print(io, "MyType(")
# #     end

# #     show(io, x.x)
# #     print(io, ", ")
# #     show(io, x.y)
# #     print(io, ")")
# # end

# https://docs.julialang.org/en/v1/manual/metaprogramming/
Base.@doc "Usage: `@Tensor A^μν_ρ^τ` where `A` is a symbol and `μ,ν,ρ,τ` are indices."
macro Tensor(expr::Expr)
    @show expr
    @show expr.args
    @show string(expr)
    symbol = string(expr)[1] # extract the tensor name
    indices = []
    if length(string(expr.args[2])) > 1 && #"name" part should only be one character - if more, then should be _
       string(expr.args[2])[nextind(string(expr.args[2]), 1)] == '_'
        push!(indices, [v => lower for v in split(string(expr.args[2])[nextind(string(expr.args[2]), 1, 2):end], "")]...)
    else
        @assert length(string(expr.args[2])) == 1 "Tensor symbol must be a single character, not $(expr.args[2])"
    end
    head = expr.args[3]
    while hasproperty(head, :args) && head.args[1] == :^
        # @show head
        indices = parse_indices(head, indices, false)

        head = head.args[3]
    end
    # @show head
    indices = parse_indices(head, indices, true)
    return quote
        Tensor(
            $(string(symbol)),
            [$([:($(esc(Symbol(index.first))) => $(esc(Symbol(index.second)))) for index in indices]...)]
        )
    end
end

function parse_indices(head, indices, leaf)
    str_head = [v for (i, v) in enumerate(string(leaf ? head : head.args[2]))]
    # @show str_head
    if length(str_head) == 1
        indices = [indices..., str_head[1] => upper]
    else
        underscore = findfirst('_' .== str_head)
        # underscore = underscore === nothing ? length(string(head)) : underscore
        if underscore != nothing
            underscoreless_str = [str_head[1:underscore-1]..., str_head[underscore+1:end]...]
            indices = [indices..., [v => i ≥ underscore ? lower : upper for (i, v) in enumerate(underscoreless_str)]...]
        else
            indices = [indices..., [v => upper for (i, v) in enumerate(str_head)]...]
        end
    end
    return indices
end

Base.@doc"Usage: `@Tensor A_bc` where `A` is a symbol and `b`/`c` are indices."
macro Tensor(expr::Symbol)
    length(string(expr)) == 1 && return Tensor(string(expr))
    index_count = length(string(expr)) - 2
    # split the symbol A_bc at the underscore
    symbol, indices = split(string(expr), "_")
    return quote
        Tensor(
            $(string(symbol)),
            [$([:($(esc(Symbol(index))) => lower) for index in indices]...)]
        )
    end
end

macro Tensor(expr::String)
    length(expr) == 1 && return Tensor(expr) #simple no index tensor
    symbol = expr[1:1]

    indexpos = 2#depreciate this
    #allow more flexible names (may have bars and dots maybe or other eg "det" for op name etc)
    @assert contains("_^([{", expr[indexpos]) "Tensor name should be one character."#throw error if not starting with following
    positions = []
    indices = []
    currpos = upper
    pause = false
    asym = false
    sym = false
    slotnum = 0 #current real index pos (for sym/asym logic)
    functionoflist = []
    funof = false
    tabstart = 0
    tabaccum = []
    tableaulist = []
    coordstring = ""
    for (i, char) in enumerate(expr)
        i == 1 && continue #skip the first character since its the symbol, and we already extracted it
        char = string(char)
        if char == "^"
            currpos = upper
        elseif char == "_"
            currpos = lower
        elseif char == "|"#pause/unpause sym/asym
            pause = (!pause)
        elseif char == "["
            tabaccum = Vector{Int64}()
            asym = true
            tabstart = slotnum + 1
        elseif char == "]"
            asym = false
            push!(tableaulist, YoungTableau(fill(1, slotnum - tabstart + 1), tabaccum))
        elseif char == "("
            tabaccum = Vector{Int64}()
            sym = true
            tabstart = slotnum + 1
        elseif char == ")"
            sym = false
            push!(tableaulist, YoungTableau([length(tabaccum)], tabaccum))
        elseif char == "{"
            funof = true
            println("starting function of")
        else
            #increase the real indexpos
            if !funof
                slotnum += 1
                push!(indices, Symbol(char))
                push!(positions, currpos)
                if !pause && asym
                    push!(tabaccum, slotnum)
                elseif !pause && sym
                    push!(tabaccum, slotnum)
                end
            elseif char == "," || char == "}"
                push!(functionoflist, :($(esc(Symbol(coordstring)))))
                println("added function of: ", coordstring)
                coordstring = ""
                if char == "}"
                    funof = false
                end
            else
                coordstring *= char
            end
        end
        # indexpos+=1
        # @show tableaulist
    end
    functions_expr = Expr(:vect, functionoflist...)
    #pair indices and pos into a vector
    indplacement = quote
        [$([:($(esc(Symbol(idx))) => $pos) for (idx, pos) in zip(indices, positions)]...)]
    end
    return :(Tensor($symbol, [$indplacement...], [$tableaulist...], $functions_expr))
end