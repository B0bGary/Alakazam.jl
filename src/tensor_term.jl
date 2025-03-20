# include("indices.jl")
# include("custom_exceptions.jl")

SymorNum = Union{SymbolicUtils.Symbolic,Number}
abstract type TensorSuperType end
abstract type TotallyAntisymmetricSuperType <: TensorSuperType end
abstract type TotallySymmetricSuperType <: TensorSuperType end

mutable struct DummyPatternTensor <: TensorSuperType
    name::String
    indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}

    function DummyPatternTensor(name::String="?", indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}=Vector{Pair{IndexSuperType,IndexPosition}}(),
        function_of::Vector{Coordinate}=Vector{Coordinate}())  #lets you have scalars with no indices
        new(name, indices, function_of, Dict{String,Number}())
    end

    function DummyPatternTensor(name::String="?", indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}=Vector{Pair{IndexSuperType,IndexPosition}}())  #lets you have scalars with no indices
        new(name, indices, [], Dict{String,Number}())
    end
    function DummyPatternTensor(name::String="?", indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}=Vector{Pair{IndexSuperType,IndexPosition}}(), function_of::Coordinate=Coordinate("", Vector{IndexSet}()))  #lets you have scalars with no indices
        new(name, indices, [function_of], Dict{String,Number}())
    end
    function DummyPatternTensor(name::String="?", indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}=Vector{Pair{IndexSuperType,IndexPosition}}(), function_of::Coordinate=Coordinate("", Vector{IndexSet}()), weights::Dict{String,Number}=Dict{String,Number}())  #lets you have scalars with no indices
        new(name, indices, [function_of],)
    end

end



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
            indices = vcat(indices,arg)
        elseif arg isa Pair{I,IndexPosition} where {I<:IndexSuperType}
            push!(indices,arg)
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


# struct SymmetricTensor <:TotallySymmetricSuperType
# end
mutable struct SymmetricTensor <: TotallySymmetricSuperType
    name::String
    indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}
    tableaux::Vector{YoungTableau}

    function SymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(),
        function_of::Vector{Coordinate}=Vector{Coordinate}()) where {I<:IndexSuperType}#lets you have scalars with no indices
        w = Dict{String,Number}()
        for coord in function_of
            w[coord.name] = 1
        end
        new(name, indices, function_of, w, [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end


    function SymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(),
        function_of::Coordinate=nothing) where {I<:IndexSuperType} #lets you have scalars with no indices
        w = Dict{String,Number}()
        if function_of != nothing
            w[function_of.name] = 1
        end
        new(name, indices, [function_of], w, [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function SymmetricTensor(name::String,
        function_of::Coordinate=nothing) #lets you have scalars with no indices
        w = Dict{String,Number}()
        if function_of != nothing
            w[function_of.name] = 1
        end
        new(name, Vector{Pair{IndexSuperType,IndexPosition}}(), [function_of], w, [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function SymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}()) where {I<:IndexSuperType}#lets you have scalars with no indices

        new(name, indices, [], Dict{String,Number}(), [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end
    function SymmetricTensor(T::DummyPatternTensor) #
        new(T.name, T.indices, T.function_of, T.weights, [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function SymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition},}(),
        function_of::Vector{Coordinate}=Vector{Coordinate}(), weights::Dict{String,Number}=Dict{String,Number}()) where {I<:IndexSuperType}#lets you have scalars with no indices
        w = Dict{String,Number}()
        for coord in function_of
            w[coord.name] = 1
        end
        new(name, indices, function_of, merge(w, weights), [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end


    function SymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(),
        function_of::Coordinate=nothing, weights::Dict{String,Number}=Dict{String,Number}()) where {I<:IndexSuperType} #lets you have scalars with no indices
        w = Dict{String,Number}()
        if function_of != nothing
            w[function_of.name] = 1
        end

        new(name, indices, [function_of], merge(w, weights), [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function SymmetricTensor(name::String,
        function_of::Coordinate=nothing, weights::Dict{String,Number}=Dict{String,Number}()) #lets you have scalars with no indices
        w = Dict{String,Number}()
        if function_of != nothing
            w[function_of.name] = 1
        end
        new(name, Vector{Pair{IndexSuperType,IndexPosition}}(), [function_of], merge(w, weights), [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function SymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(), weights::Dict{String,Number}=Dict{String,Number}()) where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices, [], weights, [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function SymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(), function_of::Vector{Coordinate}=Vector{Coordinate}(),
        weights::Dict{String,Number}=Dict{String,Number}(), tableaux::Vector{YoungTableau}=Vector{YoungTableau}()) where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices, function_of, weights, [YoungTableau([size(indices)[1]], collect(1:size(indices)[1]))])
    end

    # weights::Dict{String,Number}=Dict()

end

# struct AntisymmetricTensor <:TotallyAntisymmetricSuperType
# end
mutable struct AntisymmetricTensor <: TotallyAntisymmetricSuperType
    name::String
    indices::Vector{Pair{I,IndexPosition}} where {I<:IndexSuperType}
    function_of::Vector{Coordinate}
    weights::Dict{String,Number}
    tableaux::Vector{YoungTableau}

    function AntisymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(),
        function_of::Vector{Coordinate}=Vector{Coordinate}()) where {I<:IndexSuperType}#lets you have scalars with no indices
        w = Dict{String,Number}()
        for coord in function_of
            w[coord.name] = 1
        end
        new(name, indices, function_of, w, [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end


    function AntisymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(),
        function_of::Coordinate=nothing) where {I<:IndexSuperType} #lets you have scalars with no indices
        w = Dict{String,Number}()
        if function_of != nothing
            w[function_of.name] = 1
        end
        new(name, indices, [function_of], w, [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function AntisymmetricTensor(name::String,
        function_of::Coordinate=nothing) #lets you have scalars with no indices
        w = Dict{String,Number}()
        if function_of != nothing
            w[function_of.name] = 1
        end
        new(name, Vector{Pair{IndexSuperType,IndexPosition}}(), [function_of], w, [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function AntisymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}()) where {I<:IndexSuperType}#lets you have scalars with no indices

        new(name, indices, [], Dict{String,Number}(), [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end
    function AntisymmetricTensor(T::DummyPatternTensor) #
        new(T.name, T.indices, T.function_of, T.weights, [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function AntisymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition},}(),
        function_of::Vector{Coordinate}=Vector{Coordinate}(), weights::Dict{String,Number}=Dict{String,Number}()) where {I<:IndexSuperType}#lets you have scalars with no indices
        w = Dict{String,Number}()
        for coord in function_of
            w[coord.name] = 1
        end
        new(name, indices, function_of, merge(w, weights), [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end


    function AntisymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(),
        function_of::Coordinate=nothing, weights::Dict{String,Number}=Dict{String,Number}()) where {I<:IndexSuperType} #lets you have scalars with no indices
        w = Dict{String,Number}()
        if function_of != nothing
            w[function_of.name] = 1
        end

        new(name, indices, [function_of], merge(w, weights), [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function AntisymmetricTensor(name::String,
        function_of::Coordinate=nothing, weights::Dict{String,Number}=Dict{String,Number}()) #lets you have scalars with no indices
        w = Dict{String,Number}()
        if function_of != nothing
            w[function_of.name] = 1
        end
        new(name, Vector{Pair{IndexSuperType,IndexPosition}}(), [function_of], merge(w, weights), [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end

    function AntisymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(), weights::Dict{String,Number}=Dict{String,Number}()) where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices, [], weights, [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end

    # weights::Dict{String,Number}=Dict()
    function AntisymmetricTensor(name::String, indices::Vector{Pair{I,IndexPosition}}=Vector{Pair{IndexSuperType,IndexPosition}}(), function_of::Vector{Coordinate}=Vector{Coordinate}(),
        weights::Dict{String,Number}=Dict{String,Number}(), tableaux::Vector{YoungTableau}=Vector{YoungTableau}()) where {I<:IndexSuperType}#lets you have scalars with no indices
        new(name, indices, function_of, weights, [YoungTableau([1 for _ in 1:size(indices)[1]], collect(1:size(indices)[1]))])
    end
end

# function Tensor(name::String, indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}(),
#     function_of::Coordinate=nothing) #lets you have scalars with no indices
#     new(name, indices,[function_of])
# end

# function Tensor(name::String, indices::Vector{Pair{Index, IndexPosition}}=Vector{Pair{Index, IndexPosition}}()) #lets you have scalars with no indices
#     new(name, indices,[])
# end




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
    function TensorTerm(term::TensorSuperType)
        new([term], term.indices)
    end
end


mutable struct TensorExpression
    summands::Vector{Pair{SymorNum,TensorTerm}}
    function TensorExpression(summand::TensorTerm)
        merge_sum(new(Vector((1, summand))))
    end
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
"Usage: `@Tensor A^μν_ρ^τ` where `A` is a symbol and `μ,ν,ρ,τ` are indices."
macro Tensor(expr::Expr)
    # @show expr
    symbol = string(expr)[1]
    indices = []
    if length(string(expr.args[2])) > 1 && 
       string(expr.args[2])[nextind(string(expr.args[2]),1)] == '_'
        push!(indices, [ v => lower for v in split(string(expr.args[2])[nextind(string(expr.args[2]),1,2):end], "")]...)
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
    str_head = [ v for (i,v) in enumerate(string(leaf ? head : head.args[2])) ]
    # @show str_head
    if length(str_head) == 1
            indices = [indices..., str_head[1] => upper]
    else
        underscore = findfirst('_' .== str_head)
        # underscore = underscore === nothing ? length(string(head)) : underscore
        if underscore != nothing
            underscoreless_str = [str_head[1:underscore-1]..., str_head[underscore+1:end]...]
            indices = [indices..., [ v => i ≥ underscore ? lower : upper for (i,v) in enumerate(underscoreless_str)]...]
        else
            indices = [indices..., [ v => upper for (i,v) in enumerate(str_head)]...]
        end
    end
    return indices
end

"Usage: `@Tensor A_bc` where `A` is a symbol and `b`/`c` are indices."
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