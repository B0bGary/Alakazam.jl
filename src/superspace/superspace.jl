

mutable struct Superspace
    grassmann_basis::Vector{TensorSuperType}
    covariant_basis::Vector{CovDerivSuperType}
    dimension::Int
    function Superspace(grassmann_basis::Vector{T}, covariant_basis::Vector{D}) where {T<:TensorSuperType,D<:CovDerivSuperType}
        dimct=1
        for grass in grassmann_basis
            grading(grass) != 1 && error("Superspace basis must be Grassmann odd")
            for index in indices(grass)
                dimct *= !isLightcone(index) ? dimension(index[1]) : index[1].name=="±" || index[1].name=="∓" ? 2 : 1
            end
        end
        new(grassmann_basis, covariant_basis, dimct)
    end
    Superspace(grassmann_basis::TensorSuperType, covariant_basis::CovDerivSuperType) = Superspace([grassmann_basis], [covariant_basis])
end


"""
Number of independent components the index `i` of a Grassmann coordinate carries.
"""
function grassmann_dimension(i::Pair{I,IndexPosition}) where {I<:IndexSuperType}
    isLightcone(i) || return dimension(i[1])
    return (i[1].name == "±" || i[1].name == "∓") ? 2 : 1
end

"""
The name of the Grassmann coordinate a covariant derivative differentiates along. These are built
as `Derivative("D", [x, θ], …)`, the bosonic coordinate first and the Grassmann one last.
"""
_grassmann_coord_name(c::CovDerivSuperType) = coord(c)[end].name

"""
    grassmann_order(g) -> Int

The highest power of the Grassmann `g` that can be non-zero, namely its number of independent
components. Anticommutativity kills anything above it.
"""
grassmann_order(g::TensorSuperType) = prod(grassmann_dimension(i) for i in indices(g); init=1)

"""
    drop_nilpotent(exp, basis) -> TensorExpression

Discard every term carrying a higher power of some Grassmann in `basis` than it can hold.
"""
function drop_nilpotent(exp::TensorExpression, basis::Vector{T}) where {T<:TensorSuperType}
    out = exp
    for g in basis
        isempty(function_of(g)) &&
            throw(ArgumentError("Grassmann $g is not a function of any coordinate, so it carries no weight"))
        out = drop_weight_above(out, function_of(g)[1].name, grassmann_order(g))
    end
    return out
end
drop_nilpotent(exp::Tensorial, basis::Vector{<:TensorSuperType}) = drop_nilpotent(TensorExpression(exp), basis)
drop_nilpotent(exp, SS::Superspace) = drop_nilpotent(exp, SS.grassmann_basis)

"""
    derive(exp, basis)
    derive(exp, superspace)

Perform derive() as usualy, and then discard the terms carrying nilpotent products of grassmanns.
"""
derive(exp, basis::Vector{T}) where {T<:TensorSuperType} = drop_nilpotent(derive(exp), basis)
derive(exp, SS::Superspace) = derive(exp, SS.grassmann_basis)

"""
    grassmann_coefficient(exp, basis, weights) -> TensorExpression

The coefficient of `∏θi^{wi}` in a superfield expansion, for the Grassmanns `basis` and the powers `weights`. 
This is the naive component, the coefficient as written, in the convention that the Grassmanns stand to the left of it.
This is similar to project_to_component, which extracts the coeffs in rigourous way via integration/differentiation.
"""
function grassmann_coefficient(exp::TensorExpression, basis::Vector{T},
    weights::Vector{Int}) where {T<:TensorSuperType}
    length(basis) == length(weights) ||
        throw(ArgumentError("expected one weight per Grassmann, got $(length(weights)) for $(length(basis))"))
    gnames = String[name(g) for g in basis]
    picked = exp
    for (g, w) in zip(basis, weights)
        isempty(function_of(g)) && throw(ArgumentError("Grassmann $g is not a function of any coordinate"))
        picked = get_weighted_terms(picked, function_of(g)[1].name, w)
    end
    total = sum(weights)
    acc = CoeffPair[]
    for (c, tt) in summands(simplify(picked))
        fs = terms(tt)
        gpos = [k for k in eachindex(fs) if name(fs[k]) in gnames]
        rpos = [k for k in eachindex(fs) if !(name(fs[k]) in gnames)]
        length(gpos) == total ||
            throw(ArgumentError("term $tt carries $(length(gpos)) Grassmanns, expected $total"))
        # ordered as `basis` lists them, stably so repeats keep their relative order
        sort!(gpos; by=k -> findfirst(==(name(fs[k])), gnames), alg=Base.Sort.MergeSort)
        sgn = _graded_perm_sign(fs, vcat(gpos, rpos))
        push!(acc, CoeffPair(sgn * c, TensorTerm(TensorSuperType[fs[k] for k in rpos])))
    end
    isempty(acc) && return zero()
    return simplify(TensorExpression(acc))
end
grassmann_coefficient(exp::Tensorial, basis, weights) =
    grassmann_coefficient(TensorExpression(exp), basis, weights)

grassmann_coefficient(exp::Tensorial, basis::Superspace, weights) =
    grassmann_coefficient(TensorExpression(exp), basis.grassmann_basis, weights)

grassmann_coefficient(exp::Tensorial, basis::Superspace, weight1::Int64, weight2::Int64) =
    grassmann_coefficient(TensorExpression(exp), basis.grassmann_basis, [weight1, weight2])

integrate(Tex::Tensorial, SS::Superspace, output::Bool=true) = integrate(TensorExpression(Tex), SS.covariant_basis, output)

#shorthand
∫(Tex::TensorExpression, SS::Superspace, output::Bool=false) = integrate(Tex, SS, output)
∫(Tex::TensorExpression, covariant::Vector{D}, output::Bool=false) where D<:CovDerivSuperType = integrate(Tex, covariant, output)

"""
Perform superspace integration.
"""
function integrate(Tex::TensorExpression, covariant::Vector{D}, output::Bool=true) where D<:CovDerivSuperType
    # susymeasure=measure(Tex,covariant)
    texc=copy(Tex)
    outds = Pair{Int,DerivativeSuperType}[]
    total = sum(prod(grassmann_dimension(i) for i in indices(c)) for c in covariant)
    current = 0
    for c in eachindex(covariant)
        covar=covariant[c]
        copies=prod(grassmann_dimension(i) for i in indices(covar))
        for d in 1:copies
            current += 1
            dnew=fresh_indices(covar, texc)
            push!(outds, c => dnew)
            if output
                print("\rEvaluating superspace integral $current/$total.")
                flush(stdout)
            end

            texc=derive(dnew(texc))
            #`[D, Dbar]` integrates dθ then dθ̄, `[Dbar, D]` the other way round, and the two differ by a sign since D and D̄ anticommute.
            if length(covariant)==2 && c==1 #complex conj vars, we can drop weights
                texc=drop_weight_above(texc, _grassmann_coord_name(covariant[2]), copies) #drops terms that are  zero() by grassmann
            elseif length(covariant)==2 && c==2 #complex conj vars, we can drop weights
                texc=drop_weight_above(texc, _grassmann_coord_name(covariant[1]), copies-d) #drops terms that are  zero() by grassmann... at this stage can drop thetas
            end
        end
    end
    output && println()
    return hard_merge(simplify(_superspace_measure(outds) * texc))
end

"""
The measure completing a full superspace integral or projection. `outds` lists the derivatives applied,
each with the covariant derivative it is a copy of.
"""
function _superspace_measure(outds)
    groups = Dict{Tuple{Int,Int},Vector{IndexPair}}()     # (covariant, slot) => indices
    for (c, d) in outds, (k, p) in enumerate(op_indices(d))
        isLightcone(p) && continue
        push!(get!(groups, (c, k), IndexPair[]), op_pos(p))
    end
    measure = 1
    for v in values(groups)
        length(v) == dimension(v[1]) || continue
        measure *= EpsilonTensor(v) // factorial(length(v))
    end
    return measure
end


integrate(Tex::TensorExpression, covariant::CovDerivSuperType) = integrate(Tex, [covariant])

"""
Partially integrates over superspace to project to the desired component.
"""
function project_to_component(Tex::TensorExpression, covariant::Vector{D}, theta::Int, theta_bar::Int) where D<:CovDerivSuperType
    texc=derive(copy(Tex))
    outds = Pair{Int,DerivativeSuperType}[]
    texc=get_weighted_terms(texc, _grassmann_coord_name(covariant[1]), theta)
    length(covariant) > 1 && (texc=get_weighted_terms(texc, _grassmann_coord_name(covariant[2]), theta_bar))
    (theta==0 && theta_bar==0) && return texc
    for c in eachindex(covariant)
        copies = c==1 ? theta : theta_bar # how many times to project
        for d in 1:copies
            dnew=fresh_indices(covariant[c], texc)
            push!(outds, c => dnew)
            texc=derive(dnew(texc))
        end
    end
    for c in covariant
        texc=drop_weight_above(texc, _grassmann_coord_name(c), 0)
    end
    return hard_merge(simplify(_superspace_measure(outds) * texc))
end