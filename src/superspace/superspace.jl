

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


integrate(Tex::Tensorial, SS::Superspace, output::Bool=true) = integrate(TensorExpression(Tex), SS.covariant_basis, output)
"""
Perform superspace integration.
"""
function integrate(Tex::TensorExpression, covariant::Vector{D}, output::Bool=true) where D<:CovDerivSuperType
    # susymeasure=measure(Tex,covariant)
    texc=copy(Tex)
    outds= DerivativeSuperType[]
    total = sum(prod(dimension(i) for i in indices(c)) for c in covariant)
    current = 0
    for c in eachindex(covariant)
        covar=covariant[c]
        copies=prod(dimension(i) for i in indices(covar))
        for d in 1:copies
            current += 1
            dnew=fresh_indices(covar, texc)
            push!(outds, dnew)
            if output
                print("\rEvaluating superspace integral $current/$total.")
                flush(stdout)
            end

            texc=derive(dnew(texc))
            # @show texc
            if length(covariant)==2 && c==1 #complex conj vars, we can drop weights
                texc=drop_weight_above(texc, "Θ", copies) #drops terms that are  zero() by grassmann
            elseif length(covariant)==2 && c==2 #complex conj vars, we can drop weights
                texc=drop_weight_above(texc, "θ", copies-d) #drops terms that are  zero() by grassmann... at this stage can drop thetas
            end
        end
    end
    output && println()
    #now, epsilon measure factor
    epsinds=Dict{Int,Vector{Pair{Index,IndexPosition}}}()
    for ds in eachindex(outds)
        deriv=outds[ds]
        inds=indices(deriv)
        conj = length(covariant)==1 ? 0 : ds>length(outds)/2 ? 1 : 0
        for i in eachindex(inds)
            runninds=get(epsinds, i+conj*length(inds), Vector{Pair{Index,IndexPosition}}())
            push!(runninds, op_pos(inds[i]))
            epsinds[i+conj*length(inds)]=runninds
        end
    end
    neweps=1
    # @show epsinds
    for (k, v) in epsinds
        isLightcone(v[1])&&continue
        # @show k,v
        neweps*=EpsilonTensor(v)//factorial(length(v))
    end
    return hard_merge(simplify(neweps*texc))
end

# function integrate(Tex::TensorExpression, covariant::Vector{D}) where D<:CovDerivSuperType
#     susymeasure=measure(Tex, covariant)
#     derived=derive(susymeasure)
#     projected=drop_weight_above(derived, "θ", 0)
#     return projected
# end
integrate(Tex::TensorExpression, covariant::CovDerivSuperType) = integrate(Tex, [covariant])
# """
# Determine the superspace measure.
# """
# function measure(Tex::TensorExpression, covariant::Vector{D}) where D<:CovDerivSuperType

#     #for each deriv in covariant, determine the number of deriv copies needed (ie spinor times flav)
#     outds=[]
#     texc=copy(Tex)
#     for covar in covariant
#         copies=prod(dimension(i) for i in indices(covar))
#         for d in 1:copies
#             dnew=fresh_indices(covar, texc)
#             push!(outds, dnew)
#             texc=dnew(texc)
#         end
#     end

#     @show outds
#     @show texc
#     #then, construct a new copy with new indies not used already in the tensor expression

#     #return also the totally antisym symbols and normalisation one needs to also return

# end

"""
Partially integrates over superspace to project to the desired component.
"""
function project_to_component(Tex::TensorExpression, covariant::Vector{D}, theta::Int, theta_bar::Int) where D<:CovDerivSuperType
    # susymeasure=measure(Tex,covariant)
    texc=copy(Tex)
    outds=DerivativeSuperType[]
    total = theta+theta_bar
    current = 0
    texc=get_weighted_terms(texc, "θ", theta)
    texc=get_weighted_terms(texc, "Θ", theta_bar)
    (theta==0 && theta_bar==0) && return texc
    for c in eachindex(covariant)
        covar=covariant[c]
        copies = c==1 ? theta : theta_bar # how many times to project
        for d in 1:copies
            current += 1
            dnew=fresh_indices(covar, texc)
            push!(outds, dnew)
            # print("\rEvaluating superspace integral $current/$total.")
            # flush(stdout)
            texc=derive(dnew(texc))
            # if length(covariant)==2 && c==1 #complex conj vars, we can drop weights
            #     texc=drop_weight_above(texc, "Θ", copies) #drops terms that are  zero() by grassmann
            # elseif length(covariant)==2 && c==2 #complex conj vars, we can drop weights
            #     texc=drop_weight_above(texc, "θ", copies-d) #drops terms that are  zero() by grassmann... at this stage can drop thetas
            # end
        end
    end
    # println() 
    #now, epsilon measure factor
    epsinds=Dict{Int,Vector{Pair{Index,IndexPosition}}}()
    for ds in eachindex(outds)
        deriv=outds[ds]
        inds=indices(deriv)
        conj = length(covariant)==1 ? 0 : ds>length(outds)/2 ? 1 : 0
        for i in eachindex(inds)
            runninds=get(epsinds, i+conj*length(inds), Vector{Pair{Index,IndexPosition}}())
            push!(runninds, op_pos(inds[i]))
            epsinds[i+conj*length(inds)]=runninds
        end
    end
    texc=drop_weight_above(texc, "θ", 0)
    texc=drop_weight_above(texc, "Θ", 0)
    neweps=1
    # @show epsinds
    for (k, v) in epsinds
        isLightcone(v[1])&&continue
        # @show k,v
        neweps*=EpsilonTensor(v)//factorial(length(v))
    end

    return hard_merge(simplify(neweps*texc))
end