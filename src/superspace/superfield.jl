

"""
Generates a product of grassmann variables (I think this can be replaced with appension...?)
"""
function generate_grassmann_product(grassmann::TensorSuperType, len::Integer)
    # length exceeds number of terms allowed by antisymmetry
    (len > prod(st.range for st in get_index_sets(grassmann))) && return Alakazam. zero()
    len == 0 && return 1
    term = TensorTerm(grassmann)
    for i in 1:(len-1)
        term *= set_indices!(copy(grassmann), Pair{IndexSuperType,IndexPosition}[!isLightcone(ind[1]) ? safe_pop!(get_unused_indices(term, ind[1].index_set), ind[1].index_set) => ind[2] : ind for ind in indices(grassmann)])
    end
    return term
end



generate_superfield_naive(name::String, SS::Superspace, xc::Coordinate) = generate_superfield_naive(name, SS.grassmann_basis, xc)
generate_superfield_naive(name::String, extra_inds::Vector{Pair{I,IndexPosition}}, SS::Superspace, xc::Coordinate) where I<:IndexSuperType = generate_superfield_naive(name, extra_inds, SS.grassmann_basis, xc)
generate_superfield_naive(name::String, grassmann_basis::Vector{T}, xc::Coordinate) where T<:TensorSuperType = generate_superfield_naive(name, Vector{Pair{IndexSuperType,IndexPosition}}(), grassmann_basis, xc)

"""
Generates a superfield expansion over the specified superspace, without performing any reductions based on symmetries.
User can supply external indices that the superfield carries. One can also construct using an explicit superspace.
The Grassmanns are written before the component field.
"""
function generate_superfield_naive(name::String, extra_inds::Vector{Pair{I,IndexPosition}}, grassmann_basis::Vector{T}, xc::Coordinate) where T<:TensorSuperType where I<:IndexSuperType
    gdict=Dict()
    gdictbar=Dict()
    gs=[gdict, gdictbar]
    tex= zero()
    length(grassmann_basis) <= 2 || error("Basis should be a single real Grassmann, or two complex conjugate Grassmann.")
    for i in 1:length(grassmann_basis)
        g=grassmann_basis[i]
        spinor_dim = grassmann_dimension(indices(g)[1])

        flavour_dim=length(indices(g))==2 ? dimension(indices(g)[2][1]) : 1 #needs to work for flav or no flav
        max_order=spinor_dim*flavour_dim
        for or in 0:max_order
            base_term=generate_grassmann_product(g, or)
            # order_terms=generate_independant_structure(base_term) #vector of tuples (term,spinortab,flavtab)
            gs[i][or]=[base_term]
        end
    end
    # @show gdict
    for (o1, v1) in gdict
        for trm in v1
            spt1=indices(grassmann_basis[1])[1][1].index_set #spinor index type
            flt1 = length(indices(grassmann_basis[1]))==2 ? indices(grassmann_basis[1])[2][1].index_set : nothing #flav index type if exists
            if length(keys(gdictbar))>0
                # @show gdictbar
                for (o2, v2) in gdictbar
                    for trm2 in v2
                        spt2=indices(grassmann_basis[2])[1][1].index_set #spinor index type
                        flt2 = length(indices(grassmann_basis[2]))==2 ? indices(grassmann_basis[2])[2][1].index_set : nothing #flav index type if exists

                        newterm=Tensor("_", extra_inds)×trm×trm2 #tensor appension so no repeated inds.
                        # @show newterm
                        newterm=TensorTerm(terms(isTensorTerm(newterm))[2:end]) #hack to get ext inds unused
                        # @show free_indices_ordered(newterm)
                        freeinds = op_pos(free_indices_ordered(newterm))
                        sp1inds = [i for i in freeinds if i[1].index_set == spt1]
                        if spt2 == spt1
                            sp2inds = []
                        else
                            sp2inds = [i for i in freeinds if i[1].index_set == spt2]
                        end
                        fl1inds = flt1 === nothing ? [] :
                                  [i for i in freeinds if i[1].index_set == flt1]
                        if flt2 == flt1
                            fl2inds = []
                        else
                            fl2inds = flt2 === nothing ? [] :
                                      [i for i in freeinds if i[1].index_set == flt2]
                        end

                        combined_inds = Pair{IndexSuperType,IndexPosition}[vcat(extra_inds, sp1inds, sp2inds, fl1inds, fl2inds)...]
                        ψ=Tensor(name, combined_inds, xc, "("*string(o1)*","*string(o2)*")"=>lower)
                        tex += newterm * ψ
                    end
                end
            else
                # @show trm
                newterm=Tensor("_", extra_inds)×trm #tensor appension so no repeated inds.
                # @show newterm
                newterm=TensorTerm(terms(isTensorTerm(newterm))[2:end])
                freeinds = op_pos(free_indices_ordered(newterm))

                sp1inds = [i for i in freeinds if i[1].index_set == spt1]

                fl1inds = flt1 === nothing ? [] :
                          [i for i in freeinds if i[1].index_set == flt1]

                combined_inds = Pair{IndexSuperType,IndexPosition}[
                    vcat(extra_inds, sp1inds, fl1inds)...
                ]
                ψ = Tensor(name, combined_inds, xc, "("*string(o1)*")"=>lower)
                tex += newterm * ψ
            end
        end

    end

    return tex
end


generate_superfield(name::String, superspace::Superspace, xc::Coordinate;
    kwargs...) =
    generate_superfield(name, Pair{IndexSuperType,IndexPosition}[], superspace.grassmann_basis, xc; kwargs...)
generate_superfield(name::String,  extra_inds::Vector{Pair{I,IndexPosition}},superspace::Superspace, xc::Coordinate;
    kwargs...) where I<:IndexSuperType =
    generate_superfield(name, extra_inds, superspace.grassmann_basis, xc; kwargs...)

generate_superfield(name::String, grassmann_basis::Vector{T}, xc::Coordinate;
    kwargs...) where T<:TensorSuperType =
    generate_superfield(name, Pair{IndexSuperType,IndexPosition}[], grassmann_basis, xc; kwargs...)

function generate_superfield(name::String, extra_inds::Vector{Pair{I,IndexPosition}},
    grassmann_basis::Vector{T}, xc::Coordinate;
    S_type=nothing, F_type=nothing,
    pair_tensor=nothing) where {T<:TensorSuperType,I<:IndexSuperType}
    gdict = Dict()
    gdictbar = Dict()
    gs = [gdict, gdictbar]
    tex = Alakazam. zero()
    nextra = length(extra_inds)
    length(grassmann_basis) <= 2 || error("Basis should be a single real Grassmann, or two complex conjugate Grassmann.")
    for i in 1:length(grassmann_basis)
        g = grassmann_basis[i]
        spinor_dim = !isLightcone(indices(g)[1]) ? dimension(indices(g)[1][1]) :
                     indices(g)[1][1].name == "±" || indices(g)[1][1].name == "∓" ? 2 : 1
        flavour_dim = length(indices(g)) == 2 ? dimension(indices(g)[2][1]) : 1  # works for flav or no flav
        max_order = spinor_dim * flavour_dim
        for or in 0:max_order
            base_term = generate_grassmann_product(g, or)
            # vector of tuples (term, spinortab, flavtab)
            order_terms = generate_independant_structure(base_term;
                S_type=S_type, F_type=F_type, pair_tensor=pair_tensor)
            gs[i][or] = order_terms
        end
    end
    for (o1, v1) in gdict
        for trm in v1
            t1, sp1, fl1 = trm
            spt1 = indices(grassmann_basis[1])[1][1].index_set
            flt1 = length(indices(grassmann_basis[1])) == 2 ? indices(grassmann_basis[1])[2][1].index_set : nothing
            if length(keys(gdictbar)) > 0
                for (o2, v2) in gdictbar
                    for trm2 in v2
                        t2, sp2, fl2 = trm2
                        spt2 = indices(grassmann_basis[2])[1][1].index_set
                        flt2 = length(indices(grassmann_basis[2])) == 2 ? indices(grassmann_basis[2])[2][1].index_set : nothing
                        # tensor appension so no repeated inds; the dummy "_"
                        # reserves the extra indices, then is stripped
                        newterm = reserve_extra_indices(extra_inds, t1 × t2)
                        freeinds = op_pos(free_indices_ordered(newterm))
                        sp1inds = [i for i in freeinds if i[1].index_set == spt1]
                        sp2inds = spt2 == spt1 ? [] :
                                  [i for i in freeinds if i[1].index_set == spt2]
                        fl1inds = flt1 === nothing ? [] :
                                  [i for i in freeinds if i[1].index_set == flt1]
                        fl2inds = (flt2 == flt1 || flt2 === nothing) ? [] :
                                  [i for i in freeinds if i[1].index_set == flt2]
                        combined_inds = Pair{IndexSuperType,IndexPosition}[vcat(extra_inds, sp1inds, sp2inds, fl1inds, fl2inds)...]
                        new_tab = YoungTableau[]
                        append_offset!(new_tab, sp1; base=nextra)
                        append_offset!(new_tab, sp2; base=nextra)
                        append_offset!(new_tab, fl1; base=nextra)
                        append_offset!(new_tab, fl2; base=nextra)
                        ψ = Tensor(name, combined_inds, new_tab, xc, "(" * string(o1) * "," * string(o2) * ")" => lower)
                        tex += newterm * ψ
                    end
                end
            else
                newterm = reserve_extra_indices(extra_inds, t1)

                freeinds = op_pos(free_indices_ordered(newterm))
                sp1inds = [i for i in freeinds if i[1].index_set == spt1]
                fl1inds = flt1 === nothing ? [] :
                          [i for i in freeinds if i[1].index_set == flt1]
                combined_inds = Pair{IndexSuperType,IndexPosition}[vcat(extra_inds, sp1inds, fl1inds)...]
                new_tab = YoungTableau[]
                append_offset!(new_tab, sp1; base=nextra)
                append_offset!(new_tab, fl1; base=nextra)
                ψ = Tensor(name, combined_inds, new_tab, xc, "(" * string(o1) * ")" => lower)
                tex += newterm * ψ
            end
        end
    end
    return eliminate_epsilon(tex)
end


"""
Convert the algebra associated to indices to a symbol that will be used in rep theory calculations for Grassmann structures.
"""
group_symbol(::SLAlgebra) = :SL
group_symbol(::SUAlgebra) = :SU
group_symbol(::SOAlgebra) = :SO
group_symbol(::SpAlgebra) = :Sp
group_symbol(::LieAlgebraSuperType) = :none

"""
Checks group consistency - makes sure index range matches the group dim etc.
Usually, everything should be in the fundamental rep in these Grassmann constructions.
"""
function derived_group(index_set, override, dim::Int)
    override !== nothing && return override
    alg = index_set.algebra
    sym = group_symbol(alg)
    if sym != :none && hasproperty(alg, :N) && alg.N != dim
        error("index set '$(index_set.name)': algebra $(typeof(alg)) has N=$(alg.N) but the " *
              "index dimension is $dim. This construction only works for the fundamental rep. Override needed.")
    end
    return sym
end


"""
For "1"
"""
generate_independant_structure(term::Integer; kwargs...) = [(TensorExpression(term), nothing, nothing)]

function generate_independant_structure(term::TensorExpression; kwargs...)
    term == Alakazam. zero() && return []
    isTerm(term) == false && error("Input should be a tensor term of grassmann numbers")
    return generate_independant_structure(isTerm(term); kwargs...)
end

"""
From a product of grassmann numbers, use the algebra associated to spinor and flavour to decompose into independant structures.
This is done using the invariant tensors of each algebra (ε/ω/δ).
Reduction is only performed if it will decrease the number of indices on the coefficient field.
The reduced theta products are returned, as well as the associated tableau to be applied to the coefficient field.
"""
function generate_independant_structure(term::TensorTerm; S_type=nothing, F_type=nothing,
    pair_tensor=nothing)
    order = length(terms(term))
    order == 0 && return [(TensorExpression(1), nothing, nothing)]
    t1 = terms(term)[1]
    ninds = length(indices(t1))
    spinor_set = indices(t1)[1][1].index_set
    spinor_dim = dimension(indices(t1)[1][1])
    flav_set = ninds > 1 ? indices(t1)[2][1].index_set : nothing
    flav_dim = ninds > 1 ? dimension(indices(t1)[2][1]) : 0

    S = derived_group(spinor_set, S_type, spinor_dim)
    F = flav_set === nothing ? :none : derived_group(flav_set, F_type, flav_dim)

    # one index Pair per θ copy, opposite position, indexable by θ label
    epsspininds = op_pos([i for i in indices(term) if i[1].index_set == spinor_set])
    epsflavinds = flav_set === nothing ? nothing :
                  op_pos([i for i in indices(term) if i[1].index_set == flav_set])
    structs = independent_structures(order, spinor_dim, flav_dim, S, F)
    tex = []
    for s in structs
        nt = copy(term)
        nt = apply_structure_side(nt, s.spinor, epsspininds, spinor_set, "spinor", pair_tensor)
        nt = apply_structure_side(nt, s.flavour, epsflavinds, flav_set, "flavour", pair_tensor)
        push!(tex, (eliminate_epsilon(nt),
            residual_tableau(s.spinor),
            residual_tableau(s.flavour)))
    end
    return tex
end



"""
Applies invariant structures to a product of grassmann variables.
"""
function apply_structure_side(nt, side, side_inds, index_set, side_name::String, pair_tensor_override)
    (side === nothing || side.tableau === nothing) && return nt

    side.residual_rows === nothing && error(
        "generate_independant_structure: $(side_name) branch shape $(side.residual_shape) " *
        "($(side.kind)) is a validated composite with no elementary witness: no explicit " *
        "tensor product represents this copy. Its multiplicity is real, so it cannot be " *
        "silently dropped; handle it manually or choose different group types.")

    eps_ctor = epsilon(index_set.algebra)
    for group in side_epsilon_columns(side)
        theta_pairs = [side_inds[k] for k in group if k > 0]
        pos = theta_pairs[1][2]          # ε groups always contain >=1 θ label
        ordered = Pair{IndexSuperType,IndexPosition}[]
        for k in group
            if k > 0
                push!(ordered, side_inds[k])
            else                          # free slot: fresh index, contracted C--ε
                push!(ordered, safe_pop!(get_unused_indices(nt, index_set), index_set) => pos)
            end
        end
        nt = eps_ctor(ordered) * nt
    end

    pair_ctor = pair_tensor_override === nothing ? metric(index_set.algebra) : pair_tensor_override
    for pair in side.pair_groups
        side.is_composite && error(
            "generate_independant_structure: $(side_name) side carries a non-elementary Ω " *
            "invariant consuming θ-copies $(pair): no explicit elementary tensor product " *
            "represents it; handle this branch manually.")
        pair_ctor === nothing && error(
            "generate_independant_structure: $(side_name) side needs a " *
            (side.kind == :Sp ? "ω (symplectic)" : "δ (metric)") *
            " pair tensor on θ-copies $(pair), but metric($(typeof(index_set.algebra))) " *
            "returned nothing and no `pair_tensor` override was given.")
        nt = pair_ctor([side_inds[k] for k in pair]) * nt
    end
    return nt
end


"""
Reconstruct the coefficient tableau from the grassmann helper structure.
"""
function residual_tableau(side)
    (side === nothing || side.residual_rows === nothing || isempty(side.residual_rows)) && return nothing
    shape = [length(r) for r in side.residual_rows]
    k = sum(shape)
    return YoungTableau(shape, collect(1:k), k)
end
