@testset "Superspace" begin

    #setup
    if 1==1
        spacetime_indices = IndexSet("spacetimelc", 4,lower)
        # flavour_indices = IndexSet("flavour_indices", 2,upper)
        spinor_indices=IndexSet("spinor_indices", 2,lower)
        spinor_conj_indices=IndexSet("spinor_conj_indices", 2,lower)


        α = Index("α", spinor_indices, 1)
        β = Index("β", spinor_indices, 1)
        γ = Index("γ", spinor_indices, 1)
        χ = Index("χ", spinor_indices, 1)
        αd = Index("̇α", spinor_conj_indices, 1)
        βd = Index("̇β", spinor_conj_indices, 1)
        γd = Index("̇γ", spinor_conj_indices, 1)
        χd = Index("̇χ", spinor_conj_indices, 1)
        

        μ = Index("μ", spacetime_indices)
        ν = Index("ν", spacetime_indices)
        ρ = Index("ρ", spacetime_indices)
        τ = Index("τ", spacetime_indices)

        x = Coordinate("x", [spacetime_indices])
        θ = Coordinate("θ", [spinor_indices])
        θbar = Coordinate("Θ", [spinor_conj_indices])
        custom_sort_order(["θ","Θ","ϵ"])

        #need better way of registration
        grassman_deriv=PartialDerivative("𝛛",zero, [θ], [α => lower])
        grassman_deriv_bar=PartialDerivative("∂",zero, [θbar], [αd => lower])
        spatial_deriv=PartialDerivative("∇",zero, [x], [μ => lower])

        add_derivative_action(grassman_deriv_bar,Tensor("Θ",[βd=>lower],θbar),EpsilonTensor([βd =>lower, αd=>lower]))
        add_derivative_action(grassman_deriv,Tensor("θ",[β=>lower],θ),EpsilonTensor([β =>lower, α=>lower]))
        add_derivative_action(grassman_deriv_bar,Tensor("Θ",[βd=>upper],θbar),KroneckerDelta([αd =>lower, βd=>upper]))
        add_derivative_action(grassman_deriv,Tensor("θ",[β=>upper],θ),KroneckerDelta([α =>lower, β=>upper]))

        # Am=Tensor("A", [x], [αd=>lower])

        covar_expansion=grassman_deriv_bar+im/2*Tensor("θ",[α=>upper],θ)*Tensor("σ",[μ=>upper,α=>lower,αd=>lower])*spatial_deriv#+Am
        covar_unbar_expansion=grassman_deriv+im/2*Tensor("σ",[μ=>upper,α=>lower,αd=>lower])*Tensor("Θ",[αd=>upper],θbar)*spatial_deriv
        covar=Derivative("̄D",zero, [x,θbar], [ αd => lower])
        covar_unbar=Derivative("D",zero, [x,θ], [ α => lower])


        #test algebra of derivatives
        a1=Derivative("D",Derivative("̄D",Tensor("X", [x,θ,θbar]), [x,θbar], [γd=>lower]), [x,θ], [γ=>lower])
        @show a1
        # @show expand_derivative(a1,covar,covar_expansion)
        # @show expand_derivative(expand_derivative(a1,covar,covar_expansion),covar_unbar,covar_unbar_expansion)
        a1ex=expand_derivative(expand_derivative(expand_derivative(a1,covar,covar_expansion),covar_unbar,covar_unbar_expansion),covar,covar_expansion)
        # @show a1ex
        alex1=Alakazam.simplify(apply_derivative(a1ex))
        @show alex1
        
        a2= Derivative("̄D",Derivative("D",Tensor("X", [x,θ,θbar]), [x,θ], [γ=>lower]), [x,θbar], [γd=>lower])
        @show a2
        a2ex=expand_derivative(expand_derivative(expand_derivative(a2,covar,covar_expansion),covar_unbar,covar_unbar_expansion),covar,covar_expansion)
        @show a2ex
        alex2=apply_derivative(a2ex)
        @show alex2
        # @show Alakazam.hard_simplify(Alakazam.simplify(rename_dummies(sort(alex1+alex2))))
        # @show alex1+alex2
        @show (Alakazam.simplify(alex1+alex2))
        # @show sort(Alakazam.simplify(alex1+alex2))
        # @show rename_dummies(sort(Alakazam.simplify(alex1+alex2)))

        # acomm=Derivative("D",Derivative("̄D",Tensor("X", [x,θ,θbar]), [x,θbar], [γd=>lower]), [x,θ], [γ=>lower])+ Derivative("̄D",Derivative("D",Tensor("X", [x,θ,θbar]), [x,θ], [γ=>lower]), [x,θbar], [γd=>lower])
        # @show acomm
        # #need to add functionality to commute derivative operators...
        # acommexp=expand_derivative(expand_derivative(expand_derivative(acomm,covar,covar_expansion),covar_unbar,covar_unbar_expansion),covar,covar_expansion)
        # @show acommexp
        # acommexp=apply_derivative(acommexp)
        # @show Alakazam.simplify(apply_derivative(acommexp))

        # # X=(Tensor("x", x)+Tensor("θ",[α=>lower],θ)*Tensor("A",[α=>upper],x)+Tensor("Θ",[αd=>lower],θbar)*Tensor("B",[αd=>upper],x)
        # # +Tensor("θ",[α=>lower],θ)*Tensor("Θ",[αd=>lower],θbar)*Tensor("C",[α=>upper, αd=>upper],x)
        # # +Tensor("θ",[α=>lower],θ)*Tensor("θ",[β=>lower],θ)*Tensor("D",[α=>upper, β=>upper],x)
        # # +Tensor("Θ",[αd=>lower],θbar)*Tensor("Θ",[βd=>lower],θbar)*Tensor("E",[αd=>upper, βd=>upper],x)
        # # +Tensor("θ",[α=>lower],θ)*Tensor("θ",[β=>lower],θ)*Tensor("Θ",[αd=>lower],θbar)*Tensor("F",[α=>upper, β=>upper, αd=>upper],x) 
        # # +Tensor("θ",[α=>lower],θ)*Tensor("Θ",[αd=>lower],θbar)*Tensor("Θ",[βd=>lower],θbar)*Tensor("G",[α=>upper, βd=>upper, αd=>upper],x)
        # # +Tensor("θ",[α=>lower],θ)*Tensor("θ",[β=>lower],θ)*Tensor("Θ",[αd=>lower],θbar)*Tensor("Θ",[βd=>lower],θbar)*Tensor("H",[α=>upper, β=>upper, βd=>upper, αd=>upper],x))
        # # X=(Tensor("A", [x,θ,θbar])
        # # +sqrt(2)*Tensor("θ",[α=>lower], θ)*Tensor("ψ",[α=>upper], [x,θ,θbar])
        # # +Tensor("θ",[α=>lower], θ)*Tensor("θ",[α=>upper], θ)*Tensor("F",[x,θ,θbar]))
 
        # X=(Tensor("A", [x])+Tensor("θ",[α=>lower], θ)*Tensor("ψ",[α=>upper], [x])
        # -Tensor("θ",[α=>lower], θ)Tensor("θ",[α=>upper], θ)*Tensor("F",[x])
        # +im/2*Tensor("θ",[α=>upper], θ)Tensor("Θ",[αd=>upper], θbar)*Tensor("σ",[μ=>upper,α=>lower,αd=>lower])*Derivative("∇",Tensor("A", [x]), [x], [μ=>lower])
        # +im/2*Tensor("θ",[β=>upper], θ)Tensor("θ",[β=>lower], θ)*Tensor("Θ",[αd=>upper], θbar)*Tensor("σ",[μ=>upper,α=>lower,αd=>lower])*Derivative("∇",Tensor("ψ",[α=>upper], [x]), [x], [μ=>lower])
        # +1/4*Tensor("θ",[α=>lower], θ)*Tensor("θ",[α=>upper], θ)*Tensor("Θ",[αd=>lower],θbar)*Tensor("Θ",[αd=>upper],θbar)*Derivative("∇",Derivative("∇",Tensor("A", [x]), [x], [μ=>upper]), [x], [μ=>lower]))
 

        # @show X

        # Dbar_SUPERFIELD=Derivative("̄D",X, [x], [γd=>lower])
        # out=expand_derivative(Dbar_SUPERFIELD,covar,covar_expansion)
        # @show Alakazam.simplify(rename_dummies(apply_derivative(out)))
        # # Dchi=Alakazam.drop_weight_above(Alakazam.simplify(apply_derivative(out)),"θ",4)
        # # @show Dchi
    end
  
# b1=PartialDerivative("∇",PartialDerivative("∇",Tensor("X", [x,θ,θbar]), [x], [μ=>lower]), [x], [ν=>lower])
# @show b1
# nestb=Alakazam.to_nested(b1)
# @show nestb



    
end