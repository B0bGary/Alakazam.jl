@testset "Thesis code Test" begin
   
    

#thesis specific code to get constraints
function strip_thetas(TEx::TensorExpression)
    used_inds=union(get_dummies(TEx),get_free_indices(TEx))
    unused_flav=Vector(setdiff(flavour_indices.indices,used_inds))
    unused_spin=Vector(setdiff(spinor_indices.indices,used_inds))
    unused_conj=Vector(setdiff(spinor_conj_indices.indices,used_inds))
    new_exp=zero
    for summand in TEx.summands
        new_term=summand[1]
        tcounter=1
        tbcounter=1
        for term in summand[2].terms
            if term.name=="θ"
                # println(term)
                for indexes in eachindex(term.indices)
                    if term.indices[indexes][2]==upper
                        if term.indices[indexes][1].index_set==flavour_indices
                            dl=KroneckerDelta([term.indices[indexes], unused_flav[2*tcounter-1]=>lower])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        else
                            dl=KroneckerDelta([term.indices[indexes], unused_spin[tcounter]=>lower])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        end
                    else
                        if term.indices[indexes][1].index_set==flavour_indices
                            dl=KroneckerDelta([term.indices[indexes], unused_flav[2*tcounter-1]=>upper])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        else
                            dl=KroneckerDelta([term.indices[indexes], unused_spin[tcounter]=>upper])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        end
                    end
                    
                end
                # println(new_term)
                tcounter+=1
            elseif term.name=="Θ"
                # println(term)
                for indexes in eachindex(term.indices)
                    if term.indices[indexes][2]==upper
                        if term.indices[indexes][1].index_set==flavour_indices
                            dl=KroneckerDelta([term.indices[indexes], unused_flav[2*tbcounter]=>lower])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        else
                            dl=KroneckerDelta([term.indices[indexes], unused_conj[tbcounter]=>lower])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        end
                    else
                        if term.indices[indexes][1].index_set==flavour_indices
                            dl=KroneckerDelta([term.indices[indexes], unused_flav[2*tbcounter]=>upper])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        else
                            dl=KroneckerDelta([term.indices[indexes], unused_conj[tbcounter]=>upper])
                            #println("replacing ",term," with ",dl )
                            new_term*=dl
                        end
                    end
                    
                end
                tbcounter+=1
                # println(new_term)
            else
                new_term*=term
            end
        end
        new_exp+=new_term
    end
    return new_exp
end

function specify_indices(TEx::TensorExpression, reps::Dict{Index,Index})
    new_exp=zero
    w=copy(TEx)
    for summand in w.summands
        term=1
        for tens in summand[2].terms
            if typeof(tens)<:OperatorSuperType
                tens.operand=specify_indices(tens.operand,reps)
                tens.operand_indices=tens.operand.summands[1][2].indices
                # println(reps)
                # test=tens.operator_indices[1][1]
                # println(test )
                # println(haskey(reps,test))
                tens.operator_indices=[ haskey(reps,ind[1]) ? reps[ind[1]]=>ind[2] : ind[1]=>ind[2]  for ind in tens.operator_indices]
                # println("op should inclde", tens.operator_indices)
                tens.indices=vcat(tens.operator_indices,tens.operand_indices)
                # println("op has", tens.indices)
            else
                if !isempty(tens.indices)
                tens.indices=[ haskey(reps,ind[1]) ? reps[ind[1]]=>ind[2] : ind[1]=>ind[2]  for ind in tens.indices]
                end
            end
            if !(typeof(term)<:SymorNum) && ((typeof(term)<:TensorTerm && tens in term.terms) || (typeof(term)<:TensorExpression && tens in term.summands[1][2].terms) )
                term=zero
                break
            else
                term=term*tens
            end
        end
        new_exp+=term
    end

    return new_exp
end

function to_lc_expression(TEx::TensorExpression, reps::Dict{Index,Index})
    #replace non-dummies
    exp=specify_indices(TEx, reps)
    # println(exp)
    #cookery for dummies....
    new_exp=zero

    for summand in exp.summands
        dums=intersect(get_dummies(summand[2]),Set(spacetime_indices.indices))
        tempterm=TensorExpression(summand[2])
        # println("reset",tempterm)
        while !isempty(dums)
            dum=pop!(dums)
            tempterm=specify_indices(tempterm, Dict(dum => p))+specify_indices(tempterm, Dict(dum => m))
            # println(tempterm)
        end

        new_exp+=summand[1]*tempterm
    end
    res_ex=zero
    for summand in new_exp.summands
        flag=0
            for tens in summand[2].terms
                if typeof(tens)<:OperatorSuperType
                    if tens.operator_indices[1][1].name!=tens.operator_indices[2][1].name
                        flag=1
                     #   println("ignoring term", summand[2])
                        break
                    end
                end
            end
            if flag==0 &&term_weight(summand[2], "θ")<=2 && term_weight(summand[2], "Θ")<=2
                res_ex+=summand[1]*summand[2]
            end
    end

    # end
    out_exp=zero
    flag=0
    for summand in new_exp.summands
        for tens in summand[2].terms
            if typeof(tens)<:OperatorSuperType
                if tens.operator_indices[1][1].name!=tens.operator_indices[2][1].name
                    
                    flag=1
                    break
                end
            end
        end
        if flag==0
            # println("got here, ",summand[2] )
        out_exp+=summand[1]*summand[2]
        end
        flag=0
    end



    return out_exp
end

function to_chiral_part(TEx::TensorExpression, reps::Index)
    new_exp=zero

    for summand in TEx.summands
        term=1

        for tens in summand[2].terms
            if tens.name=="θ" || tens.name=="Θ"
                if tens.indices[1][1].name!=reps.name
                    term=zero
                    break
                end
            end
            term*=tens
        end

        new_exp+=term
    end

    return new_exp
end

#define index sets

spacetime_indices = IndexSet("spacetimelc", 2,lower)
flavour_indices = IndexSet("flavour_indices", 2,upper)
lightcone_indices = IndexSet("lightcone", 2,lower)
spinor_indices=IndexSet("spinor_indices", 2,lower)
spinor_conj_indices=IndexSet("spinor_conj_indices", 2,lower)

p = Index("+", lightcone_indices, 1)
m = Index("-", lightcone_indices, 1)

a = Index("a", flavour_indices, 1)
b = Index("b", flavour_indices,1)
c = Index("c", flavour_indices,1)
d = Index("d", flavour_indices,1)
e = Index("e", flavour_indices,1)

i = Index("i", flavour_indices,1)
j = Index("j", flavour_indices,1)
k = Index("k", flavour_indices,1)
l = Index("l", flavour_indices,1)
q = Index("q", flavour_indices,1)
n = Index("n", flavour_indices,1)
r = Index("r", flavour_indices,1)
s = Index("s", flavour_indices,1)
t = Index("t", flavour_indices,1)
u = Index("u", flavour_indices,1)
v = Index("v", flavour_indices,1)
w = Index("w", flavour_indices,1)
y = Index("y", flavour_indices,1)
z = Index("z", flavour_indices,1)

x = Coordinate("x", [spacetime_indices, spacetime_indices])
θ = Coordinate("θ", [spacetime_indices])
θbar = Coordinate("Θ", [spacetime_indices])

grassman_deriv=Derivative("𝛛",zero, [θ], [i => lower])
grassman_deriv_bar=Derivative("∂",zero, [θbar], [i => lower])
spatial_deriv=Derivative("∇",zero, [x])
minus_deriv=Derivative("δ",zero, [x])


add_derivative_action(grassman_deriv_bar,Tensor("Θ",[j=>lower],θbar),EpsilonTensor([j =>lower, i=>lower]))
add_derivative_action(grassman_deriv,Tensor("θ",[j=>lower],θ),EpsilonTensor([j =>lower, i=>lower]))
add_derivative_action(grassman_deriv_bar,Tensor("Θ",[j=>upper],θbar),KroneckerDelta([i =>lower, j=>upper]))
add_derivative_action(grassman_deriv,Tensor("θ",[j=>upper],θ),KroneckerDelta([i =>lower, j=>upper]))


covar_expansion=-1*grassman_deriv_bar+im/2*Tensor("θ",[i=>lower],θ)*spatial_deriv
covar_unbar_expansion=grassman_deriv-im/2*Tensor("Θ",[ i=>lower],θbar)*spatial_deriv
covar=Derivative("̄D",zero, [x,θbar], [ i => lower])
covar_unbar=Derivative("D",zero, [x,θ], [ i => lower])

custom_sort_order(["θ","Θ","ϵ"])

SUPERFIELD_S=(Tensor("Θ", x)
+Tensor("θ",[ i=> lower],θ)*Tensor("A",[ i=>upper],x)
-Tensor("Θ",[ i=> lower],θbar)*Tensor("̄A",[ i=>upper],x)
+Tensor("θ",[ i=> lower],θ)*Tensor("θ",[ i=> upper],θ)*Tensor("B",x)
-Tensor("Θ",[ i=> lower],θbar)*Tensor("Θ",[i=> upper],θbar)*Tensor("̄B",x)
+Tensor("θ",[ i=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("C",[  i=>upper,j=>upper],x)
+Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ i=> lower],θbar)*Tensor("S",[   i=>upper],x)
+Tensor("θ",[ i=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[ j=> upper],θbar)*Tensor("̄S",[   i=>upper],x)
+Tensor("θ",[ i=> lower],θ)*Tensor("θ",[ i=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("Θ",[ k=> upper],θbar)*Tensor("T",x)
)

SUPERFIELD_Tmmmm=(Tensor("T",[], x)
+Tensor("θ",[ q=> lower],θ)*Tensor("X",[ q=>upper],x)
-Tensor("Θ",[ q=> lower],θbar)*Tensor("̄X",[ q=>upper],x)
+Tensor("θ",[ q=> lower],θ)*Tensor("θ",[ q=> upper],θ)*Tensor("U",[],x)
-Tensor("Θ",[ q=> lower],θbar)*Tensor("Θ",[ q=> upper],θbar)*Tensor("̄U",[ ],x)
+Tensor("θ",[ q=> lower],θ)*Tensor("Θ",[ n=> lower],θbar)*Tensor("V",[  q=>upper,n=>upper],x)
+Tensor("θ",[ q=> lower],θ)*Tensor("θ",[ q=> upper],θ)*Tensor("Θ",[ r=> lower],θbar)*Tensor("S",[  r=>upper],x)
+Tensor("Θ",[ q=> lower],θbar)*Tensor("Θ",[ q=> upper],θbar)*Tensor("θ",[ r=> lower],θ)*Tensor("̄S",[  r=>upper],x)
+Tensor("θ",[ q=> lower],θ)*Tensor("θ",[ q=> upper],θ)*Tensor("Θ",[ r=> lower],θbar)*Tensor("Θ",[ r=> upper],θbar)*Tensor("Θ",[ ],x)
)

SUPERFIELD_Wmmm=(Tensor("Γ",[i=>upper], x)
+Tensor("θ",[ j=> lower],θ)*Tensor("Ω",[ i=>upper,j=>upper],x)
+Tensor("Θ",[ j=> lower],θbar)*Tensor("χ",[ i=>upper,j=>upper],x)
+Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Q",[ i=>upper],x)
+Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[j=> upper],θbar)*Tensor("R",[ i=>upper],x)
+Tensor("θ",[ k=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("Y",[  i=>upper,k=>upper,j=>upper],x)
+Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("M",[i=>upper, k=>upper],x)
+Tensor("θ",[ k=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[ j=> upper],θbar)*Tensor("N",[   i=>upper, k=>upper],x)
+Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("Θ",[ k=> upper],θbar)*Tensor("Z",[ i=> upper],x)
)



SUPERFIELD_barWmmm=(Tensor("̄Γ",[i=>upper], x)
+Tensor("θ",[ j=> lower],θ)*Tensor("̄χ",[ i=>upper,j=>upper],x)
+Tensor("Θ",[ j=> lower],θbar)*Tensor("̄Ω",[ i=>upper,j=>upper],x)
-Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("̄R",[ i=>upper],x)
-Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[j=> upper],θbar)*Tensor("̄Q",[ i=>upper],x)
+Tensor("θ",[ k=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("̄Y",[  i=>upper,j=>upper,k=>upper],x)
-Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("̄N",[i=>upper, k=>upper],x)
-Tensor("θ",[ k=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[ j=> upper],θbar)*Tensor("̄x",[   i=>upper, k=>upper],x)
+Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("Θ",[ k=> upper],θbar)*Tensor("̄Z",[ i=> upper],x)
)


SUPERFIELD_Wm=(Tensor("P",[i=>lower], x)
-im/2*Tensor("θ",[ j=> lower],θ)*Tensor("Θ",[ i=>lower,j=>upper],x)
+im/2*Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[ i=>lower,j=>upper],x)
+Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("G",[ i=>lower],x)
+Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[j=> upper],θbar)*Tensor("H",[ i=>lower],x)
+Tensor("θ",[ k=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("I",[  i=>lower,k=>upper,j=>upper],x)
+Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("J",[i=>lower, k=>upper],x)
+Tensor("θ",[ k=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[ j=> upper],θbar)*Tensor("K",[   i=>lower, k=>upper],x)
+2*im*Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("Θ",[ k=> upper],θbar)*Tensor("̄S",[ i=>lower],x)
)


SUPERFIELD_Wmbar=(Tensor("̄P",[i=>lower], x)
-im/2*Tensor("θ",[ j=> lower],θ)*Tensor("Θ",[ i=>lower,j=>upper],x)
+im/2*Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[ i=>lower,j=>upper],x)
-Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("̄H",[ i=>lower],x)
-Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[j=> upper],θbar)*Tensor("̄G",[ i=>lower],x)
+Tensor("θ",[ k=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("̄I",[  i=>lower,k=>upper,j=>upper],x)
-Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("̄K",[i=>lower, k=>upper],x)
-Tensor("θ",[ k=> lower],θ)*Tensor("Θ",[ j=> lower],θbar)*Tensor("Θ",[ j=> upper],θbar)*Tensor("̄J",[   i=>lower, k=>upper],x)
-2*im*Tensor("θ",[ j=> lower],θ)*Tensor("θ",[ j=> upper],θ)*Tensor("Θ",[ k=> lower],θbar)*Tensor("Θ",[ k=> upper],θbar)*Tensor("S",[ i=>lower],x)
)




compute_sc2=false
compute_constraints=true
constraint_type=0
# print(SUPERFIELD_Wmmm)

if compute_sc2
    sc2=SUPERFIELD_S*SUPERFIELD_Tmmmm-(SUPERFIELD_Wmbar-SUPERFIELD_Wm)*rename_dummies(SUPERFIELD_barWmmm-SUPERFIELD_Wmmm)
    sc2=drop_weight_above(drop_weight_above(sc2,"θ",3),"Θ",3)

    D1_SUPERFIELD=Derivative("̄D",sc2, [x], [ t=>lower])
    println(D1_SUPERFIELD)
    out=expand_derivative(D1_SUPERFIELD,covar,covar_expansion)
    # println(out)
    out=apply_derivative(out)
    println(out)
    out=simplify(out)
    out=drop_weight_above(drop_weight_above(out,"θ",3),"Θ",2)
    println("first deriv done")
    println(out)

    D2_SUPERFIELD=Derivative("̄D",out, [x], [ u=>lower])
    out=expand_derivative(D2_SUPERFIELD,covar,covar_expansion)
    out=apply_derivative(out)
    println(out)
    out=simplify(out)
    out=drop_weight_above(drop_weight_above(out,"θ",3),"Θ",1)

    println("2nd deriv done")
    D3_SUPERFIELD=Derivative("D",out, [x], [ v=>lower])
    out=expand_derivative(D3_SUPERFIELD,covar_unbar,covar_unbar_expansion)
    out=apply_derivative(out)
    println(out)
    out=simplify(out)
    out=drop_weight_above(drop_weight_above(out,"θ",2),"Θ",1)

    println("3rd deriv done")

    D4_SUPERFIELD=Derivative("D",out, [x], [ w=>lower])
    out=expand_derivative(D4_SUPERFIELD,covar_unbar,covar_unbar_expansion)
    out=apply_derivative(out)
    println(out)
    out=simplify(out)
    out=drop_weight_above(drop_weight_above(out,"θ",1),"Θ",1)

    println("4th deriv done")
    println(out)

    println("begin simplify")
    out=simplify(out)
    println("simplified")
    println(out)
    ttbar=get_weighted_terms(get_weighted_terms(out,"θ",0),"Θ",0)

    final_contracted=simplify(1/16*(EpsilonTensor([w =>upper, v=>upper])*(EpsilonTensor([z =>upper, y=>upper])*ttbar)))

    println("2,2 latex")
    println(to_latex(final_contracted))
end


if compute_constraints
    consr=Dict()

    if constraint_type==0
        Dbar_SUPERFIELD=Derivative("̄D",SUPERFIELD_barWmmm, [x], [ i=>lower])
        out=expand_derivative(Dbar_SUPERFIELD,covar,covar_expansion)
        out=apply_derivative(out)
        D_SUPERFIELD=Derivative("D",SUPERFIELD_Wmmm, [x], [ i=>lower])
        out2=expand_derivative(D_SUPERFIELD,covar_unbar,covar_unbar_expansion)
        out2=apply_derivative(out2)
        # println(out)
        summ=simplify(-1*out+out2) #equal to del_--S
        # summ=simplify(KroneckerDelta([z=>upper,i=>lower])*summ) #n olonger needed ? i think
        # summ=drop_weight_above(drop_weight_above(summ,"θ",3),"Θ",3)
        println(to_latex(summ))
        println(get_free_indices(summ))
        constrai=summ-apply_derivative(Derivative("δ",SUPERFIELD_S, [x])) #sets the full equation =0.

        for _i in 0:2
            for _j in 0:2
                term=get_weighted_terms(get_weighted_terms(constrai,"θ",_i),"Θ",_j)

                # grassman_deriv=Derivative("𝛛",zero, [θ], [i => lower])
                # grassman_deriv_bar=Derivative("∂",zero, [θbar], [i => lower])
                # grassman_deriv_bar.operand
                indrs=[u,v,w,y]

            ##  println(term)
                for dj=1:_j
                    term=apply_derivative(Derivative("∂",term, [θbar], [indrs[dj] => lower]))
                #  println(term)
                end
                for di=1:_i
                    term=apply_derivative(Derivative("𝛛",term, [θ], [indrs[2+di] => lower]))
                end
                term=contract_metrics(1*term,EpsilonTensor)
                term=simplify(term)
                println("constraint:"*string(_i)*","*string(_j)*"")
                println("\\begin{dmath}")
                println(to_latex(term)*"=0")
                println("\\end{dmath}")
                consr[(_i,_j)]=term
                
            end

        end


        println(consr[2,2])

    elseif constraint_type==1

        Dbar_SUPERFIELD=Derivative("̄D",SUPERFIELD_Tmmmm, [x], [ i=>lower])
        out=expand_derivative(Dbar_SUPERFIELD,covar,covar_expansion)
        summ=apply_derivative(out)
        # D_SUPERFIELD=Derivative("D",SUPERFIELD_Wmmm, [x], [ i=>lower])
        # out2=expand_derivative(D_SUPERFIELD,covar_unbar,covar_unbar_expansion)
        # out2=apply_derivative(out2)
        # println(out)
        # summ=simplify(-1*out+out2) #equal to del_--S
        summ=simplify(KroneckerDelta([z=>upper,i=>lower])*summ)
        # summ=drop_weight_above(drop_weight_above(summ,"θ",3),"Θ",3)
        # print(to_latex(summ))
        constrai=summ-apply_derivative(Derivative("δ",SUPERFIELD_Wm, [x])) #sets the full equation =0.
    
        for _i in 0:2
            for _j in 0:2
                term=get_weighted_terms(get_weighted_terms(constrai,"θ",_i),"Θ",_j)
    
                # grassman_deriv=Derivative("𝛛",zero, [θ], [i => lower])
                # grassman_deriv_bar=Derivative("∂",zero, [θbar], [i => lower])
                # grassman_deriv_bar.operand
                indrs=[u,v,w,y]
    
              ##  println(term)
                for dj=1:_j
                    term=apply_derivative(Derivative("∂",term, [θbar], [indrs[dj] => lower]))
                  #  println(term)
                end
                for di=1:_i
                    term=apply_derivative(Derivative("𝛛",term, [θ], [indrs[2+di] => lower]))
                end
                term=contract_metrics(1*term,EpsilonTensor)
                term=simplify(term)
                println("constraint:"*string(_i)*","*string(_j)*"")
                println("\\begin{dmath}")
                println(to_latex(term)*"=0")
                println("\\end{dmath}")
    
                
            end
    
        end
    
    elseif constraint_type==2


    end



end



end