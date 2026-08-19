# const SymorNum = Union{SymbolicUtils.Symbolic,Number}
# const SymorNum = Union{
#     SymbolicUtils.BasicSymbolic{SymbolicUtils.SafeReal},
#     SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal},
#     SymbolicUtils.BasicSymbolic{SymbolicUtils.TreeReal},
#     Number}
const SymorNum = Union{
    SymbolicUtils.BasicSymbolic{SymbolicUtils.SymReal},
    Number}

abstract type IndexSuperType end
#Supertype of tensor, term, expression etc to allow common functions
abstract type Tensorial end
abstract type TensorSuperType <: Tensorial end
#basic data structure. in case we want to extend later eg for ops in some way
abstract type TensorDataSuperType end


# keep these, since we need to check if metrics in some checks
abstract type MetricSuperType <: TensorSuperType end
abstract type SymMetricSuperType <: MetricSuperType end
abstract type AntisymMetricSuperType <: MetricSuperType end
abstract type HermMetricSuperType <: MetricSuperType end #hermitian metric has different indices
abstract type EpsilonTensorSuperType <: AntisymMetricSuperType end
abstract type KroneckerSuperType <: TensorSuperType end #kronecker is treated not as a metric, since different level indices


abstract type  MatrixSuperType<: TensorSuperType end
abstract type PauliSuperType <: MatrixSuperType end #new supertype for pauli matrices, since they have special properties that we want to be able to check for in simplification rules
abstract type GammaSuperType <: MatrixSuperType end #new supertype for gamma matrices, since they have special properties that we want to be able to check for in simplification rules
#made this a seperate type just to be safe
abstract type Gamma5SuperType <: MatrixSuperType end #new supertype for gamma matrices, since they have special properties that we want to be able to check for in simplification rules

abstract type DiracVectorSuperType <:MatrixSuperType end # bras and kets, don't want commute since inner prod and outer prod are different

abstract type CliffordTensorSuperType <: TensorSuperType end
abstract type GammaTensorSuperType <: CliffordTensorSuperType end
abstract type PauliTensorSuperType <: TensorSuperType end
abstract type BraSuperType <:DiracVectorSuperType end 
abstract type KetSuperType <:DiracVectorSuperType end 



abstract type OperatorSuperType <: TensorSuperType end

abstract type DifferentiableFunctionSuperType <:OperatorSuperType end #differentiable functions? for chain rule?


abstract type LeibnizOperatorSuperType <: OperatorSuperType end #things like derivatives that have the (graded) Leibniz derivation property.

abstract type DerivativeSuperType <: LeibnizOperatorSuperType end
abstract type PartialSuperType <: DerivativeSuperType end
abstract type CovDerivSuperType <: DerivativeSuperType end
abstract type GeneratorDerivSuperType <: DerivativeSuperType end

abstract type VariationalSuperType <: PartialSuperType end #we will see if partial is the best way...

abstract type LieBracketSuperType <:LeibnizOperatorSuperType end #lie brackets have a derivation property


abstract type HomomorphismSuperType <: OperatorSuperType end #these won't obey Leibniz, rather O(AB)=O(A)O(B) (eg like the det)
abstract type AntimorphismSuperType <: OperatorSuperType end #these won't obey Leibniz, rather O(AB)=O(B)O(A) (herm conj, transpose)
abstract type InvolutionSuperType <: AntimorphismSuperType end #O(O(A))=A
abstract type DaggerSuperType <: InvolutionSuperType end #incase we need dagger specific features...

abstract type CyclicSuperType <: OperatorSuperType end # O(ABC)=O(BCA)
abstract type TraceSuperType <: CyclicSuperType end # O(ABC)=O(BCA)


abstract type  LieAlgebraSuperType end
abstract type NumberSet end