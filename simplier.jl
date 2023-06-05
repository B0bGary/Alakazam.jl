temp=bt00
println(temp)
#  

# temp=eliminate_epsilon(rename_dummies(eliminate_epsilon(eliminate_kronecker(KroneckerDelta([s=>lower, D=>upper])*temp))))


#    temp=eliminate_epsilon(rename_dummies(eliminate_epsilon((EpsilonTensor([r=>lower, D=>lower])*eliminate_kronecker(temp)))))


println("\\begin{dmath}",to_latex(temp),"=0\\end{dmath}")