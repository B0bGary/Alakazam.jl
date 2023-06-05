
# tt=3
# tb=3

# t33=rename_dummies(get_weighted_terms(get_weighted_terms(res,"θ",tt),"Θ",tb))
# t33=project_component(t33,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],2,2)
# t33=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(t33)))))
# println("\$\$ \\theta=",tt,"\\quad \\overline{\\theta}=",tb,":\$\$")
# println("\\begin{dmath}",to_latex(t33),"=0\\end{dmath}")
# println(t33)


# tt=4
# tb=3

# t43=rename_dummies(get_weighted_terms(get_weighted_terms(res,"θ",tt),"Θ",tb))
# t43=project_component(t43,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],2,2)
# t43=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(t43)))))
# println("\$\$ \\theta=",tt,"\\quad \\overline{\\theta}=",tb,":\$\$")
# println("\\begin{dmath}",to_latex(t43),"=0\\end{dmath}")
# println(t43)


# tt=3
# tb=4

# t34=rename_dummies(get_weighted_terms(get_weighted_terms(res,"θ",tt),"Θ",tb))
# t34=project_component(t34,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],2,2)
# t34=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(t34)))))
# println("\$\$ \\theta=",tt,"\\quad \\overline{\\theta}=",tb,":\$\$")
# println("\\begin{dmath}",to_latex(t34),"=0\\end{dmath}")
# println(t34)


# include("4dconjsleepy.jl")
include("superspace.jl")
sav=res2
res2=res
tt=0
tb=0
bt00=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt00=project_component(bt00,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt00=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt00)))))

tt=1
tb=0
bt10=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt10=project_component(bt10,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt10=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt10)))))

tt=2
tb=0
bt20=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt20=project_component(bt20,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt20=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt20)))))

bt202=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt202=project_component(bt202,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],2,1)
bt202=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt202)))))

tt=3
tb=0
bt30=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt30=project_component(bt30,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt30=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt30)))))

tt=4
tb=0
bt40=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt40=project_component(bt40,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt40=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt40)))))


tt=0
tb=1
bt01=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt01=project_component(bt01,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt01=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt01)))))

tt=1
tb=1
bt11=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt11=project_component(bt11,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt11=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt11)))))


tt=2
tb=1
bt21=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt21=project_component(bt21,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt21=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt21)))))

bt212=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt212=project_component(bt212,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],2,1)
bt212=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt212)))))

tt=3
tb=1
bt31=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt31=project_component(bt31,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt31=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt31)))))

tt=4
tb=1
bt41=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt41=project_component(bt41,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt41=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt41)))))

tt=0
tb=2
bt02=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt02=project_component(bt02,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt02=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt02)))))


bt022=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt022=project_component(bt022,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,2)
bt022=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt022)))))


tt=1
tb=2
bt12=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt12=project_component(bt12,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt12=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt12)))))


bt122=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt122=project_component(bt122,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,2)
bt122=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt122)))))

tt=2
tb=2
bt22=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt22=project_component(bt22,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt22=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt22)))))

bt222=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt222=project_component(bt222,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,2)
bt222=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt222)))))

bt223=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt223=project_component(bt223,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],2,1)
bt223=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt223)))))

bt2224=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt2224=project_component(bt2224,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],2,2)
bt2224=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt2224)))))

tt=3
tb=2
bt32=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt32=project_component(bt32,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt32=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt32)))))

bt322=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt322=project_component(bt322,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,2)
bt322=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt322)))))



tt=4
tb=2
bt42=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt42=project_component(bt42,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt42=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt42)))))

bt422=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt422=project_component(bt422,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,2)
bt422=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt422)))))

tt=0
tb=3
bt03=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt03=project_component(bt03,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt03=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt03)))))

tt=1
tb=3
bt13=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt13=project_component(bt13,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt13=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt13)))))

tt=2
tb=3
bt23=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt23=project_component(bt23,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt23=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt23)))))

bt232=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt232=project_component(bt232,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt232=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt232)))))

tt=3
tb=3
bt33=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt33=project_component(bt33,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt33=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt33)))))

tt=4
tb=3
bt43=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt43=project_component(bt43,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt43=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt43)))))

tt=0
tb=4
bt04=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt04=project_component(bt04,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt04=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt04)))))

tt=1
tb=4
bt14=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt14=project_component(bt14,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt14=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt14)))))

tt=2
tb=4
bt24=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt24=project_component(bt24,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt24=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt24)))))

bt242=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt242=project_component(bt242,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],2,1)
bt242=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt242)))))

tt=3
tb=4
bt34=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt34=project_component(bt34,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt34=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt34)))))

tt=4
tb=4
bt44=rename_dummies(get_weighted_terms(get_weighted_terms(res2,"θ",tt),"Θ",tb))
bt44=project_component(bt44,θ,[i=>upper, α=>lower],tt,tb,θbar,[i=>lower, αdot=>upper],1,1)
bt44=rename_dummies(contract_metrics(simplify(eliminate_epsilon(eliminate_kronecker(bt44)))))