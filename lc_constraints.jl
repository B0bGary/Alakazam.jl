
t=1
tb=2

t11=get_weighted_terms(get_weighted_terms(res,"θ",t),"Θ",tb)


println("\$\$ \\theta=",t,"\\quad \\overline{\\theta}=",tb,":\$\$")
println("\$\$",to_latex(t11),"=0\$\$")
println(t11)

