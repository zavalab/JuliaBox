using JuMP
using Ipopt

NS=2; 
S=1:NS
rhs=zeros(NS)
pr=zeros(NS)
rhs[1]=500
rhs[2]=700
pr[1]=0.6
pr[2]=0.4

m=Model(solver=IpoptSolver());
@variable(m,x[S]>=0)
@variable(m,y[S]>=0)
@constraint(m,cons[s in S], x[s]+y[s]>=rhs[s])
@constraint(m,nonant,x[1]-x[2]==0)
@objective(m,Min,sum{pr[s]*(2*x[s]+3*y[s]), s in S})
solve(m)   

print(getvalue(x))
print(getobjectivevalue(m))
print(getdual(nonant))

#solve with dualized constraint
lam=getdual(nonant)
m=Model(solver=IpoptSolver());
@variable(m,x[S]>=0)
@variable(m,y[S]>=0)
@constraint(m,cons[s in S], x[s]+y[s]>=rhs[s])
@objective(m,Min,sum{pr[s]*(2*x[s]+3*y[s]), s in S}-lam*(x[1]-x[2]))
solve(m)  

print(getvalue(x))
print(getobjectivevalue(m))
