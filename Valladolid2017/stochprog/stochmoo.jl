
# illustration of multi-stakeholder optimization problem
# victor m. zavala, 2017

using JuMP 
using Ipopt 
using Gadfly

# number of objectives
nObj = 2

# number of stakeholders
nStake = 20;

W = rand(nStake,nObj)
for i = 1:size(W,1)
W[i,:] = W[i,:] / sum(W[i,:])
end

# container for solutions
x1v=zeros(nStake)
x2v=zeros(nStake)
cp=zeros(nStake)

for p in 1:nStake

# Model 
m = Model(solver=IpoptSolver(tol = 1e-4, max_iter = 100,linear_solver ="mumps",mu_strategy="monotone"))

@variable(m, 0 <= x1<= 1)    
@variable(m, 0 <= x2<= 1)  
@variable(m,cost)
@NLconstraint(m, cons, x1^2 + x2 == 1 ) 
@NLconstraint(m,costcons,cost== W[p,1]*(x1-1)^2 + W[p,2]*(x2-x1)^2)
@NLobjective(m, Min,cost)
   
# solve model and get solution
solve(m)
    
# collect solution    
    x1v[p]=getvalue(x1);
    x2v[p]=getvalue(x2);
    cp[p]=getvalue(cost)
end

# Compute compromise solution using cvar
P=1:nStake
alpha=0.1

m = Model(solver=IpoptSolver(tol = 1e-4, max_iter = 100,linear_solver ="mumps",mu_strategy="monotone"))

@variable(m, 0 <= x1<= 1)    
@variable(m, 0 <= x2<= 1)  
@variable(m, phi[P]>=0)
@variable(m, cost[P])
@variable(m, VaR)
@NLconstraint(m, cons, x1^2 + x2 == 1 ) 
@NLconstraint(m,costcons[p in P],cost[p]== W[p,1]*(x1-1)^2 + W[p,2]*(x2-x1)^2)
@NLconstraint(m,cvarcons[p in P],cost[p] - VaR <= phi[p])
@NLobjective(m, Min, VaR+(1/(alpha*nStake))*sum{phi[p],p in P})
   
# solve model and get solution
solve(m)


cx1=zeros(1)
cx2=zeros(1)
cx1=getvalue(x1)
cx2=getvalue(x2)

plot(
layer(x=[cx1],y=[cx2], Geom.point, Theme(default_color=color("red"))),
layer(x=x1v,y=x2v, Geom.point, Theme(default_color=color("green")))
)




