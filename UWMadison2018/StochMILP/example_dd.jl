using JuMP
using Ipopt

NS=2
S=1:NS
rhs=zeros(NS)
pr=zeros(NS)
rhs[1]=500
rhs[2]=700
pr[1]=0.6
pr[2]=0.4

# dual subproblems
function solve_subproblem(s,rhs,pr,lambda) 
    m=Model(solver=IpoptSolver(print_level=0));
    @variable(m,x>=0)
    @variable(m,y>=0)
    @constraint(m,cons, x+y>=rhs[s])
    @objective(m,Min,pr[s]*(2*x+3*y) + lambda[s]*x)
    solve(m)   
    return getvalue(x),getobjectivevalue(m)
end

# compute upper bound
function solve_upperbound(x,rhs,pr,S)
    m=Model(solver=IpoptSolver(print_level=0));
    @variable(m,y[S]>=0)
    @constraint(m,cons[s in S], x+y[s]>=rhs[s])
    @objective(m,Min,sum{pr[s]*(2*x+3*y[s]), s in S})
    solve(m)  
    return(getobjectivevalue(m))
end

Nit=4       # number of iterations
lam=0       # initial dual variable
x=zeros(NS) # container for first-stage
D=zeros(NS) # container for dual function values
lambda=zeros(NS) # constainer for duals
alpha=0.001 # step-size duals
zUB=0;      # Upper bound 

for it=1:Nit
    
    # get lower bound
    lambda[1]=-lam
    lambda[2]=+lam
    for s=1:NS
        (x[s],D[s])=solve_subproblem(s,rhs,pr,lambda)    
    end
    println('\n')
    @printf("lambda=%e\n", lam)    
    println(x)
    println(D)
    err=x[1]-x[2]
    @printf("Error=%e\n", err)
    zLB=sum(D)
    @printf("zLB=%e\n", zLB)
    
    # get upper bound for every x[s]
    for s=1:NS
        zUB=solve_upperbound(x[s],rhs,pr,S)    
        gap=zUB-zLB
        @printf("zUB=%e\n", zUB)
        @printf("Gap=%e\n", gap)
    end
    
    # update duals
    lam=lam-alpha*(err)
end