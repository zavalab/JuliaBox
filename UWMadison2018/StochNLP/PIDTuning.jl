# call libraries
using Ipopt
using JuMP

# sets
 N=100;
Tf=10;
 h=Tf/N;
 T=1:N;
Tm=1:N-1;

# set time vector
time=zeros(N);
for t=1:N
 time[t] = h*(t-1);
end

# data
 x0 =  0.0;
  K =  1.0;
 Kd =  0.5;
tau =  1.0;
xsp = -1.0;
  d = -1.0;

# Create JuMP model object
m = Model(solver=IpoptSolver())

# variables (states and inputs)
@variable(m,-2.5<=x[T]<=2.5)
@variable(m,-2.0<=u[T]<=2.0)
@variable(m, int[T])
@variable(m,cost[T])

# variables (controller design)
@variable(m, -10<= Kc <=10)
@variable(m,-100<=tauI<=100)
@variable(m,-100<=tauD<=1000)

# constraints
@constraint(m, eqdyn[t in Tm],(1/tau)*(x[t+1]-x[t])/h + x[t+1]== K*u[t+1]+Kd*d);
@constraint(m, eqcon[t in Tm], u[t+1] == Kc*(xsp-x[t])+ tauI*int[t+1] + tauD*(x[t+1]-x[t])/h);
@constraint(m, eqint[t in Tm], (int[t+1]-int[t])/h == (xsp-x[t+1]));
@constraint(m, eqinix,   x[1] == x0);
@constraint(m, eqinit, int[1] ==  0);
@constraint(m, eqcost[t in T], cost[t]==(10*(x[t]-xsp)^2 + 0.01*u[t]^2));

# objective function
@objective(m, Min, (1/N)*sum(cost[t] for t in T));

# solve problem
solve(m)

# display results
println(getvalue(Kc))
println(getvalue(tauI))
println(getvalue(tauD))
println(getvalue(cost))
