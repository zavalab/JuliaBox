
# optimal PID controller tuning in JuMP
# Victor M. Zavala
# UW-Madison, 2017

# call libraries
using Ipopt
using JuMP
using PyPlot

# sets
NS=3;
 N=100;
Tf=10;
 h=Tf/N;
 T=1:N;
Tm=1:N-1;
 S=1:NS;

# set time vector
time=zeros(N);
for t=1:N
 time[t] = h*(t-1);
end

# scenario data
K=zeros(NS);   # system gain
x0=zeros(NS);  # initial state
Kd=zeros(NS);  # disturbance gain
tau=zeros(NS); # time constaint
xsp=zeros(NS); # set-point
d=zeros(NS);   # disturbance 

  K[1] =  1.0;
 x0[1] =  0.0;
 Kd[1] =  0.5;
tau[1] =  1.0;
xsp[1] = -1.0;
  d[1] = -1.0;

  K[2] =  1.0;
 x0[2] =  0.0;
 Kd[2] =  0.5;
tau[2] =  1.0;
xsp[2] = -2.0;
  d[2] = -1.0;

  K[3] =  1.0;
 x0[3] =  0.0;
 Kd[3] =  0.5;
tau[3] =  1.0;
xsp[3] =  1.0;
  d[3] = -1.0;

# Create JuMP model object
m = Model(solver=IpoptSolver())

# variables (states and inputs)
@variable(m,-2.5<=x[T,S]<=2.5)
@variable(m,-2.0<=u[T,S]<=2.0)
@variable(m, int[T,S])
@variable(m,cost[T,S])

# variables (controller design) 
@variable(m, -10<= Kc <=10)
@variable(m,-100<=tauI<=100)
@variable(m,-100<=tauD<=1000)

# constraints
@constraint(m, eqdyn[t in Tm, s in S],(1/tau[s])*(x[t+1,s]-x[t,s])/h + x[t+1,s]== K[s]*u[t+1,s]+Kd[s]*d[s]);
@constraint(m, eqcon[t in Tm, s in S], u[t+1,s] == Kc*(xsp[s]-x[t,s])+ tauI*int[t+1,s] + tauD*(x[t+1,s]-x[t,s])/h);
@constraint(m, eqint[t in Tm, s in S], (int[t+1,s]-int[t,s])/h == (xsp[s]-x[t+1,s]));
@constraint(m, eqinix[s in S],   x[1,s] == x0[s]);
@constraint(m, eqinit[s in S], int[1,s] ==  0);
@constraint(m, eqcost[t in T, s in S], cost[t,s]==(10*(x[t,s]-xsp[s])^2 + 0.01*u[t,s]^2));

# objective function
@objective(m, Min, (1/(N*NS))*sum(cost[t,s] for t in T,s in S));

# solve problem 
solve(m)

# display results
println(getvalue(Kc))
println(getvalue(tauI))
println(getvalue(tauD))

# plot responses
x=zeros(NS,N)
for s in 1:NS
    for j=1:N
       x[s,j]=getvalue(getindex(m,:x)[j,s]) 
    end
end

plot(T, x[1,:]);
plot(T, x[2,:]);
plot(T, x[3,:]);
xlabel("Time")
ylabel("x(t)")
grid("on")


