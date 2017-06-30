function get_scenario_model(s)

m=Model()

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
@constraint(m, eqdyn[t in Tm],(1/tau[s])*(x[t+1]-x[t])/h + x[t+1]== K[s]*u[t+1]+Kd[s]*d[s]);
@constraint(m, eqcon[t in Tm], u[t+1] == Kc*(xsp[s]-x[t])+ tauI*int[t+1] + tauD*(x[t+1]-x[t])/h);
@constraint(m, eqint[t in Tm], (int[t+1]-int[t])/h == (xsp[s]-x[t+1]));
@constraint(m, eqinix,   x[1] == x0[s]);
@constraint(m, eqinit, int[1] ==  0);
@constraint(m, eqcost[t in T], cost[t]==(10*(x[t]-xsp[s])^2 + 0.01*u[t]^2));

# objective function
@objective(m, Min, (1/(N*NS))*sum(cost[t] for t in T));

return m

end
