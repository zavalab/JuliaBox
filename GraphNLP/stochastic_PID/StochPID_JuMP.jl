using JuMP
using MadNLP

# sets
NS = 5;      # Number of Scenarios
N  = 100;    # Number of timesteps
Tf = 10;     # Final time
h  = Tf/N;   # Time step, or Δt
T  = 1:N;    # Set of times
Tm = 1:N-1;  # Set of times minus one
S  = 1:NS;   # Set of scenarios

# set time vector
time = zeros(N);
for t = 1:N
  time[t] = h*(t-1);
end

# define parameters used within each scenario
K    = 1.0;            # gain
x0   = 0.0;            # starting point
Kd   = 0.5;            # disturbance gain
tau  = 1.0             # inverse time constant
d    = fill(-1.0,5);   # disturbance
xsp  = [-2.0, -1.5, -0.5, 0.5, 1.0] # set point

# define JuMP model
model = Model(() -> MadNLP.Optimizer())

# define model variables
@variable(model, -2.5 <= x[S, T] <= 2.5) # process variable
@variable(model, -2.0 <= u[S, T] <= 2.0) # process input
@variable(model, int[S, T]) # discretized integral at time t
@variable(model, cost[S,T]) # cost to objective function

# define first stage variables
@variable(model, -10 <= Kc <= 10, start=0)
@variable(model, -100 <= tauI <= 100, start=0)
@variable(model, -100 <= tauD <= 100, start=0)
    
# constrain model to follow first order linear system
@constraint(model, eqdyn[s in S, t in Tm], (1/tau)*(x[s,t+1]-x[s,t])/h + x[s,t+1] == K*u[s,t+1] + Kd*d[s]); 
# define controller operation
@constraint(model, eqcon[s in S, t in Tm], u[s,t+1] == Kc*(xsp[s] - x[s,t]) + tauI*int[s,t+1] + tauD*(x[s,t+1] - x[s,t])/h);
# discretize integral
@constraint(model, eqint[s in S, t in Tm], (int[s,t+1] - int[s,t])/h == xsp[s] - x[s,t+1]);
# set starting point
@constraint(model, eqinix[s in S], x[s,1] == x0);
# set initial integral term to zero
@constraint(model, eqinit[s in S], int[s,1] == 0);
# define cost for objective function
@constraint(model, eqcost[s in S, t in T], cost[s,t] == 100*(xsp[s]-x[s,t])^2 + 0.01*u[s,t]^2);

# minimize the cost; each scenario is equally weighted
@objective(model, Min, (1/(NS)) * sum(cost[s,t] for s in S, t in T))

optimize!(model)
