
using JuMP
using Ipopt
using PyPlot

# physical property data
eps= 0.025		             # pipe rugosity - [mm]
z= 0.80        			     # gas compressibility  - []
rhon=0.72         		     # density of air at normal conditions - [kg/m3]
R=8314.0       			     # universal gas constant [J/kgmol-K]
M=18.0    			         # gas molar mass [kg/kgmol]
pi=3.14         		     # pi
T=293.15      		         # reference temperature [K]
Cp=2.34        		         # heat capacity @ constant pressure [kJ/kg-K]
Cv=1.85        		         # heat capacity @ constant volume [kJ/kg-K]

#scaling factors
ffac=(1e+6*rhon)/(24*3600)                     # from scmx10-6/day to kg/s
ffac2=(3600)/(1e+4*rhon)                       # from kg/s to scmx10-4/hr
pfac=1e+5                                      # from bar to Pa
pfac2=1e-5                                     # from Pa to bar
dfac=1e-3                                      # from mm to m
lfac=1e+3;                                     # from km to m

#parameters
diameter = 0.92
pipe_length = 60000
area = (1/4)*pi*diameter*diameter
lam = (2*log10(3.7*diameter/(eps*dfac)))^(-2); 

#other calculated constants
gam = Cp/Cv       		     	# expansion coefficient [-]
nu2 = gam*z*R*T/M  			    # gas speed of sound
om = (gam-1.0)/gam 		     	# aux constant [-]
c4 = (1/ffac2)*(Cp*T)
#lumped constants for pipeline
c1 = (pfac2/ffac2)*(nu2/area)
c2 = area*(ffac2/pfac2)
c3 = area*(pfac2/ffac2)*(8*lam*nu2)/(pi*pi*diameter^5);

horizon = 24*3600;                       # horizon length
time_grid = 1:24;                        # set of temporal mesh points
x_grid = 1:10;                           # set of spatial mesh 
dx = pipe_length / (length(x_grid) - 1)  # equally spaced spatial grid points
dt = horizon / length(time_grid);        # equally spaced time grid points
NS = 3 
S = 1:NS;                                # scenario set

gas_demands = zeros(3,24)
for s in S
    gas_demands[s,:] = fill(150,length(time_grid))
end
gas_demands[1,10:15] = 200
gas_demands[2,10:15] = 500
gas_demands[3,10:15] = 600

for s in S
    step(time_grid,gas_demands[s,:],color="blue")
end
ax = gca()
ax[:set_xlim]([1,24])
xlabel("Time [h]")
ylabel("Demand requested [scm/h *10^4]")
grid("on")

m = Model(solver = IpoptSolver())
#Define inlet, outlet, and flow variables
@variable(m,pin[S,time_grid] >= 10, start = 60)
@variable(m,pout[S,time_grid] >=10, start = 60)
@variable(m,fin[S,time_grid] >=0, start = 100)
@variable(m,fout[S,time_grid] >= 0, start = 100)
@variable(m, px[S,time_grid,x_grid] >= 0, start = 60)
@variable(m, fx[S,time_grid,x_grid] >= 0, start = 10)
# Couple flow convenience variables
@constraint(m, flow_in[s = S,t = time_grid],  fx[s,t,1] == fin[s,t])
@constraint(m, flow_out[s = S,t = time_grid], fx[s,t,x_grid[end]] == fout[s,t]);

@variable(m,ddeliver[S,time_grid] >= 0, start = 100)
@variable(m, dtarget[S,time_grid] >= 0)
@constraint(m, gasLimit[s=S, t = time_grid], ddeliver[s,t] <= dtarget[s,t]);
@constraint(m,deliver[s=S, t=time_grid],fout[s,t] == ddeliver[s,t])
@constraint(m, demand[s=S, t = time_grid], dtarget[s,t] == gas_demands[s,t]);

@variable(m, dp[S,time_grid] >= 0 , start = 10)
@variable(m, 0<=pow[S,time_grid] <= 10000, start = 500);
@NLconstraint(m, powereqn[s=S,t = time_grid], pow[s,t] == c4*fin[s,t]*(((pin[s,t]+dp[s,t])/pin[s,t])^om-1));

#Auxillary variables
@variable(m, slack1[S,time_grid,x_grid] >= 0, start = 10)  #auxiliary variable 
@variable(m, slack2[S,time_grid,x_grid] >= 0, start = 10)  #auxiliary variable
@variable(m, slack3[S,time_grid,x_grid] >= 0, start = 10)  #auxiliary variable

#Auxillary constraints (let's us avoid the fractional terms)
@NLconstraint(m, slackeq1[s=S, t = time_grid, x = x_grid],  slack1[s,t,x]*px[s,t,x] - c3*fx[s,t,x]*fx[s,t,x] == 0)
@NLconstraint(m, slackeq2[s=S, t = time_grid, x = x_grid],  slack2[s,t,x]*px[s,t,x] - 2*c1*fx[s,t,x] == 0)
@NLconstraint(m, slackeq3[s=S, t = time_grid, x = x_grid],  slack3[s,t,x]*px[s,t,x]*px[s,t,x] - c1*fx[s,t,x]*fx[s,t,x] == 0)

#PDE equations
@constraint(m, mass[s=S, t = time_grid[1:end-1], x = x_grid[1:end-1]], 
(px[s,t+1,x]-px[s,t,x])/dt + c1*(fx[s,t+1,x+1]-fx[s,t+1,x])/dx == 0 )
@constraint(m, momentum[s=S, t = time_grid[1:end-1], x = x_grid[1:end-1]], 
(fx[s,t+1,x]-fx[s,t,x])/dt == -slack2[s,t+1,x]*(fx[s,t+1,x+1]-fx[s,t+1,x])/dx +
slack3[s,t+1,x]*(px[s,t+1,x+1]-px[s,t+1,x])/dx -c2*(px[s,t+1,x+1]-px[s,t+1,x])/dx - slack1[s,t+1,x])

#Boundary Conditions
@constraint(m,p_in[s=S, t = time_grid],  pin[s,t] == 54)
@constraint(m,press_in[s=S, t = time_grid],  px[s,t,1] == pin[s,t] + dp[s,t])
@constraint(m,press_out[s=S, t = time_grid], px[s,t,x_grid[end]] == pout[s,t])

#Initial Conditions
@constraint(m, mass_ss[s=S, t = 1, x = x_grid[1:end-1]],(fx[s,t,x+1] - fx[s,t,x]) == 0)
@constraint(m, momentum_ss[s=S, t = 1, x = x_grid[1:end-1]], -c2*(px[s,t,x+1] - px[s,t,x])/dx - slack1[s,t,x] == 0);

#Periodic terminal constraint
@variable(m,linepack[S,time_grid])
@constraint(m,linepack_cons[s=S], linepack[s,time_grid[end]] >= linepack[s,time_grid[1]]);
@constraint(m,linepack_def[s in S,t = time_grid],linepack[s,t] == sum( fx[s,t,x] for x in x_grid)*dx);

compressor_cost = 0.1  # compression cost in USD/kWh
@variable(m, powercost[S])
@constraint(m,boostcosteqn[s = S], powercost[s] == sum(compressor_cost*pow[s,t] for t = time_grid)*dt/3600)

gas_cost = 1000        # value of gas delivered
@variable(m, demandcost[s = S])
@constraint(m, integratedGasCost[s = S], demandcost[s] == sum(gas_cost*ddeliver[s,t] for t = time_grid));

@objective(m,Min,(1/NS)*sum(powercost[s] - demandcost[s] for s = S));

@constraint(m,nonant[s in S, t in time_grid; s>=1], pow[s,t]==pow[1,t]);

solve(m)

d_profile=zeros(NS,24)
for s=S
    for t=time_grid
d_profile[s,t] = getvalue(ddeliver[s,t])
    end
    step(time_grid,d_profile[s,:],color="red")
end
ax = gca()
ax[:set_xlim]([1,24])
xlabel("Time [h]")
ylabel("Demand delivered [scm/h *10^4]")
grid("on")

for s=S
    for t=time_grid
d_profile[s,t] = getvalue(pow[s,t])
    end
    step(time_grid,d_profile[s,:],color="blue")
end
ax = gca()
ax[:set_xlim]([1,24])
xlabel("Time [h]")
ylabel("Compressor Power [-]")
grid("on")


