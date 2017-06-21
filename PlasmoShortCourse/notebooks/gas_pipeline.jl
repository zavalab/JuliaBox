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
Tgas = 293.15      		     # reference temperature [K]
Cp = 2.34        		     # heat capacity @ constant pressure [kJ/kg-K]
Cv = 1.85        		     # heat capacity @ constant volume [kJ/kg-K]
U = 1.0*0.1     		     # pipe heat transfer coefficient [J/m2-s-K]

#scaling factors
ffac=(1e+6*rhon)/(24*3600)                     # from scmx10-6/day to kg/s
ffac2=(3600)/(1e+4*rhon)                       # from kg/s to scmx10-4/hr
pfac=1e+5                                      # from bar to Pa
pfac2=1e-5                                     # from Pa to bar
dfac=1e-3                                      # from mm to m
lfac=1e+3;                                      # from km to m

#parameters
diameter = 0.92
pipe_length = 60000
area = (1/4)*pi*diameter*diameter
lam = (2*log10(3.7*diameter/(eps*dfac)))^(-2);

#other calculated constants
gam = Cp/Cv       		     	# expansion coefficient [-]
nu2 = gam*z*R*Tgas/M  			# gas speed of sound
om = (gam-1.0)/gam 		     	# aux constant [-]
c4 = (1/ffac2)*(Cp*Tgas)
#lumped constants for pipeline
c1 = (pfac2/ffac2)*(nu2/area)
c2 = area*(ffac2/pfac2)
c3 = area*(pfac2/ffac2)*(8*lam*nu2)/(pi*pi*diameter^5);

horizon = 24*3600;
time_grid = 1:24;
x_grid = 1:10;
dx = pipe_length / (length(x_grid) - 1)  #equally spaced spatial grid points
dt = horizon / length(time_grid);         #equally spaced time grid points

m = Model(solver = IpoptSolver())
#Define inlet, outlet, and flow variables
@variable(m,pin[time_grid] >= 0, start = 60)
@variable(m,pout[time_grid] >= 0, start = 60)
@variable(m,fin[time_grid] >=0, start = 100)
@variable(m,fout[time_grid] >= 0, start = 100)
@variable(m, px[time_grid,x_grid] >= 0, start = 60)
@variable(m, fx[time_grid,x_grid] >= 0, start = 10)
# Couple flow convenience variables
@constraint(m, flow_in[t = time_grid],  fx[t,1] == fin[t])
@constraint(m, flow_out[t = time_grid], fx[t,x_grid[end]] == fout[t]);

@variable(m,fdeliver[time_grid] >= 0, start = 100)
@variable(m, fdemand[time_grid] >= 0)
@constraint(m, gasLimit[t = time_grid], fdeliver[t] <= fdemand[t]);
@constraint(m,deliver[t=time_grid],fout[t] == fdeliver[t])
@constraint(m, demand[t = time_grid], fdemand[t] == 150);

@variable(m, dp[time_grid] >= 0 , start = 10)
@variable(m, pow[time_grid] >= 0, start = 500);
@NLconstraint(m, powereqn[t = time_grid], pow[t] == c4*fin[t]*(((pin[t]+dp[t])/pin[t])^om-1));

#Auxillary variables
@variable(m, slack1[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable for friction loss term
@variable(m, slack2[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable
@variable(m, slack3[time_grid,x_grid] >= 0, start = 10)  #auxiliary variable

#Auxillary constraints (let's us avoid the fractional terms)
@NLconstraint(m, slackeq1[t = time_grid, x = x_grid],  slack1[t,x]*px[t,x] - c3*fx[t,x]*fx[t,x] == 0)
@NLconstraint(m, slackeq2[t = time_grid, x = x_grid],  slack2[t,x]*px[t,x] - 2*c1*fx[t,x] == 0)
@NLconstraint(m, slackeq3[t = time_grid, x = x_grid],  slack3[t,x]*px[t,x]*px[t,x] - c1*fx[t,x]*fx[t,x] == 0)

#PDE equations
@constraint(m, mass[t = time_grid[1:end-1], x = x_grid[1:end-1]],
    (px[t+1,x]-px[t,x])/dt + c1*(fx[t+1,x+1]-fx[t+1,x])/dx == 0 )
@constraint(m, momentum[t = time_grid[1:end-1], x = x_grid[1:end-1]],
    (fx[t+1,x]-fx[t,x])/dt == -slack2[t+1,x]*(fx[t+1,x+1]-fx[t+1,x])/dx +
                            slack3[t+1,x]*(px[t+1,x+1]-px[t+1,x])/dx -c2*(px[t+1,x+1]-px[t+1,x])/dx - slack1[t+1,x])

#Boundary Conditions
@constraint(m,p_in[t = time_grid],  pin[t] == 50)
#@constraint(m,f_in[t = time_grid],  fin[t] == 80)
@constraint(m,press_in[t = time_grid],  px[t,1] == pin[t] + dp[t])
@constraint(m,press_out[t = time_grid], px[t,x_grid[end]] == pout[t])

#Initial Conditions
@constraint(m, mass_ss[t = 1, x = x_grid[1:end-1]],(fx[t,x+1] - fx[t,x]) == 0)
@constraint(m, momentum_ss[t = 1, x = x_grid[1:end-1]], -c2*(px[t,x+1] - px[t,x])/dx - slack1[t,x] == 0);

#Periodic terminal constraint
@variable(m,linepack[time_grid])
@constraint(m,linepack_def[t = time_grid],linepack[t] == sum(fx[t,x] for x in x_grid)*dx)
@constraint(m,linepack_cons, linepack[time_grid[end]] >= linepack[time_grid[1]])

compressor_cost = 0.0
@variable(m, powercost)
@constraint(m,boostcosteqn, powercost == sum(compressor_cost*pow[t] for t = time_grid)*dt/3600)

gas_cost = -1000
@variable(m, demandcost)
@constraint(m, integratedGasCost, demandcost == sum(gas_cost*fdeliver[t] for t = time_grid));

@objective(m,Min,powercost + demandcost);
