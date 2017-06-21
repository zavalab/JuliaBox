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
#We will use the same for each pipeline
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
time_grid = 1:48;
x_grid = 1:10;
dx = pipe_length / (length(x_grid) - 1)  #equally spaced spatial grid points
dt = horizon / length(time_grid);         #equally spaced time grid points

links = 1:13
plinks = [1,13]
alinks = collect(2:12)
nodes = 1:14
supplies = [1]
demands = [1]
supply_locs = [1]
demand_locs = [14]
lstartloc = collect(1:13) #link start locations
lendloc = collect(2:14)   #link end locations


m = Model(solver = IpoptSolver())
#Define inlet, outlet, and flow variables
@variable(m,pin[links,time_grid] >= 0, start = 60)
@variable(m,pout[links,time_grid] >= 0, start = 60)
@variable(m,fin[links,time_grid] >=0, start = 100)
@variable(m,fout[links,time_grid] >= 0, start = 100)
@variable(m, px[links,time_grid,x_grid] >= 0, start = 60)
@variable(m, fx[links,time_grid,x_grid] >= 0, start = 10)
# Couple flow convenience variables
@constraint(m, flow_in[link = links,t = time_grid],  fx[link,t,1] == fin[link,t])
@constraint(m, flow_out[link = links,t = time_grid], fx[link,t,x_grid[end]] == fout[link,t]);

# supply variable
@variable(m,sgen[supplies,time_grid] >= 0, start = 100)

# Demand variables
@variable(m,fdeliver[demands,time_grid] >= 0, start = 100)
@variable(m, fdemand[demands,time_grid] >= 0)
@constraint(m, gasLimit[d = demands,t = time_grid], fdeliver[d,t] <= fdemand[d,t]);
@constraint(m, demand[d = demands,t = time_grid], fdemand[d,t] == 150);

@variable(m, dp[link = alinks,time_grid] >= 0 , start = 10)
@variable(m, pow[link = alinks,time_grid] >= 0, start = 500);
@NLconstraint(m,powereqn[link = alinks,t = time_grid], pow[link,t] == c4*fin[link,t]*(((pin[link,t]+dp[link,t])/pin[link,t])^om-1));

#Auxillary variables
@variable(m, slack1[links,time_grid,x_grid] >= 0, start = 10)  #auxiliary variable for friction loss term
@variable(m, slack2[links,time_grid,x_grid] >= 0, start = 10)  #auxiliary variable
@variable(m, slack3[links,time_grid,x_grid] >= 0, start = 10)  #auxiliary variable

#Auxillary constraints (let's us avoid the fractional terms)
@NLconstraint(m, slackeq1[link = links,t = time_grid, x = x_grid],  slack1[link,t,x]*px[link,t,x] - c3*fx[link,t,x]*fx[link,t,x] == 0)
@NLconstraint(m, slackeq2[link = links,t = time_grid, x = x_grid],  slack2[link,t,x]*px[link,t,x] - 2*c1*fx[link,t,x] == 0)
@NLconstraint(m, slackeq3[link = links,t = time_grid, x = x_grid],  slack3[link,t,x]*px[link,t,x]*px[link,t,x] - c1*fx[link,t,x]*fx[link,t,x] == 0)

#PDE equations
@constraint(m, mass[link = links,t = time_grid[1:end-1], x = x_grid[1:end-1]],
    (px[link,t+1,x]-px[link,t,x])/dt + c1*(fx[link,t+1,x+1]-fx[link,t+1,x])/dx == 0 )
@constraint(m, momentum[link = links,t = time_grid[1:end-1], x = x_grid[1:end-1]],
    (fx[link,t+1,x]-fx[link,t,x])/dt == -slack2[link,t+1,x]*(fx[link,t+1,x+1]-fx[link,t+1,x])/dx +
                            slack3[link,t+1,x]*(px[link,t+1,x+1]-px[link,t+1,x])/dx -c2*(px[link,t+1,x+1]-px[link,t+1,x])/dx - slack1[link,t+1,x])

#Boundary Conditions
#@constraint(m,f_in[t = time_grid],  fin[t] == 80)
@constraint(m,press_in1[link = plinks,t = time_grid], px[link,t,1] == pin[link,t])
@constraint(m,press_in2[link = alinks,t = time_grid], px[link,t,1] == pin[link,t] + dp[link,t])
@constraint(m,press_out[link = links,t = time_grid], px[link,t,x_grid[end]] == pout[link,t])

#Initial Conditions
@constraint(m, mass_ss[link = links,t = 1, x = x_grid[1:end-1]],(fx[link,t,x+1] - fx[link,t,x]) == 0)
@constraint(m, momentum_ss[link = links,t = 1, x = x_grid[1:end-1]], -c2*(px[link,t,x+1] - px[link,t,x])/dx - slack1[link,t,x] == 0);

# node balances
@constraint(m,nodebalance[node = nodes,t = time_grid],
sum(fout[link,t] for link in links if lendloc[link] == node)
+ sum(sgen[supply,t] for supply in supplies if supply_locs[supply] == node)
- sum(fin[link,t] for link in links if lstartloc[link] == node)
- sum(fdeliver[demand,t] for demand in demands if demand_locs[demand] == node) == 0)

#Periodic terminal constraint
@variable(m,linepack[links,time_grid])
@constraint(m,linepack_def[link = links,t = time_grid],linepack[link,t] == sum(fx[link,t,x] for x in x_grid)*dx)
@constraint(m,linepack_cons[link = links], linepack[link,time_grid[end]] >= linepack[link,time_grid[1]])

compressor_cost = 0.1
@variable(m, powercost[alinks])
@constraint(m,boostcosteqn[link = alinks], powercost[link] == sum(compressor_cost*pow[link,t] for t = time_grid)*dt/3600)

gas_cost = -1000
@variable(m, demandcost[demands])
@constraint(m, integratedGasCost[d = demands], demandcost[d] == sum(gas_cost*fdeliver[d,t] for t = time_grid));

@objective(m,Min,sum(powercost[link] for link in alinks) + sum(demandcost[d] for d in demands));
