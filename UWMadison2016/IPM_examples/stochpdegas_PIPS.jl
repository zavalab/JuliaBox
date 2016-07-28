push!(LOAD_PATH, pwd())
using Ipopt
using Plasmo
using JuMP
import MPI
# sets
TF=24*3600                              # horizon time - [s]
Nt=24                                 # number of temporal grid points
Nx=3                                 # number of spatial grid points
S=32                                 # number of scenarios

TIMEG=1:Nt                           # set of temporal grid points
TIMEGm=1:Nt-1                        # set of temporal grid points minus 1
DIS=1:Nx                             # set of spatial grid points
SCENG=1:S                            # scenario set

# links
type LinkData                        # set of links
     name::ASCIIString
     startloc::ASCIIString           # start node
     endloc::ASCIIString             # end node
     diam::Float64                   # link diameter - mm
     length::Float64                 # link length - km
     ltype::ASCIIString              # link type, passive or active
     c1                              # aux constant 
     c2                              # aux constant 
     c3                              # aux constant 
     dx                              # spatial grid spacing - [m]     
     lam                             # friction coefficient - []
     A                               # pipe transveral area - [m^2]
end
linkDict = Dict{ASCIIString, LinkData}()


# nodes
type NodeData
     name::ASCIIString
     pmin::Float64		      # min pessure - bar
     pmax::Float64  		      # max pressure - bar
end
nodeDict = Dict{ASCIIString, NodeData}()


# supply
type SupplyData                       # set of suppliers
     name::ASCIIString
     loc::ASCIIString		      # supply location
     min::Float64   		      # min supply - scmx106/day
     max::Float64   		      # max supply - scmx106/day
end
supDict = Dict{ASCIIString, SupplyData}()


# demand
type DemandData                        # set of suppliers
     name::ASCIIString
     loc::ASCIIString    	       # demand location
     d::Float64                        # base demand - scmx106/day
     stochd			       # stochastic demands - [scmx10-4/hr]
end
demDict = Dict{ASCIIString, DemandData}()


# physical data
eps= 0.025		             # pipe rugosity - [mm]
z= 0.80        			     # gas compressibility  - []
rhon=0.72         		     # density of air at normal conditions - [kg/m3]
R=8314.0       			     # universal gas constant [J/kgmol-K]
M=18.0    			     # gas molar mass [kg/kgmol]
pi=3.14         		     # pi
nu2=0               		     # gas speed of sound [m2/s2]
Tgas = 293.15      		     # reference temperature [K]
Cp = 2.34        		     # heat capacity @ constant pressure [kJ/kg-K]
Cv = 1.85        		     # hat capacity @ constant volume [kJ/kg-K]
gam = Cp/Cv       		     # expansion coefficient [-]
om = (gam-1.0)/gam 		     # aux constant [-]
U = 1.0*0.1     		     # pipe heat transfer coefficient [J/m2-s-K]
Tamb = 20+273.15   		     # soil temperature [K]
Tsup = 30+273.15   		     # supply temperature [K]

# scaling and constants
ffac = 0			     # from scmx106/day to kg/s
ffac2 = 0              		     # from kg/s to scmx10-4/hr
pfac = 0             		     # from bar to Pa
pfac2 = 0             		     # from Pa to bar
dfac = 0             		     # from mm to m
lfac = 0             		     # from km to m
c4 = 0             		     # aux constant [kW/(scmx10-4/hr)]

# cost factors
ce = 0.1		             # cost of compression [$/kWh]
cd = 1e6         		     # demand tracking cost [-]
cT = 1e6         		     # terminal constraint cost [-]
cs =   0         		     # supply cost [$/scmx10-4]

# define temporal discretization info
dtG = 0			             # temporal grid spacing - [s]

include("set_data_paper.jl")
NODE = keys(nodeDict)
LINK = keys(linkDict)
SUP =  keys(supDict)
DEM =  keys(demDict)

function createGasModel(SCEN)
	 m = Model(solver=IpoptSolver())		 
         @variable(m, nodeDict[j].pmin<=p[SCEN, j in NODE, TIMEG]<=nodeDict[j].pmax, start= 50)             # node pressure - [bar]
	 @variable(m, 0<=dp[SCEN, j = LINK,  TIMEG; linkDict[j].ltype == "a"]<=100, start= 10)             # compressor boost - [bar]
         @variable(m, 1<=fin[SCEN, LINK, TIMEG]<=500, start= 100)                                           # flow in pipe - [scmx10-4/hr]
         @variable(m, 1<=fout[SCEN, LINK, TIMEG]<=500, start= 100)                                          # flow out pipe - [scmx10-4/hr]
	 @variable(m, 0.01<=sG[SCEN, j in SUP, TIMEG]<=supDict[j].max, start = 10)                          # supply flow - [scmx10-4/hr]
	 @variable(m, dem[SCEN, DEM, TIMEG],    start=100)                                                  # demand flow - [scmx10-4/hr]
         @variable(m, 0<=pow[SCEN, j = LINK, TIMEG; linkDict[j].ltype == "a"]<=3000, start= 1000)           # compressor power - [kW]
	 @variable(m, slack[SCEN, LINK, TIMEG, DIS]>=0, start= 10)                                          # auxiliary variable

         # define spatio-temporal variables
         @variable(m, 10<=px[SCEN, LINK, TIMEG, DIS]<=100, start= 50)                                       # link pressure profile - [bar]
         @variable(m, 1<=fx[SCEN, LINK, TIMEG, DIS]<=100, start= 100)                                       # link flow profile - [scmx10-4/hr]
        
	 # compressor equations
         @NLconstraint(m, powereq[i = SCEN, j = LINK, t = TIMEG; linkDict[j].ltype == "a"], pow[i,j,t] == c4*fin[i,j,t]*(((p[i,linkDict[j].startloc,t] + dp[i,j,t])/(p[i,linkDict[j].startloc,t]))^om-1))

	 # node balance [mass]
         @constraint(m,nodemeq[k = SCEN, i in NODE, t = TIMEG], sum{       fout[k,j,t], j in LINK ; linkDict[j].endloc==i}
                                                                 - sum{       fin[k,j,t],  j in LINK ; linkDict[j].startloc==i}
                                                                 + sum{        sG[k,j,t],  j in SUP ; supDict[j].loc == i }
                                                                 - sum{       dem[k,j,t],  j in DEM ; demDict[j].loc == i }
                                                                 ==0)
         

	 # flow equations for passive and active links
         @constraint(m, flow[i = SCEN, j = LINK, t = TIMEGm, k = 1:(Nx-1)], (px[i,j,t+1,k]-px[i,j,t,k])/dtG + linkDict[j].c1*(fx[i,j,t+1,k+1]-fx[i,j,t+1,k])/(linkDict[j].dx)==0)

         # boundary conditions flow
         @constraint(m, flow_start[i = SCEN, j = LINK, t = TIMEG], fx[i,j,t,1]==fin[i,j,t])
         @constraint(m, flow_end[i = SCEN, j = LINK, t = TIMEG], fx[i,j,t,Nx]==fout[i,j,t])

         # pressure equations for passive and active links
         @constraint(m, press[i = SCEN, j = LINK, t = TIMEGm,k = 1:(Nx-1)], (fx[i,j,t+1,k]-fx[i,j,t,k])/dtG == - linkDict[j].c2*(px[i,j,t+1,k+1]-px[i,j,t+1,k])/linkDict[j].dx - slack[i,j,t+1,k])
         @NLconstraint(m, slackeq[i = SCEN, j = LINK, t = TIMEG, k = 1:Nx],  slack[i,j,t,k]*px[i,j,t,k] - linkDict[j].c3*fx[i,j,t,k]*fx[i,j,t,k] == 0);
	
         # boundary conditions pressure, passive links
         @constraint(m, presspas_start[i = SCEN, j = LINK, t = TIMEG; linkDict[j].ltype == "p"], px[i,j,t,1] ==  p[i,linkDict[j].startloc,t])
         @constraint(m,   presspas_end[i = SCEN, j = LINK, t = TIMEG; linkDict[j].ltype == "p"], px[i,j,t,Nx] == p[i,linkDict[j].endloc,t])

         # boundary conditions, active links
         @constraint(m, pressact_start[i = SCEN, j = LINK, t = TIMEG; linkDict[j].ltype == "a"], px[i,j,t,1] ==  p[i,linkDict[j].startloc,t] + dp[i,j,t])
         @constraint(m,   pressact_end[i = SCEN, j = LINK, t = TIMEG; linkDict[j].ltype == "a"], px[i,j,t,Nx] == p[i,linkDict[j].endloc,t])
	 	 
	 # fix pressure at supply nodes
         @constraint(m, suppres[i = SCEN, j in SUP, t = TIMEG], p[i,supDict[j].loc,t] == nodeDict[supDict[j].loc].pmin)
	 

	 # discharge pressure for compressors
	 @constraint(m, dispress[i in SCEN,j in LINK,t in TIMEG; linkDict[j].ltype=="a"],  p[i,linkDict[j].startloc,t]+dp[i,j,t] <= nodeDict[linkDict[j].startloc].pmax)

	 # line pack constraints
         @constraint(m, line_packT[i in SCEN],  sum{  sum{fx[i,j,Nt,k], k in DIS}*linkDict[j].dx, j in LINK} >= sum{ sum{fx[i,j,1,k],k in DIS}*linkDict[j].dx, j = LINK})
	 
         # ss constraints
         @constraint(m, flow_ss[i = SCEN, j = LINK, t =0, k = 1:(Nx-1)], (fx[i,j,t+1,k+1]-fx[i,j,t+1,k])==0)
         @constraint(m, pres_ss[i = SCEN, j = LINK, t =0, k = 1:(Nx-1)],  - linkDict[j].c2*(px[i,j, t+1,k+1]-px[i,j,t+1,k])/linkDict[j].dx - slack[i,j,t+1,k] == 0)

     	 @objective(m, Min, 1e-3*(1.0/S)*sum{ce*pow[i,j,t]*(dtG/3600),i= SCEN, j = LINK,t = TIMEG; linkDict[j].ltype == "a"}
			      +1e-3*(1.0/S)*sum{cd*(dem[i,j,t]-demDict[j].stochd[i,t])^2, i in SCEN, j in DEM, t in TIMEG})
			    
	 # non-anticipativity constraints (note this is only for t==1)
	 if (length(SCEN)> 1)
	     bs = SCEN[1]
	     @constraint(m,  nonantdq[i in SCEN,j in LINK,t in TIMEG; linkDict[j].ltype =="a" && t ==1 && i != bs],   dp[i,j,t] ==  dp[bs,j,t])
	     @constraint(m,  nonantde[i in SCEN,j in DEM, t in TIMEG;                            t ==1 && i != bs],   dem[i,j,t]== dem[bs,j,t])
	 end
	 
	 return m
end


#=
IL=createGasModel(SCENG)
solve(IL)
=#

IL = NetModel()
@variable(IL, dp[j = LINK; linkDict[j].ltype == "a"], start= 10)                       # compressor boost - [bar]
@variable(IL, dem[DEM],    start=100)                                                  # demand flow - [scmx10-4/hr]
for s in SCENG
   single_scenario = createGasModel(s:s) 
   @addNode(IL, single_scenario, "s$s")
   @constraint(IL,  nonantdq[j in LINK,t in TIMEG; linkDict[j].ltype =="a" && t ==1],   dp[j] ==  getvariable(single_scenario, :dp)[s,j,t])
   @constraint(IL,  nonantde[j in DEM, t in TIMEG;                            t ==1],   dem[j]==  getvariable(single_scenario,:dem)[s,j,t])
end
ParPipsNlp_solve(IL)


