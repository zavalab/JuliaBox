# stochastic optimal control problem PLASMO and PIPS-NLP/IPOPT
# Yankai Cao and Victor M. Zavala
# University of Wisconsin-Madison, 2016

using Ipopt
using Plasmo
using JuMP

# sets
TF=24*3600                           # horizon time - [s]
Nt=24                                # number of temporal grid points
Nx=3                                 # number of spatial grid points
S=3                                  # number of scenarios
TIMEG=1:Nt                           # set of temporal grid points
TIMEGm=1:Nt-1                        # set of temporal grid points minus 1
DIS=1:Nx                             # set of spatial grid points
SCENG=1:S                            # scenario set
dtG = 0                              # temporal grid spacing - [s]

# links
type LinkData                        # set of links
     name::String
     startloc::String           # start node
     endloc::String             # end node
     diam::Float64                   # link diameter - mm
     length::Float64                 # link length - km
     ltype::String              # link type, passive or active
     c1                              # aux constant
     c2                              # aux constant
     c3                              # aux constant
     dx                              # spatial grid spacing - [m]
     lam                             # friction coefficient - []
     A                               # pipe transveral area - [m^2]
end
linkDict = Dict{String, LinkData}()

# nodes
type NodeData
     name::String
     pmin::Float64                # min pessure - bar
     pmax::Float64                # max pressure - bar
end
nodeDict = Dict{String, NodeData}()

# supply
type SupplyData                   # set of suppliers
     name::String
     loc::String             # supply location
     min::Float64                 # min supply - scmx106/day
     max::Float64                 # max supply - scmx106/day
end
supDict = Dict{String, SupplyData}()

# demand
type DemandData                    # set of suppliers
     name::String
     loc::String              # demand location
     d::Float64                    # base demand - scmx106/day
     stochd                        # stochastic demands - [scmx10-4/hr]
end
demDict = Dict{String, DemandData}()

# physical data
eps= 0.025            # pipe rugosity - [mm]
z= 0.80               # gas compressibility  - []
rhon=0.72             # density of air at normal conditions - [kg/m3]
R=8314.0              # universal gas constant [J/kgmol-K]
M=18.0                # gas molar mass [kg/kgmol]
pi=3.14               # pi
nu2=0                 # gas speed of sound [m2/s2]
Tgas = 293.15         # reference temperature [K]
Cp = 2.34             # heat capacity @ constant pressure [kJ/kg-K]
Cv = 1.85             # heat capacity @ constant volume [kJ/kg-K]
gam = Cp/Cv           # expansion coefficient [-]
om = (gam-1.0)/gam    # aux constant [-]
U = 1.0*0.1           # pipe heat transfer coefficient [J/m2-s-K]
Tamb = 20+273.15      # soil temperature [K]
Tsup = 30+273.15      # supply temperature [K]

# scaling and constants
ffac = 0              # from scmx106/day to kg/s
ffac2 = 0             # from kg/s to scmx10-4/hr
pfac = 0              # from bar to Pa
pfac2 = 0             # from Pa to bar
dfac = 0              # from mm to m
lfac = 0              # from km to m
c4 = 0                # aux constant [kW/(scmx10-4/hr)]

# cost factors
ce = 0.1             # cost of compression [$/kWh]
cd = 1e6             # demand tracking cost [-]
cT = 1e6             # terminal constraint cost [-]
cs =   0             # supply cost [$/scmx10-4]

# set data
include("set_PDEGas_data.jl")
NODE = keys(nodeDict)
LINK = keys(linkDict)
SUP =  keys(supDict)
DEM =  keys(demDict)

# define scenario model
include("createPDEGasModel.jl")

# create two-stage graph model
gas = GraphModel()
master = Model()
master_node=add_node(gas,master)
@variable(master, dp[j = LINK; linkDict[j].ltype == "a"], start= 10)   # compressor boost - [bar]
@variable(master, dem[DEM],    start=100)                              # demand flow - [scmx10-4/hr]

# create array of children models
child_nodes = Array{Plasmo.NodeOrEdge}(S)
for s in SCENG
    gasch = createPDEGasModel(s)
    child = add_node(gas,gasch)
    child_nodes[s]=child
    @linkconstraint(gas,[j in LINK,t in TIMEG; linkDict[j].ltype =="a" && t ==1],   dp[j] == gasch[:dp][j,t])
    @linkconstraint(gas,[j in DEM, t in TIMEG;                            t ==1],  dem[j] == gasch[:dem][j,t])
end
# solve with Ipopt
gas.solver = IpoptSolver()
solve(gas)