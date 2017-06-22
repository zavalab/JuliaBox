using JuMP
using Ipopt
using PyPlot
using Plasmo

eps= 0.025  # pipe rugosity - [mm]
z= 0.80     # gas compressibility  - []
rhon=0.72   # density of air at normal conditions - [kg/m3]
R=8314.0    # universal gas constant [J/kgmol-K]
M=18.0      # gas molar mass [kg/kgmol]
pi=3.14     # pi
Tgas = 293.15    # reference temperature [K]
Cp = 2.34        # heat capacity @ constant pressure [kJ/kg-K]
Cv = 1.85        # heat capacity @ constant volume [kJ/kg-K]

#scaling factors
ffac=(1e+6*rhon)/(24*3600)                     # from scmx10-6/day to kg/s
ffac2=(3600)/(1e+4*rhon)                       # from kg/s to scmx10-4/hr
pfac=1e+5                                      # from bar to Pa
pfac2=1e-5                                     # from Pa to bar
dfac=1e-3                                      # from mm to m
lfac=1e+3;                                      # from km to m

#We will also use a global horizon variable!
global horizon = 24*3600;
include("gas_components.jl");

#pipeline data
pipe_length = 100000
diameter = 0.92
time_grid = 1:24
x_grid = 1:3
min_pressure = 1
max_pressure = 100
min_flow = 1
max_flow = 300

#compressor data
dp_min = 1
dp_max = 20
min_power = 1
max_power = 5000
comp_cost = 0.1;

pressure_lower = 50
pressure_upper = 70

#demand parameters
gas_cost = -1000

#supply parameters
supply_cost = 0
fgen_lower = 0
fgen_upper = 125

n_supplies = 1
n_demands = 1
# node_lower_pressure = 50
# node_upper_pressure = 70


pdata = PipeData(pipe_length,diameter,time_grid,x_grid,min_pressure,max_pressure,min_flow,max_flow) #pipe data
cdata = CompData(dp_min,dp_max,min_power,max_power,comp_cost) #compressor data
supply_node_data = NodeData(time_grid,n_supplies,fgen_upper,0,54,70,gas_cost,supply_cost)
demand_node_data = NodeData(time_grid,0,fgen_upper,n_demands,39,41,gas_cost,supply_cost)
junction_node_data = NodeData(time_grid,0,fgen_upper,0,34,70,gas_cost,supply_cost);

#Create the Graph model and get the graph object from it
model = GraphModel()
graph = getgraph(model)

#add 14 gas nodes
for i = 1:14
    node = add_node!(model)
    if i == 1
        setmodel!(node,gasnode(supply_node_data))
    elseif i == 14
        setmodel!(node,gasnode(demand_node_data))
    else
        setmodel!(node,gasnode(junction_node_data))
    end
end

#Add the edges between them
for j = 1:13
    n_from = getnode(model,j)
    n_to = getnode(model,j+1)
    edge = add_edge!(model,n_from,n_to)
    if j == 1 || j == 13
        pdata = PipeData(300000,diameter,time_grid,x_grid,min_pressure,max_pressure,min_flow,max_flow) #pipe data
        #pipe_model = adiabaticpassivelink(pdata)
        pipe_model = weymouthpassivelink(pdata)
        pipe_model = steadystatestart(pipe_model,pdata)
        setmodel!(edge,pipe_model)
    else
        pdata = PipeData(100000,diameter,time_grid,x_grid,min_pressure,max_pressure,min_flow,max_flow) #pipe data
        #pipe_model = adiabaticactivelink(pdata,cdata)
        pipe_model = weymouthactivelink(pdata,cdata)
        pipe_model = steadystatestart(pipe_model,pdata)
        setmodel!(edge,pipe_model)
    end
end

#Add the linking constraints

#Node Conservation at each node
for (index,node) in getnodes(model)
    @linkconstraint(model,graph,node,[t = time_grid], 0 == sum(in_edges(graph,node)[i][:fout][t] for i = 1:length(in_edges(graph,node))) -
        sum(out_edges(graph,node)[i][:fin][t] for i = 1:length(out_edges(graph,node))) + node[:total_supplied][t] - node[:total_delivered][t])
end

#Boundary conditions on each edge
for (index,edge) in getedges(model)
    @linkconstraint(model,graph,edge,[t = time_grid],edge[:pin][t] == getconnectedfrom(graph,edge)[:pressure][t])
    @linkconstraint(model,graph,edge,[t = time_grid],edge[:pout][t] == getconnectedto(graph,edge)[:pressure][t])
end

#Fix the demand at the last node
demand = getnode(graph,14)
d = getmodel(demand)
@constraint(d,[t = time_grid],d[:demand][1,t] == 41.67)

supply = getnode(graph,1)
s = getmodel(supply)
@constraint(s,[t = time_grid],s[:pressure][t] == 54)

model.solver = IpoptSolver(tol = 1E-6)
#solve(model,graph)
