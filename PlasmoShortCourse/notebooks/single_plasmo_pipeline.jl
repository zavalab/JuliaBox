using Plasmo
using Ipopt

global horizon = 48
include("gas_components.jl")

#pipeline data
pipe_length = 30000
diameter = 0.92
time_grid = 1:48
x_grid = 1:10
min_pressure = 1
max_pressure = 100
min_flow = 1
max_flow = 300

#compressor data
dp_min = 1
dp_max = 20
min_power = 1
max_power = 100000
comp_cost = 0.1
pdata = PipeData(pipe_length,diameter,time_grid,x_grid,min_pressure,max_pressure,min_flow,max_flow) #pipe data
cdata = CompData(dp_min,dp_max,min_power,max_power,comp_cost) #compressor data

supply_node_data = NodeData(time_grid,1,300,0,50,70,-1000,0.0)
demand_node_data = NodeData(time_grid,0,0,1,50,70,-1000,0.0)

#Create the Plasmo Graph Model
model = GraphModel()
graph = getgraph(model)
n1 = add_node!(model)
n2 = add_node!(model)
pipe = add_edge!(model,n1,n2)

euler_model = euleractivelink(pdata,cdata)
euler_model = steadystatestart(euler_model,pdata)

#set the models on our nodes and edge
setmodel!(n1,gasnode(supply_node_data))
setmodel!(n2,gasnode(demand_node_data))
setmodel!(pipe,euler_model)

#set our linking constraints for node conservation
@linkconstraint(model,graph,n1,[t = time_grid], 0 == sum(in_edges(graph,n1)[i][:fout][t] for i = 1:length(in_edges(graph,n1))) -
sum(out_edges(graph,n1)[i][:fin][t] for i = 1:length(out_edges(graph,n1))) + n1[:total_supplied][t] - n1[:total_delivered][t])

@linkconstraint(model,graph,n2,[t = time_grid], 0 == sum(in_edges(graph,n2)[i][:fout][t] for i = 1:length(in_edges(graph,n2))) -
sum(out_edges(graph,n2)[i][:fin][t] for i = 1:length(out_edges(graph,n2))) + n2[:total_supplied][t] - n2[:total_delivered][t])

#set linking constraints for our pipe boundary condition
@linkconstraint(model,graph,pipe,[t = time_grid],pipe[:pin][t] == getconnectedfrom(graph,pipe)[:pressure][t])
@linkconstraint(model,graph,pipe,[t = time_grid],pipe[:pout][t] == getconnectedto(graph,pipe)[:pressure][t])

#fix the demand
d = getmodel(n2)
@constraint(d,[t = time_grid],n2[:demand][1,t] == 30)

model.solver = IpoptSolver()
solve(model,graph)

# flat = create_flat_graph_model(model)
# flat.solver = IpoptSolver()
