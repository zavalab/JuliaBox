using Ipopt
include("13_pipe_gas.jl")
include("../parse_input_data.jl")
include("../write_outputs.jl")
include("plot_13_pipe.jl")

time_grid = 1:24  #24 hour horizon
x_grid = 1:3      #3 grid points in every pipeline
const horizon = 24*3600  #the horizon is in seconds

#Read data files for nodes, links, supplies, and demands
nodeinfo = readdlm("DAT_13pipe_GAS/nodeinfo.tab")
linkinfo = readdlm("DAT_13pipe_GAS/linkinfo.tab")
supplyinfo = readdlm("DAT_13pipe_GAS/supinfo.tab")
demandinfo = readdlm("DAT_13pipe_GAS/demandinfo.tab")

#Parse the data into dictionaries
node_data = parsedlm(nodeinfo)
link_data = parsedlm(linkinfo)
supply_data = parsedlm(supplyinfo)
demand_data = parsedlm(demandinfo)

#create the network using our 13 pipe function
gas_network,supply_nodes,demand_nodes = create_13_pipe_network(node_data,link_data,supply_data,demand_data,time_grid,x_grid,horizon)

#generate an array of some constant demands
demand = fill(10*ffac*ffac2,1,48)
#Create the flattened model of the gas network
m = generatemodel!(gas_network)

# NOTE: We could have fixed demands and pressures before generating the flat model, but we can keep the graph more general if we just fix data afterwards
fix_demands!(m,demand_nodes,demand) #fix all the demand node demands using our fix_demands function
fix_supply_pressure!(m,supply_nodes,54) #fix all the supply pressures using our function

#Set the solver and solve the nonlinear program
m.solver = IpoptSolver(linear_solver = "ma57", tol = 1E-6)
solve(m)

#write the results
fpath = "./results/gas_delivered.tab"
write_gasnet_data(gas_network,length(x_grid),length(time_grid))
#write_gas_delivered(demand_nodes,fpath)
plot_13_pipe(press_profile,length(x_grid),length(time_grid))
