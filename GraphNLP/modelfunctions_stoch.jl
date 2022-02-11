#############################################################################################
# Define constants used within constraints

z = 0.8                              # gas compressibility  - []
rhon = 0.72                          # density of natural gas at normal conditions - [kg/m3]
R = 8314.0                           # universal gas constant [J/kgmol-K]
M = 18.0                             # gas molar mass of natural gas [kg/kgmol]
Tgas = 293.15                        # reference temperature [K]
Cp = 2.34                            # heat capacity @ constant pressure [kJ/kg-K]
Cv = 1.85                            # heat capacity @ constant volume [kJ/kg-K]

#scaling factors for optimization
ffac=(1e+6*rhon)/(24*3600)            # from MMSCM/day to kg/s
ffac2=(3600)/(1e+4*rhon)              # from kg/s to scmx10^-4/hr
pfac=1e+5                             # from bar to Pa
pfac2=1e-5                            # from Pa to bar
dfac=1e-3                             # from mm to m
lfac=1e+3                             # from km to m
gam = Cp/Cv                           # expansion coefficient [-]
n_poly = gam
nu2 = gam*z*R*Tgas/M                  # gas speed of sound
om = (gam-1.0)/gam                    # aux constant [-]
c4 = (1/ffac2)*(Cp*Tgas)              #[kW/(scmx10-4/hr)]

#############################################################################################


# Define a function for creating a subgraph for each junction
function create_junction_model(data,nt)
    # Create junction OptiGraph
    graph = OptiGraph()
    # Create nodes for each time point
    @optinode(graph,time_nodes[1:nt])
    # Unpack data; if the junction demands or supplies any gas, one of these values will be > 0
    n_demands = length(data[:demand_values])
    n_supplies = length(data[:supplies])

    # Iterate over all nodes to define variables and constraints on each node
    for (i,node) in enumerate(time_nodes)
        @variable(node, data[:pmin] <= pressure <= data[:pmax], start = 60)
        @variable(node, 0 <= fgen[1:n_supplies] <= 200, start = 30)
        @variable(node, fdeliver[1:n_demands] >= 0, start = 23)
        @variable(node, fdemand[1:n_demands] >= 0, start = 23)
        @variable(node, foverdemand[1:n_demands] >= 0, start=0)

        @constraint(node,[d = 1:n_demands], (fdeliver[d])*.1 <= (fdemand[d]+foverdemand[d])*.1)

        @expression(node, total_supplied, sum(fgen[s] for s = 1:n_supplies))
        @expression(node, total_delivered,sum(fdeliver[d] for d = 1:n_demands))

        # If there is a demand on this node, define an objective function
        if n_demands > 0
            @objective(node,Min,-sum(1000*(fdeliver[d]-1.25*foverdemand[d]) for d = 1:n_demands))  #Changed from 2 to 1.25
        end
    end
    return graph
end

# Define function for creating a subgraph for each pipeline
function create_pipeline_model(data,nt,nx)
    # Create pipeline OptiGraph
    graph = OptiGraph()
    #Create grid of modelnodes
    @optinode(graph,grid[1:nt,1:nx])
    #unpack data
    c1 = data[:c1]; c2 = data[:c2]; c3 = data[:c3]
    dx = data[:pipe_length] / (nx - 1)


    # Create variables on each node in the grid
    for node in grid
        @variable(node, 1 <= px <= 100,start = 50)
        @variable(node, 0 <= fx <= 100, start = 30)
        @variable(node,slack >= 0, start = 0.1)
        @NLconstraint(node, (slack*px - c3*fx*fx) == 0)
    end

    #Setup boundary variables
    @expression(graph,fin[t=1:nt],grid[:,1][t][:fx])
    @expression(graph,fout[t=1:nt],grid[:,end][t][:fx])
    @expression(graph,pin[t=1:nt],grid[:,1][t][:px])
    @expression(graph,pout[t=1:nt],grid[:,end][t][:px])
    @expression(graph,linepack[t = 1:nt],1/c1*sum(grid[t,x][:px]*dx for x in 1:nx-1))

    #Finite differencing.  Backward difference in time from t, Forward difference in space from x.
    @linkconstraint(graph, press[t = 2:nt, x = 1:nx - 1], ((grid[t,x][:px] - grid[t-1,x][:px])/dt
    + c1*(grid[t,x+1][:fx] - grid[t,x][:fx])/dx)*10 == 0 )  #~.1-.01
    @linkconstraint(graph, flow[t = 2:nt, x = 1:nx - 1],  ((grid[t,x][:fx] - grid[t-1,x][:fx])/dt)*10 ==
    (-c2*(grid[t,x+1][:px] - grid[t,x][:px])/dx - grid[t,x][:slack])*10) # ~.1-.01

    #initial steady state
    @linkconstraint(graph, ssflow[x = 1:nx-1], (grid[1,x+1][:fx] - grid[1,x][:fx])/10 == 0)   #~10
    @linkconstraint(graph, sspress[x = 1:nx-1], (-c2*(grid[1,x+1][:px] - grid[1,x][:px])/dx - grid[1,x][:slack])*10 == 0) #~ .1

    # Refill pipeline linepack at the end of time horizon
    @linkconstraint(graph,.001*linepack[end] >= .001*linepack[1])  #~100?
    return graph
end

# Define a function to create a subgraph for each compressor
function create_compressor_model(data,nt)
    # Create compressor OptiGraph
    graph = OptiGraph()
    # Create nodes at each time point
    @optinode(graph,time_nodes[1:nt])

    # Add variables, constraints, and objectives to each node
    for node in time_nodes
        @variable(node, 1 <= psuction <= 100, start = 60)
        @variable(node, 1 <= pdischarge <= 100, start = 60)
        @variable(node, 0 <= power <= 1000, start = 200)
        @variable(node, flow >= 0, start = 30)
        @variable(node, 1 <= boost <= 30,start = 1)
        @constraint(node, (pdischarge) == (psuction + boost))
        @NLconstraint(node, (power)*.01 == (c4*flow*((pdischarge/psuction)^om-1))*.01 )
        @objective(node, Min, 0.1*power*(dt/3600.0))
    end
    @expression(graph,fin[t=1:nt],time_nodes[t][:flow])
    @expression(graph,fout[t=1:nt],time_nodes[t][:flow])

    return graph
end
