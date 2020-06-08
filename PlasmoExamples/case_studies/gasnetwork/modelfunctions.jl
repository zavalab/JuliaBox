using Plasmo

z = 0.8                      # gas compressibility  - []
rhon = 0.72         		     # density of natural gas at normal conditions - [kg/m3]
R = 8314.0       			     # universal gas constant [J/kgmol-K]
M = 18.0    			         # gas molar mass of natural gas [kg/kgmol]
Tgas = 293.15      		     # reference temperature [K]
Cp = 2.34        		     # heat capacity @ constant pressure [kJ/kg-K]
Cv = 1.85        		     # heat capacity @ constant volume [kJ/kg-K]

#scaling factors for optimization
ffac=(1e+6*rhon)/(24*3600)                     # from MMSCM/day to kg/s
ffac2=(3600)/(1e+4*rhon)                       # from kg/s to scmx10^-4/hr
pfac=1e+5                                      # from bar to Pa
pfac2=1e-5                                     # from Pa to bar
dfac=1e-3                                      # from mm to m
lfac=1e+3                                      # from km to m
gam = Cp/Cv       		     	# expansion coefficient [-]
n_poly = gam
nu2 = gam*z*R*Tgas/M  			# gas speed of sound
om = (gam-1.0)/gam 		     	# aux constant [-]
c4 = (1/ffac2)*(Cp*Tgas)        #[kW/(scmx10-4/hr)]

function create_junction_model(data,nt)
    graph = ModelGraph()
    @optinode(graph,time_nodes[1:nt])
    n_demands = length(data[:demand_values])
    n_supplies = length(data[:supplies])
    for (i,node) in enumerate(time_nodes)
        @variable(node, data[:pmin] <= pressure <= data[:pmax], start = 60)
        @variable(node, 0 <= fgen[1:n_supplies] <= 200, start = 10)
        @variable(node, fdeliver[1:n_demands] >= 0, start = 10)
        @variable(node, fdemand[1:n_demands] >= 0, start = 10)

        @constraint(node,[d = 1:n_demands],fdeliver[d] <= fdemand[d])

        @expression(node, total_supplied, sum(fgen[s] for s = 1:n_supplies))
        @expression(node, total_delivered,sum(fdeliver[d] for d = 1:n_demands))
        @expression(node, total_delivercost,sum(1000*fdeliver[d] for d = 1:n_demands))

        if n_demands > 0
            @objective(node,Min,total_delivercost)
        end
    end
    return graph
end

function create_pipeline_model(data,nt,nx)
    #unpack data
    c1 = data[:c1]; c2 = data[:c2]; c3 = data[:c3]
    dx = data[:pipe_length] / (nx - 1)

    #Create pipeline modelgraph
    graph = ModelGraph()
    #Create grid of modelnodes
    @optinode(graph,grid[1:nt,1:nx])

    #Create variables on each node in the grid
    for node in grid
        @variable(node, 1 <= px <= 100,start = 10)
        @variable(node, 0 <= fx <= 100, start = 10)
        @variable(node,slack >= 0, start = 1)
        @NLnodeconstraint(node, slack*px - c3*fx*fx == 0)
    end

    #Setup boundary variables
    @expression(graph,fin[t=1:nt],grid[:,1][t][:fx])
    @expression(graph,fout[t=1:nt],grid[:,end][t][:fx])
    @expression(graph,pin[t=1:nt],grid[:,1][t][:px])
    @expression(graph,pout[t=1:nt],grid[:,end][t][:px])
    @expression(graph,linepack[t = 1:nt],1/c1*sum(grid[t,x][:px]*dx for x in 1:nx-1))

    #Finite differencing.  Backward difference in time from t, Forward difference in space from x.
    @linkconstraint(graph, press[t = 2:nt, x = 1:nx - 1], (grid[t,x][:px] - grid[t-1,x][:px])/dt
    + c1*(grid[t,x+1][:fx] - grid[t,x][:fx])/dx == 0 )
    @linkconstraint(graph, flow[t = 2:nt, x = 1:nx - 1],  (grid[t,x][:fx] - grid[t-1,x][:fx])/dt ==
    -c2*(grid[t,x+1][:px] - grid[t,x][:px])/dx - grid[t,x][:slack])

    #initial steady state
    @linkconstraint(graph, ssflow[x = 1:nx-1], grid[1,x+1][:fx] - grid[1,x][:fx] == 0)
    @linkconstraint(graph, sspress[x = 1:nx-1], -c2*(grid[1,x+1][:px] - grid[1,x][:px])/dx - grid[1,x][:slack] == 0)

    #Refill pipeline linepack
    @linkconstraint(graph,linepack[end] >= linepack[1])
    return graph
end

function create_compressor_model(data,nt)
    graph = ModelGraph()
    @optinode(graph,time_nodes[1:nt])
    for node in time_nodes
        @variable(node, 1 <= psuction <= 100, start = 60)   #This could also be inforced at the suction node
        @variable(node, 1 <= pdischarge <= 100, start = 60)
        @variable(node, 0 <= power <= 1000, start = 500)
        @variable(node, flow >= 0, start = 10)
        @variable(node, 1 <= boost <= 30,start = 1)
        @constraint(node, pdischarge == psuction + boost)
        @NLnodeconstraint(node, power == c4*flow*((pdischarge/psuction)^om-1) )
        @objective(node, Min, 0.1*power*(dt/3600.0))
    end
    @expression(graph,fin[t=1:nt],time_nodes[t][:flow])
    @expression(graph,fout[t=1:nt],time_nodes[t][:flow])

    return graph
end
